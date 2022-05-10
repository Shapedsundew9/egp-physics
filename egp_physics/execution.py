"""The GC Execution environment."""

from logging import getLogger, NullHandler, DEBUG, captureWarnings
from .ep_type import asstr

_logger = getLogger(__name__)
_logger.addHandler(NullHandler())
captureWarnings(True)
_LOG_DEBUG = _logger.isEnabledFor(DEBUG)

# Count the number of references to each imported object
# If it reaches 0 the import can be removed.
# Keys = import name: values = reference count
import_references = {}

# _GMS is a parameter of pGC codons which are dynamically created.
_GMS = None
def set_gms(gms):
    """pGC's may need to find GC's in the GMS."""
    global _GMS
    _GMS = gms

def name_func(ref):
    ref_str = ref.to_bytes(8,'big', signed=True).hex()
    return f'ref_{ref_str}'

def callable_exists(ref):
    return name_func(ref) in globals()

def write_arg(iab, c):
    return "(" + ", ".join([str(c[arg[1]]) if arg[0] == 'C' else arg[0].lower() + "[" + str(arg[1]) + "]" for arg in iab]) + ",)"

def callable_string(gc):
    # TODO: Add a depth parameter that generates functions that have less codons as single functions.
    string = "# ref: " + str(gc['ref']) + "\n"
    string += "# i = (" + ", ".join((asstr(i) for i in gc['igraph'].input_if())) + ")\n"
    string += "def " + name_func(gc['ref'])
    string += "(i):\n" if len(gc['inputs']) else "():\n"
    if _LOG_DEBUG and len(gc['inputs']): string += f"\t_logger.debug(f\"{gc['ref']}: i = {{i}}\")\n"
    graph = gc['igraph'].app_graph
    if gc.get('meta_data', None) is None or not 'function' in gc['meta_data']:
        c = graph['C'] if 'C' in graph else tuple()
        if gc.get('gca_ref', None) is not None:
            string += "\ta = " + name_func(gc['gca_ref'])
            string += "(" + write_arg(graph['A'], c) + ")\n" if 'A' in gc['graph'] else "()\n"
            if _LOG_DEBUG: string += f"\t_logger.debug(f\"{gc['ref']}: a = {{a}}\")\n"
        if gc.get('gcb_ref', None) is not None:
            string += "\tb = " + name_func(gc['gcb_ref'])
            string += "(" + write_arg(graph['B'], c) + ")\n" if 'B' in gc['graph'] else "()\n"
            if _LOG_DEBUG: string += f"\t_logger.debug(f\"{gc['ref']}: b = {{b}}\")\n"
        retval = write_arg(graph['O'], c) if 'O' in graph else 'None'
        if _LOG_DEBUG: string += f"\t_logger.debug(f\"{gc['ref']}: return {retval} = {{{retval}}}\")\n"
        string += "\treturn " + retval + "\n\n\n"
    else:
        format_dict = {'c' + str(i): v for i, v in enumerate(graph['C'])} if 'C' in graph else {}
        format_dict.update({'i' + str(i): 'i[{}]'.format(i) for i in range(len(gc['inputs']))})
        code = gc['meta_data']['function']['python3']['0']
        if 'code' in code: string += "\t" + code['code'].format(**format_dict) + "\n"
        formated_inline = code['inline'].format(**format_dict)
        if _LOG_DEBUG: string += f"\t_logger.debug(f\"{gc['ref']}: {formated_inline} = {{{formated_inline}}}\")\n"
        string += "\treturn (" + formated_inline + ",)\n\n\n"
    if _LOG_DEBUG: _logger.debug(f"Callable string created:\n{string}")
    return string

def create_callable(gc):
    """Create a callable function from a gGC in the global namespace.

    Since functions are added to the global namespace multiple GP instances
    may create the same function.
    """
    global import_references
    global_name = name_func(gc['ref'])
    if global_name not in globals():

        # Import imports into the global namespace
        if gc.get('meta_data', None) is not None and 'function' in gc['meta_data']:
            python = gc['meta_data']['function']['python3']['0']
            if 'imports' in python:
                if _LOG_DEBUG:
                    _logger.debug(f"Importing: {python}")
                for impt in python['imports']:
                    if impt['name'] not in globals():
                        string = "from {module} import {object} as {name}\n".format(**impt)
                        if _LOG_DEBUG: _logger.debug(f"New imports executable: {string}")
                        exec(string, globals())
                        import_references[impt['name']] = 1
                    else:
                        import_references[impt['name']] += 1

        exec(callable_string(gc), globals())
        return exec_wrapper(globals()[global_name])
    if _LOG_DEBUG:
        _logger.warning(f'Function {global_name}() already exists!')
        _logger.debug(f'gc creating existing exec function is {gc}.')
        assert 'exec' in gc
        assert gc['exec'] is not None
        assert False
    return gc['exec']

def remove_callable(gc):
    # NOTE: You cannot unload an imported module.
    # This builds up cruft. How much of an issue is this?
    # Added a cruft counter. Maybe we do a full shutdown and self restart
    # periodically to remove cruft? Generally I do not approve as we should
    # tidy up after ourselves but in this case we are forced too.
    gc['exec'] = None
    name = name_func(gc['ref'])
    if name in globals():
        if _LOG_DEBUG:
            _logger.debug(f'Deleting {name} from GP execution environment.')
        del globals()[name]
    if 'meta_data' in gc and isinstance(gc['meta_data'], dict) and 'function' in gc['meta_data']:
        python = gc['meta_data']['function']['python3']['0']
        if 'imports' in python:
            for impt in python['imports']:
                import_references[impt['name']] -= 1
            cruft = [name for name, count in import_references.items() if not count]
            if cruft:
                _logger.warning(f"{len(cruft)} imports are no longer used: {cruft}")


def exec_wrapper(func):
    """Wrap func with a try...catch.

    Args
    ----
    func (callable): Function to wrap with a try..catch

    Returns
    -------
    (callable): Wrapped func
    """
    def wrapped(*args, **kwargs):
        try:
            retval = func(*args, **kwargs)
        except Exception as e:
            _logger.debug(f'Exception occured in execution wrapper for {func.__name__}: {e}')
            retval = (None,)
        return retval
    return wrapped


def lazy_exec(inputs, gc):
    """Create a default executable.

    When executed lazy_exec() replaces itself as the GC callable
    with the actual callable wrapped to catch any exceptions.
    """
    # FIXME: This is flawed because GCA/GCB may not exist. Needs to be recursive
    # but that requires access to the gene pool it is in.
    gc['exec'] = exec_wrapper(create_callable(gc))
    return gc['exec'](inputs)