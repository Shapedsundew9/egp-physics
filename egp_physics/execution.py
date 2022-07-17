"""The GC Execution environment."""

from logging import getLogger, NullHandler, DEBUG, captureWarnings
from .ep_type import asstr
from .gc_graph import const_idx, ref_idx


_logger = getLogger(__name__)
_logger.addHandler(NullHandler())
captureWarnings(True)
_LOG_DEBUG = _logger.isEnabledFor(DEBUG)
_OVER_MAX = 1 << 64
_MASK = _OVER_MAX - 1

# Count the number of references to each imported object
# If it reaches 0 the import can be removed.
# Keys = import name: values = reference count
import_references = {}

class magic_mapping_list(list):
    """Limited use class used for mapping function argument templates to default values."""

    def __init__(self, name, template=True):
        """Vanilla empty list with a default for non-existant elements.

        Args
        ----
        name (str): See __getitem__()
        template (bool): Make default returned string a template (enclose in {})
        """
        super().__init__()
        self.name = name
        self.template = template

    def __repr__(self):
        """Default representation is just the template text."""
        return '{' + self.name + '}' if self.template else self.name

    def __getitem__(self, idx):
        """Return a default value if index does not exist.

        Otherwise works like a vanilla list.

        Args
        ----
        idx (int): Index into list.

        Returns
        -------
        element or (str): Default value is a string constructed from name.
        """
        try:
            retval = super().__getitem__(idx)
        except IndexError:
            retval = self.name + '[' + str(idx) + ']'
            if self.template:
                retval = '{' + retval + '}'
        return retval


class magic_mapping_dict(dict):
    """Limited use class used for mapping function argument templates to default values."""

    def __init__(self, template=True):
        """Vanilla dict with default for non-existent keys.

        Args
        ----
        template (bool): Make default returned string a template (enclose in {})
        """
        self.template = template

    def __missing__(self, key):
        """Return a default magic_mapping_list if key does not exist.

        Otherwise works like a vanilla dict.

        Args
        ----
        key (hashable): Key into dict.

        Returns
        -------
        element or (magic_mapping_list): Default magic_mapping_list is initialised with the key.
        """
        return magic_mapping_list(key, self.template)


# GMS is a parameter of pGC codons which are dynamically created.
_GMS = None
def set_gms(gms):
    """pGC's may need to find GC's in the GMS."""
    global _GMS
    _GMS = gms


def write_args(arg_list):
    """Create an argument string.

    Args
    ----
    (list(str)): From arg_list() function.

    Returns
    -------
    (str): Executable tuple definition e.g.:
        'tuple()' # for no arguments
        '(arg1, arg2, arg3)'
        '(arg1,)'
    """
    if arg_list:
        if len(arg_list) > 1:
            return '(' + ','.join(arg_list) + ')'
        return '(' + arg_list[0] + ',)'
    return 'tuple()'


def arg_list(iab, c, ref_str, template=False):
    """Create an argument list

    All functions take a single tuple argument.
    The arguments can be templated using {} string fomratting to allow them
    to be substituted when inlining.

    Args
    ----
    iab (graph[iab]): I, A or B row definition froma GC application graph
    c (graph['C']): Row C definition or None if row C is not defined
    ref_str (str): String representation of the GC reference.
    template (bool): If True wrap variable names in {}

    Returns
    -------
    (list(str)): Argument string e.g. ['{a_023c9f0d[0]}', '{i_894dea21[4]}', 'constant_string',]
    """
    args = []
    if iab is not None:
        for arg in iab:
            row = arg[ref_idx.ROW]
            if row == 'C':
                args.append(str(c[arg[const_idx.VALUE]]))
            else:
                args.append('{' + row.lower() + '_' + ref_str + '[' + str(arg[ref_idx.INDEX]) + ']}')
            if template:
                args[-1] = '{{' + args[-1] + '}}'
    return args


def code_body(gc, gpc, max_depth, marker):
    """Create a string representing the code body of the callable for gc.

    The code body naming convention is each rows output is labelled with the
    lower case row letter followed by f"_{gc['ref']:016x}".

    Args
    ----
    gc (gGC): The gGC for which to define the callable
    gpc (GPC): The Gene Pool Cache from which to retrieve any sub-GC's
    max_depth (int): The maximum number code lines in the function
    marker (str): Prefix for function names.

    Returns
    -------
    (str) The code body string for gc.
    (bool) True if GCA was inlined
    (bool) True if GCB was inlined
    """
    # TODO: Optimisations
    #   1. Map inputs and outputs to save copying tuples
    #   2. Do not create 1 element tuples within a code body
    #   3. Do not pass empty tuples as zero input inputs

    graph = gc['igraph'].app_graph
    ref_str = f"{(_OVER_MAX + gc['ref']) & _MASK:016x}"
    space = max_depth - 1 # For the return statement
    cb_dict = {'A': '', 'B': '', 'O': ''}
    referenced = []

    _logger.debug(f'ref = {ref_str}')
    _logger.debug(f'graph = {graph}')

    # Codon special case
    if gc.get('meta_data', None) is not None and 'function' in gc['meta_data']:
        format_dict = {'c' + str(i): v for i, v in enumerate(graph['C'])} if 'C' in graph else {}
        format_dict.update({'i' + str(i): f'{{i_{ref_str}[{i}]}}' for i in range(len(gc['inputs']))})
        code = gc['meta_data']['function']['python3']['0']
        if 'code' in code: cb_dict['A'] = "\t" + code['code'].format_map(format_dict) + "\n"
        formatted_inline = f'\t{{o_{ref_str}}} = (' + code['inline'].format_map(format_dict) + ',)\n'
        cb_dict['A'] += formatted_inline
        #if _LOG_DEBUG: cb_dict['O'] += f"\t_logger.debug(f\"{ref_str}: {formatted_inline} = {{{{{formatted_inline}}}}}\")\n"
    else:
        # Get the code body of GCx. Create them if need be.
        cb_dict['O'] = f"\t{{o_{ref_str}}} = {write_args(arg_list(graph.get('O'), graph.get('C'), ref_str))}\n"
        for gcx_ref_key, row in filter(lambda x: gc[x[0]] is not None, (('gca_ref', 'A'), ('gcb_ref', 'B'))):
            # If GCA exists check there is enough room without exceeding max_depth
            # not forgetting GCB will need at least 1 line if it exists.
            gcx_ref = gc[gcx_ref_key]
            gcx = gpc[gcx_ref]
            gcx_ref_str = f'{(_OVER_MAX + gcx_ref) & _MASK:016x}'
            gcx_lc = gcx['cb'].count('\n')
            gcx_graph = gcx['igraph'].app_graph
            gcx_i = gcx['num_inputs'] > 0
            gcx_o = gcx['num_outputs'] > 0
            gcx_lines = sum(x in gcx_graph for x in 'ABC') + gcx_i + gcx_o
            _logger.debug(f"row = {row}")
            _logger.debug(f"gc['gc{row}_ref'] = {gcx_ref_str}")
            if space > gcx_lines: # Can inline at least the sub-GC
                if gcx_i:
                    cb_dict[row] += f"\t{{{'i_' + gcx_ref_str}}} = {write_args(arg_list(graph.get(row), graph.get('C'), ref_str))}\n"
                if gcx_lc < space: # Can inline the sub-GC's inlined code body
                    space -= gcx_lc
                    callable_imports(gcx) # gcx may have imports inlined.
                    cb_dict[row] += gcx['cb'].split('\treturn')[0] # Remove output row
                else: # Just the sub-GC code body
                    space -= gcx_lines
                    for sgcx_ref, srow in filter(lambda x: gcx[x[0]] is not None, (('gca_ref', 'A'), ('gcb_ref', 'B'))):
                        sgcx_ref_str = f'{(_OVER_MAX + gcx[sgcx_ref]) & _MASK:016x}'
                        cb_dict[row] += f"\t{{{srow.lower()}_{gcx_ref_str}}} = {marker}{sgcx_ref_str}"
                        cb_dict[row] += f"({write_args(arg_list(gcx_graph.get(srow), gcx_graph.get('C'), gcx_ref_str))})\n"
                        referenced.append(gcx[sgcx_ref])
                    cb_dict[row] += f"\t{{{'o_' + gcx_ref_str}}} = {write_args(arg_list(gcx_graph.get(srow), gcx_graph.get('C'), gcx_ref_str))}\n"
                if gcx_o:
                    cb_dict[row] += f"\t{{{row.lower() + '_' + ref_str}}} = {{{'o_' + gcx_ref_str}}}\n"

            else: # Direct call - take 1 line
                space -= 1
                cb_dict[row] += f"\t{{{row.lower()}_{ref_str}}} = {marker}{gcx_ref_str}"
                cb_dict[row] += f"({write_args(arg_list(graph.get(row), graph.get('C'), ref_str))})\n"
                referenced.append(gcx_ref)

            _logger.debug(f'cb_dict[{row}]: {cb_dict[row]}')

    cb_dict['O'] += f"\treturn o_{ref_str}\n"
    cb = (cb_dict['A'] + cb_dict['B'] + cb_dict['O'])
    if _LOG_DEBUG:
        _logger.debug(f'GC graph rows: {list(graph.keys())}')
        _logger.debug(f'Referenced: {[f"{(_OVER_MAX + r) & _MASK:016x}" for r in referenced]}')
        _logger.debug(f'cb: {cb}')

    return cb, referenced


def callable_string(gc, marker):
    """Create a string that can be exec()'s into a callable function.

    Callables are optimised to reduce call depth. sub-GC's are inlined up to
    a maximum function size of max_depth codons.

    Args
    ----
    gc (gGC): The gGC for which to define the callable
    marker (str): Prefix for function names.

    Returns
    -------
    (str) The string defining the callable function for gc
    """
    # TODO: There are various levels of debug to put in here.
    ref_str = f"{(_OVER_MAX + gc['ref']) & _MASK:016x}"
    string = "# ref: " + str(gc['ref']) + "\n"
    string += "# i = (" + ", ".join((asstr(i) for i in gc['igraph'].input_if())) + ")\n"
    string += "def " + marker + ref_str
    string += f'(i_{ref_str}):\n' if len(gc['inputs']) else "():\n"
    if _LOG_DEBUG and len(gc['inputs']):
        string += f"\t_logger.debug(f\"{ref_str}: i_{ref_str} = {{i_{ref_str}}}\")\n"
    graph = gc['igraph'].app_graph
    if gc.get('meta_data', None) is None or not 'function' in gc['meta_data']:
        string += gc['cb'].format_map(magic_mapping_dict(False))
    else:
        format_dict = {'c' + str(i): v for i, v in enumerate(graph['C'])} if 'C' in graph else {}
        format_dict.update({f'i{i}': f'i_{ref_str}[{i}]' for i in range(len(gc['inputs']))})
        code = gc['meta_data']['function']['python3']['0']
        if 'code' in code: string += "\t" + code['code'].format_map(format_dict) + "\n"
        formatted_inline = code['inline'].format_map(format_dict)
        if _LOG_DEBUG:
            string += f"\t_logger.debug(f\"{ref_str}: {formatted_inline} = {{{formatted_inline}}}\")\n"
        string += "\treturn (" + formatted_inline + ",)\n\n\n"
    if _LOG_DEBUG: _logger.debug(f"Callable string created:\n{string}")
    return string


def create_callable(gc, gpc, max_depth=20, marker='ref_'):
    """Create a callable function from a gGC in the global namespace.

    Since functions are added to the global namespace multiple GP instances
    may create the same function.

    Callables are optimised to reduce call depth. sub-GC's are inlined up to
    a maximum function size of max_depth code lines.

    Any imports specified are added to the global imports.

    Args
    ----
    gc (gGC): The gGC for which to define the callable
    gpc (GPC): The Gene Pool Cache from which to retrieve any sub-GC's
    max_depth (int): The maximum number of codons in the function
    marker (str): Prefix for function names.

    Returns
    -------
    (callable) The callable function for gc
    """
    # The GC depth could be very deep and cause recursion issues so this is a looping implementation.
    # Walk the GC graph to look for GC's that do not yet have a code body 'cb' defined.
    # If a cb is defined then all sub-GC's have cb's defined and callables as necessary.
    if gc['cb'] is None:
        path = [gc]
        inspect = [gc]
        while inspect:
            g = inspect.pop()
            for ref in ('gca_ref', 'gcb_ref'):
                gcx_ref = g[ref]
                if gcx_ref is not None:   # No path
                    gcx = gpc[gcx_ref]
                    if gcx['cb'] is None:  # Path not yet processed
                        path.append(gcx)
                        inspect.append(gcx)

        # Path contains all the executables that need to be created in reverse order
        # Because GC's can be inlined there is a difference between having a defined code block and
        # having a defined callable.
        for sgc in reversed(path):
            sgc['cb'], referenced = code_body(sgc, gpc, max_depth, marker)
            #if sgc.get('meta_data', None) is not None and 'function' in sgc['meta_data']:
            #    sgc['callable'] = create_exec(sgc, marker)
            #else:
            for ref in referenced:
                rgc = gpc[ref]
                if rgc['callable'] is None:
                    rgc['callable'] = create_exec(rgc, marker)

    # Make the top level executable
    # Note that the executable is not wrapped at this point as it is likely to become a sub-GC
    # in the future and the wrapper will just be overhead.
    if gc['callable'] is None:
        gc['callable'] = create_exec(gc, marker)

    # Wrap the returned function
    return exec_wrapper(gc['callable'])


def callable_imports(gc):
    """Import imports into the global namespace.

    Args
    ----
    gc (gGC): The gGC for which to define the callable
    """
    global import_references
    if gc.get('meta_data', None) is not None and 'function' in gc['meta_data']:
        python = gc['meta_data']['function']['python3']['0']
        if 'imports' in python:
            if _LOG_DEBUG:
                _logger.debug(f"Importing: {python}")
            for impt in python['imports']:
                if impt['name'] not in globals():
                    string = "from {module} import {object} as {name}\n".format_map(impt)
                    if _LOG_DEBUG: _logger.debug(f"New imports executable: {string}")
                    exec(string, globals())
                    import_references[impt['name']] = 1
                else:
                    import_references[impt['name']] += 1


def create_exec(gc, marker):
    """Create an executable from a callable string in the global namespace.

    Args
    ----
    gc (gGC): The gGC for which to define the callable
    marker (str): Prefix for function names.

    Returns
    -------
    (callable or None) The string defining the callable function for gc
    """
    global_name = f"{marker}{(_OVER_MAX + gc['ref']) & _MASK:016x}"
    if global_name not in globals():
        callable_imports(gc)
        exec(callable_string(gc, marker), globals())
        return globals()[global_name]


def remove_callable(gc, marker='ref_'):
    # NOTE: You cannot unload an imported module.
    # This builds up cruft. How much of an issue is this?
    # Added a cruft counter. Maybe we do a full shutdown and self restart
    # periodically to remove cruft? Generally I do not approve as we should
    # tidy up after ourselves but in this case we are forced too.
    global import_references
    gc['callable'] = None
    name = f"{marker}{(_OVER_MAX + gc['ref']) & _MASK:016x}"
    if name in globals():
        if _LOG_DEBUG:
            _logger.debug(f'Deleting {name} from GP execution environment.')
        del globals()[name]
    if 'meta_data' in gc and isinstance(gc['meta_data'], dict) and 'function' in gc['meta_data']:
        python = gc['meta_data']['function']['python3']['0']
        if 'imports' in python:
            print(import_references)
            print(python['imports'])
            for impt in python['imports']:
                import_references[impt['name']] -= 1
            cruft = [name for name, count in import_references.items() if not count]
            if cruft:
                _logger.warning(f"{len(cruft)} imports are no longer used: {cruft}")


def exec_wrapper(func):
    """Wrap func with a try...catch.

    This should only be used for the top level GC as it is overhead.

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

            # Catches issues with destroying GC exec functions e.g. 'name {marker}{gc['ref']:08x} is not defined'
            # TODO: Are there legitimate cases?
            assert str(e)[-14:] != 'is not defined'
            retval = None
        return retval
    return wrapped