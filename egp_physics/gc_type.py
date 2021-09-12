"""Defines the *GC types.

GC types come in several flavours related to their stage in the
GC lifecycle.

TODO: A diagram

GC types are principly dictionaries with some validation to support
development and encapsulation of the storage (which may need to be
imlemented differently in the future to reduce resident RAM costs)

The overall GC definition is stored in a cerberus schema and includes
some meta data to differenciate the types.

Type instances are 'only' self consistent. For population consistency
checks see the relevant collections e.g. genomic_library, gene_pool
"""

from copy import deepcopy

from .ep_type import asint, vtype
from .gc_graph import gc_graph
from .generic_validator import SCHEMA, generic_validator, random_reference

# GC types
__GC = '_'
_EGC = 'e'
_PGC = 'p'
_AGC = 'a'
_XGC = 'x'
_MGC = 'm'


# GitHub Markdown Emoji's
_DEFINABLE_MD = ":large_blue_circle:"
_REQUIRED_MD = ":black_circle:"


def _get_schema(t):
    """Logic for extracting the specific GC type schema from the general GC schema.

    Sets the nullable field which may differ between GC types.

    Args
    ----
    t (str): A valid GC type letter
    t_key_set (iterable(str)): Fields

    Returns
    -------
    (dict): A cerberus schema for GC type t
    """
    schema = {k: deepcopy(v) for k, v in filter(lambda x: t in x[1]['meta']['types'], SCHEMA.items())}
    for v in schema.values():
        v['required'] = True
        v['nullable'] = t in v['meta'].get('nullable_types', tuple())
        if t not in v['meta'].get('default_types', tuple()):
            if 'default' in v:
                del v['default']
            if 'default_setter' in v:
                del v['default_setter']
        if t not in v['meta'].get('check_types', tuple()) and 'check_with' in v:
            del v['check_with']
    return schema


def interface_definition(xputs, vt=vtype.TYPE_OBJECT):
    """Create an interface definition from xputs.

    Used to define the inputs or outputs of a GC from an iterable
    of types, objects or EP definitions.

    Args
    ----
    xputs (iterable(object)): Object is of the type defined by vt.
    vt (vtype): The interpretation of the object. See definition of vtype.

    Returns
    -------
    tuple(list(ep_type_int), list(ep_type_int), list(int)): A list of the xputs as EP type in value
        format, a list of the EP types in xputs in value format in ascending order and a list of
        indexes into it in the order of xputs.
    """
    xput_eps = tuple((asint(x, vt) for x in xputs))
    xput_types = sorted(set(xput_eps))
    return xput_eps, xput_types, [xput_types.index(x) for x in xput_eps]


class _GC(dict):
    """Abstract base class for genetic code (GC) types.

    Validity is in the context of a steady state GC i.e. an INVALID field will not cause
    a validate() failure but the GC will not be stable (valid).
    """

    validator = generic_validator(SCHEMA)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setdefault('_ref', random_reference())

    def validate(self):
        """Validate all required key:value pairs are correct."""
        if not self.validator.validate(self):
            raise ValueError(f"Validation FAILED with:\n{self.validator.error_str()}")


class eGC(_GC):
    """Embryonic GC type.

    Embryonic GC's have the minimal set of fields necessary to form a GC but the values
    are not necessarily valid. A Minimal GC (mGC) has the same fields as an
    embryonic GC but provides the valid minimal field values guarantee.
    """

    validator = generic_validator(_get_schema(_EGC), allow_unknown=True)

    def __init__(self, gc={}, inputs=tuple(), outputs=tuple(), vt=vtype.OBJECT, sv=True):
        """Construct.

        Args
        ----
        gc (a _GC dervived object): GC to ensure is eGC compliant.
        inputs (iterable(object)): GC inputs. Object is of the type defined by vt.
        outputs (iterable(object)): GC outputs. Object is of the type defined by vt.
        vt (vtype): The interpretation of the object. See vtype definition.
        sv (bool): Suppress validation. If True the eGC will not be validated on construction.
        """
        # TODO: Consider lazy loading fields
        super().__init__(gc)
        graph_inputs, self['input_types'], self['inputs'] = interface_definition(inputs, vt)
        graph_outputs, self['output_types'], self['outputs'] = interface_definition(outputs, vt)
        self.setdefault('_gca', None)
        self.setdefault('_gcb', None)
        if '_graph' not in self:
            _graph = gc_graph()
            _graph.add_inputs(graph_inputs)
            _graph.add_outputs(graph_outputs)
            self['_graph'] = _graph
        if not sv:
            self.validate()


class mGC(_GC):
    """Minimal GC type.

    Minimal GC's have the minimal set of fields with valid values necessary to form a GC.
    A Minimal GC (mGC) has the same fields as an embryonic GC but provides the valid values guarantee.
    """

    validator = generic_validator(_get_schema(_MGC), allow_unknown=True)

    def __init__(self, gc={}, _graph=gc_graph(), _gca=None, _gcb=None, sv=True):
        """Construct.

        Args
        ----
        gc (a _GC dervived object): GC to ensure is mGC compliant.
        _gca (int or None): gca reference.
        _gcb (int or None): gcb reference.
        sv (bool): Suppress validation. If True the eGC will not be validated on construction.
        """
        super().__init__(gc)
        self.setdefault('_graph', _graph)
        self.setdefault('_gca', _gca)
        self.setdefault('_gcb', _gcb)
        if 'inputs' not in self:
            inputs = self['_graph'].input_if()
            outputs = self['_graph'].output_if()
            _, self['input_types'], self['inputs'] = interface_definition(inputs, vtype.EP_TYPE_INT)
            _, self['output_types'], self['outputs'] = interface_definition(outputs, vtype.EP_TYPE_INT)
        if not sv:
            self.validate()


def _md_string(gct, key):
    """Generate the GitHub MD string relevant to the key for gct."""
    if key not in gct.validator.schema:
        return ""
    if 'default' in gct.validator.schema[key] or 'default_setter' in gct.validator.schema[key]:
        return _DEFINABLE_MD
    return _REQUIRED_MD


def md_table():
    """Create a GitHub Markdown table showing the requirements of each field for each GC type."""
    gcts = (
        _GC(sv=True),
        eGC(sv=True),
        mGC(sv=True)
    )
    with open('gc_type_table.md', 'w') as file_ptr:
        file_ptr.write("GC Type Field Requirements\n")
        file_ptr.write("==========================\n\n")
        file_ptr.write(_DEFINABLE_MD + ': Defined if not set, ' + _REQUIRED_MD + ': Required.\n\n')
        file_ptr.write('| Field | ' + ' | '.join((x.__class__.__qualname__ for x in gcts)) + ' |\n')
        file_ptr.write('| --- | ' + ' | '.join(('---' for _ in gcts)) + ' |\n')
        for k in sorted(SCHEMA.keys()):
            file_ptr.write(f'| {k} | ' + ' | '.join((_md_string(gct, k) for gct in gcts)) + ' |\n')
