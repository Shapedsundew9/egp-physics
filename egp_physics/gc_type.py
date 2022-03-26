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

from copy import copy, deepcopy
from hashlib import blake2b
from logging import DEBUG, NullHandler, getLogger

from .ep_type import asint, vtype
from .gc_graph import gc_graph
from .execution import create_callable, exec_wrapper
from .generic_validator import SCHEMA, generic_validator, random_reference

_logger = getLogger(__name__)
_logger.addHandler(NullHandler())
_LOG_DEBUG = _logger.isEnabledFor(DEBUG)

# GC types
__GC = '_'
_EGC = 'e'
_PGC = 'p'
_AGC = 'a'
_XGC = 'x'
_MGC = 'm'
_NGC = 'n'
_GGC = 'g'


# GitHub Markdown Emoji's
_DEFINABLE_MD = ":large_blue_circle:"
_REQUIRED_MD = ":black_circle:"


# Evolve a pGC after this many 'uses'.
# MUST be a power of 2
M_CONSTANT = 1 << 4
M_MASK = M_CONSTANT - 1
NUM_PGC_LAYERS = 16
# With M_CONSTANT = 16 & NUM_PGC_LAYERS = 16 it will take 16**16 (== 2**64 == 18446744073709551616)
# population individual evolutions to require a 17th layer (and that is assuming all PGC's are
# children of the one in the 16th layer). Thats about 5.8 billion evolutions per second for
# 100 years. A million super fast cores doing 5.8 million per second...only an outside chance
# of hitting the limit if Erasmus becomes a global phenomenon and is not rewritten! Sensibly future proofed.


# FIXME: This is duplicated in egp_physics.gc_type. Consider creating a seperate module of
# field definitions.
# PROPERTIES must define the bit position of all the properties listed in
# the "properties" field of the entry_format.json definition.
PROPERTIES = {
    "extended": 1 << 0,
    "constant": 1 << 1,
    "conditional": 1 << 2,
    "deterministic": 1 << 3,
    "memory_modify": 1 << 4,
    "object_modify": 1 << 5,
    "physical": 1 << 6,
    "arithmetic": 1 << 16,
    "logical": 1 << 17,
    "bitwise": 1 << 18,
    "boolean": 1 << 19,
    "sequence": 1 << 20
}
PHYSICAL_PROPERTY = PROPERTIES['physical']
LAYER_COLUMNS = (
    "evolvability",
    "fitness",
    "e_count",
    "f_count",
    "if"
)
LAYER_COLUMNS_RESET = {
    "e_count": 1,
    "f_count": 1
}


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
    return xput_eps, xput_types, bytes([xput_types.index(x) for x in xput_eps])


def unordered_interface_hash(input_eps, output_eps):
    """Create a 64-bit hash of the population interface definition.

    The interface hash is order agnostic i.e.

    (float, int, float) has the same hash as (float, float, int) has
    the same hash as (int, float, float).

    Args
    ----
    input_eps (iterable(int)): Iterable of input EP types.
    output_eps (iterable(int)): Iterable of output EP types.

    Returns
    -------
    (int): 64 bit hash as a signed 64 bit int.
    """
    h = blake2b(digest_size=8)
    for i in sorted(input_eps):
        h.update(i.to_bytes(2, 'big'))
    for o in sorted(output_eps):
        h.update(o.to_bytes(2, 'big'))
    a = int.from_bytes(h.digest(), 'big')
    return (0x7FFFFFFFFFFFFFFF & a) - (a & (1 << 63))


def ordered_interface_hash(input_types, output_types, inputs, outputs):
    """Create a 64-bit hash of the population interface definition.

    The interface hash is specific to the order and type in the inputs
    and outputs. This is important in population individuals.

    Args
    ----
    input_types (iterable(int)): Iterable of input EP types in ascending order.
    output_types (iterable(int)): Iterable of output EP types in ascending order.
    inputs (bytes): Indices into input_types for the input parameters.
    outputs (bytes): Indices into output_types for the input parameters.

    Returns
    -------
    (int): 64 bit hash as a signed 64 bit int.
    """
    h = blake2b(digest_size=8)
    for i in input_types:
        h.update(i.to_bytes(2, 'big'))
    for o in sorted(output_types):
        h.update(o.to_bytes(2, 'big'))
    h.update(inputs)
    h.update(outputs)
    a = int.from_bytes(h.digest(), 'big')
    return (0x7FFFFFFFFFFFFFFF & a) - (a & (1 << 63))


def ref_from_sig(sig):
    """Create a reference from a signature.

    Args
    ----
    sig (bytes): GC signature bytes object

    Returns
    -------
    (int): A 64 bit reference.
    """
    a = int.from_bytes(sig[:8], 'little')
    return (0x7FFFFFFFFFFFFFFF & a) - (a & (1 << 63))


def is_pgc(gc):
    """Determine if a GC is a PGC.

    Args
    ----
    gc(dict-like): A GC dict-like object.

    Returns
    -------
    (bool): True if gc is a pGC else False
    """
    if _LOG_DEBUG:
        # More juicy test for consistency
        # TODO: More conditions can be added
        it = gc.get('input_types', [])
        i = gc.get('inputs', [])
        ot = gc.get('output_types', [])
        o = gc.get('outputs', [])
        pgc_inputs = bool(it) and it[0] == asint('egp_physics.gc_type_gGC') and len(i) == 1
        pgc_outputs = bool(ot) and ot[0] == asint('egp_physics.gc_type_gGC') and len(o) == 1
        check = (pgc_inputs and pgc_outputs) == (gc.get('pgc_fitness', None) is not None)
        if not check:
            ValueError(f"PGC is not a PGC!: {gc['ref']}\n\t{pgc_inputs}, {pgc_outputs}, {gc.get('pgc_fitness', None)},"
                          f" {(pgc_inputs and pgc_outputs)}, {(gc.get('pgc_fitness', None) is not None)}")
    return gc.get('pgc_fitness', None) is not None


class _GC(dict):
    """Abstract base class for genetic code (GC) types.

    Validity is in the context of a steady state GC i.e. an INVALID field will not cause
    a validate() failure but the GC will not be stable (valid).

    All derived classes mutate the the supplied dict-like object to make
    it a valid _GC type. The dict is NOT copied.
    """

    validator = generic_validator(SCHEMA)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if 'ref' not in self:
            self['ref'] = self._ref_from_sig('signature', random_ref=True)

    def _ref_from_sig(self, s, random_ref=False):
        """Make a reference from the first 8 bytes of the signature if it exists."""
        if s not in self:
            if random_ref:
                return random_reference()
            return None
        elif self[s] is not None:
            return ref_from_sig(self[s])
        return None

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

    def __init__(self, gc={}, inputs=None, outputs=None, vt=vtype.OBJECT, sv=True):
        """Construct.

        NOTE: gc will be modified

        Args
        ----
        gc (a dict-like object): GC to ensure is eGC compliant.
        inputs (iterable(object)): GC inputs. Object is of the type defined by vt.
        outputs (iterable(object)): GC outputs. Object is of the type defined by vt.
        vt (vtype): The interpretation of the object. See vtype definition.
        sv (bool): Suppress validation. If True the eGC will not be validated on construction.
        """
        super().__init__(gc)

        if inputs is not None:
            graph_inputs, self['input_types'], self['inputs'] = interface_definition(inputs, vt)
        else:
            graph_inputs = []
        if outputs is not None:
            graph_outputs, self['output_types'], self['outputs'] = interface_definition(outputs, vt)
        else:
            graph_outputs = []
        self.setdefault('gca_ref', self._ref_from_sig('gca'))
        self.setdefault('gcb_ref', self._ref_from_sig('gcb'))
        self.setdefault('generation', 0)
        self['modified'] = True
        if 'igraph' not in self:
            if 'graph' in self:
                igraph = gc_graph(self['graph'])
                graph_inputs = igraph.input_if()
                graph_outputs = igraph.output_if()
            else:
                igraph = gc_graph()
                igraph.add_inputs(graph_inputs)
                igraph.add_outputs(graph_outputs)
                self['graph'] = igraph.application_graph()
            self['igraph'] = igraph
        elif 'graph' not in self:
            self['graph'] = self['igraph'].application_graph()

        if not sv:
            self.validate()


class mGC(_GC):
    """Minimal GC type.

    Minimal GC's have the minimal set of fields with valid values necessary to form a GC.
    A Minimal GC (mGC) has the same fields as an embryonic GC but provides the valid values guarantee.
    """

    validator = generic_validator(_get_schema(_MGC), allow_unknown=True)

    def __init__(self, gc={}, igraph=gc_graph(), sv=True):
        """Construct.

        gc combined with igraph must be in a steady state.
        if gc['igraph'] exists igraph will be ignored.

        Args
        ----
        gc (a _GC dervived object): GC to ensure is mGC compliant.
        gca_ref (int or None): gca reference.
        gcb_ref (int or None): gcb reference.
        sv (bool): Suppress validation. If True the mGC will not be validated on construction.
        """
        super().__init__(gc)
        self.setdefault('igraph', igraph)
        self.setdefault('gca_ref', self._ref_from_sig('gca'))
        self.setdefault('gcb_ref', self._ref_from_sig('gcb'))
        self.setdefault('generation', 0)
        if 'inputs' not in self:
            inputs = self['igraph'].input_if()
            outputs = self['igraph'].output_if()
            _, self['input_types'], self['inputs'] = interface_definition(inputs, vtype.EP_TYPE_INT)
            _, self['output_types'], self['outputs'] = interface_definition(outputs, vtype.EP_TYPE_INT)
        if not sv:
            self.validate()


class gGC(_GC):
    """Gene Pool GC type.

    Gene pool GC types hold a lot of transient data.
    """

    validator = generic_validator(_get_schema(_GGC), allow_unknown=True)
    higher_layer_cols = tuple((col for col in filter(lambda x: x[0] == '_', validator.schema.keys())))

    def __init__(self, gc={}, modified=True, population=None, sv=True):
        """Construct.

        Ensure all fields are defined as required for the Gene Pool.

        Args
        ----
        gc (a _GC dervived object): GC to ensure is eGC compliant.
        sv (bool): Suppress validation. If True the eGC will not be validated on construction.
        """
        # TODO: Consider lazy loading fields
        if isinstance(gc, gGC):
            self = gc
        else:
            super().__init__(gc)
            self.setdefault('modified', modified)
            self.setdefault('population', population)
            self.setdefault('pgc_ref', self._ref_from_sig('pgc'))
            self.setdefault('ancestor_a_ref', self._ref_from_sig('ancestor_a'))
            self.setdefault('gca_ref', self._ref_from_sig('gca'))
            self.setdefault('gcb_ref', self._ref_from_sig('gcb'))
            self.setdefault('igraph', gc_graph(self.get('graph', {})))
            self.setdefault('evolved', [True])
            self.setdefault('generation', 0)
            self.setdefault('offspring_count', 0)
            if 'inputs' not in self:
                inputs = self['igraph'].input_if()
                outputs = self['igraph'].output_if()
                _, self['input_types'], self['inputs'] = interface_definition(inputs, vtype.EP_TYPE_INT)
                _, self['output_types'], self['outputs'] = interface_definition(outputs, vtype.EP_TYPE_INT)

            # Every GC must have the callable created but only individuals get a wrapped version
            self.setdefault('exec', create_callable(self))
            if population:
                self['exec'] = exec_wrapper(self['exec'])
            for col in filter(lambda x: x[1:] in gc.keys(), gGC.higher_layer_cols):
                gc[col] = copy(gc[col[1:]])

            # PGCs have special fields in the Gene Pool
            if is_pgc(self) and 'pgc_f_valid' not in self:
                if _LOG_DEBUG:
                    _logger.debug(f"{self['pgc_fitness']}")
                self['pgc_delta_fitness'] = [0.0] * NUM_PGC_LAYERS
                self['pgc_previous_fitness'] = copy(self['pgc_fitness'])
                self['pgc_f_valid'] = [f > 0.0 for f in self['pgc_fitness']]
            else:
                self.setdefault('fitness', 0.0)
                self.setdefault('survivability', 0.0)

            if _LOG_DEBUG:
                if self['ref'] == self['gca_ref']:
                    raise ValueError('GC ref == GCA ref. A GC cannot self reference.')
                if self['ref'] == self['gcb_ref']:
                    raise ValueError('GC ref == GCB ref. A GC cannot self reference.')

            # TODO: This should be logger based validation.
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
        mGC(sv=True),
        gGC(sv=True)
    )
    with open('gc_type_table.md', 'w') as file_ptr:
        file_ptr.write("GC Type Field Requirements\n")
        file_ptr.write("==========================\n\n")
        file_ptr.write(_DEFINABLE_MD + ': Defined if not set, ' + _REQUIRED_MD + ': Required.\n\n')
        file_ptr.write('| Field | ' + ' | '.join((x.__class__.__qualname__ for x in gcts)) + ' |\n')
        file_ptr.write('| --- | ' + ' | '.join(('---' for _ in gcts)) + ' |\n')
        for k in sorted(SCHEMA.keys()):
            file_ptr.write(f'| {k} | ' + ' | '.join((_md_string(gct, k) for gct in gcts)) + ' |\n')
