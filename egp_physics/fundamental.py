"""The operations that can be performed on a GC."""
from copy import copy, deepcopy
from logging import DEBUG, Logger, NullHandler, getLogger
from random import randint
from typing import LiteralString, Callable, Literal

from egp_types.eGC import eGC
from egp_types.ep_type import interface_definition, vtype
from egp_types.gc_graph import gc_graph
from egp_types.reference import ref_str
from egp_types.aGC import aGC
from egp_types.dGC import dGC
from egp_types.egp_typing import EndPointType
from egp_types.end_point import dst_end_point
from egp_stores.gene_pool import gene_pool

from .egp_typing import Work


# Logging
_logger: Logger = getLogger(__name__)
_logger.addHandler(NullHandler())
_LOG_DEBUG: bool = _logger.isEnabledFor(DEBUG)

# Filler
_EMPTY_GC_GRAPH = gc_graph()

# Steady state exception filters.
_LMT: LiteralString = ' AND NOT ({exclude_column} = ANY({exclusions})) ORDER BY RANDOM() LIMIT 1'

# TODO: Replace with a localisation hash?
_IT: LiteralString = '{input_types}'
_OT: LiteralString = '{output_types}'
_ITS: LiteralString = '{itypes}::SMALLINT[]'
_OTS: LiteralString = '{otypes}::SMALLINT[]'
_IDX: LiteralString = '{inputs} = {iidx}'
_ODX:  LiteralString = '{outputs} = {oidx}'

_MATCH_TYPE_0_SQL: LiteralString = 'WHERE ' + _IT + ' = ' + _ITS + ' AND ' + _IDX + ' AND ' + _OT + ' = ' + _OTS + ' AND ' + _ODX + _LMT
_MATCH_TYPE_1_SQL: LiteralString = 'WHERE ' + _IT + ' = ' + _ITS + ' AND ' + _OT + ' = ' + _OTS + ' AND ' + _ODX + _LMT
_MATCH_TYPE_2_SQL: LiteralString = 'WHERE ' + _IT + ' = ' + _ITS + ' AND ' + _IDX + ' AND ' + _OT + ' = ' + _OTS + _LMT
_MATCH_TYPE_3_SQL: LiteralString = 'WHERE ' + _IT + ' = ' + _ITS + ' AND ' + _OT + ' = ' + _OTS + _LMT
_MATCH_TYPE_4_SQL: LiteralString = 'WHERE ' + _IT + ' <@ ' + _ITS + ' AND ' + _OT + ' = ' + _OTS + _LMT
_MATCH_TYPE_5_SQL: LiteralString = 'WHERE ' + _IT + ' <@ ' + _ITS + ' AND ' + _OT + ' @> ' + _OTS + _LMT
_MATCH_TYPE_6_SQL: LiteralString = 'WHERE ' + _IT + ' <@ ' + _ITS + ' AND ' + _OT + ' && ' + _OTS + _LMT
_MATCH_TYPE_7_SQL: LiteralString = 'WHERE ' + _IT + ' && ' + _ITS + ' AND ' + _OT + ' && ' + _OTS + _LMT
_MATCH_TYPE_8_SQL: LiteralString = 'WHERE ' + _OT + ' && ' + _OTS + ' ' + _LMT
_MATCH_TYPE_9_SQL: LiteralString = 'WHERE ' + _IT + ' && ' + _ITS + ' ' + _LMT
# Catch for when xtypes is an empty set.
_MATCH_TYPE_10_SQL: LiteralString = 'WHERE ' + _OT + ' = ' + _OTS + ' ' + _LMT
_MATCH_TYPE_11_SQL: LiteralString = 'WHERE ' + _IT + ' = ' + _ITS + ' ' + _LMT


_MATCH_TYPES_SQL: tuple[LiteralString, ...] = (
    _MATCH_TYPE_0_SQL,
    _MATCH_TYPE_1_SQL,
    _MATCH_TYPE_2_SQL,
    _MATCH_TYPE_3_SQL,
    _MATCH_TYPE_4_SQL,
    _MATCH_TYPE_5_SQL,
    _MATCH_TYPE_6_SQL,
    _MATCH_TYPE_7_SQL,
    _MATCH_TYPE_8_SQL,
    _MATCH_TYPE_9_SQL,
    _MATCH_TYPE_10_SQL,
    _MATCH_TYPE_11_SQL
)
_NUM_MATCH_TYPES: int = len(_MATCH_TYPES_SQL)


# PGC Constants
RANDOM_PGC_SIGNATURE: bytes = b'\x00' * 32
_PGC_PARENTAL_PROTECTION_FACTOR: float = 0.75
_POPULATION_PARENTAL_PROTECTION_FACTOR: float = 0.75


def default_dict_gc(ref: Callable[[], int]) -> dGC:
    """Create a default dGC with ref generated from ref()."""
    return {
        'gc_graph': _EMPTY_GC_GRAPH,
        'ancestor_a_ref': 0,
        'ancestor_b_ref': 0,
        'ref': ref(),
        'gca_ref': 0,
        'gcb_ref': 0
    }


def clone(gc_to_clone: aGC, ref: Callable[[], int], copy_graph: bool = True) -> dGC:
    """Create a minimal close of gc.

    Args
    ----
    gc_to_clone: GC to clone

    Returns
    -------
    Minimal clone of gc as a dict.
    """
    return {
        'ref': ref(),
        'ancestor_a_ref': gc_to_clone['ref'],
        'ancestor_b_ref': 0,
        # If gc is a codon then it does not have a GCA
        'gca_ref': gc_to_clone['gca_ref'] if gc_to_clone['gca_ref'] is not None else gc_to_clone['ref'],
        'gcb_ref': gc_to_clone['gcb_ref'],
        # More efficient than reconstructing
        'gc_graph': deepcopy(gc_to_clone['gc_graph']) if copy_graph else gc_graph()
    }


def proximity_select(gms: gene_pool, xputs: dict) -> aGC:
    """Select a genetic code to at least partially connect inputs to outputs.

    The selection is weighted random based on the suitability of the candidates
    based on 'physics' alone.

    For performance (and the fact it is unknown in terms of success at the time of
    development) the selection process follows 2 steps.

        a) The type of match is randomly selected from the list defined by
            _MATCH_TYPES_SQL using the proximity_weights data from the meta table.
        b) From the candidates found by a) randomly select one*.

    In the event no candidates are found for type of match N p_count will be incremented

    for match N and type of match N+1 will be
    attempted. If no matches are found for match type 3 then None is returned.

    TODO: *The performance of this function needs careful benchmarking. The higher index
    match types could return a lot of results for the random selection step.
    TODO: Consider match types where the order of inputs/outputs does not matter.

    Args
    ----
    xputs (dict): The inputs & outputs required of the candidates.
        {
            'input_types': list(int) - list of ascending gc_types in the inputs.
            'inputs': list(int) - Ordered indexes into list above defining the GC input interface.
            'output_types': list(int) - list of ascending gc_types in the inputs.
            'outputs': list(int) - Ordered indexes into list above defining the GC output interface.
            'exclude_column': (str) - Column with values to exclude
            'exclusions': list(object) - Exclude rows with these values in 'exclude_column' column.
        }
        Other key:value pairs will be ignored.

    Returns
    -------
    (int, dict): (match_type, agc) or None
    """
    # TODO: Lots of short queries is inefficient. Ideas:
    #   a) Specific query support from GPC (interface hash)
    #   b) https://stackoverflow.com/questions/42089781/sql-if-select-returns-nothing-then-do-another-select ?
    #   c) Cache general queries (but this means missing out on new options) and randomly select from a list of candidates.
    #   d) Batch queries (but this is architecturally tricky)
    #   e) Optimize DB for these queries.
    #   f) Cache queries at the DB, in the parent process & in the sub-process?
    #   g) This should first search the GP and fallback to the GL
    match_type = randint(0, _NUM_MATCH_TYPES - 1)
    agc = tuple(gms.select(_MATCH_TYPES_SQL[match_type], literals=xputs))
    while not agc and match_type < _NUM_MATCH_TYPES - 1:
        if _LOG_DEBUG:
            _logger.debug(f'Proximity selection match_type {match_type} found no candidates.')
        match_type += 1
        agc = tuple(gms.select(_MATCH_TYPES_SQL[match_type], literals=xputs))
    if agc:
        if _LOG_DEBUG:
            _logger.debug(f'Proximity selection match_type {match_type} found a candidate.')
            _logger.debug(f"Candidate: {agc[0]}")
        return agc[0]
    return None  # FIXME: Cannot be None


def steady_state_exception(gms: gene_pool, fgc: aGC) -> Work:
    """Define what GC must be inserted to complete or partially complete the fgc graph.

    fgc is analysed to determine what end point destinations are unconnected and the highest row
    on which one or more of the unconnected destination endpoints resides.

    Candidates to react (plug the gap) are found from the gms (genetic material store)
    The source of genetic material should be the most local scope.

    If no candidates are found something is very wrong!

    Args
    ----
    gms: The gene_pool.
    fgc: aGC with incomplete graph.

    Returns
    -------
    Insertion work: (target_gc, insert_gc, 'A', 'B' or 'O')
    """
    if _LOG_DEBUG:
        _logger.debug(f"Steady state exception thrown for GC ref {ref_str(fgc.get('ref', 0))}.")
    fgc_gcg: gc_graph = fgc['gc_graph']
    unconnected_eps: tuple[dst_end_point, ...] = tuple(fgc_gcg.i_graph.dst_unref_filter())


    # Find viable source types above the highest row.
    filter_func = fgc_gcg.rows_filter(fgc_gcg.src_rows[above_row], fgc_gcg.src_filter())
    inputs = list((ep[ep_idx.TYPE] for ep in filter(filter_func, fgc_gcg.graph.values())))

    xputs = {
        'exclude_column': 'signature',
        'exclusions': list()
    }
    _, xputs['itypes'], xputs['iidx'] = interface_definition(inputs, vtype.EP_TYPE_INT)
    _, xputs['otypes'], xputs['oidx'] = interface_definition(outputs, vtype.EP_TYPE_INT)

    # Find a gc based on the criteria
    insert_gc = proximity_select(gms, xputs)
    if insert_gc is None:
        _logger.warning("Steady state exception failed to find a candidate in the GMS.")
        return None

    return (fgc, eGC(insert_gc), above_row)


def xGC_inherit(child, parent, pgc):
    """Inherit generic properties from parent to child.

    child is modified.
    parent is modified.

    Args
    ----
    child (xGC): Offspring of parent and pgc.
    parent (xGC): Parent of child.
    pgc (pGC): pGC that operated on parent to product child.
    """
    # TODO: What about survivability? Treat like the above/ something like it?

    child['population_uid'] = parent['population_uid']
    child['ancestor_a_ref'] = parent['ref']
    child['pgc_ref'] = pgc['ref']
    child['generation'] = parent['generation'] + 1
    child['effective_pgc_refs'] = copy(parent['effective_pgc_refs'])
    child['effective_pgc_fitness'] = copy(parent['effective_pgc_fitness'])

    parent['offspring_count'] += 1
