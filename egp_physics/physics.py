"""The operations that can be performed on a GC."""
from collections.abc import Iterable
from copy import copy, deepcopy
from logging import DEBUG, Logger, NullHandler, getLogger
from pprint import pformat
from random import randint
from typing import Literal, LiteralString, Callable, cast

from egp_execution.execution import create_callable
from egp_stores.gene_pool_cache import gene_pool_cache
from egp_stores.gene_pool import gene_pool
from egp_types.eGC import eGC
from egp_types.xGC import xGC, pGC, gGC
from egp_types.ep_type import interface_definition, vtype
from egp_types.gc_graph import DST_EP, SRC_EP, gc_graph
from egp_types.gc_type_tools import M_MASK, NUM_PGC_LAYERS, is_pgc
from egp_types.internal_graph import internal_graph
from egp_types.reference import ref_str
from egp_types.aGC import aGC
from egp_types.dGC import dGC
from egp_types.egp_typing import Row, VALID_ROW_SOURCES

from numpy import array, float32, isfinite
from numpy.random import choice as weighted_choice
from .egp_typing import WorkStack, Work, NewGCDef, InsertRow

_logger: Logger = getLogger(__name__)
_logger.addHandler(NullHandler())
_LOG_DEBUG: bool = _logger.isEnabledFor(DEBUG)

# Filler
_EMPTY_GC_GRAPH = gc_graph()
_SLASH_N = '\n'

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
    """Create a default DictGC with ref generated from ref()."""
    return {
        'gc_graph': _EMPTY_GC_GRAPH,
        'ancestor_a_ref': 0,
        'ancestor_b_ref': 0,
        'ref': ref(),
        'gca_ref': 0,
        'gcb_ref': 0
    }


def _insert_graph_case_0(tig: internal_graph, iig: internal_graph, rig: internal_graph) -> None:
    """Insert igc into tgc case 0."""
    _logger.debug("Case 0: Stack")
    rig.update(iig.copy_row('I', True))
    rig.update(iig.insert_row_as('A'))
    rig.update(tig.insert_row_as('B'))
    rig.update(tig.copy_row('O', True))
    rig.update(rig.direct_connect('I', 'A'))
    rig.update(rig.direct_connect('B', 'O'))


def _insert_graph_case_1(tig: internal_graph, iig: internal_graph, rig: internal_graph) -> None:
    """Insert igc into tgc case 1."""
    _logger.debug("Case 1: No row A or B")
    rig.update(iig.insert_row_as('A'))
    rig.update(tig.copy_row('O'))


def _insert_graph_case_2(tig: internal_graph, iig: internal_graph, rig: internal_graph) -> None:
    """Insert igc into tgc case 2."""
    _logger.debug("Case 2: No row B and insert above A")
    rig.update(iig.insert_row_as('A'))
    rig.update(tig.move_row('A', 'B'))
    rig.update(tig.copy_row('O'))
    rig.redirect_refs('O', DST_EP, 'A', 'B')


def _insert_graph_case_3(tig: internal_graph, iig: internal_graph, rig: internal_graph) -> None:
    """Insert igc into tgc case 3."""
    _logger.debug("Case 3: No row B and insert below A")
    rig.update(tig.copy_rows(('A', 'O')))
    rig.update(iig.insert_row_as('B'))


def _insert_graph_case_4(tig: internal_graph, iig: internal_graph, rig: internal_graph, fig: internal_graph) -> None:
    """Insert igc into tgc case 4."""
    _logger.debug("Case 4: Has rows A & B and insert above A")
    fig.update(tig.copy_rows(('I', 'C'), True))
    fig.update(iig.insert_row_as('A'))
    fig.update(tig.move_row('A', 'B'))
    fig.update(fig.direct_connect('B', 'O'))
    fig.update(fig.append_connect('A', 'O'))
    rig.update(rig.direct_connect('I', 'A'))
    rig.update(fig.move_row_cls('O', DST_EP, 'A', SRC_EP, True))
    rig.update(tig.copy_rows(('B', 'O')))


def _insert_graph_case_5(tig: internal_graph, iig: internal_graph, rig: internal_graph, fig: internal_graph) -> None:
    """Insert igc into tgc case 5."""
    _logger.debug("Case 5: Has rows A & B and insert above B")
    fig.update(tig.copy_rows(('I', 'C'), True))
    fig.update(tig.copy_rows_dst_eps(('A',)))
    fig.update(tig.copy_rows_src_eps(('A',), True))
    fig.update(iig.insert_row_as('B'))
    fig.update(fig.direct_connect('A', 'O'))
    fig.update(fig.append_connect('B', 'O'))
    rig.update(rig.direct_connect('I', 'A'))
    rig.update(fig.move_row_cls('O', DST_EP, 'A', SRC_EP, True))
    rig.update(tig.copy_rows(('B', 'O')))


def _insert_graph_case_6(tig: internal_graph, iig: internal_graph, rig: internal_graph, fig: internal_graph) -> None:
    """Insert igc into tgc case 6."""
    _logger.debug("Case 6: Has rows A & B and insert above O")
    fig.update(tig.copy_rows(('I', 'C'), True))
    fig.update(tig.copy_rows_dst_eps(('A', 'B')))
    fig.update(tig.copy_rows_src_eps(('A', 'B'), True))
    fig.update(fig.direct_connect('A', 'O'))
    fig.update(fig.append_connect('B', 'O'))
    rig.update(rig.direct_connect('I', 'A'))
    rig.update(fig.move_row_cls('O', DST_EP, 'A', SRC_EP, True))
    rig.update(iig.insert_row_as('B'))
    rig.update(tig.copy_row('O', True))


def _insert_graph_case_7(tig: internal_graph, iig: internal_graph, fig: internal_graph) -> None:
    """Insert igc into tgc case 7."""
    fig.update(tig.move_row_cls('A', DST_EP, 'I', SRC_EP, True))
    fig.update(iig.insert_row_as('A'))
    fig.update(tig.move_row('A', 'B', True))
    fig.update(fig.direct_connect('I', 'B'))
    fig.update(tig.move_row_cls('A', SRC_EP, 'O', DST_EP, True))
    fig.update(fig.direct_connect('B', 'O'))


def _insert_graph_case_8(tig: internal_graph, iig: internal_graph, fig: internal_graph) -> None:
    """Insert igc into tgc case 8."""
    fig.update(tig.move_row_cls('A', DST_EP, 'I', SRC_EP, True))
    fig.update(tig.copy_row('A', True))
    fig.update(fig.direct_connect('I', 'A'))
    fig.update(iig.insert_row_as('B'))
    fig.update(tig.move_row_cls('A', SRC_EP, 'O', DST_EP, True))
    fig.update(fig.direct_connect('A', 'O'))


def _insert_graph_case_9(tig: internal_graph, iig: internal_graph, fig: internal_graph) -> None:
    """Insert igc into tgc case 9."""
    fig.update(tig.move_row_cls('B', DST_EP, 'I', SRC_EP, True))
    fig.update(iig.insert_row_as('A'))
    fig.update(tig.copy_row('B', True))
    fig.update(fig.direct_connect('I', 'B'))
    fig.update(tig.move_row_cls('B', SRC_EP, 'O', DST_EP, True))
    fig.update(fig.direct_connect('B', 'O'))


def _insert_graph_case_10(tig: internal_graph, iig: internal_graph, fig: internal_graph) -> None:
    """Insert igc into tgc case 10."""
    fig.update(tig.move_row_cls('B', DST_EP, 'I', SRC_EP, True))
    fig.update(tig.move_row('B', 'A', True))
    fig.update(fig.direct_connect('I', 'A'))
    fig.update(iig.insert_row_as('B'))
    fig.update(tig.move_row_cls('B', SRC_EP, 'O', DST_EP, True))
    fig.update(fig.direct_connect('B', 'O'))


def _insert_graph_case_11(tig: internal_graph, iig: internal_graph, rig: internal_graph) -> None:
    """Insert igc into tgc case 11."""
    _logger.debug("Case 11: Inverse Stack")
    rig.update(tig.copy_row('I', True))
    rig.update(tig.insert_row_as('A'))
    rig.update(iig.insert_row_as('B'))
    rig.update(iig.copy_row('O', True))
    rig.update(rig.direct_connect('I', 'A'))
    rig.update(rig.direct_connect('B', 'O'))


def _insert_graph(tgcg: gc_graph, igcg: gc_graph, above_row: InsertRow) -> tuple[gc_graph, gc_graph]:
    """Insert igc_gcg graph into tgc_gcg graph above row above_row.

    See https://docs.google.com/spreadsheets/d/1YQjrM91e5x30VUIRzipNYX3W7yiFlg6fy9wKbMTx1iY/edit?usp=sharing
    for more details (ask permission to view).
    graph is not modified.

    Args
    ----
    tgcg: Internal gc_graph format gc_graph to insert into.
    igcg: Internal gc_graph format gc_graph to insert.
    above_row: A valid row. Note that 'C' and 'U' are no-ops.

    Returns
    -------
    rig, fig: Result internal graph and formed internal graph
    """
    tig: internal_graph = tgcg.i_graph
    iig: internal_graph = igcg.i_graph
    fig: internal_graph = internal_graph()
    rig: internal_graph = internal_graph()
    if above_row != 'I':
        if tgcg.has_f:
            rig = deepcopy(tgcg.i_graph)
        else:
            rig.update(tig.copy_rows_src_eps(('I', 'C'), True))

    # TODO: There are opportunities to reduce overhead by making some internal_graph manipulation functions
    # act on self rather than returning a dictionary to update (into self)
    if above_row == 'I':
        _insert_graph_case_0(tig, iig, rig)
    elif above_row == 'Z':
        _insert_graph_case_11(tig, iig, rig)
    elif not tgcg.has_a:
        _insert_graph_case_1(tig, iig, rig)
    elif not tgcg.has_f:
        if not tgcg.has_b:
            if above_row == 'A':
                _insert_graph_case_2(tig, iig, rig)
            else:
                _insert_graph_case_3(tig, iig, rig)
        else:
            if above_row == 'A':
                _insert_graph_case_4(tig, iig, rig, fig)
            elif above_row == 'B':
                _insert_graph_case_5(tig, iig, rig, fig)
            else:
                _insert_graph_case_6(tig, iig, rig, fig)
    else:
        if above_row == 'A':
            _insert_graph_case_7(tig, iig, fig)
        elif above_row == 'O':
            _insert_graph_case_8(tig, iig, fig)
        elif above_row == 'B':
            _insert_graph_case_9(tig, iig, fig)
        elif above_row == 'P':
            _insert_graph_case_10(tig, iig, fig)

    # Pre-completed logging
    if _LOG_DEBUG:
        _logger.debug(f"tgc ({type(tgcg)}):\n{pformat(tgcg)}")
        _logger.debug(f"igc ({type(tgcg)}):\n{pformat(igcg)}")
        _logger.debug(f"Pre-completed rig ({type(rig)}):\n{pformat(rig)}")
        _logger.debug(f"Pre-completed fig ({type(fig)}):\n{pformat(fig)}")

    # Tidy up dangling references and normalize (make new connections)
    rig.complete_references()
    rgc_gc_graph = gc_graph(i_graph=rig)
    rgc_gc_graph.normalize()
    if fig:
        fig.complete_references()
        fgc_gc_graph = gc_graph(i_graph=fig)
        fgc_gc_graph.normalize()
    else:
        fgc_gc_graph = gc_graph()

    # Post-completed logging
    if _LOG_DEBUG:
        _logger.debug(f"Completed rig:\n{pformat(rig)}")
        _logger.debug(f"Completed fig:\n{pformat(fig)}")

    return rgc_gc_graph, fgc_gc_graph


def stablize(gms: gene_pool, tgc: aGC) -> NewGCDef:
    """If tgc is not stable force a steady state exceptiion"""
    if tgc['gc_graph'].is_stable():
        if _LOG_DEBUG:
            assert tgc['gc_graph'].validate()
            _logger.debug('Target GC is stable.')
        return (tgc, {})
    if _LOG_DEBUG:
        _logger.debug('Target GC is unstable.')
    return _recursive_insert_gc(gms, [steady_state_exception(gms, tgc)])


def _insert_gc(gms: gene_pool, tgc: aGC, igc: aGC, above_row: InsertRow) -> NewGCDef:
    """Insert igc into tgc above row 'above_row'.

    Args
    ----
    gms: A source of genetic material & references. Needs to be in the context of the sub-process.
    tgc: GC to insert insert_gc into.
    igc: GC to insert into target_gc.
    above_row: Insert above this row.

    Returns
    -------
    (rgc, {ref: fgc}): First fgc is the igc followed by fgc's
    created to stabilise rgc.
    """
    if _LOG_DEBUG:
        _logger.debug('Inserting into Target GC.')
    return _recursive_insert_gc(gms, [(tgc, igc, above_row)])


def _insert_gc_case_0(tgc: aGC, igc: aGC, rgc: aGC) -> None:
    """Insert igc data into tgc case 0."""
    _logger.debug("Case 0: Stack")
    rgc['gca_ref'] = igc['ref']
    rgc['gcb_ref'] = tgc['ref']


def _insert_gc_case_1(igc: aGC, rgc: aGC) -> None:
    """Insert igc data into tgc case 1."""
    _logger.debug("Case 1")
    rgc['gca_ref'] = igc['ref']


def _insert_gc_case_2(tgc: aGC, igc: aGC, rgc: aGC) -> None:
    """Insert igc data into tgc case 2."""
    _logger.debug("Case 2")
    rgc['gca_ref'] = igc['ref']
    if tgc['gca_ref'] is not None:
        rgc['gcb_ref'] = tgc['gca_ref']
    else:
        rgc['gcb_ref'] = tgc['ref']


def _insert_gc_case_3(tgc: aGC, igc: aGC, rgc: aGC) -> None:
    """Insert igc data into tgc case 3."""
    _logger.debug("Case 3")
    if tgc['gca_ref'] is not None:
        rgc['gca_ref'] = tgc['gca_ref']
    else:
        rgc['gca_ref'] = tgc['ref']
    rgc['gcb_ref'] = igc['ref']


def _insert_gc_case_4(tgc: aGC, igc: aGC, rgc: aGC, fgc: aGC) -> None:
    """Insert igc data into tgc case 4."""
    _logger.debug("Case 4")
    fgc['gca_ref'] = igc['ref']
    fgc['gcb_ref'] = tgc['gca_ref']
    rgc['gca_ref'] = fgc['ref']
    rgc['gcb_ref'] = tgc['gcb_ref']


def _insert_gc_case_5(tgc: aGC, igc: aGC, rgc: aGC, fgc: aGC) -> None:
    """Insert igc data into tgc case 5."""
    _logger.debug("Case 5")
    fgc['gca_ref'] = tgc['gca_ref']
    fgc['gcb_ref'] = igc['ref']
    rgc['gca_ref'] = fgc['ref']
    rgc['gcb_ref'] = tgc['gcb_ref']


def _insert_gc_case_6(tgc: aGC, igc: aGC, rgc: aGC, fgc: aGC) -> None:
    """Insert igc data into tgc case 6."""
    _logger.debug("Case 6")
    fgc['gca_ref'] = tgc['gca_ref']
    fgc['gcb_ref'] = tgc['gcb_ref']
    rgc['gca_ref'] = fgc['ref']
    rgc['gcb_ref'] = igc['ref']


def _insert_gc_case_7(tgc: aGC, igc: aGC, rgc: aGC, fgc: aGC) -> None:
    """Insert igc data into tgc case 7."""
    _logger.debug("Case 7")
    rgc['gca_ref'] = fgc['ref']
    rgc['gcb_ref'] = tgc['gcb_ref']
    fgc['gca_ref'] = igc['ref']
    fgc['gcb_ref'] = tgc['gca_ref']


def _insert_gc_case_8(tgc: aGC, igc: aGC, rgc: aGC, fgc: aGC) -> None:
    """Insert igc data into tgc case 8."""
    _logger.debug("Case 8")
    rgc['gca_ref'] = fgc['ref']
    rgc['gcb_ref'] = tgc['gcb_ref']
    fgc['gca_ref'] = tgc['gca_ref']
    fgc['gcb_ref'] = igc['ref']


def _insert_gc_case_9(tgc: aGC, igc: aGC, rgc: aGC, fgc: aGC) -> None:
    """Insert igc data into tgc case 9."""
    _logger.debug("Case 9")
    rgc['gca_ref'] = tgc['gca_ref']
    rgc['gcb_ref'] = fgc['ref']
    fgc['gca_ref'] = igc['ref']
    fgc['gcb_ref'] = tgc['gcb_ref']


def _insert_gc_case_10(tgc: aGC, igc: aGC, rgc: aGC, fgc: aGC) -> None:
    """Insert igc data into tgc case 10."""
    _logger.debug("Case 10")
    rgc['gca_ref'] = tgc['gca_ref']
    rgc['gcb_ref'] = fgc['ref']
    fgc['gca_ref'] = tgc['gcb_ref']
    fgc['gcb_ref'] = igc['ref']


def _insert_gc_case_11(tgc: aGC, igc: aGC, rgc: aGC) -> None:
    """Insert igc data into tgc case 11."""
    _logger.debug("Case 11: Inverse Stack")
    rgc['gca_ref'] = tgc['ref']
    rgc['gcb_ref'] = igc['ref']


def _recursive_insert_gc(gms: gene_pool, work_stack: WorkStack) -> NewGCDef:
    """Recursively insert GC's until the target GC is stable.

    A work stack is used to avoid actual recursion.

    While there is work to do:
        Pop work off the work stack.
        If the target_gc has a B row then to insert a GC will
        require two new ones to be made: One to combine the top two
        GC's and a second to bind that GC to the bottom one (fgc & tgc respectively).
        fgc is not restricted to have the same inputs and outputs as the
        target_gc unlike tgc which must. If either fgc or tgc is not
        in a steady state then a steady state exception is thrown the handler of
        which will return insertion work to complete the graph. That work will be
        added to the top of the work stack. Steady state GC's are added to the
        return value, fgc to the back and tgc to the front, thus when the function
        returns the correct tgc will be at the head of the fgc_list.

        Args
        ----
        work_stack: List of (Target GC, Insert GC, above row) worl to do

        Returns
        -------
        A new GC definition structure.
    """
    fgc_dict: dict[int, aGC] = {}
    new_tgc_flag: bool = False
    new_tgc: dGC = default_dict_gc(int)
    while work_stack:
        if _LOG_DEBUG:
            _logger.debug(f"Work stack depth: {len(work_stack)}")

        fgc: dGC = default_dict_gc(gms.next_reference)
        rgc: dGC = default_dict_gc(gms.next_reference)
        work: Work = work_stack.pop(0)
        tgc = work[0]
        igc = work[1]
        above_row = work[2]

        if _LOG_DEBUG:
            _logger.debug(f"Work: Target={ref_str(tgc.get('ref', 0))}, "
                          f"Insert={ref_str(igc.get('ref', 0))}, Above Row={above_row}")

        # Insert into the graph
        tgcg: gc_graph = tgc['gc_graph']
        igcg: gc_graph = igc['gc_graph']
        fgcg: gc_graph
        rgcg: gc_graph
        rgcg, fgcg = _insert_graph(tgcg, igcg, above_row)
        rgcg_steady: bool = rgcg.normalize()
        fgcg_steady: bool = False
        if fgcg:
            fgcg_steady = fgcg.normalize()
            if _LOG_DEBUG:
                _logger.debug(f"Normalized fgc:\n{pformat(fgcg)}")
            fgc['graph'] = fgcg.app_graph
            fgc['gc_graph'] = fgcg
            fgc['ancestor_a_ref'] = igc['ref']
            fgc['ancestor_b_ref'] = tgc['ref']
        if _LOG_DEBUG:
            _logger.debug(f"Normalized rgc:\n{pformat(rgcg)}")
        rgc['graph'] = rgcg.app_graph
        rgc['gc_graph'] = rgcg
        rgc['ancestor_a_ref'] = tgc['ref']
        rgc['ancestor_b_ref'] = igc['ref']

        # Insert into the GC
        # The insert_gc is always referenced in the tree of the final rgc
        fgc_dict[igc['ref']] = igc
        fgc_dict[tgc['ref']] = tgc
        if above_row == 'I':
            _insert_gc_case_0(tgc, igc, rgc)
        elif above_row == 'Z':
            _insert_gc_case_11(tgc, igc, rgc)
        elif not tgcg.has_a:
            _insert_gc_case_1(igc, rgc)
        elif not tgcg.has_f:
            if not tgcg.has_b:
                if above_row == 'A':
                    _insert_gc_case_2(tgc, igc, rgc)
                else:
                    _insert_gc_case_3(tgc, igc, rgc)
            else:
                if above_row == 'A':
                    _insert_gc_case_4(tgc, igc, rgc, fgc)
                elif above_row == 'B':
                    _insert_gc_case_5(tgc, igc, rgc, fgc)
                else:
                    _insert_gc_case_6(tgc, igc, rgc, fgc)
        else:

            if above_row == 'A':
                _insert_gc_case_7(tgc, igc, rgc, fgc)
            elif above_row == 'O':
                _insert_gc_case_8(tgc, igc, rgc, fgc)
            elif above_row == 'B':
                _insert_gc_case_9(tgc, igc, rgc, fgc)
            elif above_row == 'P':
                _insert_gc_case_10(tgc, igc, rgc, fgc)

        # rgc['ref'] must be new & replace any previous mentions
        # of target_gc['ref'] in fgc_dict[*][...ref fields...]
        # and the new_tgc if it is defined.
        #
        # In the case where target_gc is unstable it is not in fgc_dict
        # but will appear
        # TODO: There must be a more efficient way of doing this
        # FIXME: I think this breaks if tgc == igc
        # IDEA: Could we use a data class for references?
        new_ref: int = rgc['ref']
        old_ref: int = tgc['ref']
        if _LOG_DEBUG:
            _logger.debug(f"Replacing {ref_str(old_ref)} with {ref_str(new_ref)}.")
        for nfgc in fgc_dict.values():
            for ref in ('gca_ref', 'gcb_ref', 'ancestor_a_ref', 'ancestor_b_ref'):
                if nfgc[ref] == old_ref:
                    nfgc[ref] = new_ref
        if new_tgc_flag:
            for ref in ('gca_ref', 'gcb_ref', 'ancestor_a_ref', 'ancestor_b_ref'):
                if new_tgc[ref] == old_ref:
                    new_tgc[ref] = new_ref

        # Check we have valid graphs
        # FGC is added to the work stack ahead of RGC
        # which ensures the new TGC is correctly set
        # (i.e. does not end up being the RGC of an unstable FGC)
        if fgcg:
            if not fgcg_steady:
                if _LOG_DEBUG:
                    _logger.debug(f"FGC ref {ref_str(fgc['ref'])} is unstable.")
                work_stack.insert(0, steady_state_exception(gms, fgc))
            else:
                if _LOG_DEBUG:
                    assert fgcg.validate(), "FGC validation failed."
                    _logger.debug(f"FGC ref {ref_str(fgc['ref'])} added to fgc_dict.")
                fgc_dict[fgc['ref']] = fgc

        if not rgcg_steady:
            if _LOG_DEBUG:
                _logger.debug(f"RGC ref {ref_str(rgc['ref'])} is unstable.")
            work_stack.insert(0, steady_state_exception(gms, rgc))
        elif not new_tgc_flag:
            if _LOG_DEBUG:
                assert rgcg.validate(), "RGC validation failed."
                _logger.debug(f"Resultant GC defined:\n{rgc}")
            new_tgc_flag = True
            new_tgc = rgc
        else:
            if _LOG_DEBUG:
                assert rgcg.validate(), "RGC validation failed."
                _logger.debug(f"RGC ref {ref_str(rgc['ref'])} added to fgc_dict.")
            fgc_dict[rgc['ref']] = rgc

        if _LOG_DEBUG:
            _logger.debug(f"fgc_dict: {[ref_str(x) for x in fgc_dict]}")

    if _LOG_DEBUG:
        _logger.debug(f"fgc_dict details:\n{_SLASH_N.join(ref_str(k) + ':' + _SLASH_N + str(v) for k, v in fgc_dict.items())}")
        # TODO: target_gc & new_tgc interface must be the same. Validate.

    return (new_tgc, fgc_dict)


def gc_insert(gms: gene_pool, tgc: aGC, igc: aGC, above_row: Literal['I', 'A', 'B', 'O']) -> xGC:
    """Insert insert_gc into target_gc above row 'above_row'.

    If insert_gc is None then the target_gc is assessed for stability. If it
    is stable then [target_gc] will be returned otherwise a steady state exception
    is 'thrown'.

    NOTE: If a steady state exception occurs for which a candidate cannot
    be found in the GMS this function returns None.

    Args
    ----
    gms: A source of genetic material & references. Needs to be in the context of the sub-process.
    tgc: The target GC.
    igc: GC to insert into tgc.
    above_row: One of 'I', 'A', 'B' or 'O'.

    Returns
    -------
    The GC created as a result of the insertion.
    """
    new_gc_definition: NewGCDef = _insert_gc(gms, tgc, igc, above_row)
    gms.pool[new_gc_definition[0]['ref']] = new_gc_definition[0]
    gms.pool.update(new_gc_definition[1])
    return gms.pool[new_gc_definition[0]['ref']]


def gc_stack(gms: gene_pool, bottom_gc: aGC, top_gc: aGC, invert=False) -> xGC:
    """Stack two GC's.

    top_gc is stacked on top of bottom_gc to create a new gc.
    If the new gc is invalid it is repaired with a steady state exception.

    Args
    ----
    gms: A source of genetic material.
    top_gc: GC to stack on top
    bottom_gc: GC to put on the bottom.
    invert: If True the top_gc is stacked below the bottom_gc.

    Returns
    -------
    rgc: Resultant GC with a valid graph
    """
    new_gc_definition: NewGCDef = _insert_gc(gms, bottom_gc, top_gc, ('I', 'Z')[invert])
    gms.pool[new_gc_definition[0]['ref']] = new_gc_definition[0]
    gms.pool.update(new_gc_definition[1])
    return gms.pool[new_gc_definition[0]['ref']]


def _clone(gc_to_clone: aGC, ref: Callable[[], int], copy_graph: bool = True) -> dGC:
    """Create a minimal clone of gc.

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
        'gca_ref': gc_to_clone['gca_ref'] if gc_to_clone['gca_ref'] else gc_to_clone['ref'],
        'gcb_ref': gc_to_clone['gcb_ref'],
        # More efficient than reconstructing
        'gc_graph': deepcopy(gc_to_clone['gc_graph']) if copy_graph else gc_graph()
    }


def gc_remove(gms: gene_pool, tgc: aGC, row: Row):
    """Remove row A, B, C, P or O from tgc['graph'] to create rgc.

    The subsequent invalid graph is normalised and used to create rgc.
    Removing a row is likely to result in an invalid graph which is
    repaired using recursive steady state exceptions.

    In the case row = 'P' or 'O' this is a single to remove row F
    and remove the path through 'P' or 'O'. NB: GC will still have a valid
    row 'O'.

    Args
    ----
    gms: A source of genetic material.
    tgc: Target xGC to modify.
    row: Either 'A', 'B', 'C', 'P' or 'O'.

    Returns
    -------
    rgc: Resultant minimal GC with a valid graph
    """
    rgc: dGC = _clone(tgc, gms.next_reference, False)
    if _LOG_DEBUG:
        _logger.debug(f"Minimally cloned {ref_str(tgc['ref'])} to {ref_str(rgc['ref'])}")
        _logger.debug(f'Removing row {row}.')
    rgc_igraph: internal_graph = internal_graph()
    rgc_igraph.update(tgc['gc_graph'].i_graph.copy_row('I'))
    if row != 'C':
        rgc_igraph.update(tgc['gc_graph'].i_graph.copy_row('C'))
    if row == 'A':
        rgc['gca_ref'] = tgc['gcb_ref']
        rgc['gcb_ref'] = 0
        if  tgc['gc_graph'].has_b:
            rgc_igraph.update(tgc['gc_graph'].i_graph.move_row('B', 'A', has_f=tgc['gc_graph'].has_f))
    elif row == 'B':
        rgc['gca_ref'] = tgc['gca_ref']
        rgc_igraph.update(tgc['gc_graph'].i_graph.copy_row('A'))
        rgc['gcb_ref'] = 0
    elif row == 'P':
        rgc_igraph.update(tgc['gc_graph'].i_graph.copy_row('O'))
    elif row == 'O':
        rgc_igraph.update(tgc['gc_graph'].i_graph.move_row('P', 'O'))
    rgc['gc_graph'].normalize()
    return _pgc_epilogue(gms, rgc)


def _pgc_epilogue(gms: gene_pool, agc: aGC) -> xGC:
    """Stabilises aGC and updates the gms.

    pGC's may return a single aGC or None if a valid aGC could
    not be created. Since pGC's may be stacked each pGC must
    be able to handle a None input.

    fGC's created as a side effect of of pGC activity are
    added to the GP local cache.

    Args
    ----
    gms: A source of genetic material.
    agc: Unstable GC

    Returns
    -------
    rgc: Resultant gGC
    """
    if _LOG_DEBUG:
        _logger.debug(f'pGC epilogue with agc = {agc}')
    new_gc_definition: NewGCDef = stablize(gms, agc)
    gms.pool[new_gc_definition[0]['ref']] = new_gc_definition[0]
    gms.pool.update(new_gc_definition[1])
    return gms.pool[new_gc_definition[0]['ref']]


def gc_remove_all_connections(gms: gene_pool, tgc: xGC) -> xGC:
    """Remove all the connections in gc's graph.

    The subsequent invalid graph is normalised and used to create rgc.
    rgc is very likely to have an invalid graph which is
    repaired using recursive steady state exceptions.

    NOTE: If a steady state exception occurs for which a candidate cannot
    be found in the GMS this function returns None.

    Args
    ----
    gms: A source of genetic material.
    tgc: Target xGC to modify.

    Returns
    -------
    rgc: Resultant xGC with a valid graph
    """
    dgc: dGC = _clone(tgc, gms.next_reference)
    if _LOG_DEBUG:
        _logger.debug(f"Minimally cloned {ref_str(tgc['ref'])} to {ref_str(dgc['ref'])}")
    dgc['gc_graph'].remove_all_connections()
    dgc['gc_graph'].normalize()
    return _pgc_epilogue(gms, dgc)


def gc_add_input(gms: gene_pool, tgc: xGC) -> xGC:
    """Add a random input to the GC.

    The subsequent graph is normalised and used to create rgc.
    If rgc has an invalid graph it is
    repaired using recursive steady state exceptions.

    NOTE: If a steady state exception occurs for which a candidate cannot
    be found in the GMS this function returns None.

    Args
    ----
    gms: A source of genetic material.
    tgc: Target xGC to modify.

    Returns
    -------
    rgc: Resultant xGC with a valid graph or None
    """
    dgc: dGC = _clone(tgc, gms.next_reference)
    if _LOG_DEBUG:
        _logger.debug(f"Minimally cloned {ref_str(tgc['ref'])} to {ref_str(dgc['ref'])}")
    dgc['gc_graph'].add_input()
    dgc['gc_graph'].normalize()
    return _pgc_epilogue(gms, dgc)


def gc_remove_input(gms: gene_pool, tgc: xGC) -> xGC:
    """Remove a random input from the GC.

    The subsequent graph is normalised and used to create rgc.
    If rgc has an invalid graph it is
    repaired using recursive steady state exceptions.

    NOTE: If a steady state exception occurs for which a candidate cannot
    be found in the GMS this function returns None.

    Args
    ----
    gms: A source of genetic material.
    tgc: Target xGC to modify.

    Returns
    -------
    rgc: Resultant minimal GC with a valid graph or None
    """
    dgc: dGC = _clone(tgc, gms.next_reference)
    if _LOG_DEBUG:
        _logger.debug(f"Minimally cloned {ref_str(tgc['ref'])} to {ref_str(dgc['ref'])}")
    dgc['gc_graph'].remove_input()
    dgc['gc_graph'].normalize()
    return _pgc_epilogue(gms, dgc)


def gc_add_output(gms: gene_pool, tgc: xGC) -> xGC:
    """Add a random output to the GC.

    The subsequent graph is normalised and used to create rgc.
    If rgc has an invalid graph it is
    repaired using recursive steady state exceptions.

    NOTE: If a steady state exception occurs for which a candidate cannot
    be found in the GMS this function returns None.

    Args
    ----
    gms: A source of genetic material.
    tgc: Target xGC to modify.

    Returns
    -------
    rgc: Resultant minimal GC with a valid graph or None
    """
    dgc: dGC = _clone(tgc, gms.next_reference)
    if _LOG_DEBUG:
        _logger.debug(f"Minimally cloned {ref_str(tgc['ref'])} to {ref_str(dgc['ref'])}")
    dgc['gc_graph'].add_output()
    dgc['gc_graph'].normalize()
    return _pgc_epilogue(gms, dgc)


def gc_remove_output(gms: gene_pool, tgc: xGC) -> xGC:
    """Add a random input from the GC.

    The subsequent graph is normalised and used to create rgc.
    If rgc has an invalid graph it is
    repaired using recursive steady state exceptions.

    NOTE: If a steady state exception occurs for which a candidate cannot
    be found in the GMS this function returns None.

    Args
    ----
    gms: A source of genetic material.
    tgc: Target xGC to modify.

    Returns
    -------
    rgc: Resultant minimal GC with a valid graph or None
    """
    dgc: dGC = _clone(tgc, gms.next_reference)
    if _LOG_DEBUG:
        _logger.debug(f"Minimally cloned {ref_str(tgc['ref'])} to {ref_str(dgc['ref'])}")
    dgc['gc_graph'].add_output()
    dgc['gc_graph'].normalize()
    return _pgc_epilogue(gms, dgc)


def gc_remove_constant(gms: gene_pool, tgc: xGC) -> xGC:
    """Remove a random constant from the GC.

    The subsequent graph is normalised and used to create rgc.
    If rgc has an invalid graph it is
    repaired using recursive steady state exceptions.

    NOTE: If a steady state exception occurs for which a candidate cannot
    be found in the GMS this function returns None.

    Args
    ----
    gms: A source of genetic material.
    tgc: Target xGC to modify.

    Returns
    -------
    rgc: Resultant minimal GC with a valid graph or None
    """
    dgc: dGC = _clone(tgc, gms.next_reference)
    if _LOG_DEBUG:
        _logger.debug(f"Minimally cloned {ref_str(tgc['ref'])} to {ref_str(dgc['ref'])}")
    dgc['gc_graph'].remove_constant()
    dgc['gc_graph'].normalize()
    return _pgc_epilogue(gms, dgc)


def interface_proximity_select(gms: gene_pool, xputs: dict[str, bytes | list | str]) -> aGC:
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

    TODO: The occurance rate of this function needs to be carefully monitored. It probably
    should push back some cost into the system. Mutations should evolve to be stable.
    This is a key metric to the maturity of the system.    
    TODO: *The performance of this function needs careful benchmarking. The higher index
    match types could return a lot of results for the random selection step.
    TODO: Consider match types where the order of inputs/outputs does not matter.
    TODO: To find a 'match' we can not only loosen the criteria we can also
    reach out further to the microbiome (which can reach out to the biome and so on)

    Args
    ----
    xputs: The inputs & outputs required of the candidates.
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
    aGC
    """
    # TODO: Lots of short queries is inefficient. Ideas:
    #   a) Specific query support from GPC (interface hash)
    #   b) https://stackoverflow.com/questions/42089781/sql-if-select-returns-nothing-then-do-another-select ?
    #   c) Cache general queries (but this means missing out on new options) and randomly select from a list of candidates.
    #   d) Batch queries (but this is architecturally tricky)
    #   e) Optimize DB for these queries.
    #   f) Cache queries at the DB, in the parent process & in the sub-process?
    #   g) This should first search the GP and fallback to the GL
    # TODO: Radical idea: Above is for truely random GC selection. Rare events. In general 'nearby' GCs
    #   will be close relations in the gene pool. 
    match_type: int = randint(0, _NUM_MATCH_TYPES - 1)
    agc = tuple(gms.select(_MATCH_TYPES_SQL[match_type], literals=xputs))
    while not agc and match_type < _NUM_MATCH_TYPES - 1:
        if _LOG_DEBUG:
            _logger.debug(f'Proximity selection match_type {match_type} found no candidates.')
        match_type += 1
        agc = tuple(gms.select(_MATCH_TYPES_SQL[match_type], literals=xputs))
    if _LOG_DEBUG:
        _logger.debug(f'Proximity selection match_type {match_type} found a candidate.')
        _logger.debug(f"Candidate: {agc[0]}")
    return agc[0]


def steady_state_exception(gms: gene_pool, fgc: aGC) -> Work:
    """Define what GC must be inserted to complete or partially complete the fgc graph.

    fgc is analysed to determine what end point destinations are unconnected and the highest row
    on which one or more of the unconnected destination endpoints resides.

    Candidates to react (plug the gap) are found from the gms (genetic material source) which is
    either a gene_pool or genomic_library instance. NB: The source of genetic material should be
    the most local scope.

    If no candidates are found in the GMS then None is returned.

    Args
    ----
    gms: A source of genetic material
    fgc (fGC): fGC with incomplete graph.

    Returns
    -------
    A workstack: (target_gc, insert_gc, 'A', 'B' or 'O')
    """
    if _LOG_DEBUG:
        _logger.debug(f"Steady state exception thrown for GC ref {ref_str(fgc['ref'])}.")
    fgc_graph: gc_graph = fgc['gc_graph']

    # Find unconnected destination endpoints. Determine highest row & endpoint types.
    above_row: str = 'Z'
    outputs: list[int] = []
    for ep in fgc_graph.i_graph.dst_unref_filter():
        if ep.row < above_row:
           _above_row = ep.row
        outputs.append(ep.typ)

    # Find viable source types above the highest row.
    inputs: list[int] = [ep.typ for ep in fgc_graph.i_graph.src_rows_filter(VALID_ROW_SOURCES[fgc['gc_graph'].has_f][cast(Row, above_row)])]

    xputs: dict[str, bytes | list | str] = {
        'exclude_column': 'signature',
        'exclusions': list()
    }
    _, xputs['itypes'], xputs['iidx'] = interface_definition(inputs, vtype.EP_TYPE_INT)
    _, xputs['otypes'], xputs['oidx'] = interface_definition(outputs, vtype.EP_TYPE_INT)

    # Find a gc based on the criteria
    insert_gc: aGC = interface_proximity_select(gms, xputs)
    return (fgc, insert_gc, cast(InsertRow, above_row))


def create_SMS(gms: gene_pool, pgc: pGC, ggc: gGC):
    """Create a Super Mutation Sequence (SMS).

    pgc has positively mutated gc's parent to create ggc.
    The chains of mutations that led to this positive effect are called SMS's.
    This function creates them as stand alone pGC's in the GP.

    Modifies gp & ggc.

    Args
    ----
    gp (gene_pool): The gene_pool that contains pGC and its creators.
    pgc (pGC): A physical GC that positively changed ggc's parent to create ggc.
    ggc (gGC): A target population individual.
    """
    # Removed: This should be handled through primatives
    # i.e. Make these mutation pGC's
    # 1. Get xGC ancestor pGC
    # 2. Get xGC pGC
    # 3. Stack 1. on 2.
    # 4. Get delta fitness for xGC and ancestor(s)
    # 5. if 4. > 0.0 then 3. else no-op
    #  etc. These pGC's will then get fitter or die out. 


def pGC_fitness(gms: gene_pool, pgc: pGC, ggc: gGC, delta_fitness: float) -> int:
    """Update the fitness of the pGC and the pGC's that created it.

    pgc is modified.
    -1.0 <= delta_fitness <= 1.0

    pGC's are checked for evolution

    Args
    ----
    gms: The gene_pool that contains pGC and its creators.
    pgc: A physical GC.
    ggc: The gGC created by pgc changing fitness from its parent by delta_fitness
    delta_fitness: The change in fitness of the GC pGC mutated.

    Returns
    -------
    The number of pGC evolutions that occured as a result of the fitness update.
    """
    depth = 0
    _pGC_fitness(pgc, delta_fitness, depth)
    delta_fitness = pgc['pgc_delta_fitness'][depth]
    evolved: bool = evolve_physical(gms, pgc, depth)
    evolutions = int(evolved)
    pgc_creator: xGC | None = gms.pool.get(pgc['pgc_ref'])
    while evolved and pgc_creator is not None:
        depth += 1
        pgc_creator = cast(pGC, pgc_creator)
        _pGC_evolvability(pgc_creator, delta_fitness, depth)
        _pGC_fitness(pgc_creator, delta_fitness, depth)
        delta_fitness = pgc_creator['pgc_delta_fitness'][depth]
        evolved = evolve_physical(gms, pgc_creator, depth)
        evolutions += evolved
        pgc = pgc_creator
        pgc_creator = gms.pool.get(pgc_creator['pgc_ref'], None)
    return evolutions


def _pGC_fitness(pgc: pGC, delta_fitness: float, depth: int) -> float:
    """Update the fitness of the pGC.

    pgc is modified.
    -1.0 <= delta_fitness <= 1.0 or None which maps to -1.0

    A pGC that has a net neutral affect on target GC's has a fitness of 0.5
    i.e. it increases fitness by as much as it decreases it.
    This can only be true if compute resource impact is part of fitness.

    pGC's with fitness below 0.5 are not necessarily useless. They may be effective
    as a sub-GC of an SMS or be highly evolvable.

    Args
    ----
    pgc: pGC to update.
    xgc: gGC pgc mutated.
    delta_fitness: The change in fitness of the GC pGC mutated.
    depth: The layer in the environment pgc is at.

    Returns
    -------
    Mapped delta_fitness
    """
    old_count: int = pgc['pgc_f_count'][depth]
    pgc['pgc_f_count'][depth] += 1
    pgc['pgc_fitness'][depth] = (old_count * pgc['pgc_fitness'][depth] + (delta_fitness / 2 + 0.5)) / pgc['pgc_f_count'][depth]
    return delta_fitness


def _pGC_evolvability(pgc: pGC, delta_fitness: float, depth: int) -> None:
    """Update the evolvability of a PGC.

    pgc is modified.
    -1.0 <= delta_fitness <= 1.0

    Args
    ----
    pgc (pGC): PGC to update.
    delta_fitness (float): Difference in fitness between this GC & its offspring.
    depth (int): The layer in the environment pgc is at.
    """
    increase: float = 0.0 if delta_fitness < 0 else delta_fitness
    old_count: int = pgc['pgc_e_count'][depth]
    pgc['pgc_e_count'][depth] += 1
    pgc['pgc_evolvability'][depth] = (old_count * pgc['pgc_evolvability'][depth] + increase) / pgc['pgc_e_count'][depth]


def population_GC_evolvability(ggc: gGC, delta_fitness: float) -> None:
    """Update the evolvability of a population GC.

    xgc is modified.
    -1.0 <= delta_fitness <= 1.0

    Args
    ----
    xgc (pGC): xGC to update.
    delta_fitness (float): Difference in fitness between this GC & its offspring.
    """
    increase: float = 0.0 if delta_fitness < 0 else delta_fitness
    old_count: int = ggc['e_count']
    ggc['e_count'] += 1
    ggc['evolvability'] = (old_count * ggc['evolvability'] + increase) / ggc['e_count']


def evolve_physical(gms: gene_pool, pgc: pGC, depth: int) -> bool:
    """Evolve the pgc as needed.

    pgc is checked to see if it meets evolution criteria. If it does
    it is evolved & the gene pool updated with its offspring.

    pGC's evolve when they have been 'used' M_CONSTANT times (which must
    be a positive integer power of 2). The number of uses is persisted in
    all scopes.

    Args
    ----
    gp (gene_pool): The gene pool containing pgc.
    pgc (pGC): The pgc to evolve as necessary.
    depth (int): The layer in the environment pgc is at.

    Returns
    -------
    (bool): True if the pGC was evolved else False
    """
    if not (pgc['pgc_f_count'][depth] & M_MASK):
        pgc['pgc_delta_fitness'][depth] = 0.0
        ppgc: list[pGC] = select_pGC(gms, [pgc], depth + 1)
        wrapped_ppgc_callable = create_callable(ppgc, gms.pool)
        result = wrapped_ppgc_callable((pgc,))
        if result is None:
            # pGC went pop - should not happen very often
            _logger.warning(f"ppGC {ref_str(pgc['ref'])} threw an exception when called.")
            offspring = None
        else:
            offspring = result[0]

        if offspring is not None:
            if _LOG_DEBUG:
                assert isinstance(offspring, xGC)
            pGC_inherit(offspring, pgc, ppgc)
        return True
    return False


def select_pGC(gms: gene_pool, xgc_refs: Iterable[int], depth: int = 0) -> list[pGC]:
    """Select a pgc to evolve xgc.

    The Gene Pool Cache is refreshed on a periodic basis by the expiry of the sub-process
    and a write back to the Gene Pool. At this point at least some pGC's that are not direct decendants
    of this Gene Pool Cache will be bought in.

    Selection moves 

    A pGC is selected to act on each xGC.

    pGC selection is achieved as follows:

    If xgc['next_pgc_ref'] is set then that is the pGC used else
    a choice is made from effective_pgc_refs, pGC's with positive
    fitness (>0.5) and pGC's with negative fitness with relative
    weights of 4, 2 and 1 respectively.
    In the event that any of these categories has no
    pGC's in the weight is reduced to 0.
    NOTE: There can never be no pGC's in any category.

    Args
    ----
    gp: The gene pool containing pgc.
    xgc_refs: Iterable of GC references to find a pGC for.
    depth: The layer in the environment to find a pGC.

    Returns
    -------
    The pGCs to evolve xgcs.
    """
    all_pgcs = positive_pgcs = negative_pgcs = None
    positive_category_weight = 2
    negative_category_weight = 1
    matched_pgcs = []
    for xgc_ref in xgc_refs:
        xgc = gms[xgc_ref]

        # Intergenerational pGC selection
        next_pgc_ref = xgc['next_pgc_ref']
        if next_pgc_ref is not None:
            xgc['next_pgc_ref'] = None
            matched_pgcs.append(gms.pool[next_pgc_ref])
        else:
            # Selection category selection
            effective_category_weight = 0 if xgc['effective_pgc_refs'] is None else 4
            _weights = (effective_category_weight, positive_category_weight, negative_category_weight)
            category_weights = array(_weights, dtype=float32)
            normalised_category_weights = category_weights / category_weights.sum()
            category = weighted_choice(3, p=normalised_category_weights)

            # Only do this once if it is needed as it is expensive
            if category > 0:
                if all_pgcs is None:
                    # TODO: Need a better data structure
                    # NOTE: *_weights have no bias so the probability of a pGC with fitness
                    # 0.501 being selected relative to a pGC of fitness 0.600 is 1% - that make sense
                    # at the time of writing.
                    all_pgcs = tuple(gc for gc in gms.pool.values() if is_pgc(gc))
                    if _LOG_DEBUG:
                        assert all_pgcs, 'There are no viable pGCs in the GP!'
                        _logger.debug(f'{len(all_pgcs)} pGCs in the local GP cache.')

                # If we end up with a category 3 we have to loop otherwise we can break out
                while True:
                    # pGC's that have had a net positive affect on target fitness
                    # Only do this once if it is needed as it is expensive
                    if positive_pgcs is None and category == 1:
                        positive_pgcs = tuple(gc for gc in all_pgcs if gc['pgc_fitness'][depth] > 0.5)
                        redo = False
                        if positive_pgcs:
                            positive_weights = array([pgc['pgc_fitness'][depth] for pgc in positive_pgcs], dtype=float32) - 0.5
                            positive_normalised_weights = positive_weights / positive_weights.sum()
                        else:
                            positive_category_weight = 0
                            category = 2
                        if _LOG_DEBUG:
                            _logger.debug(f'{len(positive_pgcs)} positive pGCs in the local GP cache.')
                            assert positive_pgcs and all(isfinite(positive_weights)), "Not all positive pGCs have a finite weight!"

                    if category == 1:
                        matched_pgcs.append(weighted_choice(positive_pgcs, p=positive_normalised_weights))
                        break

                    # pGC's that have had a net negative affect on target fitness
                    # Only do this once if it is needed as it is expensive
                    # Category must == 2 at this point
                    if negative_pgcs is None:
                        negative_pgcs = tuple(gc for gc in all_pgcs if gc['pgc_fitness'][depth] <= 0.5)
                        if negative_pgcs:
                            negative_weights = array([pgc['pgc_fitness'][depth] for pgc in negative_pgcs], dtype=float32)
                            negative_normalised_weights = negative_weights / negative_weights.sum()
                        else:
                            negative_category_weight = 0
                            category = 1
                        if _LOG_DEBUG:
                            _logger.debug(f'{len(negative_pgcs)} negative pGCs in the local GP cache.')
                            assert negative_pgcs and all(isfinite(negative_weights)), "Not all negative pGCs have a finite weight!"

                    if category == 2:
                        matched_pgcs.append(weighted_choice(negative_pgcs, p=negative_normalised_weights))
                        break

                    # An effective category == 3 (would have broken out of the while loop before now if it wasn't)
                    if _LOG_DEBUG:
                        _logger.debug(f'Category 3 selection event.')
                        assert category == 1, "In a category 3 event the next iteration must have category == 1!"
                        assert positive_pgcs is None, "Category 3 pGC selection can only be reached from a first iteration category 2 selection!"
                        assert not negative_pgcs, "Category 3 pGC selection can only be reached if there are no negative pGCs!"

            else:  # Category == 0
                fitness = array(xgc['effective_pgc_fitness'], dtype=float32)
                normalised_weights = fitness / fitness.sum()
                matched_pgcs.append(gms[weighted_choice(xgc['effective_pgc_refs'], p=normalised_weights)])

    return matched_pgcs


def pGC_inherit(child, parent, pgc):
    """Inherit PGC only properties from parent to child.

    child is modified.
    parent is modified.

    Parental protection works differently for a PGC than a population individual because a
    PGC child must survive before being characterized where as a population child is
    characterized immediately.

    Args
    ----
    child (xGC): Offspring of parent and pgc.
    parent (xGC): Parent of child.
    pgc (pGC): pGC that operated on parent to product child.
    """
    # TODO: A better data structure would be quicker
    child['pgc_fitness'] = [f * _PGC_PARENTAL_PROTECTION_FACTOR for f in parent['pgc_fitness']]
    child['pgc_f_count'] = [2] * NUM_PGC_LAYERS
    child['pgc_evolvability'] = [f * _PGC_PARENTAL_PROTECTION_FACTOR for f in parent['pgc_evolvability']]
    child['pgc_e_count'] = [2] * NUM_PGC_LAYERS

    child['_pgc_fitness'] = [0.0] * NUM_PGC_LAYERS
    child['_pgc_f_count'] = [0] * NUM_PGC_LAYERS
    child['_pgc_evolvability'] = [0.0] * NUM_PGC_LAYERS
    child['_pgc_e_count'] = [0] * NUM_PGC_LAYERS

    xGC_inherit(child, parent, pgc)


def population_GC_inherit(child, parent, pgc):
    """Inherit population properties from parent to child.

    child is modified.
    parent is modified.

    The child must have been characterized i.e. fitness and survivability defined.

    Parental protection for a population individual can only be of benefit.

    Args
    ----
    child (xGC): Offspring of parent and pgc.
    parent (xGC): Parent of child.
    pgc (pGC): pGC that operated on parent to product child.
    """
    if _LOG_DEBUG:
        if not all((field in child for field in ('fitness', 'survivability'))):
            raise ValueError('Child GC has not been characterized.')

    # There is no way of characterising first
    if parent['e_count'] == 1:
        child['evolvability'] = 1.0
        child['e_count'] = 1
    else:
        child['evolvability'] = parent['evolvability']
        child['e_count'] = max((2, parent['e_count'] >> 1))

    inherited_survivability = parent['survivability'] * _POPULATION_PARENTAL_PROTECTION_FACTOR
    inherited_fitness = parent['fitness'] * _POPULATION_PARENTAL_PROTECTION_FACTOR
    child['survivability'] = max((child['survivability'], inherited_survivability))
    child['fitness'] = max((child['fitness'], inherited_fitness))
    xGC_inherit(child, parent, pgc)


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
