"""The operations that can be performed on a GC."""
from copy import deepcopy
from logging import DEBUG, Logger, NullHandler, getLogger
from pprint import pformat
from typing import Literal

from egp_stores.gene_pool import gene_pool
from egp_types.xGC import xGC
from egp_types.gc_graph import DST_EP, SRC_EP, gc_graph
from egp_types.internal_graph import internal_graph
from egp_types.reference import ref_str
from egp_types.aGC import aGC
from egp_types.dGC import dGC
from .egp_typing import WorkStack, Work, NewGCDef, InsertRow

from .fundamental import default_dict_gc, steady_state_exception


_logger: Logger = getLogger(__name__)
_logger.addHandler(NullHandler())
_LOG_DEBUG: bool = _logger.isEnabledFor(DEBUG)
_SLASH_N = '\n'


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
