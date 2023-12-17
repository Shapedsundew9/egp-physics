"""The operations that can be performed on a GC."""
from copy import deepcopy
from logging import DEBUG, Logger, NullHandler, getLogger
from typing import Literal, LiteralString, cast, Callable
from random import randint

from egp_stores.gene_pool import gene_pool
from egp_types.xGC import xGC
from egp_types.gc_graph import DST_EP, SRC_EP, gc_graph
from egp_types.internal_graph import internal_graph
from egp_types.reference import ref_str
from egp_types.ep_type import interface_definition, vtype
from egp_types.aGC import aGC
from egp_types.dGC import dGC
from egp_types.egp_typing import Row, VALID_ROW_SOURCES
from .egp_typing import WorkStack, Work, NewGCDef, InsertRow


_logger: Logger = getLogger(__name__)
_logger.addHandler(NullHandler())
_LOG_DEBUG: bool = _logger.isEnabledFor(DEBUG)
_SLASH_N = "\n"


# Filler
_EMPTY_GC_GRAPH = gc_graph()


# Steady state exception filters.
_LMT: LiteralString = (
    " AND NOT ({exclude_column} = ANY({exclusions})) ORDER BY RANDOM() LIMIT 1"
)

# TODO: Replace with a localisation hash?
_IT: LiteralString = "{input_types}"
_OT: LiteralString = "{output_types}"
_ITS: LiteralString = "{itypes}::SMALLINT[]"
_OTS: LiteralString = "{otypes}::SMALLINT[]"
_IDX: LiteralString = "{inputs} = {iidx}"
_ODX: LiteralString = "{outputs} = {oidx}"

_MATCH_TYPE_0_SQL: LiteralString = (
    "WHERE "
    + _IT
    + " = "
    + _ITS
    + " AND "
    + _IDX
    + " AND "
    + _OT
    + " = "
    + _OTS
    + " AND "
    + _ODX
    + _LMT
)
_MATCH_TYPE_1_SQL: LiteralString = (
    "WHERE " + _IT + " = " + _ITS + " AND " + _OT + " = " + _OTS + " AND " + _ODX + _LMT
)
_MATCH_TYPE_2_SQL: LiteralString = (
    "WHERE " + _IT + " = " + _ITS + " AND " + _IDX + " AND " + _OT + " = " + _OTS + _LMT
)
_MATCH_TYPE_3_SQL: LiteralString = (
    "WHERE " + _IT + " = " + _ITS + " AND " + _OT + " = " + _OTS + _LMT
)
_MATCH_TYPE_4_SQL: LiteralString = (
    "WHERE " + _IT + " <@ " + _ITS + " AND " + _OT + " = " + _OTS + _LMT
)
_MATCH_TYPE_5_SQL: LiteralString = (
    "WHERE " + _IT + " <@ " + _ITS + " AND " + _OT + " @> " + _OTS + _LMT
)
_MATCH_TYPE_6_SQL: LiteralString = (
    "WHERE " + _IT + " <@ " + _ITS + " AND " + _OT + " && " + _OTS + _LMT
)
_MATCH_TYPE_7_SQL: LiteralString = (
    "WHERE " + _IT + " && " + _ITS + " AND " + _OT + " && " + _OTS + _LMT
)
_MATCH_TYPE_8_SQL: LiteralString = "WHERE " + _OT + " && " + _OTS + " " + _LMT
_MATCH_TYPE_9_SQL: LiteralString = "WHERE " + _IT + " && " + _ITS + " " + _LMT
# Catch for when xtypes is an empty set.
_MATCH_TYPE_10_SQL: LiteralString = "WHERE " + _OT + " = " + _OTS + " " + _LMT
_MATCH_TYPE_11_SQL: LiteralString = "WHERE " + _IT + " = " + _ITS + " " + _LMT


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
    _MATCH_TYPE_11_SQL,
)
_NUM_MATCH_TYPES: int = len(_MATCH_TYPES_SQL)


def default_dict_gc(ref: Callable[[], int]) -> dGC:
    """Create a default dGC with ref generated from ref()."""
    return {
        "gc_graph": _EMPTY_GC_GRAPH,
        "ancestor_a_ref": 0,
        "ancestor_b_ref": 0,
        "ref": ref(),
        "gca_ref": 0,
        "gcb_ref": 0,
    }


def _insert_graph_case_0(
    tig: internal_graph, iig: internal_graph, rig: internal_graph
) -> None:
    """Insert igc into tgc case 0."""
    _logger.debug("Case 0: Stack")
    rig.update(iig.copy_row("I", True))
    rig.update(iig.as_row("A"))
    rig.update(tig.as_row("B"))
    rig.update(tig.copy_row("O", True))
    rig.update(rig.direct_connect("I", "A"))
    rig.update(rig.direct_connect("B", "O"))


def _insert_graph_case_1(
    tig: internal_graph, iig: internal_graph, rig: internal_graph
) -> None:
    """Insert igc into tgc case 1."""
    _logger.debug("Case 1: No row A or B")
    rig.update(iig.as_row("A"))
    rig.update(tig.copy_row("O"))


def _insert_graph_case_2(
    tig: internal_graph, iig: internal_graph, rig: internal_graph
) -> None:
    """Insert igc into tgc case 2."""
    _logger.debug("Case 2: No row B and insert above A")
    rig.update(iig.as_row("A"))
    rig.update(tig.move_row("A", "B"))
    rig.update(tig.copy_row("O"))
    rig.redirect_refs("O", DST_EP, "A", "B")


def _insert_graph_case_3(
    tig: internal_graph, iig: internal_graph, rig: internal_graph
) -> None:
    """Insert igc into tgc case 3."""
    _logger.debug("Case 3: No row B and insert below A")
    rig.update(tig.copy_rows(("A", "O")))
    rig.update(iig.as_row("B"))


def _insert_graph_case_4(
    tig: internal_graph, iig: internal_graph, _: internal_graph, fig: internal_graph
) -> None:
    """Insert igc into tgc case 4."""
    _logger.debug("Case 4: Has rows A & B and insert above A")
    fig.update(tig.embed("A", "B"))
    fig.update(iig.as_row("A"))


def _insert_graph_case_5(
    tig: internal_graph, iig: internal_graph, _: internal_graph, fig: internal_graph
) -> None:
    """Insert igc into tgc case 5."""
    _logger.debug("Case 5: Has rows A & B and insert above B")
    fig.update(tig.embed("A", "A"))
    fig.update(iig.as_row("B"))


def _insert_graph_case_6(
    tig: internal_graph, iig: internal_graph, _: internal_graph, fig: internal_graph
) -> None:
    """Insert igc into tgc case 6."""
    _logger.debug("Case 6: Has rows A & B and insert above O")
    fig.update(tig.embed("B", "A"))
    fig.update(iig.as_row("B"))


def _insert_graph_case_7(
    tig: internal_graph, iig: internal_graph, fig: internal_graph
) -> None:
    """Insert igc into tgc case 7."""
    fig.update(tig.embed("A", "B"))
    fig.update(iig.as_row("A"))


def _insert_graph_case_8(
    tig: internal_graph, iig: internal_graph, fig: internal_graph
) -> None:
    """Insert igc into tgc case 8."""
    fig.update(tig.embed("A", "A"))
    fig.update(iig.as_row("B"))


def _insert_graph_case_9(
    tig: internal_graph, iig: internal_graph, fig: internal_graph
) -> None:
    """Insert igc into tgc case 9."""
    fig.update(tig.embed("B", "B"))
    fig.update(iig.as_row("A"))


def _insert_graph_case_10(
    tig: internal_graph, iig: internal_graph, fig: internal_graph
) -> None:
    """Insert igc into tgc case 10."""
    fig.update(tig.embed("B", "A"))
    fig.update(iig.as_row("B"))


def _insert_graph_case_11(
    tig: internal_graph, iig: internal_graph, rig: internal_graph
) -> None:
    """Insert igc into tgc case 11."""
    _logger.debug("Case 11: Inverse Stack")
    rig.update(tig.copy_row("I", True))
    rig.update(tig.as_row("A"))
    rig.update(iig.as_row("B"))
    rig.update(iig.copy_row("O", True))
    rig.update(rig.direct_connect("I", "A"))
    rig.update(rig.direct_connect("B", "O"))


def _insert_graph(
    tgcg: gc_graph, igcg: gc_graph, above_row: InsertRow
) -> tuple[gc_graph, gc_graph]:
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
    if above_row != "I" and above_row != "Z":
        if tgcg.has_row("F"):
            rig = deepcopy(tig)
        else:
            rig.update(tig.copy_rows_src_eps(("I", "C"), True))

    # TODO: There are opportunities to reduce overhead by making some internal_graph manipulation functions
    # act on self rather than returning a dictionary to update (into self)
    if above_row == "I":
        _insert_graph_case_0(tig, iig, rig)
    elif above_row == "Z":
        _insert_graph_case_11(tig, iig, rig)
    elif not tgcg.has_row("A"):
        _insert_graph_case_1(tig, iig, rig)
    elif not tgcg.has_row("F"):
        if not tgcg.has_row("B"):
            if above_row == "A":
                _insert_graph_case_2(tig, iig, rig)
            else:
                _insert_graph_case_3(tig, iig, rig)
        else:
            rig = deepcopy(tig)
            if above_row == "A":
                _insert_graph_case_4(tig, iig, rig, fig)
            elif above_row == "B":
                _insert_graph_case_5(tig, iig, rig, fig)
            else:
                _insert_graph_case_6(tig, iig, rig, fig)
    else:
        if above_row == "A":
            _insert_graph_case_7(tig, iig, fig)
        elif above_row == "O":
            _insert_graph_case_8(tig, iig, fig)
        elif above_row == "B":
            _insert_graph_case_9(tig, iig, fig)
        elif above_row == "P":
            _insert_graph_case_10(tig, iig, fig)

    # Pre-completed logging
    if _LOG_DEBUG:
        _logger.debug(f"tgc ({type(tgcg)}):\n{tgcg.i_graph}")
        _logger.debug(f"igc ({type(tgcg)}):\n{igcg.i_graph}")
        _logger.debug(f"Pre-completed rig ({type(rig)}):\n{rig}")
        _logger.debug(f"Pre-completed fig ({type(fig)}):\n{fig}")

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
        _logger.debug(f"Completed rig:\n{rig}")
        _logger.debug(f"Completed fig:\n{fig}")

    return rgc_gc_graph, fgc_gc_graph


def _insert_gc(gms: gene_pool, tgc: aGC, igc: aGC, above_row: InsertRow, stabilize: bool = True) -> NewGCDef:
    """Insert igc into tgc above row 'above_row'.

    Args
    ----
    gms: A source of genetic material & references. Needs to be in the context of the sub-process.
    tgc: GC to insert insert_gc into.
    igc: GC to insert into target_gc.
    above_row: Insert above this row.
    stablize: If True the graph will be stablized before returning.

    Returns
    -------
    (rgc, {ref: fgc}): First fgc is the igc followed by fgc's
    created to stabilise rgc.
    """
    if _LOG_DEBUG:
        _logger.debug("Inserting into Target GC.")
    return _recursive_insert_gc(gms, [(tgc, igc, above_row)], stabilize)


def _insert_gc_case_0(tgc: aGC, igc: aGC, rgc: aGC) -> None:
    """Insert igc data into tgc case 0."""
    _logger.debug("Case 0: Stack")
    rgc["gca_ref"] = igc["ref"]
    rgc["gcb_ref"] = tgc["ref"]


def _insert_gc_case_1(igc: aGC, rgc: aGC) -> None:
    """Insert igc data into tgc case 1."""
    _logger.debug("Case 1")
    rgc["gca_ref"] = igc["ref"]


def _insert_gc_case_2(tgc: aGC, igc: aGC, rgc: aGC) -> None:
    """Insert igc data into tgc case 2."""
    _logger.debug("Case 2")
    rgc["gca_ref"] = igc["ref"]
    if tgc["gca_ref"] is not None:
        rgc["gcb_ref"] = tgc["gca_ref"]
    else:
        rgc["gcb_ref"] = tgc["ref"]


def _insert_gc_case_3(tgc: aGC, igc: aGC, rgc: aGC) -> None:
    """Insert igc data into tgc case 3."""
    _logger.debug("Case 3")
    if tgc["gca_ref"] is not None:
        rgc["gca_ref"] = tgc["gca_ref"]
    else:
        rgc["gca_ref"] = tgc["ref"]
    rgc["gcb_ref"] = igc["ref"]


def _insert_gc_case_4(tgc: aGC, igc: aGC, rgc: aGC, fgc: aGC) -> None:
    """Insert igc data into tgc case 4."""
    _logger.debug("Case 4")
    fgc["gca_ref"] = igc["ref"]
    fgc["gcb_ref"] = tgc["gca_ref"]
    rgc["gca_ref"] = fgc["ref"]
    rgc["gcb_ref"] = tgc["gcb_ref"]


def _insert_gc_case_5(tgc: aGC, igc: aGC, rgc: aGC, fgc: aGC) -> None:
    """Insert igc data into tgc case 5."""
    _logger.debug("Case 5")
    fgc["gca_ref"] = tgc["gca_ref"]
    fgc["gcb_ref"] = igc["ref"]
    rgc["gca_ref"] = fgc["ref"]
    rgc["gcb_ref"] = tgc["gcb_ref"]


def _insert_gc_case_6(tgc: aGC, igc: aGC, rgc: aGC, fgc: aGC) -> None:
    """Insert igc data into tgc case 6."""
    _logger.debug("Case 6")
    fgc["gca_ref"] = tgc["gcb_ref"]
    fgc["gcb_ref"] = igc["ref"]
    rgc["gca_ref"] = tgc["gca_ref"]
    rgc["gcb_ref"] = fgc["ref"]


def _insert_gc_case_7(tgc: aGC, igc: aGC, rgc: aGC, fgc: aGC) -> None:
    """Insert igc data into tgc case 7."""
    _logger.debug("Case 7")
    rgc["gca_ref"] = fgc["ref"]
    rgc["gcb_ref"] = tgc["gcb_ref"]
    fgc["gca_ref"] = igc["ref"]
    fgc["gcb_ref"] = tgc["gca_ref"]


def _insert_gc_case_8(tgc: aGC, igc: aGC, rgc: aGC, fgc: aGC) -> None:
    """Insert igc data into tgc case 8."""
    _logger.debug("Case 8")
    rgc["gca_ref"] = fgc["ref"]
    rgc["gcb_ref"] = tgc["gcb_ref"]
    fgc["gca_ref"] = tgc["gca_ref"]
    fgc["gcb_ref"] = igc["ref"]


def _insert_gc_case_9(tgc: aGC, igc: aGC, rgc: aGC, fgc: aGC) -> None:
    """Insert igc data into tgc case 9."""
    _logger.debug("Case 9")
    rgc["gca_ref"] = tgc["gca_ref"]
    rgc["gcb_ref"] = fgc["ref"]
    fgc["gca_ref"] = igc["ref"]
    fgc["gcb_ref"] = tgc["gcb_ref"]


def _insert_gc_case_10(tgc: aGC, igc: aGC, rgc: aGC, fgc: aGC) -> None:
    """Insert igc data into tgc case 10."""
    _logger.debug("Case 10")
    rgc["gca_ref"] = tgc["gca_ref"]
    rgc["gcb_ref"] = fgc["ref"]
    fgc["gca_ref"] = tgc["gcb_ref"]
    fgc["gcb_ref"] = igc["ref"]


def _insert_gc_case_11(tgc: aGC, igc: aGC, rgc: aGC) -> None:
    """Insert igc data into tgc case 11."""
    _logger.debug("Case 11: Inverse Stack")
    rgc["gca_ref"] = tgc["ref"]
    rgc["gcb_ref"] = igc["ref"]


def _recursive_insert_gc(gms: gene_pool, work_stack: WorkStack, stablize: bool = True) -> NewGCDef:
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
            _logger.debug(
                f"Work: Target={ref_str(tgc.get('ref', 0))}, "
                f"Insert={ref_str(igc.get('ref', 0))}, Above Row={above_row}"
            )

        # Insert into the graph
        tgcg: gc_graph = tgc["gc_graph"]
        igcg: gc_graph = igc["gc_graph"]
        fgcg: gc_graph
        rgcg: gc_graph
        rgcg, fgcg = _insert_graph(tgcg, igcg, above_row)
        rgcg_steady: bool = rgcg.normalize()
        fgcg_steady: bool = False
        if fgcg:
            fgcg_steady = fgcg.normalize()
            if _LOG_DEBUG:
                _logger.debug(f"Normalized fgc:\n{fgcg.i_graph}")
            fgc["graph"] = fgcg.app_graph
            fgc["gc_graph"] = fgcg
            fgc["ancestor_a_ref"] = igc["ref"]
            fgc["ancestor_b_ref"] = tgc["ref"]
        if _LOG_DEBUG:
            _logger.debug(f"Normalized rgc:\n{rgcg.i_graph}")
        rgc["graph"] = rgcg.app_graph
        rgc["gc_graph"] = rgcg
        rgc["ancestor_a_ref"] = tgc["ref"]
        rgc["ancestor_b_ref"] = igc["ref"]

        # Insert into the GC
        # The insert_gc is always referenced in the tree of the final rgc
        fgc_dict[igc["ref"]] = igc
        fgc_dict[tgc["ref"]] = tgc
        if above_row == "I":
            _insert_gc_case_0(tgc, igc, rgc)
        elif above_row == "Z":
            _insert_gc_case_11(tgc, igc, rgc)
        elif not tgcg.has_row("A"):
            _insert_gc_case_1(igc, rgc)
        elif not tgcg.has_row("F"):
            if not tgcg.has_row("B"):
                if above_row == "A":
                    _insert_gc_case_2(tgc, igc, rgc)
                else:
                    _insert_gc_case_3(tgc, igc, rgc)
            else:
                if above_row == "A":
                    _insert_gc_case_4(tgc, igc, rgc, fgc)
                elif above_row == "B":
                    _insert_gc_case_5(tgc, igc, rgc, fgc)
                else:
                    _insert_gc_case_6(tgc, igc, rgc, fgc)
        else:
            if above_row == "A":
                _insert_gc_case_7(tgc, igc, rgc, fgc)
            elif above_row == "O":
                _insert_gc_case_8(tgc, igc, rgc, fgc)
            elif above_row == "B":
                _insert_gc_case_9(tgc, igc, rgc, fgc)
            elif above_row == "P":
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
        new_ref: int = rgc["ref"]
        old_ref: int = tgc["ref"]
        if _LOG_DEBUG:
            _logger.debug(f"Replacing {ref_str(old_ref)} with {ref_str(new_ref)}.")
        for nfgc in fgc_dict.values():
            for ref in ("gca_ref", "gcb_ref", "ancestor_a_ref", "ancestor_b_ref"):
                if nfgc[ref] == old_ref:
                    nfgc[ref] = new_ref
        if new_tgc_flag:
            for ref in ("gca_ref", "gcb_ref", "ancestor_a_ref", "ancestor_b_ref"):
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
                fgc_dict[fgc["ref"]] = fgc

        # We may not want to immediately stabilise the graph
        if not rgcg_steady and stablize:
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
            # Order of addition important.
            fgc_dict[rgc["ref"]] = rgc

        if _LOG_DEBUG:
            _logger.debug(f"fgc_dict: {[ref_str(x) for x in fgc_dict]}")

    if _LOG_DEBUG:
        _logger.debug(
            f"fgc_dict details:\n{_SLASH_N.join(ref_str(k) + ':' + _SLASH_N + str(v) for k, v in fgc_dict.items())}"
        )
        # TODO: target_gc & new_tgc interface must be the same. Validate.

    return (new_tgc, fgc_dict)


def gc_insert(
    gms: gene_pool, tgc: aGC, igc: aGC, above_row: Literal["I", "A", "B", "O"]
) -> xGC:
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
    gms.pool[new_gc_definition[0]["ref"]] = new_gc_definition[0]
    gms.pool.update(new_gc_definition[1])
    return gms.pool[new_gc_definition[0]["ref"]]


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
    new_gc_definition: NewGCDef = _insert_gc(gms, bottom_gc, top_gc, ("I", "Z")[invert])
    gms.pool[new_gc_definition[0]["ref"]] = new_gc_definition[0]
    gms.pool.update(new_gc_definition[1])
    return gms.pool[new_gc_definition[0]["ref"]]


def interface_proximity_select(
    gms: gene_pool, xputs: dict[str, bytes | list | str]
) -> aGC:
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
            _logger.debug(
                f"Proximity selection match_type {match_type} found no candidates."
            )
        match_type += 1
        agc = tuple(gms.select(_MATCH_TYPES_SQL[match_type], literals=xputs))
    if _LOG_DEBUG:
        _logger.debug(f"Proximity selection match_type {match_type} found a candidate.")
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
    """
    if _LOG_DEBUG:
        _logger.debug(
            f"Steady state exception thrown for GC ref {ref_str(fgc['ref'])}."
        )
    fgc_graph: gc_graph = fgc["gc_graph"]

    # Find unconnected destination endpoints. Determine highest row & endpoint types.
    above_row: str = "Z"
    outputs: list[int] = []
    for ep in fgc_graph.i_graph.dst_unref_filter():
        if ep.row < above_row:
            above_row = ep.row
        outputs.append(ep.typ)

    # Find viable source types above the highest row.
    inputs: list[int] = [
        ep.typ
        for ep in fgc_graph.i_graph.src_rows_filter(
            VALID_ROW_SOURCES[fgc["gc_graph"].has_row("F")][cast(Row, above_row)]
        )
    ]

    xputs: dict[str, bytes | list | str] = {
        "exclude_column": "signature",
        "exclusions": list(),
    }
    _, xputs["itypes"], xputs["iidx"] = interface_definition(inputs, vtype.EP_TYPE_INT)
    _, xputs["otypes"], xputs["oidx"] = interface_definition(outputs, vtype.EP_TYPE_INT)

    # Find a gc based on the criteria
    insert_gc: aGC = interface_proximity_select(gms, xputs)
    return (fgc, insert_gc, cast(InsertRow, above_row))
