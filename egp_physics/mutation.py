"""General mutations."""
from __future__ import annotations
from logging import DEBUG, Logger, NullHandler, getLogger
from typing import TYPE_CHECKING

from egp_types.xGC import xGC
from egp_types.internal_graph import internal_graph
from egp_types.reference import ref_str
from egp_types.aGC import aGC
from egp_types.dGC import dGC
from egp_types.egp_typing import Row

from .fundamental import clone, pgc_epilogue


_logger: Logger = getLogger(__name__)
_logger.addHandler(NullHandler())
_LOG_DEBUG: bool = _logger.isEnabledFor(DEBUG)


# Circular import & runtime import avoidance
if TYPE_CHECKING:
    from egp_stores.gene_pool import gene_pool


def gc_remove(gms: gene_pool, tgc: aGC, row: Row) -> xGC:
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
    rgc: dGC = clone(tgc, gms.next_reference, False)
    if _LOG_DEBUG:
        _logger.debug(
            f"Minimally cloned {ref_str(tgc['ref'])} to {ref_str(rgc['ref'])}"
        )
        _logger.debug(f"Removing row {row}.")
    rgc_igraph: internal_graph = internal_graph()
    rgc_igraph.update(tgc["gc_graph"].i_graph.copy_row("I"))
    if row != "C":
        rgc_igraph.update(tgc["gc_graph"].i_graph.copy_row("C"))
    if row == "A":
        rgc["gca_ref"] = tgc["gcb_ref"]
        rgc["gcb_ref"] = 0
        if tgc["gc_graph"].has_row("B"):
            rgc_igraph.update(
                tgc["gc_graph"].i_graph.move_row("B", "A", has_f=tgc["gc_graph"].has_row("F"))
            )
    elif row == "B":
        rgc["gca_ref"] = tgc["gca_ref"]
        rgc_igraph.update(tgc["gc_graph"].i_graph.copy_row("A"))
        rgc["gcb_ref"] = 0
    elif row == "P":
        rgc_igraph.update(tgc["gc_graph"].i_graph.copy_row("O"))
    elif row == "O":
        rgc_igraph.update(tgc["gc_graph"].i_graph.move_row("P", "O"))
    rgc["gc_graph"].normalize()
    return pgc_epilogue(gms, rgc)


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
    dgc: dGC = clone(tgc, gms.next_reference)
    if _LOG_DEBUG:
        _logger.debug(
            f"Minimally cloned {ref_str(tgc['ref'])} to {ref_str(dgc['ref'])}"
        )
    dgc["gc_graph"].remove_all_connections()
    dgc["gc_graph"].normalize()
    return pgc_epilogue(gms, dgc)


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
    dgc: dGC = clone(tgc, gms.next_reference)
    if _LOG_DEBUG:
        _logger.debug(
            f"Minimally cloned {ref_str(tgc['ref'])} to {ref_str(dgc['ref'])}"
        )
    dgc["gc_graph"].add_input()
    dgc["gc_graph"].normalize()
    return pgc_epilogue(gms, dgc)


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
    dgc: dGC = clone(tgc, gms.next_reference)
    if _LOG_DEBUG:
        _logger.debug(
            f"Minimally cloned {ref_str(tgc['ref'])} to {ref_str(dgc['ref'])}"
        )
    dgc["gc_graph"].remove_input()
    dgc["gc_graph"].normalize()
    return pgc_epilogue(gms, dgc)


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
    dgc: dGC = clone(tgc, gms.next_reference)
    if _LOG_DEBUG:
        _logger.debug(
            f"Minimally cloned {ref_str(tgc['ref'])} to {ref_str(dgc['ref'])}"
        )
    dgc["gc_graph"].add_output()
    dgc["gc_graph"].normalize()
    return pgc_epilogue(gms, dgc)


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
    dgc: dGC = clone(tgc, gms.next_reference)
    if _LOG_DEBUG:
        _logger.debug(
            f"Minimally cloned {ref_str(tgc['ref'])} to {ref_str(dgc['ref'])}"
        )
    dgc["gc_graph"].remove_output()
    dgc["gc_graph"].normalize()
    return pgc_epilogue(gms, dgc)


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
    dgc: dGC = clone(tgc, gms.next_reference)
    if _LOG_DEBUG:
        _logger.debug(
            f"Minimally cloned {ref_str(tgc['ref'])} to {ref_str(dgc['ref'])}"
        )
    dgc["gc_graph"].remove_constant()
    dgc["gc_graph"].normalize()
    return pgc_epilogue(gms, dgc)
