"""The operation that can be performed on a GC dictionary."""
from copy import copy, deepcopy
from logging import DEBUG, NullHandler, getLogger
from pprint import pformat
from random import randint, choice
from numpy import array, float32

from .ep_type import vtype
from .gc_graph import (DST_EP, SRC_EP, ep_idx, gc_graph, hash_ep, hash_ref,
                       ref_idx)
from .gc_type import eGC, interface_definition, mGC, M_MASK, PHYSICAL_PROPERTY, LAYER_COLUMNS, LAYER_COLUMNS_RESET
from .utils.reference import random_reference

_logger = getLogger(__name__)
_logger.addHandler(NullHandler())
_LOG_DEBUG = _logger.isEnabledFor(DEBUG)


# Steady state exception filters.
_EXCLUSION_LIMIT =  ' AND NOT ({exclude_column} = ANY({exclusions})) ORDER BY RANDOM() LIMIT 1'

# FIXME: Replace with a localisation hash
_MATCH_TYPE_0_SQL = ('WHERE {input_types} = {itypes} AND {inputs} = {iidx} AND {output_types} = {otypes} AND {outputs} = {oidx}'
                     + _EXCLUSION_LIMIT)
_MATCH_TYPE_1_SQL = 'WHERE {input_types} = {itypes} AND {output_types} = {otypes} AND {outputs} = {oidx}' + _EXCLUSION_LIMIT
_MATCH_TYPE_2_SQL = 'WHERE {input_types} = {itypes} AND {inputs} = {iidx} AND {output_types} = {otypes}' + _EXCLUSION_LIMIT
_MATCH_TYPE_3_SQL = 'WHERE {input_types} = {itypes} AND {output_types} = {otypes}' + _EXCLUSION_LIMIT
_MATCH_TYPE_4_SQL = 'WHERE {input_types} <@ {itypes} AND {output_types} = {otypes}' + _EXCLUSION_LIMIT
_MATCH_TYPE_5_SQL = 'WHERE {input_types} <@ {itypes} AND {output_types} @> {otypes}' + _EXCLUSION_LIMIT
_MATCH_TYPE_6_SQL = 'WHERE {input_types} <@ {itypes} AND {output_types} && {otypes}' + _EXCLUSION_LIMIT
_MATCH_TYPE_7_SQL = 'WHERE {input_types} && {itypes} AND {output_types} && {otypes}' + _EXCLUSION_LIMIT
_MATCH_TYPE_8_SQL = 'WHERE {output_types} && {otypes} ' + _EXCLUSION_LIMIT
_MATCH_TYPE_9_SQL = 'WHERE {input_types} && {itypes} ' + _EXCLUSION_LIMIT
_MATCH_TYPES_SQL = (
    _MATCH_TYPE_0_SQL,
    _MATCH_TYPE_1_SQL,
    _MATCH_TYPE_2_SQL,
    _MATCH_TYPE_3_SQL,
    _MATCH_TYPE_4_SQL,
    _MATCH_TYPE_5_SQL,
    _MATCH_TYPE_6_SQL,
    _MATCH_TYPE_7_SQL,
    _MATCH_TYPE_8_SQL,
    _MATCH_TYPE_9_SQL
)
_NUM_MATCH_TYPES = len(_MATCH_TYPES_SQL)


# Initial GC query
_INITIAL_GC_SQL = ('WHERE {input_types} = {itypes} AND {inputs} = {iidx} '
                   'AND {output_types} = {otypes} AND {outputs} = {oidx} '
                   'AND NOT ({exclude_column} = ANY({exclusions})) ORDER BY RANDOM() LIMIT {limit}')

RANDOM_PGC_SIGNATURE = b'00000000000000000000000000000000'

def _copy_row(igc, rows, ep_type=None):
    """Copy the internal format definition of a row.

    If ep_type is None all endpoints regardless of end point type are copied.

    Args
    ----
    igc (dict): Internal gc_graph format gc_graph.
    rows (str): Valid row letters as a string e.g. 'IC'
    ep_type (bool): SRC_EP or DST_EP

    Returns
    -------
    (dict): An internal graph format dictionary containing the row.
    """
    if ep_type is not None:
        def filter_func(x): return x[1][ep_idx.EP_TYPE] == ep_type and x[1][ep_idx.ROW] in rows
    else:
        def filter_func(x): return x[1][ep_idx.ROW] in rows
    return {k: deepcopy(ep) for k, ep in filter(filter_func, igc.items())}


def _copy_clean_row(igc, rows, ep_type=None):
    """Copy the internal format definition of a row removing references.

    If ep_type is None all endpoints regardless of end point type are copied.

    Args
    ----
    igc (dict): Internal gc_graph format gc_graph.
    rows (str): Valid row letters as a string e.g. 'IC'
    ep_type (bool): SRC_EP or DST_EP

    Returns
    -------
    (dict): An internal graph format dictionary containing the clean row.
    """
    if ep_type is not None:
        def filter_func(x): return x[1][ep_idx.EP_TYPE] == ep_type and x[1][ep_idx.ROW] in rows
    else:
        def filter_func(x): return x[1][ep_idx.ROW] in rows
    copied_row = {k: copy(ep) for k, ep in filter(filter_func, igc.items())}
    for ep in copied_row.values():
        ep[ep_idx.REFERENCED_BY] = []
    return copied_row


def _move_row(igc, src_row, src_ep_type, dst_row, dst_ep_type, clean=False):
    """Move a row definition to a different row.

    The endpoints moved are filtered by src_row & ep_type.
    If ep_type is None all endpoints regardless of end point type are moved.
    The moved rows are cleaned of references.

    Args
    ----
    igc (dict): Internal gc_graph format gc_graph.
    src_row (str): A valid row letter
    src_ep_type (bool or None): SRC_EP or DST_EP or None
    dst_row (str): A valid row letter
    dst_ep_type (bool or None): SRC_EP or DST_EP or None
    clean (bool): Remove references in dst_row if True

    Returns
    -------
    (dict): A gc_graph internal format containing the destination row endpoints.
    """
    if src_ep_type is not None:
        def filter_func(x): return x[ep_idx.EP_TYPE] == src_ep_type and x[ep_idx.ROW] == src_row
    else:
        def filter_func(x): return x[ep_idx.ROW] == src_row
    dst_eps = [deepcopy(ep) for ep in filter(filter_func, igc.values())]
    if _LOG_DEBUG:
        _logger.debug("Moving {} to row {} ep_type {}".format(dst_eps, dst_row, dst_ep_type))
    for ep in dst_eps:
        ep[ep_idx.ROW] = dst_row
        if clean:
            ep[ep_idx.REFERENCED_BY] = []
    if dst_ep_type is not None:
        for ep in dst_eps:
            ep[ep_idx.EP_TYPE] = dst_ep_type
    return {hash_ep(ep): ep for ep in dst_eps}


def _direct_connect(igc, src_row, dst_row):
    """Create dst_row and directly connect it to src_row.

    A direct connection means that dst_row has the same number, gc type and
    order of destination end points as src_row has source endpoints.
    src_row SRC_EPs should exist in fgc (else no-op).
    dst_row DST_EPs will be created in fgc.
    If any dst_row DST_EPs exist in fgc behaviour is undefined.

    Args
    ----
    igc (dict): Internal gc_graph format dict to be updated.
    src_row (str): 'I', 'C', 'A', or 'B'
    dst_row (str): Valid destination for src_row.

    Returns
    -------
    (dict): A gc_graph internal format containing the destination row endpoints.
    """
    connected_row = {}
    def filter_func(x): return x[ep_idx.EP_TYPE] and x[ep_idx.ROW] == src_row
    for src_ep in tuple(filter(filter_func, igc.values())):
        dst_ep = [DST_EP, dst_row, src_ep[ep_idx.INDEX], src_ep[ep_idx.TYPE], [[src_row, src_ep[ep_idx.INDEX]]]]
        connected_row[hash_ep(dst_ep)] = dst_ep
    return connected_row


def _append_connect(igc, src_row, dst_row):
    """Append endpoints to dst_row and directly connect them to src_row.

    A direct connection means that dst_row has the same number, gc type and
    order of destination end points as src_row has source endpoints.
    src_row SRC_EPs should exist in fgc (else no-op).
    dst_row DST_EPs should exist in fgc.

    Args
    ----
    igc (dict): Internal gc_graph format dict to be updated.
    src_row (str): 'I', 'C', 'A', or 'B'
    dst_row (str): Valid destination for src_row.

    Returns
    -------
    (dict): A gc_graph internal format containing the destination row endpoints.
    """
    connected_row = {}

    # Find the next endpoint index in the destination row
    def filter_func(x): return not x[ep_idx.EP_TYPE] and x[ep_idx.ROW] == dst_row
    next_idx = max([dst_ep[ep_idx.INDEX] for dst_ep in filter(filter_func, igc.values())]) + 1

    # Append a destination endpoint for every source endpoint
    def filter_func(x): return x[ep_idx.EP_TYPE] and x[ep_idx.ROW] == src_row
    for src_ep in tuple(filter(filter_func, igc.values())):
        dst_ep = [DST_EP, dst_row, next_idx, src_ep[ep_idx.TYPE], [[src_row, src_ep[ep_idx.INDEX]]]]
        connected_row[hash_ep(dst_ep)] = dst_ep
        next_idx += 1
    return connected_row


def _redirect_refs(igc, row, ep_type, old_ref_row, new_ref_row):
    """Redirects references on row from old_ref_row to new_ref_row.

    Args
    ----
    igc (dict): Internal gc_graph format dict to be updated.
    row (str): A valid row
    ep_type (bool): The old_ref_row ep_type.
    old_ref_row (str): A valid row
    new_ref_row (str): A valid row

    Returns
    -------
    (dict): Modified igc
    """
    def filter_func(x): return x[ep_idx.EP_TYPE] == ep_type and x[ep_idx.ROW] == row
    for ep in filter(filter_func, igc.values()):
        for ref in ep[ep_idx.REFERENCED_BY]:
            if ref[ref_idx.ROW] == old_ref_row:
                ref[ref_idx.ROW] = new_ref_row
    return igc


def _insert_as(igc, row):
    """Create an internal gc_format dict with igc as row.

    Args
    ----
    igc (dict): Internal gc_graph format gc_graph.
    row (str): 'A' or 'B'

    Returns
    -------
    (dict): Internal gc_format dict with igc as row.
    """
    ret_val = {}
    for ep in filter(lambda x: x[ep_idx.ROW] in ('I', 'O'), igc.values()):
        cep = copy(ep)
        cep[ep_idx.ROW] = row
        cep[ep_idx.EP_TYPE] = not ep[ep_idx.EP_TYPE]
        cep[ep_idx.REFERENCED_BY] = []
        ret_val[hash_ep(cep)] = cep
    return ret_val


def _complete_references(igc):
    """Complete any incomplete references in igc.

    An incomplete reference is when a destination references
    a source but the source does not reference the destination.
    This is usually a side effect of insertion.
    igc is modified.

    Args
    ----
    igc (dict): Internal gc_graph format gc_graph.
    """
    for dst_ep in filter(lambda x: x[ep_idx.EP_TYPE] == DST_EP and x[ep_idx.REFERENCED_BY], igc.values()):
        dst_ref = [dst_ep[ep_idx.ROW], dst_ep[ep_idx.INDEX]]
        for ref in dst_ep[ep_idx.REFERENCED_BY]:
            src_ep = igc[hash_ref(ref, SRC_EP)]
            src_refs = src_ep[ep_idx.REFERENCED_BY]
            if dst_ref not in src_refs:
                src_refs.append(dst_ref)


def _insert(igc_gcg, tgc_gcg, above_row):  # noqa: C901
    """Insert igc into the internal graph above row above_row.

    See https://docs.google.com/spreadsheets/d/1YQjrM91e5x30VUIRzipNYX3W7yiFlg6fy9wKbMTx1iY/edit?usp=sharing
    for more details (ask permission to view).
    graph is not modified.

    Args
    ----
    igc (gc_graph): Internal gc_graph format gc_graph to insert.
    tgc (gc_graph): Internal gc_graph format gc_graph to insert into.
    above_row (str): 'A', 'B' or 'O'

    Returns
    -------
    (gc_graph, gc_graph): rgc, fgc
    """
    tgc = tgc_gcg.graph
    igc = igc_gcg.graph
    rgc = _copy_clean_row(tgc, 'IC')
    fgc = {}
    if not tgc_gcg.has_a():
        if _LOG_DEBUG:
            _logger.debug("Case 1: No row A or B")
        rgc.update(_insert_as(igc, 'A'))
        rgc.update(_copy_row(tgc, 'O'))
    elif not tgc_gcg.has_b():
        if above_row == 'A':
            if _LOG_DEBUG:
                _logger.debug("Case 2: No row B and insert above A")
            rgc.update(_insert_as(igc, 'A'))
            rgc.update(_move_row(tgc, 'A', None, 'B', None))
            rgc.update(_redirect_refs(_copy_row(tgc, 'O'), 'O', DST_EP, 'A', 'B'))
        else:
            if _LOG_DEBUG:
                _logger.debug("Case 3: No row B and insert below A")
            rgc.update(_copy_row(tgc, 'A'))
            rgc.update(_insert_as(igc, 'B'))
            rgc.update(_copy_row(tgc, 'O'))
    else:
        if above_row == 'A':
            if _LOG_DEBUG:
                _logger.debug("Case 4: Has rows A & B and insert above A")
            fgc.update(_copy_clean_row(tgc, 'IC'))
            fgc.update(_insert_as(igc, 'A'))
            fgc.update(_move_row(tgc, 'A', None, 'B', None))
            fgc.update(_direct_connect(fgc, 'B', 'O'))
            fgc.update(_append_connect(fgc, 'A', 'O'))
            rgc.update(_direct_connect(rgc, 'I', 'A'))
            rgc.update(_move_row(fgc, 'O', None, 'A', SRC_EP, True))
            rgc.update(_copy_row(tgc, 'B'))
            rgc.update(_copy_row(tgc, 'O'))
        elif above_row == 'B':
            if _LOG_DEBUG:
                _logger.debug("Case 5: Has rows A & B and insert above B")
            fgc.update(_copy_clean_row(tgc, 'IC'))
            fgc.update(_copy_row(tgc, 'A'))
            fgc.update(_insert_as(igc, 'B'))
            fgc.update(_direct_connect(fgc, 'A', 'O'))
            fgc.update(_append_connect(fgc, 'B', 'O'))
            rgc.update(_direct_connect(rgc, 'I', 'A'))
            rgc.update(_move_row(fgc, 'O', None, 'A', SRC_EP, True))
            rgc.update(_copy_row(tgc, 'B'))
            rgc.update(_copy_row(tgc, 'O'))
        else:
            if _LOG_DEBUG:
                _logger.debug("Case 6: Has rows A & B and insert above O")
            rgc.update(_direct_connect(rgc, 'I', 'A'))
            rgc.update(_move_row(tgc, 'O', DST_EP, 'A', SRC_EP, True))
            rgc.update(_insert_as(igc, 'B'))
            rgc.update(_direct_connect(rgc, 'A', 'O'))

    # Case 1 is special because rgc is invalid by definition. In this case a
    # gc_graph normalization is forced to try and avoid the inevitable steady
    # state exception.
    if _LOG_DEBUG:
        _logger.debug("Pre-completed rgc:\n{}".format(pformat(rgc)))
    if not tgc_gcg.has_a():
        rgc_graph = gc_graph()
        rgc_graph.inject_graph(rgc)
        rgc_graph.normalize()
    else:
        _complete_references(rgc)
        rgc_graph = gc_graph()
        rgc_graph.inject_graph(rgc)
    if _LOG_DEBUG:
        _logger.debug("Completed rgc:\n{}".format(pformat(rgc)))

    if fgc:
        if _LOG_DEBUG:
            _logger.debug("Pre-completed fgc:\n{}".format(pformat(fgc)))
        _complete_references(fgc)
        fgc_graph = gc_graph()
        fgc_graph.inject_graph(fgc)
        if _LOG_DEBUG:
            _logger.debug("Completed fgc:\n{}".format(pformat(fgc)))
    else:
        fgc_graph = 0
    return rgc_graph, fgc_graph


def gc_insert(gms, target_gc, insert_gc=None, above_row=None):  # noqa: C901
    """Insert insert_gc into target_gc above row 'above_row'.

    If insert_gc is None then the target_gc is assessed for stability. If it
    is stable then [target_gc] will be returned otherwise a steady state exception
    is 'thrown'.

    NOTE: If a steady state exception occurs for which a candidate cannot
    be found in the GMS this function returns None.

    A work stack is used to avoid recursion.

    Insertion work is pushed onto the work_stack.
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

        If tgc has row F no insertion occurs.

    Args
    ----
    gp (gene_pool or genomic_library): A source of genetic material.
    target_gc (eGC): eGC to insert insert_gc into.
    insert_gc (eGC): eGC to insert into target_gc.
    above_row (string): One of 'A', 'B' or 'O'.

    Returns
    -------
    (rgc, {ref: fgc}): List of fGC's. First FGC is the insert_gc followed by fgc's
    created to stabilise rgc.
    """
    if target_gc['igraph'].has_f():
        return []
    if above_row is None:
        above_row = 'ABO'[randint(0, 2)]

    rgc_graph = target_gc['igraph']
    rgc = {
        'graph': rgc_graph.app_graph,
        'igraph': rgc_graph,
        'ref': random_reference(),
        'gca_ref': target_gc['gca_ref'],
        'gcb_ref': target_gc['gcb_ref']
    }

    # If there is no gc_insert return the target if it is stable or
    # throw a steady state exception.
    if insert_gc is None:
        if rgc_graph.is_stable():
            if _LOG_DEBUG: _logger.debug('Target GC is stable & nothing to insert.')
            return [target_gc]
        if _LOG_DEBUG: _logger.debug('Target GC is unstable & nothing to insert.')
        work_stack = [steady_state_exception(gms, rgc)]
    else:
        if _LOG_DEBUG: _logger.debug('Inserting into Target GC.')
        work_stack = [(rgc, insert_gc, above_row)]

    fgc_dict = {}
    new_tgc = None
    while work_stack and work_stack[0] is not None:
        if _LOG_DEBUG:
            _logger.debug("Work stack depth: {}".format(len(work_stack)))
        fgc = {}
        rgc = {}
        target_gc, insert_gc, above_row = work_stack.pop(0)
        if _LOG_DEBUG:
            _logger.debug("Work: Target={}, Insert={}, Above Row={}".format(target_gc['ref'], insert_gc['ref'], above_row))
        # TODO: Get rid of None (make it None)

        # Insert into the graph
        tgc_graph = target_gc['igraph'] if 'igraph' in target_gc else gc_graph(target_gc['graph'])
        igc_graph = insert_gc['igraph'] if 'igraph' in insert_gc else gc_graph(insert_gc['graph'])
        rgc_graph, fgc_graph = _insert(igc_graph, tgc_graph, above_row)
        if fgc_graph:
            fgc_steady = fgc_graph.normalize()
            if _LOG_DEBUG:
                _logger.debug("Normalized fgc:\n{}".format(pformat(fgc_graph)))
            fgc['graph'] = fgc_graph.app_graph
            fgc['igraph'] = fgc_graph
        rgc_steady = rgc_graph.normalize()
        if _LOG_DEBUG:
            _logger.debug("Normalized rgc:\n{}".format(pformat(rgc_graph)))
        rgc['graph'] = rgc_graph.app_graph
        rgc['igraph'] = rgc_graph

        # Insert into the GC
        # The insert_gc is always referenced in the tree of the final rgc
        fgc_dict[insert_gc['ref']] = insert_gc
        if not tgc_graph.has_a():  # Case 1
            if _LOG_DEBUG:
                _logger.debug("Case 1")
            rgc['gca_ref'] = insert_gc['ref']
            rgc['gcb_ref'] = None
        elif not tgc_graph.has_b():
            if above_row == 'A':  # Case 2
                if _LOG_DEBUG:
                    _logger.debug("Case 2")
                rgc['gca_ref'] = insert_gc['ref']
                if target_gc['gca_ref'] is not None:
                    rgc['gcb_ref'] = target_gc['gca_ref']
                else:
                    rgc['gcb_ref'] = target_gc['ref']
                    fgc_dict[target_gc['ref']] = target_gc
            else:  # Case 3
                if _LOG_DEBUG:
                    _logger.debug("Case 3")
                if target_gc['gca_ref'] is not None:
                    rgc['gca_ref'] = target_gc['gca_ref']
                else:
                    rgc['gca_ref'] = target_gc['ref']
                    fgc_dict[target_gc['ref']] = target_gc
                rgc['gcb_ref'] = insert_gc['ref']
        else:  # Has row A & row B
            if above_row == 'A':  # Case 4
                if _LOG_DEBUG:
                    _logger.debug("Case 4")
                fgc['gca_ref'] = insert_gc['ref']
                fgc['gcb_ref'] = target_gc['gca_ref']
                fgc['ref'] = random_reference()
                rgc['gca_ref'] = fgc['ref']
                rgc['gcb_ref'] = target_gc['gcb_ref']
            elif above_row == 'B':  # Case 5
                if _LOG_DEBUG:
                    _logger.debug("Case 5")
                fgc['gca_ref'] = target_gc['gca_ref']
                fgc['gcb_ref'] = insert_gc['ref']
                fgc['ref'] = random_reference()
                rgc['gca_ref'] = fgc['ref']
                rgc['gcb_ref'] = target_gc['gcb_ref']
            else:  # Case 6
                if _LOG_DEBUG:
                    _logger.debug("Case 6")
                rgc['gca_ref'] = target_gc['ref']
                rgc['gcb_ref'] = insert_gc['ref']
                fgc_dict[target_gc['ref']] = target_gc
        rgc['ref'] = target_gc['ref']

        # Check we have valid graphs
        # FGC is added to the work stack ahead of RGC
        # which ensures the new TGC is correctly set
        # (i.e. does not end up being the RGC of an unstable FGC)
        if fgc_graph:
            if not fgc_steady:
                if _LOG_DEBUG: _logger.debug(f"FGC ref {fgc['ref']} is unstable.")
                work_stack.insert(0, steady_state_exception(gms, fgc))
            else:
                if _LOG_DEBUG: _logger.debug(f"FGC ref {fgc['ref']} added to fgc_dict.")
                fgc_dict[fgc['ref']] = mGC(gc=fgc)

        if not rgc_steady:
            if _LOG_DEBUG: _logger.debug(f"RGC ref {rgc['ref']} is unstable.")
            work_stack.insert(0, steady_state_exception(gms, rgc))
        else:
            if new_tgc is None:
                if _LOG_DEBUG: _logger.debug(f"Resultant GC defined {rgc['ref']}:\n {pformat(rgc)}")
                new_tgc = rgc
            else:
                if _LOG_DEBUG: _logger.debug(f"RGC ref {rgc['ref']} added to fgc_dict.")
                fgc_dict[rgc['ref']] = mGC(gc=rgc)

        if _LOG_DEBUG:
            _logger.debug("fgc_dict: {}".format(list(fgc_dict.keys())))

    if _LOG_DEBUG:
        _logger.debug("fgc_dict details:\n{}".format(pformat(fgc_dict)))
    return (None, None) if work_stack else (new_tgc, fgc_dict)


# stablise() is more understanable context than gc_insert() in some use cases.
stablise = gc_insert


def gc_stack(gms, top_gc, bottom_gc):
    """Stack two GC's.

    top_gc is stacked on top of bottom_gc to create a new gc.
    See gc_graph.stack() for the definition of stacking.

    If the new gc is invalid it is repaired with a steady state exception.

    Args
    ----
    NOTE: If a steady state exception occurs for which a candidate cannot
    be found in the GMS this function returns None.

    Args
    ----
    gms (genetic_material_store): A source of genetic material.
    top_gc (mGC): GC to stack on top
    bottom_gc (mGC): GC to put on the bottom.

    Returns
    -------
    rgc (mGC): Resultant minimal GC with a valid graph or None
    """
    _gc = {'gca_ref': top_gc['ref'], 'gcb_ref': bottom_gc['ref']}
    igraph = top_gc.stack(bottom_gc['igraph'])
    if igraph is None: return None
    rgc = mGC(_gc, igraph=igraph)
    if not rgc['igraph'].normalize():
        rgc, _ = gc_insert(gms, *steady_state_exception(rgc, gms))
    return rgc

def gc_remove(gms, tgc, abpo=None):
    """Remove row A, B, P or O from tgc['graph'] to create rgc.

    If the row removed is A or B then GCA or GCB is set to None.
    If the row is P then both rows F and P are removed.
    If the row is O then both rows F and O are removed, P is copied to O
    and P is removed.
    The subsequent invalid graph is normalised and used to create rgc.
    Removing a row is likely to result in an invalid graph which is
    repaired using recursive steady state exceptions.

    NOTE: If a steady state exception occurs for which a candidate cannot
    be found in the GMS this function returns None.

    Args
    ----
    gms (gene_pool or genomic_library): A source of genetic material.
    tgc (xgc): Target xGC to modify.
    abpo (str): Either 'A', 'B', 'P' or 'O'. If None a row is removed at random.

    Returns
    -------
    rgc (mGC): Resultant minimal GC with a valid graph or None
    """
    rgc = {}
    if abpo is None:
        abpo = choice('ABPO')
    rgc_graph = deepcopy(tgc['igraph'])
    rgc_graph.remove_rows(abpo)
    if abpo == 'A':
        rgc['gca_ref'] = tgc['gcb_ref']
        rgc['gcb_ref'] = None
        rgc_graph.move_refs('B', 'A')
    elif abpo == 'B':
        rgc['gca_ref'] = tgc['gca_ref']
        rgc['gcb_ref'] = None
    elif abpo == 'P':
        rgc_graph.remove_rows('F')
    elif abpo == 'O':
        rgc_graph.remove_rows('F')
        rgc_graph.graph.update(_move_row(rgc_graph.graph, 'P', None, 'O', None))
        rgc_graph.move_refs('P', 'O')
    rgc['ref'] = random_reference()
    if not rgc_graph.normalize():
        rgc, _ = gc_insert(gms, *steady_state_exception(rgc, gms))
        return rgc
    return mGC(rgc)


def proximity_select(gms, xputs):
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
    #   a) Cache queries (but this means missing out on new options)
    #   b) Batch queries (but this is architecturally tricky)
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
    return None


def steady_state_exception(gms, fgc):
    """Define what GC must be inserted to complete or partially complete the fgc graph.

    fgc is analysed to determine what end point destinations are unconnected and the highest row
    on which one or more of the unconnected destination endpoints resides.

    Candidates to react (plug the gap) are found from the gms (genetic material source) which is
    either a gene_pool or genomic_library instance. NB: The source of genetic material should be
    the most local scope.

    If no candidates are found in the GMS then None is returned.

    Args
    ----
    gms (object): Either a gene_pool or genomic_library instance.
    fgc (fGC): fGC with incomplete graph.

    Returns
    -------
    (fGC, fGC, str): (target_gc, insert_gc, 'A', 'B' or 'O') or None
    """
    if _LOG_DEBUG: _logger.debug(f"Steady state exception thrown for GC ref {fgc['ref']}.")
    fgc_graph = fgc['igraph']

    # Find unconnected destination endpoints. Determine highest row & endpoint types.
    dst_list = list(filter(fgc_graph.unreferenced_filter(fgc_graph.dst_filter()), fgc_graph.graph.values()))
    above_row = 'Z'
    outputs = []
    for ep in dst_list:
        if ep[ep_idx.ROW] < above_row:
            above_row = ep[ep_idx.ROW]
        outputs.append(ep[ep_idx.TYPE])

    # Find viable source types above the highest row.
    filter_func = fgc_graph.rows_filter(fgc_graph.src_rows[above_row], fgc_graph.src_filter())
    inputs = list((ep[ep_idx.TYPE] for ep in filter(filter_func, fgc_graph.graph.values())))

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


def pGC_fitness(gp, pgc, delta_fitness):
    """Update the fitness of the pGC and the pGC's that created it.

    pgc is modified.
    -1.0 <= delta_fitness <= 1.0

    pGC's are checked for evolution

    Args
    ----
    gp (gene_pool): The gene_pool that contains pGC and its creators.
    pgc (pGC): A physical GC.
    delta_fitness (float): The change in fitness of the GC pGC mutated.
    """
    depth = 1
    _pGC_fitness(pgc, delta_fitness, depth)
    delta_fitness = pgc['delta_fitness'][depth]
    evolved = evolve_physical(gp, pgc, depth)
    pgc = gp.pool.get(pgc['pgc'], None)
    while evolved and pgc is not None:
        depth += 1
        xGC_evolvability(pgc, delta_fitness, depth)
        _pGC_fitness(pgc, delta_fitness, depth)
        delta_fitness = pgc['delta_fitness'][depth]
        evolved = evolve_physical(gp, pgc, depth)
        pgc = gp.pool.get(pgc['pgc'], None)


def _pGC_fitness(pgc, delta_fitness, depth):
    """Update the fitness of the pGC.

    pgc is modified.
    -1.0 <= delta_fitness <= 1.0

    Args
    ----
    pgc (pGC): pGC to update.
    delta_fitness (float): The change in fitness of the GC pGC mutated.
    depth (int): The layer in the environment pgc is at.
    """
    pgc['if'] = 0.0 if delta_fitness < 0 else delta_fitness
    old_count = pgc['f_count'][depth]
    pgc['delta_fitness'][depth] += pgc['if'] - pgc['fitness'][depth]
    pgc['f_count'][depth] += 1
    pgc['fitness'][depth] = (old_count * pgc['fitness'][depth] + pgc['if']) / pgc['f_count'][depth]


def xGC_evolvability(xgc, delta_fitness, depth):
    """Update the evolvability of any type of GC.

    xgc is modified.
    -1.0 <= delta_fitness <= 1.0

    Args
    ----
    xgc (pGC): xGC to update.
    delta_fitness (float): Difference in fitness between this GC & its offspring.
    depth (int): The layer in the environment pgc is at.
    """
    increase = 0.0 if delta_fitness < 0 else delta_fitness
    old_count = xgc['e_count'][depth]
    xgc['e_count'][depth] += 1
    xgc['evolvability'][depth] = (old_count * xgc['evolvability'][depth] + increase) / xgc['e_count'][depth]


def evolve_physical(gp, pgc, depth):
    """Evolve the pgc as needed.

    pgc is checked to see if it meets evolution criteria. If it does
    it is evolved & the gene pool updated with it offspring.

    PGC's evolve when they have been 'used' M_CONSTANT times (which must
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
    if not (pgc['f_count'][depth] & M_MASK):
        pgc['delta_fitness'][depth] = 0.0
        gp.layer_evolutions[depth] += 1
        ppgc = select_pGC(gp, pgc, depth)
        offspring = ppgc.exec((pgc,))
        xGC_inherit(offspring, pgc, ppgc)
        return True
    return False


def cull_physical(pgc_iter, layer):
    """Remove one pGC in layer .

    Note that the pGC is not deleted if it is present in other layers,
    its f_count for the layer are just set to 0.

    Args
    ----
    pgc_iter iter(dict): pGC's to consider. Dict must have:
        'f_count'
        'fitness'
    depth (int): The layer in the environment to select a victim from.
    """
    func = lambda x: x['f_count'][layer]
    filtered_pool = tuple(filter(func, pgc_iter))
    weights = array([1.0 - i['fitness'][layer] for i in filtered_pool], float32)
    weights /= sum(weights)
    victim = choice(filtered_pool, None, False, weights)
    victim['f_count'][layer] = 0 #### Layer valid!


def create_layer(gp):
    """Create a new physical GC layer.

    Creation means the depth of the environment increases by adding a
    new top layer of pGC's.

    NOTE: Only pGC's can exist in layers > 0
    """
    for pgc in filter(lambda x: x['properties'] & PHYSICAL_PROPERTY, gp.pool.values()):
        for col in LAYER_COLUMNS:
            pgc[col].append(pgc[col][-1])
        for col, value in LAYER_COLUMNS_RESET.items():
            pgc[col][-1] = value if pgc[col][-2] else 0
    gp.num_layers += 1
    gp.layer_evolutions.append(0)


def select_pGC(gp, xgc, depth):
    """Select a pgc to evolve xgc.

    A pgc is found from layer depth + 1.
    If layer depth + 1 does not exist then it is created.

    Args
    ----
    gp (gene_pool): The gene pool containing xgc.
    xgc (xGC): The GC to find a pgc for.
    depth (int): The layer in the environment xgc is at.

    Returns
    -------
    pgc (pGC): A pGC to evolve xgc.
    """
    if depth == gp.num_layers:
        create_layer(gp)

    # TODO: More sophisticated selection algorithm
    # This one picks a pGC randomly weighted by fitness.
    next_layer = depth + 1
    func = lambda x: x['f_count'][next_layer]

    # OPTIMIZATION: Weights & filtered_pool could be cached for a depth?
    filtered_pool = tuple(filter(func, gp.pool.values()))
    weights = array([1.0 - i['if'][next_layer] for i in filtered_pool], float32)
    weights /= sum(weights)
    return choice(filtered_pool, None, False, weights)


def xGC_inherit(child, parent, pgc):
    """Inherit some properties from parent to child.

    child is modified.
    parent is modified.

    Args
    ----
    child (xGC): Offspring of parent and pgc.
    parent (xGC): Parent of child.
    pgc (pGC): pGC that operated on parent to product child.
    """
    child['fitness'] = copy(parent['fitness'])
    for f_count in parent['f_count']:
        child['f_count'] = 1 if f_count else 0
    child['evolvability'] = copy(parent['evolvability'])
    for f_count in parent['e_count']:
        child['e_count'] = 1 if f_count else 0

    child['population'] = parent['population']
    child['ancestor_a_ref'] = parent['ref']
    child['pgc_ref'] = pgc['ref']
    child['generation'] = parent['generation'] + 1

    parent['offspring_count'] += 1
