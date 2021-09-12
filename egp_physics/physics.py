"""The operation that can be performed on a GC dictionary."""
from copy import copy, deepcopy
from logging import DEBUG, NullHandler, getLogger
from pprint import pformat
from random import randint

from .ep_type import vtype
from .gc_graph import (DST_EP, SRC_EP, ep_idx, gc_graph, hash_ep, hash_ref,
                       ref_idx)
from .gc_type import eGC, interface_definition, mGC
from .utils.reference import random_reference

_logger = getLogger(__name__)
_logger.addHandler(NullHandler())
_LOGIT = _logger.getEffectiveLevel() == DEBUG


# Steady state exception filters.
_MATCH_TYPE_0_SQL = ('WHERE "input_types"={input_types} AND "inputs"={inputs} '
                     'AND "output_types"={output_types} AND "outputs"={outputs} '
                     'AND {exclude_column} NOT IN {exclusions} ORDER BY RANDOM() LIMIT 1')
_MATCH_TYPE_1_SQL = ('WHERE "input_types"={input_types} '
                     'AND "output_types"={output_types} AND "outputs"={outputs} '
                     'AND {exclude_column} NOT IN {exclusions} ORDER BY RANDOM() LIMIT 1')
_MATCH_TYPE_2_SQL = ('WHERE "input_types"={input_types} AND "inputs"={inputs} '
                     'AND "output_types"={output_types} '
                     'AND {exclude_column} NOT IN {exclusions} ORDER BY RANDOM() LIMIT 1')
_MATCH_TYPE_3_SQL = ('WHERE "input_types"={input_types} AND "output_types"={output_types}'
                     'AND {exclude_column} NOT IN {exclusions} ORDER BY RANDOM() LIMIT 1')
_MATCH_TYPE_4_SQL = ('WHERE "input_types"<@{input_types} AND "output_types"={output_types}'
                     'AND {exclude_column} NOT IN {exclusions} ORDER BY RANDOM() LIMIT 1')
_MATCH_TYPE_5_SQL = ('WHERE "input_types"<@{input_types} AND "output_types"@>{output_types}'
                     'AND {exclude_column} NOT IN {exclusions} ORDER BY RANDOM() LIMIT 1')
_MATCH_TYPE_6_SQL = ('WHERE "input_types"<@{input_types} AND "output_types"&&{output_types}'
                     'AND {exclude_column} NOT IN {exclusions} ORDER BY RANDOM() LIMIT 1')
_MATCH_TYPE_7_SQL = ('WHERE "input_types"&&{input_types} AND "output_types"&&{output_types}'
                     'AND {exclude_column} NOT IN {exclusions} ORDER BY RANDOM() LIMIT 1')
_MATCH_TYPE_8_SQL = ('WHERE "output_types"&&{output_types} '
                     'AND {exclude_column} NOT IN {exclusions} ORDER BY RANDOM() LIMIT 1')
_MATCH_TYPE_9_SQL = ('WHERE "input_types"&&{input_types} '
                     'AND {exclude_column} NOT IN {exclusions} ORDER BY RANDOM() LIMIT 1')
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
_INITIAL_GC_SQL = ('WHERE "input_types"={input_types} AND "inputs"={inputs} '
                   'AND "output_types"={output_types} AND "outputs"={outputs} '
                   'AND {exclude_column} NOT IN {exclusions} ORDER BY RANDOM() LIMIT {limit}')


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
    if _LOGIT:
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
        if _LOGIT:
            _logger.debug("Case 1: No row A or B")
        rgc.update(_insert_as(igc, 'A'))
        rgc.update(_copy_row(tgc, 'O'))
    elif not tgc_gcg.has_b():
        if above_row == 'A':
            if _LOGIT:
                _logger.debug("Case 2: No row B and insert above A")
            rgc.update(_insert_as(igc, 'A'))
            rgc.update(_move_row(tgc, 'A', None, 'B', None))
            rgc.update(_redirect_refs(_copy_row(tgc, 'O'), 'O', DST_EP, 'A', 'B'))
        else:
            if _LOGIT:
                _logger.debug("Case 3: No row B and insert below A")
            rgc.update(_copy_row(tgc, 'A'))
            rgc.update(_insert_as(igc, 'B'))
            rgc.update(_copy_row(tgc, 'O'))
    else:
        if above_row == 'A':
            if _LOGIT:
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
            if _LOGIT:
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
            if _LOGIT:
                _logger.debug("Case 6: Has rows A & B and insert above O")
            rgc.update(_direct_connect(rgc, 'I', 'A'))
            rgc.update(_move_row(tgc, 'O', DST_EP, 'A', SRC_EP, True))
            rgc.update(_insert_as(igc, 'B'))
            rgc.update(_direct_connect(rgc, 'A', 'O'))

    # Case 1 is special because rgc is invalid by definition. In this case a
    # gc_graph normalization is forced to try and avoid the inevitable steady
    # state exception.
    if _LOGIT:
        _logger.debug("Pre-completed rgc:\n{}".format(pformat(rgc)))
    if not tgc_gcg.has_a():
        rgc_graph = gc_graph()
        rgc_graph.inject_graph(rgc)
        rgc_graph.normalize()
    else:
        _complete_references(rgc)
        rgc_graph = gc_graph()
        rgc_graph.inject_graph(rgc)
    if _LOGIT:
        _logger.debug("Completed rgc:\n{}".format(pformat(rgc)))

    if fgc:
        if _LOGIT:
            _logger.debug("Pre-completed fgc:\n{}".format(pformat(fgc)))
        _complete_references(fgc)
        fgc_graph = gc_graph()
        fgc_graph.inject_graph(fgc)
        if _LOGIT:
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
    ([fgc]): List of fGC's. Element 0 may replace target_gc.
    Subsequent fGC's are children of element 0.
    """
    if target_gc['_graph'].has_f():
        return []
    if above_row is None:
        above_row = 'ABO'[randint(0, 2)]

    # If there is no gc_insert return the target if it is stable or
    # throw a steady state exception.
    if insert_gc is None:
        rgc_graph = target_gc['_graph']
        if rgc_graph.is_stable():
            return [target_gc]
        rgc = {'graph': rgc_graph.app_graph, '_graph': rgc_graph, '_ref': random_reference()}
        work_stack = [steady_state_exception(gms, rgc)]
    else:
        work_stack = [(target_gc, insert_gc, above_row)]

    fgc_list = []
    while work_stack and work_stack[0] is not None:
        if _LOGIT:
            _logger.debug("Work stack depth: {}".format(len(work_stack)))
        fgc = {}
        rgc = {}
        target_gc, insert_gc, above_row = work_stack.pop(0)
        if _LOGIT:
            _logger.debug("Work: Target={}, Insert={}, Above Row={}".format(target_gc['_ref'], insert_gc['_ref'], above_row))
        # TODO: Get rid of None (make it None)

        # Insert into the graph
        tgc_graph = target_gc['_graph'] if '_graph' in target_gc else gc_graph(target_gc['graph'])
        igc_graph = insert_gc['_graph'] if '_graph' in insert_gc else gc_graph(insert_gc['graph'])
        rgc_graph, fgc_graph = _insert(igc_graph, tgc_graph, above_row)
        if fgc_graph:
            fgc_steady = fgc_graph.normalize()
            if _LOGIT:
                _logger.debug("Normalized fgc:\n{}".format(pformat(fgc_graph)))
            fgc['graph'] = fgc_graph.app_graph
            fgc['_graph'] = fgc_graph
        rgc_steady = rgc_graph.normalize()
        if _LOGIT:
            _logger.debug("Normalized rgc:\n{}".format(pformat(rgc_graph)))
        rgc['graph'] = rgc_graph.app_graph
        rgc['_graph'] = rgc_graph

        # Insert into the GC
        if not tgc_graph.has_a():  # Case 1
            if _LOGIT:
                _logger.debug("Case 1")
            rgc['_gca'] = insert_gc['_ref']
            rgc['_gcb'] = None
        elif not tgc_graph.has_b():
            if above_row == 'A':  # Case 2
                if _LOGIT:
                    _logger.debug("Case 2")
                rgc['_gca'] = insert_gc['_ref']
                rgc['_gcb'] = target_gc['_gca'] if target_gc['_gca'] is not None else target_gc['_ref']
            else:  # Case 3
                if _LOGIT:
                    _logger.debug("Case 3")
                rgc['_gca'] = target_gc['_gca'] if target_gc['_gca'] is not None else target_gc['_ref']
                rgc['_gcb'] = insert_gc['_ref']
        else:  # Has row A & row B
            if above_row == 'A':  # Case 4
                if _LOGIT:
                    _logger.debug("Case 4")
                fgc['_gca'] = insert_gc['_ref']
                fgc['_gcb'] = target_gc['_gca']
                fgc['_ref'] = random_reference()
                rgc['_gca'] = fgc['_ref']
                rgc['_gcb'] = target_gc['_gcb']
            elif above_row == 'B':  # Case 5
                if _LOGIT:
                    _logger.debug("Case 5")
                fgc['_gca'] = target_gc['_gca']
                fgc['_gcb'] = insert_gc['_ref']
                fgc['_ref'] = random_reference()
                rgc['_gca'] = fgc['_ref']
                rgc['_gcb'] = target_gc['_gcb']
            else:  # Case 6
                if _LOGIT:
                    _logger.debug("Case 6")
                rgc['_gca'] = target_gc['_ref']
                rgc['_gcb'] = insert_gc['_ref']
        rgc['_ref'] = random_reference()

        # Check we have valid graphs
        if fgc_graph:
            if not fgc_steady:
                work_stack.insert(0, steady_state_exception(gms, fgc))
            else:
                fgc_list.append(mGC(gc=fgc))

        if not rgc_steady:
            work_stack.insert(0, steady_state_exception(gms, rgc))
        else:
            fgc_list.insert(0, mGC(gc=rgc))

    if _LOGIT:
        _logger.debug("fgc_list: {}".format(pformat(fgc_list)))
    return None if work_stack else fgc_list


def gc_remove(gms, tgc, abpo):
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
    gp (gene_pool or genomic_library): A source of genetic material.
    tgc (xgc): Target xGC to modify.
    abpo (str): Either 'A', 'B', 'P' or 'O'.

    Returns
    -------
    rgc (pgc): Resultant partial GC with a valid graph or None
    """
    rgc = {}
    rgc_graph = deepcopy(tgc['_graph'])
    rgc_graph.remove_rows(abpo)
    if abpo == 'A':
        rgc['_gca'] = tgc['_gcb']
        rgc['_gcb'] = None
        rgc_graph.move_refs('B', 'A')
    elif abpo == 'B':
        rgc['_gca'] = tgc['_gca']
        rgc['_gcb'] = None
    elif abpo == 'P':
        rgc_graph.remove_rows('F')
    elif abpo == 'O':
        rgc_graph.remove_rows('F')
        rgc_graph.graph.update(_move_row(rgc_graph.graph, 'P', None, 'O', None))
        rgc_graph.move_refs('P', 'O')
    rgc['_ref'] = random_reference()
    if not rgc_graph.normalize():
        rgc = gc_insert(gms, *steady_state_exception(rgc, gms))[0]
    return rgc


def mutate_graph(tgc, n=1):
    """Disconnect random connections in tgc['graph'].

    This is a random physical operation (like a cosmic ray strike) that
    disrupts the graph of the GC by deleting connections and then restoring a steady
    state.

    n is the number of connections to delete. 0 is effectively an no-op.
    If n > number of connections in the graph all connections in the graph will be
    mutated.

    NOTE: Connections may be reformed the same.

    Args
    ----
        tgc (xgc): GC's graph to mutate.
        n (int): Number of connections to delete.

    Returns
    -------
        (xgc): A new GC in a steady state.
    """
    graph = tgc['_graph']
    graph.random_remove_connection(n)
    assert graph.normalize()


def proximity_select(xputs, gms):
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
    match_type = randint(0, _NUM_MATCH_TYPES - 1)
    agc = gms.select(_MATCH_TYPES_SQL[match_type], literals=xputs, container='pkdict')
    while not agc and match_type < _NUM_MATCH_TYPES:
        if _LOGIT:
            _logger.debug(f'Proximity selection match_type {match_type} found no candidates.')
        match_type += 1
        agc = gms.select(_MATCH_TYPES_SQL[match_type], literals=xputs, container='pkdict')
    return None if not agc else agc[0]


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
    fgc_graph = fgc['_graph']

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
        'exclusions': tuple()
    }
    _, xputs['input_types'], xputs['inputs'] = interface_definition(inputs, vtype.EP_TYPE_INT)
    _, xputs['output_types'], xputs['outputs'] = interface_definition(outputs, vtype.EP_TYPE_INT)

    # Find a gc based on the criteria
    insert_gc = gms.proximity_select(xputs)
    if insert_gc is None:
        _logger.warning("Steady state exception failed to find a candidate in the GMS.")
        return None

    return (fgc, eGC(insert_gc), above_row)


def initial_gc(inputs, outputs, gms, exclusions=tuple(), num=1, vt=vtype.OBJECT):
    """Create num valid GC's with the specified inputs & outputs.

    GC's matching the criteria in the GMS will be given preference over
    creating new GC's. If there are more GC's in the GMS that meet the
    criteria than num then the returned GC's will be randomly selected.
    If there are less then valid GC's with the correct inputs and outputs
    will be created randomly.

    If the input & output specification is such that no valid GC can be found
    or created None is returned.

    Args
    ----
    inputs (iterable(object or type)): Input objects of vt type.
    outputs (iteratble(object or type)): Output objects of vt type.
    gms (object): Genetic Material Source: A gene_pool or genomic_library object.
    num (int): The number of GC's to create
    vt (vtype): See vtype definition.

    Returns
    -------
    (list(eGC)) or None
    """
    xputs = {
        'exclude_column': 'signature',
        'exclusions': exclusions,
        'limit': num
    }
    _, xputs['input_types'], xputs['inputs'] = interface_definition(inputs, vt)
    _, xputs['output_types'], xputs['outputs'] = interface_definition(outputs, vt)

    gcs = [eGC(gc=gc) for gc in gms.select(_INITIAL_GC_SQL, literals=xputs, container='dict')]
    gcs.extend((eGC(inputs=inputs, outputs=outputs, vt=vt) for _ in range(num - len(gcs))))
    return gcs
