"""The operations that can be performed on a GC."""
from __future__ import annotations
from copy import copy, deepcopy
from logging import DEBUG, Logger, NullHandler, getLogger
from random import randint
from typing import LiteralString, Callable, cast, TYPE_CHECKING

from egp_execution.execution import create_callable
from egp_types.xGC import xGC, pGC, gGC
from egp_types.ep_type import interface_definition, vtype
from egp_types.gc_graph import gc_graph
from egp_types.gc_type_tools import M_MASK, NUM_PGC_LAYERS
from egp_types.reference import ref_str
from egp_types.aGC import aGC
from egp_types.dGC import dGC
from egp_types.egp_typing import Row, VALID_ROW_SOURCES

from .egp_typing import Work, NewGCDef, InsertRow
from .pgc_bhpt import pgc_bhpt
from .insertion import _recursive_insert_gc


# Circular import & runtime import avoidance
if TYPE_CHECKING:
    from egp_stores.gene_pool import gene_pool


_logger: Logger = getLogger(__name__)
_logger.addHandler(NullHandler())
_LOG_DEBUG: bool = _logger.isEnabledFor(DEBUG)


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


# PGC Constants
RANDOM_PGC_SIGNATURE: bytes = b"\x00" * 32
_PGC_PARENTAL_PROTECTION_FACTOR: float = 0.75
_POPULATION_PARENTAL_PROTECTION_FACTOR: float = 0.75


# Binary History Probability Table setup for pGC selection
# Each pGC layer has its own BHPT.
_PGC_BHPT: list[pgc_bhpt] = []


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


def stablize(gms: gene_pool, tgc: aGC) -> NewGCDef:
    """If tgc is not stable force a steady state exceptiion"""
    if tgc["gc_graph"].is_stable():
        if _LOG_DEBUG:
            assert tgc["gc_graph"].validate()
            _logger.debug("Target GC is stable.")
        return (tgc, {})
    if _LOG_DEBUG:
        _logger.debug("Target GC is unstable.")
    return _recursive_insert_gc(gms, [steady_state_exception(gms, tgc)])


def clone(gc_to_clone: aGC, ref: Callable[[], int], copy_graph: bool = True) -> dGC:
    """Create a minimal clone of gc.

    Args
    ----
    gc_to_clone: GC to clone

    Returns
    -------
    Minimal clone of gc as a dict.
    """
    return {
        "ref": ref(),
        "ancestor_a_ref": gc_to_clone["ref"],
        "ancestor_b_ref": 0,
        # If gc is a codon then it does not have a GCA
        "gca_ref": gc_to_clone["gca_ref"]
        if gc_to_clone["gca_ref"]
        else gc_to_clone["ref"],
        "gcb_ref": gc_to_clone["gcb_ref"],
        # More efficient than reconstructing
        "gc_graph": deepcopy(gc_to_clone["gc_graph"]) if copy_graph else gc_graph(),
    }


def pgc_epilogue(gms: gene_pool, agc: aGC) -> xGC:
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
        _logger.debug(f"pGC epilogue with agc = {agc}")
    new_gc_definition: NewGCDef = stablize(gms, agc)
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

    Args
    ----
    gms: A source of genetic material
    fgc (fGC): fGC with incomplete graph.

    Returns
    -------
    A workstack: (target_gc, insert_gc, 'A', 'B' or 'O')
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
            _above_row = ep.row
        outputs.append(ep.typ)

    # Find viable source types above the highest row.
    inputs: list[int] = [
        ep.typ
        for ep in fgc_graph.i_graph.src_rows_filter(
            VALID_ROW_SOURCES[fgc["gc_graph"].has_f][cast(Row, above_row)]
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
    delta_fitness = pgc["pgc_delta_fitness"][depth]
    evolved: bool = evolve_physical(gms, pgc, depth)
    evolutions = int(evolved)
    pgc_creator: xGC | None = gms.pool.get(pgc["pgc_ref"])
    while evolved and pgc_creator is not None:
        depth += 1
        pgc_creator = cast(pGC, pgc_creator)
        _pGC_evolvability(pgc_creator, delta_fitness, depth)
        _pGC_fitness(pgc_creator, delta_fitness, depth)
        delta_fitness = pgc_creator["pgc_delta_fitness"][depth]
        evolved = evolve_physical(gms, pgc_creator, depth)
        evolutions += evolved
        pgc = pgc_creator
        pgc_creator = gms.pool.get(pgc_creator["pgc_ref"], None)
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
    old_count: int = pgc["pgc_f_count"][depth]
    pgc["pgc_f_count"][depth] += 1
    pgc["pgc_fitness"][depth] = (
        old_count * pgc["pgc_fitness"][depth] + (delta_fitness / 2 + 0.5)
    ) / pgc["pgc_f_count"][depth]
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
    old_count: int = pgc["pgc_e_count"][depth]
    pgc["pgc_e_count"][depth] += 1
    pgc["pgc_evolvability"][depth] = (
        old_count * pgc["pgc_evolvability"][depth] + increase
    ) / pgc["pgc_e_count"][depth]


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
    old_count: int = ggc["e_count"]
    ggc["e_count"] += 1
    ggc["evolvability"] = (old_count * ggc["evolvability"] + increase) / ggc["e_count"]


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
    if not (pgc["pgc_f_count"][depth] & M_MASK):
        pgc["pgc_delta_fitness"][depth] = 0.0
        ppgc: pGC = select_pGC(gms, pgc, depth + 1)
        wrapped_ppgc_callable = create_callable(ppgc, gms.pool)
        result = wrapped_ppgc_callable((pgc,))
        if result is None:
            # pGC went pop - should not happen very often
            _logger.warning(
                f"ppGC {ref_str(pgc['ref'])} threw an exception when called."
            )
            return False
        offspring = result[0]
        if _LOG_DEBUG:
            assert isinstance(offspring, xGC)
        pGC_inherit(offspring, pgc, ppgc)
        return True
    return False


def select_pGC(gms: gene_pool, xgc: xGC, depth: int = 0) -> pGC:
    """Select a pgc to evolve xgc.

    The Gene Pool Cache is refreshed on a periodic basis by the expiry of the sub-process
    and a write back to the Gene Pool. At this point at least some pGC's that are not direct decendants
    of this Gene Pool Cache will be bought in.

    A pGC is selected to act on each xGC.

    pGC selection makes use of temporal locality & fitness. The pGC's that have been used most recently and
    are most successful have the highest probability of being selected.

    Args
    ----
    gp: The gene pool containing pgc.
    xgc_refs: Iterable of GC references to find a pGC for.
    depth: The xGC layer for which to find a pGC. Used as the ref_cache() UID.

    Returns
    -------
    The pGCs to evolve xgcs.
    """
    if not _PGC_BHPT:
        _PGC_BHPT.extend(
            (
                pgc_bhpt(2**13, gms.pool, 0),
                pgc_bhpt(2**13, gms.pool, 1),
                pgc_bhpt(2**13, gms.pool, 2),
                pgc_bhpt(2**13, gms.pool, 3),
                pgc_bhpt(2**10, gms.pool, 4),
                pgc_bhpt(2**10, gms.pool, 5),
                pgc_bhpt(2**10, gms.pool, 6),
                pgc_bhpt(2**10, gms.pool, 7),
                pgc_bhpt(2**7, gms.pool, 8),
                pgc_bhpt(2**7, gms.pool, 9),
                pgc_bhpt(2**7, gms.pool, 10),
                pgc_bhpt(2**7, gms.pool, 11),
                pgc_bhpt(2**4, gms.pool, 12),
                pgc_bhpt(2**4, gms.pool, 13),
                pgc_bhpt(2**4, gms.pool, 14),
                pgc_bhpt(2**4, gms.pool, 15),
            )
        )
    return cast(pGC, gms.pool[_PGC_BHPT[depth].get()])


def pGC_inherit(child, parent, pgc) -> None:
    """Inherit pGC only properties from parent to child.

    child is modified.
    parent is modified.

    Parental protection works differently for a PGC than a population individual because a
    pGC child must survive before being characterized where as a population child is
    characterized immediately.

    Args
    ----
    child (xGC): Offspring of parent and pgc.
    parent (xGC): Parent of child.
    pgc (pGC): pGC that operated on parent to product child.
    """
    # TODO: A better data structure would be quicker
    child["pgc_fitness"] = [
        f * _PGC_PARENTAL_PROTECTION_FACTOR for f in parent["pgc_fitness"]
    ]
    child["pgc_f_count"] = [2] * NUM_PGC_LAYERS
    child["pgc_evolvability"] = [
        f * _PGC_PARENTAL_PROTECTION_FACTOR for f in parent["pgc_evolvability"]
    ]
    child["pgc_e_count"] = [2] * NUM_PGC_LAYERS

    child["_pgc_fitness"] = [0.0] * NUM_PGC_LAYERS
    child["_pgc_f_count"] = [0] * NUM_PGC_LAYERS
    child["_pgc_evolvability"] = [0.0] * NUM_PGC_LAYERS
    child["_pgc_e_count"] = [0] * NUM_PGC_LAYERS

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
        if not all((field in child for field in ("fitness", "survivability"))):
            raise ValueError("Child GC has not been characterized.")

    # There is no way of characterising first
    if parent["e_count"] == 1:
        child["evolvability"] = 1.0
        child["e_count"] = 1
    else:
        child["evolvability"] = parent["evolvability"]
        child["e_count"] = max((2, parent["e_count"] >> 1))

    inherited_survivability = (
        parent["survivability"] * _POPULATION_PARENTAL_PROTECTION_FACTOR
    )
    inherited_fitness = parent["fitness"] * _POPULATION_PARENTAL_PROTECTION_FACTOR
    child["survivability"] = max((child["survivability"], inherited_survivability))
    child["fitness"] = max((child["fitness"], inherited_fitness))
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

    child["population_uid"] = parent["population_uid"]
    child["ancestor_a_ref"] = parent["ref"]
    child["pgc_ref"] = pgc["ref"]
    child["generation"] = parent["generation"] + 1
    child["effective_pgc_refs"] = copy(parent["effective_pgc_refs"])
    child["effective_pgc_fitness"] = copy(parent["effective_pgc_fitness"])

    parent["offspring_count"] += 1
