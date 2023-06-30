"""The operations that can be performed on a GC."""
from collections.abc import Iterable
from copy import copy
from logging import DEBUG, Logger, NullHandler, getLogger
from random import randint
from typing import Union

from egp_execution.execution import create_callable
from egp_stores.gene_pool_cache import gene_pool_cache
from egp_types.eGC import eGC
from egp_types.xGC import gGC, pGC, xGC
from egp_types.ep_type import interface_definition, vtype
from egp_types.gc_graph import gc_graph
from egp_types.gc_type_tools import M_MASK, NUM_PGC_LAYERS, is_pgc
from egp_types.reference import ref_str

from .insertion import _insert_gc

from numpy import array, float32, isfinite
from numpy.random import choice as weighted_choice

_logger: Logger = getLogger(__name__)
_logger.addHandler(NullHandler())
_LOG_DEBUG: bool = _logger.isEnabledFor(DEBUG)


# PGC Constants
RANDOM_PGC_SIGNATURE: bytes = b'\x00' * 32
_PGC_PARENTAL_PROTECTION_FACTOR: float = 0.75
_POPULATION_PARENTAL_PROTECTION_FACTOR: float = 0.75


def _pgc_epilogue(gms, xgc):
    """Stabilises xGC and updates the gms.

    pGC's may return a single gGC or None if a valid gGC could
    not be created. Since pGC's may be stacked each pGC must
    be able to handle a None input.

    fGC's created as a side effect of of pGC activity are
    added to the GP local cache.

    Args
    ----
    gms (gene_pool or genomic_library): A source of genetic material.
    xgc (xGC): Unstable GC

    Returns
    -------
    rgc (gGC): Resultant gGC or None
    """
    if _LOG_DEBUG:
        _logger.debug(f'PGC epilogue with xgc = {xgc}')
    if xgc is not None:
        rgc, fgcs = _insert_gc(gms, xgc)
        if rgc is not None:
            # TODO: Yuk - need to de-mush physics & GP. gGC is a GP concept not a GMS one
            # 7-May-2022: Hmmm! But GP is a GMS and should fallback to GL when looking for a GC
            # In fact gms in the parameters should be GP?
            ggcs = gGC((rgc, *fgcs.values()))
            return ggcs[0]
    return None


def create_SMS(gp, pgc, ggc):
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
    parent = gp.pool.get(ggc['ancestor_a_ref'])
    lineage = [ggc, parent]
    assert pgc['ref'] == ggc['pgc_ref'], 'pGC providied did not create ggc!'
    parent['effective_pgc_refs'].append(pgc['ref'])
    sms = pgc

    # Consecutive increases in fitness
    increase = 0.0
    while lineage[-1] is not None and lineage[-1]['fitness'] < lineage[-2]['fitness']:
        sms = gc_stack(gp, lineage[-1]['pgc_ref'], sms)
        increase += lineage[-2] - lineage[-1]
        parent['effective_pgc_refs'].append(sms['ref'])
        lineage.append(gp.pool.get(lineage[-1]['ancestor_a_ref']))
        if _LOG_DEBUG:
            assert is_pgc(sms), 'Super Mutation Sequence is not a pGC!'

    # Consecutive reductions in fitness
    while lineage[-1] is not None and lineage[-1]['fitness'] >= lineage[-2]['fitness']:

        # Even though this SMS may be net negative there is a possibility of mutation
        # to extract a net positive that is worth keeping.
        sms = gc_stack(gp, lineage[-1]['pgc_ref'], sms)

        # If the total increase in this SMS is still net >0.0 add it as an effective pGC
        increase += lineage[-2] - lineage[-1]
        if increase > 0.0:
            parent['effective_pgc_refs'].append(sms['ref'])
        lineage.append(gp.pool.get(lineage[-1]['ancestor_a_ref']))
        if _LOG_DEBUG:
            assert is_pgc(sms), 'Super Mutation Sequence is not a pGC!'

    terminal = lineage.pop()
    if terminal is None:
        assert lineage[-1]['generation'] == 0, 'SMS lineage terminated but oldest ancestor is not generation 0!'

    # If increase is still >0.0 then we have an SMS chain
    if increase > 0.0 and terminal is not None and terminal['sms_ref'] is not None:
        sms = gc_stack(gp, terminal['sms_ref'], sms)
        parent['effective_pgc_refs'].append(sms['ref'])
        if _LOG_DEBUG:
            assert is_pgc(sms), 'Super Mutation Sequence is not a pGC!'

    # Record the positive SMS
    parent['sms_ref'] = sms['ref']


# TODO: use pGC_fitness to call create_SMS when delta_fitness is positive
def pGC_fitness(gp: gene_pool_cache, pgc: gGC, ggc: gGC, delta_fitness: Union[float, None]) -> int:
    """Update the fitness of the pGC and the pGC's that created it.

    pgc is modified.
    -1.0 <= delta_fitness <= 1.0

    pGC's are checked for evolution

    Args
    ----
    gp: The gene_pool that contains pGC and its creators.
    pgc: A physical GC.
    ggc: The gGC created by pgc changing fitness from its parent by delta_fitness: May be None.
    delta_fitness: The change in fitness of the GC pGC mutated. May be None.

    Returns
    -------
    The number of pGC evolutions that occured as a result of the fitness update.
    """
    depth = 0
    _pGC_fitness(pgc, ggc, delta_fitness, depth)
    delta_fitness = pgc['pgc_delta_fitness'][depth]
    evolved = evolve_physical(gp, pgc, depth)
    evolutions = int(evolved)
    pgc_creator = gp.pool.get(pgc['pgc_ref'], None)
    while evolved and pgc_creator is not None:
        depth += 1
        _pGC_evolvability(pgc_creator, delta_fitness, depth)
        _pGC_fitness(pgc_creator, pgc, delta_fitness, depth)
        delta_fitness = pgc_creator['pgc_delta_fitness'][depth]
        evolved = evolve_physical(gp, pgc_creator, depth)
        evolutions += evolved
        pgc = pgc_creator
        pgc_creator = gp.pool.get(pgc_creator['pgc_ref'], None)
    return evolutions


def _pGC_fitness(pgc: gGC, xgc: gGC, delta_fitness: Union[float, None], depth: int) -> float:
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
    old_count = pgc['pgc_f_count'][depth]
    pgc['pgc_f_count'][depth] += 1

    if delta_fitness is None:
        delta_fitness = -1.0

    pgc['pgc_fitness'][depth] = (old_count * pgc['pgc_fitness'][depth] + (delta_fitness / 2 + 0.5)) / pgc['pgc_f_count'][depth]

    return delta_fitness


def _pGC_evolvability(pgc, delta_fitness, depth):
    """Update the evolvability of a PGC.

    pgc is modified.
    -1.0 <= delta_fitness <= 1.0

    Args
    ----
    pgc (pGC): PGC to update.
    delta_fitness (float): Difference in fitness between this GC & its offspring.
    depth (int): The layer in the environment pgc is at.
    """
    increase = 0.0 if delta_fitness < 0 else delta_fitness
    old_count = pgc['pgc_e_count'][depth]
    pgc['pgc_e_count'][depth] += 1
    pgc['pgc_evolvability'][depth] = (old_count * pgc['pgc_evolvability'][depth] + increase) / pgc['pgc_e_count'][depth]


def evolve_physical(gp, pgc, depth):
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
        gp.layer_evolutions[depth] += 1  # FIXME: Ugh!

        ppgc = select_pGC(gp, pgc, depth + 1)
        wrapped_ppgc_callable = create_callable(ppgc, gp.pool)
        result = wrapped_ppgc_callable((pgc,))
        if result is None:
            # pGC went pop - should not happen very often
            _logger.warning(f"ppGC {ref_str(pgc['ref'])} threw an exception when called.")
            offspring = None
        else:
            offspring = result[0]

        if offspring is not None:
            if _LOG_DEBUG:
                assert isinstance(offspring, _gGC)
            offspring['pgc_ref'] = ppgc['ref']
            pGC_inherit(offspring, gp[offspring['ancestor_a_ref']], gp[offspring['ancestor_b_ref']], ppgc)
        return True
    return False


def select_pGC(gp: gene_pool_cache, xgc_refs: Iterable[int], depth: int = 0) -> list[_gGC]:
    """Select a pgc to evolve xgc.

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
        xgc = gp[xgc_ref]

        # Intergenerational pGC selection
        next_pgc_ref = xgc['next_pgc_ref']
        if next_pgc_ref is not None:
            xgc['next_pgc_ref'] = None
            matched_pgcs.append(gp.pool[next_pgc_ref])
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
                    all_pgcs = tuple(gc for gc in gp.pool.values() if is_pgc(gc))
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
                matched_pgcs.append(gp[weighted_choice(xgc['effective_pgc_refs'], p=normalised_weights)])

    return matched_pgcs


def pGC_inherit(descendant: pGC, ancestor_a: xGC, ancestor_b: xGC, pgc: pGC):
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
    descendant['pgc_fitness'] = [f * _PGC_PARENTAL_PROTECTION_FACTOR for f in ancestor_a['pgc_fitness']]
    descendant['pgc_f_count'] = [2] * NUM_PGC_LAYERS
    descendant['pgc_evolvability'] = [f * _PGC_PARENTAL_PROTECTION_FACTOR for f in ancestor_a['pgc_evolvability']]
    descendant['pgc_e_count'] = [2] * NUM_PGC_LAYERS

    descendant['_pgc_fitness'] = [0.0] * NUM_PGC_LAYERS
    descendant['_pgc_f_count'] = [0] * NUM_PGC_LAYERS
    descendant['_pgc_evolvability'] = [0.0] * NUM_PGC_LAYERS
    descendant['_pgc_e_count'] = [0] * NUM_PGC_LAYERS

    xGC_inherit(descendant, ancestor_a, ancestor_b, pgc)


def population_GC_inherit(descedant: gGC, ancestor_a: xGC, ancestor_b: xGC, pgc: pGC):
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
        if not all((field in descedant for field in ('fitness', 'survivability'))):
            raise ValueError('Child GC has not been characterized.')

    # There is no way of characterising first
    if ancestor_a['e_count'] == 1:
        descedant['evolvability'] = 1.0
        descedant['e_count'] = 1
    else:
        descedant['evolvability'] = ancestor_a['evolvability']
        descedant['e_count'] = max((2, ancestor_a['e_count'] >> 1))

    inherited_survivability = ancestor_a['survivability'] * _POPULATION_PARENTAL_PROTECTION_FACTOR
    inherited_fitness = ancestor_a['fitness'] * _POPULATION_PARENTAL_PROTECTION_FACTOR
    descedant['survivability'] = max((descedant['survivability'], inherited_survivability))
    descedant['fitness'] = max((descedant['fitness'], inherited_fitness))
    xGC_inherit(descedant, ancestor_a, ancestor_b, pgc)


def xGC_inherit(descendant: xGC, ancestor_a: xGC, ancestor_b: xGC, pgc: pGC):
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

    descendant['population_uid'] = ancestor_a['population_uid']
    descendant['ancestor_a_ref'] = ancestor_a['ref']
    descendant['pgc_ref'] = pgc['ref']
    descendant['generation'] = ancestor_a['generation'] + 1
    descendant['effective_pgc_refs'] = copy(ancestor_a['effective_pgc_refs'])
    descendant['effective_pgc_fitness'] = copy(ancestor_a['effective_pgc_fitness'])

    ancestor_a['offspring_count'] += 1
