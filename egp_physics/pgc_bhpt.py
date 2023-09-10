"""pGC Binary History Probability Tables.

This module contains the pgc_bhpts class and supporting functions.
Each layer of pGCs has a BHPT with a history length of 64 and a consideration length of 64.
They are initialised with a single static (never updated) entry representing the Gene Pool Cache (GPC).
When the GPC entry is selected a pGC is pulled randomly weighted by fitness from the GPC
The new pGC is entered into the BHPT with its history initialised by its fitness
then updated by its success in evolving a target GC.
When a pGC is selected from the BHPT it is used to evolve a target GC and its success is used to update
its history in the BHPT.
"""

from numpy import int8, where, double, zeros, int64
from numpy.typing import NDArray
from egp_physics.binary_history_probability_table import (
    binary_history_probability_table,
    default_state_weights,
)
from egp_stores.gene_pool_cache import gene_pool_cache
from egp_physics.pgc import pGC


# Binary History Probability Table setup for pGC selection
# Each pGC layer has its own BHPT. The next deeper layer tends quickly
# toward being called * times less than the previous so having a large
# number of entries is not useful.
_BHPT_H_LENGTH: int = 64
_BHPT_C_LENGTH: int = 64
_BHPT_MWSP: bool = True
_BHPT_DEFER: bool = False
_PGC_BHPT_GP_ENTRY_HISTORY: NDArray[int8] = zeros(_BHPT_C_LENGTH, dtype=int8)
_PGC_BHPT_GP_ENTRY_HISTORY[_BHPT_C_LENGTH // 4 :] = 1
_PGC_BHPT_WEIGHTS: NDArray[double] = default_state_weights(_BHPT_C_LENGTH)
_PGC_BHPT_WEIGHTS_SUM: double = _PGC_BHPT_WEIGHTS.sum()
_PGC_FITNESS_MAPPING_TO_HISTORY: NDArray[int8] = zeros(
    (128, _BHPT_C_LENGTH), dtype=int8
)


# Create a look up table mapping pGC fitness to BHPT history (consideration length)
# A granularity of 128 is accurate enough as the fitness from the GPC
# is for all time rather than the local environment.
for i in range(128):
    fitness: double = _PGC_BHPT_WEIGHTS_SUM * double(i) / 127.0
    for c in range(_BHPT_C_LENGTH):
        if fitness > _PGC_BHPT_WEIGHTS[c]:
            fitness -= _PGC_BHPT_WEIGHTS[c]
            _PGC_FITNESS_MAPPING_TO_HISTORY[i][c] = int8(1)


class pgc_bhpt(binary_history_probability_table):
    def __init__(self, size: int, gpc: gene_pool_cache, depth: int) -> None:
        """Creates a pGC BHPT.

        Args
        ----
        size: The size of the BHPT in entries (I).
        gpc: The Gene Pool Cache.
        depth: The depth of the pGC in the layers.
        """
        super().__init__(size, _BHPT_H_LENGTH, _BHPT_C_LENGTH, _BHPT_MWSP, _BHPT_DEFER)
        self._gpc: gene_pool_cache = gpc
        self._depth: int = depth
        self._refs: NDArray[int64] = zeros(size, dtype=int64)
        super()[0] = _PGC_BHPT_GP_ENTRY_HISTORY

    def __contains__(self, ref: int64) -> bool:
        """Returns True if the pGC is in the BHPT.

        Args
        ----
        ref: The GC reference to check.

        Returns
        -------
        True if the GC reference is in the BHPT.
        """
        return ref in self._refs

    def __setitem__(self, ref: int64, state: bool) -> None:
        """Sets the pGC history in the BHPT."""
        # This is about the same speed as maintaining an LRU cache unless self._refs size >> 8k
        index: int = where(self._refs == ref)[0][0]
        return super().__setitem__(index, state)

    def _add(self, pgc: pGC) -> int64:
        """Adds a pGC to the BHPT.

        Args
        ----
        pgc: The pGC to add.

        Returns
        -------
        The pgc reference.
        """
        idx: int = self.insert()
        ref: int64 = pgc["ref"]
        self._refs[idx] = ref
        super()[idx] = _PGC_FITNESS_MAPPING_TO_HISTORY[
            int(pgc["pgc_fitness"][self._depth] * 127)
        ]
        return ref

    def get(self) -> int64:
        """Returns a GC reference from the BHPT.

        Returns
        -------
        The GC reference.
        """
        idx: int = super().get()
        if idx:
            return self._refs[idx]
        return self._add(self._gpc.random_pgc())
