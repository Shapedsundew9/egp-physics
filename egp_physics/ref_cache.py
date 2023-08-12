"""Reference cache functions.

GC reference caches provide a temporal locality selection mechanism considering performance history.

A reference cache starts empty.
Each time the cache is updated with a reference a hit or miss boolean is passed. A hit (True) means the referenced
GC was useful (e.g. increased fitness). A miss (False) means it was not.
The boolean is recorded as a 1 (True) or 0 (False) in a history bit-string for that reference where the MSB is the
most recent update.

If the cache has reached capacity then the least performant GC is replaced (or randomly from references with all
zeros history).

Selecting a reference from the cache is done with a weighted random selection using the hit/miss history for each
reference as a weight. The history bit string is MSB aligned so behaves like a 2's complement number e.g. 

   Update in order of occurance = Bit string history = Weight
a) Miss, Miss, Miss, Hit        = 1000               = 8
b) Miss, Miss, Hit, Miss        = 0100               = 4
c) Hit                          = 1xxx               = 8
d) Hit, Hit, Hit                = 111x               = 14
e) Hit, Miss, Miss, Miss        = 0001               = 1
f) Miss, Miss, Miss, Miss       = 0000               = 0

where x = no history = 0. NOTE: that a limited length recent history is maintained (HISTORY_DEPTH)

If all references in the cache have 0 weight then the selection function returns None.
"""
from functools import lru_cache
from logging import Logger, NullHandler, getLogger, DEBUG
from numpy import int64, int32, float32, zeros, where, argmin
from numpy.typing import NDArray
from typing import Generator
from random import choice
from numpy.random import choice as np_choice


_logger: Logger = getLogger(__name__)
_logger.addHandler(NullHandler())
_LOG_DEBUG: bool = _logger.isEnabledFor(DEBUG)


# Cache store UID's - causes an exception if too many are used
MAX_NUM_CACHE_STORES = 128
HISTORY_DEPTH = 49  # FYI: 64 bit float mantissa = 52 bits.
HISTORY_MSB = HISTORY_DEPTH - 1
CACHE_ENTRY_LIMIT: int = 2 ** 13
assert CACHE_ENTRY_LIMIT * (2 ** HISTORY_DEPTH - 1) < 2 ** 63 - 1, 'Risk of probability normalisation denominator overflowing.'


class ref_cache():

    def __init__(self, uid: int) -> None:
        """Creates a new reference cache store."""
        self._uid: int = uid
        self.fitness: NDArray[float32] = zeros((CACHE_ENTRY_LIMIT,), dtype=float32) 
        self._history: NDArray[int64] = zeros((CACHE_ENTRY_LIMIT,), dtype=int64)
        self._refs: NDArray[int64] = zeros((CACHE_ENTRY_LIMIT,), dtype=int64)
        self._idx_dict: dict[int64, int] = {}
        self.stats: dict[bool, NDArray[int32]] = {True: zeros((CACHE_ENTRY_LIMIT,), dtype=int32), False: zeros((CACHE_ENTRY_LIMIT,), dtype=int32)}

        # Seed the cache with a single entry representing a pass-through to GP.
        self._set_pass_through()
        _logger.info(f'Reference cache store UID: {uid} created.')

    def __setitem__(self, ref: int64, hit: bool) -> None:
        """Update cache with a hit (or miss) for ref."""
        # TODO: power of 2 drop off may be to extreme, Consider implementing
        # a scheme where when history is updated the relative weight (in a new array)
        # is updated by an easily defined function.
        if (idx := self._idx_dict.get(ref)) is not None:
            self.stats[hit][idx] += 1
            # If there has been a hit we can calculate the actual fitness.
            if self.stats[True][idx]:
                self.fitness[idx] = self.stats[True][idx] / (self.stats[True][idx] + self.stats[False][idx])
        elif len(self._idx_dict) < CACHE_ENTRY_LIMIT:
            new_idx: int = len(self._idx_dict)
            self._idx_dict[ref] = new_idx 
            self._refs[new_idx] = ref
            self.stats[hit][new_idx] = 1
            # In the event of a miss on the first update the fitness is set to the inverse of the cache size.
            self.fitness[new_idx] = 1.0 if hit else 1.0 / CACHE_ENTRY_LIMIT
        else:
            victims: NDArray[int64] = self._refs[where(self._history == 0)]
            victim: int64 = choice(victims) if len(victims) else self._refs[argmin(self._history)]
            victim_idx: int = self._idx_dict[victim]
            del self._idx_dict[victim]
            self._idx_dict[ref] = victim_idx
            self.stats[hit][victim_idx] = 1
            # In the event of a miss on the first update the fitness is set to the inverse of the cache size.
            self.fitness[victim_idx] = 1.0 if hit else 1.0 / CACHE_ENTRY_LIMIT

    def _set_pass_through(self) -> None:
        """Set the pass through reference (0) history."""
        self[int64(0)] = True

    def select(self) -> int64:
        """History weighted random selection from the cache"""
        history_sum: int64 = self._history.sum()
        return np_choice(self._refs, p = self._history / history_sum)

    def invalidate(self) -> None:
        """Invalidate all cache entries."""
        # Must not destroy the containers though.
        self._history.fill(0)
        self._refs.fill(0)
        self._idx_dict.clear()
        self.stats[True] = 0
        self.stats[False] = 0
        self._set_pass_through()
        _logger.info(f'Reference cache store UID: {self._uid} invalidated.')


class ref_cache_store():

    def __init__(self) -> None:
        self._caches: dict[int, ref_cache] = {}

    def __getitem__(self, uid: int) -> ref_cache:
        """Creates a new reference cache store if it does not exist."""
        if uid not in self._caches:
            assert len(self._caches) < MAX_NUM_CACHE_STORES, "Reference cache limit has been reached."
            self._caches[uid] = ref_cache(uid)
        return self._caches[uid]

    def invalidate_all_caches(self) -> None:
        """Invalidate all caches."""
        for cache in self._caches.values():
            cache.invalidate()

    def __delitem__(self, uid: int) -> None:
        """Delete a cache."""
        del self._caches[uid]


# Global cache store
global_ref_cache_store: ref_cache_store = ref_cache_store()