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
from numpy import int64, zeros, where, argmin
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
_CACHES: dict[int, tuple[NDArray[int64], NDArray[int64], dict[int64, int], dict[bool, int]]] = {}
assert CACHE_ENTRY_LIMIT * (2 ** HISTORY_DEPTH - 1) < 2 ** 63 - 1, 'Risk of probability normalisation denominator overflowing.'


def _ref_cache_store_generator() -> Generator[tuple[NDArray[int64], NDArray[int64], dict[int64, int], dict[bool, int]], None, None]:
    """Creates a generator of reference cache stores."""
    for _ in range(MAX_NUM_CACHE_STORES):
        yield (zeros((CACHE_ENTRY_LIMIT,), dtype=int64), zeros((CACHE_ENTRY_LIMIT,), dtype=int64), {}, {True: 0, False: 0})


_RCSG: Generator[tuple[NDArray[int64], NDArray[int64], dict[int64, int], dict[bool, int]], None, None] = _ref_cache_store_generator()
def get_cache(uid: int) -> tuple[NDArray[int64], NDArray[int64], dict[int64, int], dict[bool, int]]:
    """Creates a new reference cache store.
    
    The LRU cache wrapper ensures the same cache is returned for the same uid.

    Returns
    -------
    (Hit history, ref, ref to idx dict)
    """
    _logger.info(f'Reference cache store UID: {uid} created.')
    if uid not in _CACHES:
        assert len(_CACHES) < MAX_NUM_CACHE_STORES, "Reference cache limit has been reached."
        _CACHES[uid] = next(_RCSG)
    return _CACHES[uid]


# Alias
create_cache = get_cache


def update_cache(uid: int, ref: int64, hit: bool) -> None:
    """Update cache 'uid' with a hit (or miss) for ref."""
    history, refs, idx_dict, stats = get_cache(uid)
    if (idx := idx_dict.get(ref)) is not None:
        history[idx] = (int(hit) << HISTORY_MSB) + (history[idx] >> 1)
    elif len(idx_dict) < CACHE_ENTRY_LIMIT:
        new_idx: int = len(idx_dict)
        idx_dict[ref] = new_idx 
        history[new_idx] = int(hit) << HISTORY_MSB
        refs[new_idx] = ref
    else:
        victims: NDArray[int64] = refs[where(history == 0)]
        victim: int64 = choice(victims) if len(victims) else refs[argmin(history)]
        victim_idx: int = idx_dict[victim]
        del idx_dict[victim]
        idx_dict[ref] = victim_idx
        history[victim_idx] = int(hit) << HISTORY_MSB
    stats[hit] += 1

def select_from_cache(uid: int) -> int64 | None:
    """History weighted random selection from the cache 'uid'"""
    history, refs, _, __ = get_cache(uid)
    history_sum: int64 = history.sum()
    return np_choice(refs, p = history / history_sum) if history_sum else None


def cache_metrics(uid: int) -> tuple[NDArray[int64], NDArray[int64], int, int]:
    """Get some statistics about the cache.
    
    Returns
    -------
    (Current histories, References, Total number of hits ever, Total number of misses ever)
    """
    history, refs, _, stats = get_cache(uid)
    return history, refs, stats[True], stats[False]


def invalidate_cache(uid: int) -> None:
    """Invalidate all cache entries."""
    # Must not destroy the containers though.
    history, refs, idx_dict, stats = get_cache(uid)
    history.fill(0)
    refs.fill(0)
    idx_dict.clear()
    stats[True] = 0
    stats[False] = 0


def invalidate_all_caches() -> None:
    """Invalidate all caches."""
    for uid in range(MAX_NUM_CACHE_STORES):
        invalidate_cache(uid)
