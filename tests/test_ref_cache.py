"""Tests for the reference cache implementation"""
import pytest
from egp_physics.ref_cache import update_cache, select_from_cache, cache_metrics, invalidate_cache, create_cache, get_cache, invalidate_all_caches
from egp_physics.ref_cache import HISTORY_MSB, CACHE_ENTRY_LIMIT
from numpy.random import randint, shuffle
from numpy.random import choice as np_choice
from numpy.typing import NDArray
from numpy import int64, arange, iinfo, logical_and, unpackbits
from collections import Counter
from tqdm import trange
from typing import Literal


INT64_MIN = iinfo(int64).min
INT64_MAX = iinfo(int64).max


# Cache select statistics
# See test_select_from_cache()
# Mininum and maximum count for each reference observed in 65536 runs.
STATS: dict[int64 | int, dict[Literal['min', 'max'], int]] = {
    48: {'min': 64817, 'max': 66352},
    47: {'min': 31988, 'max': 33396},
    46: {'min': 15839, 'max': 16983},
    45: {'min': 7816, 'max': 8542},
    44: {'min': 3809, 'max': 4356},
    43: {'min': 1858, 'max': 2243},
    42: {'min': 893, 'max': 1162},
    41: {'min': 410, 'max': 612},
    40: {'min': 189, 'max': 322},
    39: {'min': 83, 'max': 177},
    38: {'min': 33, 'max': 100},
    37: {'min': 11, 'max': 59},
    36: {'min': 3, 'max': 33},
    35: {'min': 1, 'max': 23},
    34: {'min': 1, 'max': 14},
    33: {'min': 1, 'max': 10},
    32: {'min': 1, 'max': 7},
    31: {'min': 1, 'max': 5},
    30: {'min': 1, 'max': 5},
    29: {'min': 1, 'max': 4},
    28: {'min': 1, 'max': 3},
    27: {'min': 1, 'max': 3},
    26: {'min': 1, 'max': 2},
    25: {'min': 1, 'max': 2},
    24: {'min': 1, 'max': 2},
    23: {'min': 1, 'max': 1},
    22: {'min': 1, 'max': 1},
    21: {'min': 1, 'max': 1},
    20: {'min': 1, 'max': 1},
    19: {'min': 1, 'max': 1},
    18: {'min': 1, 'max': 1},
    17: {'min': 1, 'max': 1},
    16: {'min': 1, 'max': 1},
    15: {'min': 1, 'max': 1},
    14: {'min': 1, 'max': 1},
    13: {'min': 1, 'max': 1},
    12: {'min': 1, 'max': 1}
}


def test_create_cache() -> None:
    """Simple instanciation"""
    invalidate_all_caches()
    create_cache(1)
    update_cache(1, int64(1), True)
    history, refs, hits, misses = cache_metrics(1)
    assert history[0] == 1 << HISTORY_MSB
    assert refs[0] == 1
    assert hits == 1
    assert misses == 0


def test_multiple_updates() -> None:
    """Check that the cache is updated correctly when the same reference is updated multiple times."""
    invalidate_all_caches()
    create_cache(1)
    update_cache(1, int64(1), True)
    update_cache(1, int64(1), False)
    history, refs, hits, misses = cache_metrics(1)
    assert history[0] == 1 << (HISTORY_MSB - 1)
    assert refs[0] == 1
    assert hits == 1
    assert misses == 1


@pytest.mark.parametrize("case", range(128))
def test_new_updates(case: int) -> None:
    "Generate a random number of hits & misses and check the cache statistics."
    if not case:
        invalidate_all_caches()
    create_cache(case)
    size: int = randint(CACHE_ENTRY_LIMIT)
    for ref in arange(1, size + 1, dtype=int64):
        update_cache(case, ref, bool(randint(2)))
    history, refs, hits, misses = cache_metrics(case)
    assert (refs != 0).sum() == size
    assert hits + misses == size
    assert hits == (history > 0).sum()
    assert misses == logical_and(history == 0, refs != 0).sum()


def test_histories() -> None:
    """Check that the cache is updated correctly when the same reference is updated multiple times."""
    invalidate_all_caches()
    create_cache(1)
    for _ in range(HISTORY_MSB + 1):
        refs: NDArray[int64] = arange(1, CACHE_ENTRY_LIMIT + 1, dtype=int64)
        shuffle(refs)
        for ref in refs:
            update_cache(1, ref, bool(randint(2)))
    history, refs, hits, misses = cache_metrics(1)
    assert unpackbits(history.view('uint8')).sum() == hits
    assert hits + misses == CACHE_ENTRY_LIMIT * (HISTORY_MSB + 1)


def test_history_overflow() -> None:
    """Check that the cache is updated correctly when the same reference is updated multiple times."""
    invalidate_all_caches()
    create_cache(1)
    loop_data = {
        0: [0, 0],
        1: [0, 0]
    }
    for loop, data in loop_data.items():
        for _ in range(HISTORY_MSB + 1):
            refs: NDArray[int64] = arange(1, CACHE_ENTRY_LIMIT + 1, dtype=int64)
            shuffle(refs)
            for ref in refs:
                update_cache(1, ref, bool(randint(2)))
        history, _, hits, __ = cache_metrics(1)
        data[0] = unpackbits(history.view('uint8')).sum()
        data[1] = CACHE_ENTRY_LIMIT * (HISTORY_MSB + 1) - data[0]
    _, __, hits, misses = cache_metrics(1)
    assert loop_data[1][0] + loop_data[0][0] == hits
    assert loop_data[1][1] + loop_data[0][1] == misses


def test_select_from_cache() -> None:
    """Check that selection from cache statistics are as expected.
    
    Create a cache with HISTORY_MSB + 1 entries with the reference the same as the index.
    Each entry having the number of hits equal to its reference number follolwed by misses
    to equal HISTORY_MSB + 1 updates.

    If HISTORY_MSB + 1 is 8, the cache will look like this:
    ref: 0 history: 0b00000000
    ref: 1 history: 0b00000001
    ref: 2 history: 0b00000011
    ref: 3 history: 0b00000111
    ref: 4 history: 0b00001111
    ref: 5 history: 0b00011111
    ref: 6 history: 0b00111111
    ref: 7 history: 0b01111111
    ref: 8 history: 0b11111111
    """
    invalidate_all_caches()
    create_cache(1)
    for ref in arange(0, HISTORY_MSB + 1, dtype=int64):
        for update in arange(0, HISTORY_MSB + 1, dtype=int64):
            update_cache(1, ref, update < ref)
    for flag in range(2):
        for ref, count in Counter((select_from_cache(1) for _ in range(2**17))).items():
            if count < STATS.get(ref, {'min': 0})['min'] or count > STATS.get(ref, {'max': 0})['max']:  # type: ignore
                assert not flag, '1 in 2**16 chance of being out of bounds occurred twice in a row...more than suspicious.'
                break
        return