"""Tests for the reference cache implementation"""
import pytest
from numpy import array, int32, zeros, float64
from numpy.typing import NDArray
from egp_physics.binary_history_probability_table import binary_history_probability_table as bhpt


# size, h_length, c_length, mwsp, defer, exception
_TEST_INIT_EXCEPTIONS: tuple[tuple[int, int, int, int, bool, bool], ...] = (
    (0, 1, 1, 0, False, True),
    (1, 0, 1, 0, False, True),
    (1, 1, 0, 0, False, True),
    (1, 1, 2, 0, False, True),
    (1, 1, 1, 1, False, True),
    (1, 1, 1, 1, False, True),
    (128, 80, 32, 40, False, True),
    (128, 80, 32, -8, False, False),
    (128, 80, 32, 16, False, False)
)

@pytest.mark.parametrize('size, h_length, c_length, mwsp, defer, exception', _TEST_INIT_EXCEPTIONS)
def test_bhpt_init(size: int, h_length: int, c_length: int, mwsp: int, defer: bool, exception:bool) -> None:
    """Tests the BHPT initialization."""
    if exception:
        with pytest.raises(ValueError):
            bhpt(size, h_length, c_length, mwsp, defer)
    else:
        bhpt(size, h_length, c_length, mwsp, defer)


def test_bhpt_repr() -> None:
    """Tests the BHPT string representation."""
    assert repr(bhpt(128, 64, 64, 63)) == "binary_history_probability_table(size=128, h_length=64,c_length=64, mwsp=63, defer=False)"
    assert repr(bhpt(128, 64, 64, 63, True)) == "binary_history_probability_table(size=128, h_length=64,c_length=64, mwsp=63, defer=True)"

@pytest.mark.parametrize('defer', (False, True))
def test_bhpt_stats(defer) -> None:
    """Tests the BHPT default statistics using an 8 state buffer.
    
    1. Create a BHPT with 256 entries and a history length of 8. The consideration length is also 8.
    2. Each entry history is set to a unique 8 bit state.
    3. The get() method probability distribution of indices should tend toward the default state weights.
    4. Defer should make no difference.
    """
    test_bhpt = bhpt(256, 8, 8, -1, defer)
    for i in range(256):
        for s in f'{i:08b}':
            test_bhpt[i] = s == '1'

    # Check the probabilities
    frequencies: NDArray[int32] = zeros(256, dtype=int32)
    frequencies[test_bhpt.get_many(100000)] += 1
    print(frequencies)
    weights: NDArray[float64] = frequencies / frequencies.max()
