"""Tests for the reference cache implementation"""
import pytest
from numpy import array, int32, zeros, float64, histogram, arange
from numpy.typing import NDArray
from egp_physics.binary_history_probability_table import binary_history_probability_table as bhpt


# Correct default state weights for a 56 bit consideration length
_DEFAULT_STATE_WEIGHTS: NDArray[float64] = array([
    1.0, 1.5874010519681994, 2.5198420997897464, 4.0, 6.3496042078727974, 10.079368399158986, 16.0, 25.398416831491197,
    40.317473596635935, 64.0, 101.59366732596479, 161.26989438654374, 256.0, 406.3746693038589, 645.0795775461753, 1024.0,
    1625.4986772154357, 2580.318310184701, 4096.0, 6501.994708861743, 10321.273240738805, 16384.0, 26007.97883544697,
    41285.09296295522, 65536.0, 104031.915341788, 165140.37185182067, 262144.0, 416127.661367152, 660561.4874072827,
    1048576.0, 1664510.645468608, 2642245.949629131, 4194304.0, 6658042.581874432, 10568983.798516523, 16777216.0,
    26632170.32749773, 42275935.19406609, 67108864.0, 106528681.30999091, 169103740.77626437, 268435456.0, 426114725.23996365,
    676414963.1050575, 1073741824.0, 1704458900.9598546, 2705659852.42023, 4294967296.0, 6817835603.839402, 10822639409.680946,
    17179869184.0, 27271342415.35761, 43290557638.723785, 68719476736.0, 109085369661.43044
], dtype=float64)

# size, h_length, c_length, mwsp, defer, exception
_TEST_INIT_EXCEPTIONS: tuple[tuple[int, int, int, bool, bool, bool], ...] = (
    (0, 1, 1, True, False, True),
    (1, 0, 1, True, False, True),
    (1, 1, 0, True, False, True),
    (1, 1, 1, True, False, False),
    (1, 1, 1, False, False, False),
    (128, 80, 32, False, False, False),
    (128, 80, 32, True, False, False),
    (128, 80, 80, True, False, False)
)

@pytest.mark.parametrize('size, h_length, c_length, mwsp, defer, exception', _TEST_INIT_EXCEPTIONS)
def test_bhpt_init(size: int, h_length: int, c_length: int, mwsp: bool, defer: bool, exception:bool) -> None:
    """Tests the BHPT initialization."""
    if exception:
        with pytest.raises(ValueError):
            bhpt(size, h_length, c_length, mwsp, defer)
    else:
        bhpt(size, h_length, c_length, mwsp, defer)


def test_bhpt_repr() -> None:
    """Tests the BHPT string representation."""
    assert repr(bhpt(128, 64, 64, True)) == "binary_history_probability_table(size=128, h_length=64,c_length=64, mwsp=True, defer=False)"
    assert repr(bhpt(12, 6, 6, False, True)) == "binary_history_probability_table(size=12, h_length=6,c_length=6, mwsp=False, defer=True)"

@pytest.mark.parametrize('defer', (False, True))
def test_bhpt_stats(defer) -> None:
    """Tests the BHPT default statistics using an 8 state buffer.
    
    1. Create a BHPT with 256 entries and a history length of 8. The consideration length is also 8.
    2. Each entry history is set to a unique 8 bit state.
    3. The get() method probability distribution of indices should tend toward the default state weights.
    4. Defer should make no difference.
    """
    test_bhpt = bhpt(256, 8, 8, False, defer)
    for i in range(256):
        for s in f'{i:08b}':
            test_bhpt[i] = s == '1'

    # Check the probabilities
    frequencies: NDArray[int32] = zeros(256, dtype=int32)
    frequencies, _ = histogram(test_bhpt.get_many(2**23), bins=arange(257))
    weights: NDArray[float64] = frequencies / frequencies.max()
    print(weights)
    print(test_bhpt._default_state_weights())
    expected_frequencies = (test_bhpt._h_table * test_bhpt._default_state_weights()).sum(axis=1)
    expected_weights: NDArray[float64] = expected_frequencies / expected_frequencies.max()
    print(weights / expected_weights)
    # plot


def test_bhpt_default_state_weights() -> None:
    """Tests the BHPT default state weights."""
    test_bhbpt = bhpt(1, 64, 56, False)
    assert (test_bhbpt._default_state_weights() == _DEFAULT_STATE_WEIGHTS).all()