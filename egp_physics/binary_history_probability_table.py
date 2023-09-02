"""Binary History Probability Table, BHPT."""

from logging import Logger, NullHandler, getLogger, DEBUG
from numpy import int8, int32, float64, zeros, argmin, array, log2
from numpy.typing import NDArray
from numpy.random import choice as np_choice
from typing import Iterable


_logger: Logger = getLogger(__name__)
_logger.addHandler(NullHandler())
_LOG_DEBUG: bool = _logger.isEnabledFor(DEBUG)


def default_state_weights(length: int) -> NDArray[float64]:
    """Returns the default state weights for a length of length."""
    return array([2**i for i in range(length)], dtype=float64)


class binary_history_probability_table():

    def __init__(self, size: int = 128, h_length: int = 64, c_length: int = 64, mwsp: bool = False, defer: bool = False) -> None:
        """Creates a Binary History Probability Table, BHPT.

        See [Binary History Probability Table](../docs/binary_history_probability_table.md) for more details.

        Args
        ----
        size: The size of the BHPT in entries (I).
        h_length: The length of the history in states/bits (L)
        c_length: The number of states/bits to consider for the prediction (N).
        mwsp: The minimal weight state position.
        """
        if size < 1:
            raise ValueError("Size (size or I) must be > 0.")
        if h_length < 1:
            raise ValueError("History length (h_length or L) must be > 0.")
        if h_length > 2**31-1:
            raise ValueError("History length (h_length or L) must be <= 2**31-1.")
        if c_length < 1:
            raise ValueError("Consideration length (c_length or N) must be > 0.")
        if c_length > h_length:
            raise ValueError("Consideration length (c_length or N) must be <= history length (h_length or L).")
        if int(log2(size) + c_length * 2 / 3) + 1 > 56:
            _logger.warning("BHPT size * consideration length may exceed 2**56-1 reducing or eliminating the influence of the oldest states.") 
        self._h_length: int = h_length
        self._c_length: int = c_length
        self._mwsp: bool = mwsp
        self._h_table: NDArray[int8] = zeros((size, h_length), dtype=int8)
        self._weights: NDArray[float64] = zeros(size, dtype=float64)
        self._valid: NDArray[int8] = zeros(size, dtype=int8)
        self._state_weights: NDArray[float64] = self._default_state_weights()
        self._probabilities: NDArray[float64] = zeros(size, dtype=float64)
        self._defer: bool = defer
        self._modified: bool = True  # Forces calculation of probabilities on first get()

    def __repr__(self) -> str:
        """Returns the string representation of the BHPT."""
        return f"binary_history_probability_table(size={self._h_table.shape[0]}, h_length={self._h_length}," + \
            f"c_length={self._c_length}, mwsp={self._mwsp}, defer={self._defer})"
    
    def __getitem__(self, index: int) -> NDArray[int8]:
        """Returns the history at the given index."""
        return self._h_table[index]
    
    def __setitem__(self, index: int, state: bool | list[int | bool] | NDArray[int8]) -> None:
        """Sets the most recent state history at the given index.
        
        If state is a list or NDArray, the most recent state is the last element in the list.
        """
        if isinstance(state, bool):
            self._valid[index] = 1
            self._h_table[index][1:] = self._h_table[index][0:-1]
            self._h_table[index][0] = int8(state)
        else:
            if isinstance(state, list):
                state = array(state, dtype=int8)
            if state.shape[0] > self._h_length:
                state = state[-self._h_length:]
            if not len(state):
                return # Nothing to do
            self._valid[index] = 1
            self._h_table[index][self._h_length - state.shape[0]:] = state[::-1]
        self._modified = True
        if not self._defer:
            c_table_index: NDArray[int8] = self._h_table[index][:self._c_length]
            if self._mwsp:
                c_table_index[self._c_length - 1] = 1
            self._weights[index] = (self._state_weights * c_table_index).sum()

    def _default_state_weights(self) -> NDArray[float64]:
        """Returns the default state weights."""
        return default_state_weights(self._c_length)
    
    def _update_probabilities(self) -> None:
        """Updates the probabilities based on the current state of the BHPT."""
        if self._defer:
            c_table: NDArray[int8] = self._h_table[:, :self._c_length]
            if self._mwsp:
                c_table[::][self._c_length - 1] = self._valid
            self._weights = (self._state_weights * c_table).sum(axis=1)
        self._probabilities = self._weights / self._weights.sum()

    def get(self) -> int:
        """Returns a random index based on weighted history probability"""
        if self._modified:
            self._update_probabilities()
            self._modified = False
        return np_choice(self._h_table.shape[0], p=self._probabilities)

    def get_many(self, num: int) -> NDArray[int32]:
        """Returns num random indices based on weighted history probability"""
        if self._modified:
            self._update_probabilities()
            self._modified = False
        return np_choice(self._h_table.shape[0], num, p=self._probabilities)

    def insert(self, auto_remove=True) -> int:
        """Returns the index of an empty entry in the BHPT.

        Args
        ----
        auto_remove: If true, removes the least weighted entry in the BHPT if there are no empty entries.
            In the event more than one entry have the same least weight, one of them is chosen at random.

        Returns
        -------
        The index of an empty entry in the BHPT.
        """
        num_empty: int = self._h_table.shape[0] - self._valid.sum() 
        if num_empty > 0:
            index = int(argmin(self._valid))
            self._valid[index] = 1
            return index
        if auto_remove:
            index = int(argmin(self._weights))
            self.remove(index)
            return index
        raise ValueError("BHPT is full and auto_remove is False.")
    
    def remove(self, index: int) -> None:
        """Removes the entry at the given index in the BHPT."""
        self._valid[index] = 0
        self._weights[index] = 0.0
        self._h_table[index] = 0
        self._modified = True

    def set_mwsp(self, mwsp: bool) -> None:
        """Sets the minimal weight state position (mwsp) to the given index."""
        self._mwsp = mwsp
        self._modified = True

    def set_defer(self, defer: bool) -> None:
        """Sets the defer flag to the given value.

        Args
        ----
        defer: The value to set the defer flag to.
        """
        self._defer = defer
        self._modified = True
