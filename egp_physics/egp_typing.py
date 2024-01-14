"""Common Erasmus GP Types."""
from typing import Literal
from egp_types.aGC import aGC

# Rows to insert above
InsertRow = Literal["A", "B", "F", "O", "P", "I", "Z"]

# A new GC definition
NewGCDef = tuple[aGC, dict[int, aGC]]
