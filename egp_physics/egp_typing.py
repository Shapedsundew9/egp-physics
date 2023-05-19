"""Common Erasmus GP Types."""
from typing import Literal
from egp_types.aGC import aGC

# Rows to insert above
InsertRow = Literal['A', 'B', 'F', 'O', 'P', 'I']

# Insertion work
# (Target GC, Insert GC, above row)
Work = tuple[aGC, aGC, InsertRow]
WorkStack = list[Work]

# A new GC definition
NewGCDef = tuple[aGC, dict[int, aGC]]
