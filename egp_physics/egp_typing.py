"""Common Erasmus GP Types."""
from egp_types.aGC import aGC

# Insertion work
Work = tuple[aGC, aGC, str]
WorkStack = list[Work]

# A new GC definition
NewGCDef = tuple[aGC, dict[int, aGC]]
