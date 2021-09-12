"""Validation of gc_type."""

from egp_physics.gc_type import md_table, eGC, mGC

md_table()


def test_instanciate_eGC_n0():
    """Validate instanciation of an eGC().

    An eGC has required fields so will raise a ValueError in instanciated
    without them defined.
    """
    assert eGC(sv=False)


def test_instanciate_mGC_n0():
    """Validate instanciation of an mGC().

    An mGC has required fields so will raise a ValueError in instanciated
    without them defined.
    """
    assert mGC(sv=False)
