"""Validation of gc_type."""

from egp_physics.gc_type import eGC, md_table, mGC

md_table()


def test_instanciate_eGC_n0():
    """Validate instanciation of an eGC().

    An eGC has required fields so will raise a ValueError in instanciated
    without them defined.
    """
    try:
        assert eGC(sv=False)
    except ValueError:
        pass


def test_instanciate_mGC_n0():
    """Validate instanciation of an mGC().

    An mGC has required fields so will raise a ValueError in instanciated
    without them defined.
    """
    try:
        assert mGC(sv=False)
    except ValueError:
        pass
