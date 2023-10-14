"""Test GC insertion operations.

See https://docs.google.com/spreadsheets/d/1YQjrM91e5x30VUIRzipNYX3W7yiFlg6fy9wKbMTx1iY/edit?usp=sharing
FIXME: These operations need test cases using the full stack to prove TGC == RGC after insertion.
"""

from collections import Counter
from hashlib import md5
from json import load, dump
from logging import DEBUG, NullHandler, getLogger, Logger
from os.path import dirname, join
from pprint import pformat
from random import choice, randint
from statistics import stdev
from copy import deepcopy
from egp_physics.fundamental import _insert_gc

# Load the results file.
_RESULTS_FILE = "test_physics_results.json"
with open(join(dirname(__file__), "data/", _RESULTS_FILE), "r") as file_ptr:
    results = load(file_ptr)


# Statistics are based on this number of iterations
STATS_N = 1000


# Logging
_logger: Logger = getLogger(__name__)
_logger.addHandler(NullHandler())
_LOG_DEBUG: bool = _logger.isEnabledFor(DEBUG)


# Set to True to generate test_physics_results.json
_GENERATE_RESULTS_FILE = True
_results = deepcopy(results)


def _write_results():
    """Dump the results to a local file.

    The dumped file can be used as the results file.
    """
    with open(_RESULTS_FILE, "w") as file_ptr:
        dump(_results, file_ptr, indent=4, sort_keys=True)


def test_basic_insert_1_simple():
    """Test case #1 of GC insertion."""
    _logger.info("Test case: test_basic_insert_1_simple")
    results_key = "basic_insert_1"
    tgc = eGC({"O": [["I", 0, 2]], "U": [["I", 1, 2]]}, sv=False)
    igc = eGC({"O": [["I", 0, 2]], "U": [["I", 1, 2]]}, sv=False)

    if _GENERATE_RESULTS_FILE:
        _results[results_key] = {}
        for _ in range(STATS_N):
            graph = _insert_gc(None, tgc, igc, "A")[0]["graph"]
            code = md5(bytearray(pformat(graph), encoding="ascii")).hexdigest()
            _results[results_key].setdefault(code, graph)
        _write_results()
    else:
        for above_row in "ABO":
            graph = _insert_gc(None, tgc, igc, above_row)[0]["graph"]
            code = md5(bytearray(pformat(graph), encoding="ascii")).hexdigest()
            assert code in results[results_key].keys()


def test_basic_insert_2_simple():
    """Test case #2 of GC insertion.

    Two codons with compatible input and output types are stacked.
    The compatible types mean a steady state exception is avoided
    (hence it is a 'basic' insert). Connectivity within the constraints
    of types and insertion location is random and so several variants
    may be created and one of which is correct.
    """
    _logger.info("Test case: test_basic_insert_2_simple")
    results_key = "basic_insert_2"
    tgc = mGC(
        igraph=gc_graph({"A": [["I", 0, 2], ["I", 1, 2]], "O": [["A", 0, 2]]}), sv=False
    )
    igc = eGC({"O": [["I", 0, 2]], "U": [["I", 1, 2]]}, sv=False)

    if _GENERATE_RESULTS_FILE:
        _results[results_key] = {}
        for _ in range(STATS_N):
            graph = _insert_gc(None, tgc, igc, "A")[0]["graph"]
            code = md5(bytearray(pformat(graph), encoding="ascii")).hexdigest()
            _results[results_key].setdefault(code, graph)
        _write_results()
    else:
        graph = _insert_gc(None, tgc, igc, "A")[0]["graph"]
        code = md5(bytearray(pformat(graph), encoding="ascii")).hexdigest()
        assert code in results[results_key].keys()


def test_basic_insert_3_simple():
    """Test case #3 of GC insertion."""
    _logger.info("Test case: test_basic_insert_3_simple")
    results_key = "basic_insert_3"
    tgc = mGC(
        igraph=gc_graph({"A": [["I", 0, 2], ["I", 1, 2]], "O": [["A", 0, 2]]}), sv=False
    )
    igc = eGC({"O": [["I", 0, 2]], "U": [["I", 1, 2]]}, sv=False)

    if _GENERATE_RESULTS_FILE:
        _results[results_key] = {}
        for _ in range(STATS_N):
            graph = _insert_gc(None, tgc, igc, "B")[0]["graph"]
            code = md5(bytearray(pformat(graph), encoding="ascii")).hexdigest()
            _results[results_key].setdefault(code, graph)
        _write_results()
    else:
        for above_row in "BO":
            graph = _insert_gc(None, tgc, igc, above_row)[0]["graph"]
            code = md5(bytearray(pformat(graph), encoding="ascii")).hexdigest()
            assert code in results[results_key].keys()


def test_basic_insert_4_simple():
    """Test case #4 of GC insertion.

    NOTE: rgc for case 5 is the same as in case 4 when tgc & igc are the same.
    The difference comes in fgc.
    """
    _logger.info("Test case: test_basic_insert_4_simple")
    graph = {
        "A": [["I", 1, 2], ["I", 1, 2]],
        "B": [["I", 0, 2], ["I", 1, 2]],
        "O": [["B", 0, 2]],
        "U": [["A", 0, 2]],
    }
    results_key = "basic_insert_4"
    tgc = mGC(igraph=gc_graph(graph), sv=False)
    igc = eGC({"O": [["I", 0, 2]], "U": [["I", 1, 2]]}, sv=False)

    if _GENERATE_RESULTS_FILE:
        _results[results_key] = {}
        for _ in range(STATS_N):
            rgc, fgc_dict = _insert_gc(None, tgc, igc, "A")
            for gcg in (rgc, *tuple(fgc_dict.values())):
                graph = gcg["graph"]
                code = md5(bytearray(pformat(graph), encoding="ascii")).hexdigest()
                _results[results_key].setdefault(code, graph)
        _write_results()
    else:
        graph = _insert_gc(None, tgc, igc, "A")[0]["graph"]
        code = md5(bytearray(pformat(graph), encoding="ascii")).hexdigest()
        assert code in results[results_key].keys()


def test_basic_insert_5_simple():
    """Test case #5 of GC insertion.

    NOTE: rgc for case 5 is the same as in case 4 when tgc & igc are the same.
    The difference comes in fgc.
    """
    _logger.info("Test case: test_basic_insert_5_simple")
    graph = {
        "A": [["I", 1, 2], ["I", 1, 2]],
        "B": [["I", 0, 2], ["I", 1, 2]],
        "O": [["B", 0, 2]],
        "U": [["A", 0, 2]],
    }
    results_key = "basic_insert_5"
    tgc = mGC(igraph=gc_graph(graph), sv=False)
    igc = eGC({"O": [["I", 0, 2]], "U": [["I", 1, 2]]}, sv=False)

    if _GENERATE_RESULTS_FILE:
        _results[results_key] = {}
        for _ in range(STATS_N):
            graph = _insert_gc(None, tgc, igc, "B")[0]["graph"]
            code = md5(bytearray(pformat(graph), encoding="ascii")).hexdigest()
            _results[results_key].setdefault(code, graph)
        _write_results()
    else:
        graph = _insert_gc(None, tgc, igc, "B")[0]["graph"]
        code = md5(bytearray(pformat(graph), encoding="ascii")).hexdigest()
        assert code in results[results_key].keys()


def test_basic_insert_6_simple():
    """Test case #6 of GC insertion."""
    _logger.info("Test case: test_basic_insert_6_simple")
    graph = {
        "A": [["I", 1, 2], ["I", 1, 2]],
        "B": [["I", 0, 2], ["I", 1, 2]],
        "O": [["B", 0, 2]],
        "U": [["A", 0, 2]],
    }
    results_key = "basic_insert_6"
    tgc = mGC(igraph=gc_graph(graph), sv=False)
    igc = eGC({"O": [["I", 0, 2]], "U": [["I", 1, 2]]}, sv=False)

    if _GENERATE_RESULTS_FILE:
        _results[results_key] = {}
        for _ in range(STATS_N):
            graph = _insert_gc(None, tgc, igc, "O")[0]["graph"]
            code = md5(bytearray(pformat(graph), encoding="ascii")).hexdigest()
            _results[results_key].setdefault(code, graph)
        _write_results()
    else:
        graph = _insert_gc(None, tgc, igc, "O")[0]["graph"]
        code = md5(bytearray(pformat(graph), encoding="ascii")).hexdigest()
        assert code in results[results_key].keys()


# TODO: Add reference sanity tests. Found a bug where circular references occured.


def test_basic_insert_1_stats():
    """Check the statistics of case #1 in the basic scenario.

    Case #1 has 12 equally probable possiblities.
    """
    _logger.info("Test case: test_basic_insert_1_stats")
    tgc = eGC({"O": [["I", 0, 2]], "U": [["I", 1, 2]]}, sv=False)
    igc = eGC({"O": [["I", 0, 2]], "U": [["I", 1, 2]]}, sv=False)

    def f():
        return md5(
            bytearray(
                pformat(_insert_gc(None, tgc, igc, "A")[0]["graph"]), encoding="ascii"
            )
        ).hexdigest()

    def func(x):
        return [f() for _ in range(x)]

    unlikely = 0
    for _ in range(3):
        counts = Counter(func(int(STATS_N * 3)))
        for checksum in counts:
            assert checksum in results["basic_insert_1"].keys()
        if stdev(counts.values()) > 35.115:
            _logger.debug(f"Standard deviation = {stdev(counts.values())}")
            unlikely += 1

    if unlikely == 1:
        _logger.info("Suspicious: Random connection probability > 1 in 2149")
    elif unlikely == 2:
        _logger.warn("Very suspicious: Random connection probability > 1 in 4618201")
    elif unlikely == 3:
        _logger.error(
            "Something is wrong: Random connection probability > 1 in 9924513949"
        )
    assert unlikely < 3


def test_basic_insert_2_stats():
    """Check the statistics of case #2 in the basic scenario.

    The connectivity of the inserted GC is uniform random. To verify
    randomness is just that this test checks that the standard deviation
    of the allowable variants is at least moderately likely. ;)

    NB: The maths on this may be shakey!!!

    In this scenario there are 4 possible results.
    Empirical analysis shows the standard deviation of the occurances of the
    four results in 1000 runs of 1000 insertions is min:0.816, max:35.355,
    avg:14.370, stdev:5.927.

    3.5-sigma is 1 in 2149 = 35.115
    So if 1000 inserts have a stdev > 35.115 three times in a row that is
    about a 1 in 10 billion shot.
    """
    _logger.info("Test case: test_basic_insert_2_stats")
    tgc = mGC(
        igraph=gc_graph({"A": [["I", 0, 2], ["I", 1, 2]], "O": [["A", 0, 2]]}), sv=False
    )
    igc = eGC({"O": [["I", 0, 2]], "U": [["I", 1, 2]]}, sv=False)

    def f():
        return md5(
            bytearray(
                pformat(_insert_gc(None, tgc, igc, "A")[0]["graph"]), encoding="ascii"
            )
        ).hexdigest()

    def func(x):
        return [f() for _ in range(x)]

    unlikely = 0
    for _ in range(3):
        counts = Counter(func(STATS_N))
        for checksum in counts:
            assert checksum in results["basic_insert_2"].keys()
        if stdev(counts.values()) > 35.115:
            _logger.debug(f"Standard deviation = {stdev(counts.values())}")
            unlikely += 1

    if unlikely == 1:
        _logger.info("Suspicious: Random connection probability > 1 in 2149")
    elif unlikely == 2:
        _logger.warn("Very suspicious: Random connection probability > 1 in 4618201")
    elif unlikely == 3:
        _logger.error(
            "Something is wrong: Random connection probability > 1 in 9924513949"
        )
    assert unlikely < 3


def test_basic_insert_3_stats():
    """Check the statistics of case #3 in the basic scenario.

    Case #3 has 9 equally probable possiblities.
    """
    _logger.info("Test case: test_basic_insert_3_stats")
    tgc = mGC(
        igraph=gc_graph({"A": [["I", 0, 2], ["I", 1, 2]], "O": [["A", 0, 2]]}), sv=False
    )
    igc = eGC({"O": [["I", 0, 2]], "U": [["I", 1, 2]]}, sv=False)

    def f():
        return md5(
            bytearray(
                pformat(_insert_gc(None, tgc, igc, "B")[0]["graph"]), encoding="ascii"
            )
        ).hexdigest()

    def func(x):
        return [f() for _ in range(x)]

    unlikely = 0
    for _ in range(3):
        counts = Counter(func(int(STATS_N * 2.25)))
        for checksum in counts:
            assert checksum in results["basic_insert_3"].keys()
        if stdev(counts.values()) > 35.115:
            _logger.debug(f"Standard deviation = {stdev(counts.values())}")
            unlikely += 1

    if unlikely == 1:
        _logger.info("Suspicious: Random connection probability > 1 in 2149")
    elif unlikely == 2:
        _logger.warn("Very suspicious: Random connection probability > 1 in 4618201")
    elif unlikely == 3:
        _logger.error(
            "Something is wrong: Random connection probability > 1 in 9924513949"
        )
    assert unlikely < 3


def test_basic_insert_4_stats():
    """Check the statistics of case #4 in the basic scenario.

    Case #4 FGC should have the same statistics as case #2.
    Only one RGC result is possible.
    """
    _logger.info("Test case: test_basic_insert_4_stats")
    graph = {
        "A": [["I", 1, 2], ["I", 1, 2]],
        "B": [["I", 0, 2], ["I", 1, 2]],
        "O": [["B", 0, 2]],
        "U": [["A", 0, 2]],
    }
    tgc = mGC(igraph=gc_graph(graph), sv=False)
    igc = eGC({"O": [["I", 0, 2]], "U": [["I", 1, 2]]}, sv=False)

    unlikely = 0
    for _ in range(3):
        rgc_codes = []
        fgc_codes = []
        for _ in range(STATS_N):
            result = _insert_gc(None, tgc, igc, "A")
            # TODO: Make code generation a function and reuse.
            rgc_codes.append(
                md5(
                    bytearray(pformat(result[0]["graph"]), encoding="ascii")
                ).hexdigest()
            )
            if _LOG_DEBUG:
                _logger.debug(f"RGC checksum {rgc_codes[-1]} for {result[0]['graph']}")
            for fgc in tuple(result[1].values())[1:]:
                fgc_codes.append(
                    md5(bytearray(pformat(fgc["graph"]), encoding="ascii")).hexdigest()
                )
                if _LOG_DEBUG:
                    _logger.debug(f"FGC checksum {fgc_codes[-1]} for {fgc['graph']}")
        rgc_counts = Counter(rgc_codes)
        fgc_counts = Counter(fgc_codes)

        for checksum in fgc_counts:
            assert checksum in results["basic_insert_4"].keys()
        for checksum in rgc_counts:
            assert checksum in results["basic_insert_4"].keys()
        if stdev(fgc_counts.values()) > 35.115:
            _logger.debug(f"Standard deviation = {stdev(fgc_counts.values())}")
            unlikely += 1

    if unlikely == 1:
        _logger.info("Suspicious: Random connection probability > 1 in 2149")
    elif unlikely == 2:
        _logger.warn("Very suspicious: Random connection probability > 1 in 4618201")
    elif unlikely == 3:
        _logger.error(
            "Something is wrong: Random connection probability > 1 in 9924513949"
        )
    assert unlikely < 3


def test_basic_insert_5_stats():
    """Check the statistics of case #5 in the basic scenario.

    Case #5 FGC should have the same statistics as case #3.
    Only one RGC result is possible.
    """
    _logger.info("Test case: test_basic_insert_5_stats")
    graph = {
        "A": [["I", 1, 2], ["I", 1, 2]],
        "B": [["I", 0, 2], ["I", 1, 2]],
        "O": [["B", 0, 2]],
        "U": [["A", 0, 2]],
    }
    tgc = mGC(igraph=gc_graph(graph), sv=False)
    igc = eGC({"O": [["I", 0, 2]], "U": [["I", 1, 2]]}, sv=False)

    unlikely = 0
    for _ in range(3):
        rgc_codes = []
        fgc_codes = []
        for _ in range(STATS_N):
            result = _insert_gc(None, tgc, igc, "B")
            # TODO: Make code generation a function and reuse.
            rgc_codes.append(
                md5(
                    bytearray(pformat(result[0]["graph"]), encoding="ascii")
                ).hexdigest()
            )
            if _LOG_DEBUG:
                _logger.debug(f"RGC checksum {rgc_codes[-1]} for {result[0]['graph']}")
            for fgc in tuple(result[1].values())[1:]:
                fgc_codes.append(
                    md5(bytearray(pformat(fgc["graph"]), encoding="ascii")).hexdigest()
                )
                if _LOG_DEBUG:
                    _logger.debug(f"FGC checksum {fgc_codes[-1]} for {fgc['graph']}")
        rgc_counts = Counter(rgc_codes)
        fgc_counts = Counter(fgc_codes)

        for checksum in fgc_counts:
            assert checksum in results["basic_insert_5"].keys()
        for checksum in rgc_counts:
            assert checksum in results["basic_insert_5"].keys()
        if stdev(fgc_counts.values()) > 35.115:
            _logger.debug(f"Standard deviation = {stdev(fgc_counts.values())}")
            unlikely += 1

    if unlikely == 1:
        _logger.info("Suspicious: Random connection probability > 1 in 2149")
    elif unlikely == 2:
        _logger.warn("Very suspicious: Random connection probability > 1 in 4618201")
    elif unlikely == 3:
        _logger.error(
            "Something is wrong: Random connection probability > 1 in 9924513949"
        )
    assert unlikely < 3


def test_basic_insert_6_stats():
    """Check the statistics of case #6 in the basic scenario.

    Case #6 FGC should have the same statistics as case #3.
    Only one RGC result is possible.
    """
    _logger.info("Test case: test_basic_insert_6_stats")
    graph = {
        "A": [["I", 1, 2], ["I", 1, 2]],
        "B": [["I", 0, 2], ["I", 1, 2]],
        "O": [["B", 0, 2]],
        "U": [["A", 0, 2]],
    }
    tgc = mGC(igraph=gc_graph(graph), sv=False)
    igc = eGC({"O": [["I", 0, 2]], "U": [["I", 1, 2]]}, sv=False)

    def f():
        return md5(
            bytearray(
                pformat(_insert_gc(None, tgc, igc, "O")[0]["graph"]), encoding="ascii"
            )
        ).hexdigest()

    def func(x):
        return [f() for _ in range(x)]

    unlikely = 0
    for _ in range(3):
        counts = Counter(func(STATS_N))
        for checksum in counts:
            assert checksum in results["basic_insert_6"].keys()
        if stdev(counts.values()) > 35.115:
            _logger.debug(f"Standard deviation = {stdev(counts.values())}")
            unlikely += 1

    if unlikely == 1:
        _logger.info("Suspicious: Random connection probability > 1 in 2149")
    elif unlikely == 2:
        _logger.warn("Very suspicious: Random connection probability > 1 in 4618201")
    elif unlikely == 3:
        _logger.error(
            "Something is wrong: Random connection probability > 1 in 9924513949"
        )
    assert unlikely < 3


def test_random_homogeneous_insertion():
    """Randomly insert eGC's which all have the same type inputs and outputs.

    Create a pool of 10 GC's with random numbers of inputs and outputs but all have
    the same EP type.

    Randomly select GCs from the pool to insert into each other adding the resultant
    GC back to the pool.

    At the end there is a mush of GC's. They should be valid GC's & no steady state
    exception should occur.
    """
    _logger.info("Test case: test_random_homogeneous_insertion")
    gc_list = [
        mGC(eGC(inputs=[2] * randint(1, 8), outputs=[2] * randint(1, 8)), sv=False)
        for _ in range(10)
    ]
    for _ in range(1000):
        tgc = choice(gc_list)
        igc = choice(gc_list)
        gc_list.append(_insert_gc(None, tgc, igc, choice("ABO"))[0])
        assert gc_list[-1]["igraph"].validate()
