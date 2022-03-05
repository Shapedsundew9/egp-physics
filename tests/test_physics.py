"""Test GC operations.

This test module assumes it has access to a postgresql instance as configured in
data/test_glib_config.json. The user requires database CREATE & DELETE rights.
"""

from collections import Counter
from hashlib import md5
from json import load
from logging import DEBUG, INFO, WARN, ERROR, FATAL, NullHandler, getLogger
from os.path import dirname, join
from pprint import pformat
from random import choice, randint
from statistics import stdev

from egp_physics.gc_graph import gc_graph
from egp_physics.gc_type import eGC, mGC
from egp_physics.physics import stablize

# Load the results file.
with open(join(dirname(__file__), "data/test_physics_results.json"), "r") as file_ptr:
    results = load(file_ptr)


STATS_N = 1000


# Logging
_logger = getLogger(__name__)
_logger.addHandler(NullHandler())
_LOG_DEBUG = _logger.isEnabledFor(DEBUG)
_LOG_INFO = _logger.isEnabledFor(INFO)
_LOG_WARN = _logger.isEnabledFor(WARN)
_LOG_ERROR = _logger.isEnabledFor(ERROR)
_LOG_FATAL = _logger.isEnabledFor(FATAL)


def test_basic_insert_1_simple():
    """Test case #1 of GC insertion."""
    tgc = eGC(inputs=(2, 2), outputs=(2,), sv=False)
    igc = eGC(inputs=(2, 2), outputs=(2,), sv=False)
    graph = stablize(None, tgc, igc, 'A')[0]['graph']
    code = md5(bytearray(pformat(graph), encoding='ascii')).hexdigest()
    """
    # Debug code
    graphs = {}
    for _ in range(100):
        graph = gc_insert(None, tgc, igc, 'A')[0]['graph']
        code = md5(bytearray(pformat(graph), encoding='ascii')).hexdigest()
        graphs[code] = graph
    print(pformat(graphs))
    """
    assert code in results['basic_insert_1'].keys()


def test_basic_insert_2_simple():
    """Test case #2 of GC insertion.

    Two codons with compatible input and output types are stacked.
    The compatible types mean a steady state exception is avoided
    (hence it is a 'basic' insert). Connectivity within the constraints
    of types and insertion location is random and so several variants
    may be created and one of which is correct.
    """
    tgc = mGC(igraph=gc_graph({'A': [['I', 0, 2], ['I', 1, 2]], 'O': [['A', 0, 2]]}), sv=False)
    igc = eGC(inputs=(2, 2), outputs=(2,), sv=False)
    graph = stablize(None, tgc, igc, 'A')[0]['graph']
    code = md5(bytearray(pformat(graph), encoding='ascii')).hexdigest()

    """
    # Debug code
    graphs = {}
    for _ in range(100):
        graph = gc_insert(None, tgc, igc, 'A')[0]['graph']
        code = md5(bytearray(pformat(graph), encoding='ascii')).hexdigest()
        graphs[code] = graph
    print(pformat(graphs))
    """
    assert code in results['basic_insert_2'].keys()


def test_basic_insert_3_simple():
    """Test case #3 of GC insertion."""
    tgc = mGC(igraph=gc_graph({'A': [['I', 0, 2], ['I', 1, 2]], 'O': [['A', 0, 2]]}), sv=False)
    igc = eGC(inputs=(2, 2), outputs=(2,), sv=False)
    graph = stablize(None, tgc, igc, 'B')[0]['graph']
    code = md5(bytearray(pformat(graph), encoding='ascii')).hexdigest()

    """
    # Debug code
    graphs = {}
    for _ in range(100):
        graph = gc_insert(None, tgc, igc, 'B')[0]['graph']
        code = md5(bytearray(pformat(graph), encoding='ascii')).hexdigest()
        graphs[code] = graph
    print(pformat(graphs))
    """
    assert code in results['basic_insert_3'].keys()


def test_basic_insert_4_simple():
    """Test case #4 of GC insertion."""
    graph = {
        'A': [['I', 1, 2], ['I', 1, 2]],
        'B': [['I', 0, 2], ['I', 1, 2]],
        'O': [['B', 0, 2]],
        'U': [['A', 0, 2]]
    }
    tgc = mGC(igraph=gc_graph(graph), sv=False)
    igc = eGC(inputs=(2, 2), outputs=(2,), sv=False)
    graph = stablize(None, tgc, igc, 'A')[0]['graph']
    code = md5(bytearray(pformat(graph), encoding='ascii')).hexdigest()

    """
    # Debug code
    graphs = {}
    # rgc result is the sole possible result for the simple test
    graph = gc_insert(None, tgc, igc, 'A')[0]['graph']
    code = md5(bytearray(pformat(graph), encoding='ascii')).hexdigest()
    graphs[code] = graph
    for _ in range(100):
        # fgc results are for use in test_basic_insert_4_stats test
        graph = gc_insert(None, tgc, igc, 'A')[1]['graph']
        code = md5(bytearray(pformat(graph), encoding='ascii')).hexdigest()
        graphs[code] = graph
    print(pformat(graphs))
    """
    assert code in results['basic_insert_4'].keys()


def test_basic_insert_5_simple():
    """Test case #5 of GC insertion."""
    graph = {
        'A': [['I', 1, 2], ['I', 1, 2]],
        'B': [['I', 0, 2], ['I', 1, 2]],
        'O': [['B', 0, 2]],
        'U': [['A', 0, 2]]
    }
    tgc = mGC(igraph=gc_graph(graph), sv=False)
    igc = eGC(inputs=(2, 2), outputs=(2,), sv=False)
    graph = stablize(None, tgc, igc, 'B')[0]['graph']
    code = md5(bytearray(pformat(graph), encoding='ascii')).hexdigest()

    """
    # Debug code
    graphs = {}
    # rgc result is the sole possible result for the simple test
    graph = gc_insert(None, tgc, igc, 'B')[0]['graph']
    code = md5(bytearray(pformat(graph), encoding='ascii')).hexdigest()
    graphs[code] = graph
    for _ in range(100):
        # fgc results are for use in test_basic_insert_4_stats test
        graph = gc_insert(None, tgc, igc, 'B')[1]['graph']
        code = md5(bytearray(pformat(graph), encoding='ascii')).hexdigest()
        graphs[code] = graph
    print(pformat(graphs))
    """
    assert code in results['basic_insert_5'].keys()


def test_basic_insert_6_simple():
    """Test case #6 of GC insertion."""
    graph = {
        'A': [['I', 1, 2], ['I', 1, 2]],
        'B': [['I', 0, 2], ['I', 1, 2]],
        'O': [['B', 0, 2]],
        'U': [['A', 0, 2]]
    }
    tgc = mGC(igraph=gc_graph(graph), sv=False)
    igc = eGC(inputs=(2, 2), outputs=(2,), sv=False)
    graph = stablize(None, tgc, igc, 'O')[0]['graph']
    code = md5(bytearray(pformat(graph), encoding='ascii')).hexdigest()

    """
    # Debug code
    graphs = {}
    for _ in range(100):
        graph = gc_insert(None, tgc, igc, 'O')[0]['graph']
        code = md5(bytearray(pformat(graph), encoding='ascii')).hexdigest()
        graphs[code] = graph
    print(pformat(graphs))
    """
    assert code in results['basic_insert_6'].keys()


def test_basic_insert_1_stats():
    """Check the statistics of case #1 in the basic scenario.

    Case #1 has 12 equally probable possiblities.
    """
    tgc = eGC(inputs=(2, 2), outputs=(2,), sv=False)
    igc = eGC(inputs=(2, 2), outputs=(2,), sv=False)
    def f(): return md5(bytearray(pformat(stablize(None, tgc, igc, 'A')[0]['graph']), encoding='ascii')).hexdigest()
    def func(x): return [f() for _ in range(x)]
    unlikely = 0
    for _ in range(3):
        counts = Counter(func(int(STATS_N * 3)))
        for checksum in counts:
            assert checksum in results['basic_insert_1'].keys()
        if stdev(counts.values()) > 35.115:
            _logger.debug(f"Standard deviation = {stdev(counts.values())}")
            unlikely += 1

    if unlikely == 1:
        _logger.info("Suspicious: Random connection probability > 1 in 2149")
    elif unlikely == 2:
        _logger.warn("Very suspicious: Random connection probability > 1 in 4618201")
    elif unlikely == 3:
        _logger.error("Something is wrong: Random connection probability > 1 in 9924513949")
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
    tgc = mGC(igraph=gc_graph({'A': [['I', 0, 2], ['I', 1, 2]], 'O': [['A', 0, 2]]}), sv=False)
    igc = eGC(inputs=(2, 2), outputs=(2,), sv=False)
    def f(): return md5(bytearray(pformat(stablize(None, tgc, igc, 'A')[0]['graph']), encoding='ascii')).hexdigest()
    def func(x): return [f() for _ in range(x)]
    unlikely = 0
    for _ in range(3):
        counts = Counter(func(STATS_N))
        for checksum in counts:
            assert checksum in results['basic_insert_2'].keys()
        if stdev(counts.values()) > 35.115:
            _logger.debug(f"Standard deviation = {stdev(counts.values())}")
            unlikely += 1

    if unlikely == 1:
        _logger.info("Suspicious: Random connection probability > 1 in 2149")
    elif unlikely == 2:
        _logger.warn("Very suspicious: Random connection probability > 1 in 4618201")
    elif unlikely == 3:
        _logger.error("Something is wrong: Random connection probability > 1 in 9924513949")
    assert unlikely < 3


def test_basic_insert_3_stats():
    """Check the statistics of case #3 in the basic scenario.

    Case #3 has 9 equally probable possiblities.
    """
    tgc = mGC(igraph=gc_graph({'A': [['I', 0, 2], ['I', 1, 2]], 'O': [['A', 0, 2]]}), sv=False)
    igc = eGC(inputs=(2, 2), outputs=(2,), sv=False)
    def f(): return md5(bytearray(pformat(stablize(None, tgc, igc, 'B')[0]['graph']), encoding='ascii')).hexdigest()
    def func(x): return [f() for _ in range(x)]
    unlikely = 0
    for _ in range(3):
        counts = Counter(func(int(STATS_N * 2.25)))
        for checksum in counts:
            assert checksum in results['basic_insert_3'].keys()
        if stdev(counts.values()) > 35.115:
            _logger.debug(f"Standard deviation = {stdev(counts.values())}")
            unlikely += 1

    if unlikely == 1:
        _logger.info("Suspicious: Random connection probability > 1 in 2149")
    elif unlikely == 2:
        _logger.warn("Very suspicious: Random connection probability > 1 in 4618201")
    elif unlikely == 3:
        _logger.error("Something is wrong: Random connection probability > 1 in 9924513949")
    assert unlikely < 3


def test_basic_insert_4_stats():
    """Check the statistics of case #4 in the basic scenario.

    Case #4 FGC should have the same statistics as case #2.
    Only one RGC result is possible.
    """
    graph = {
        'A': [['I', 1, 2], ['I', 1, 2]],
        'B': [['I', 0, 2], ['I', 1, 2]],
        'O': [['B', 0, 2]],
        'U': [['A', 0, 2]]
    }
    tgc = mGC(igraph=gc_graph(graph), sv=False)
    igc = eGC(inputs=(2, 2), outputs=(2,), sv=False)

    unlikely = 0
    for _ in range(3):
        rgc_codes = []
        fgc_codes = []
        for _ in range(STATS_N):
            result = stablize(None, tgc, igc, 'A')
            #TODO: Make code generation a function and reuse.
            rgc_codes.append(md5(bytearray(pformat(result[0]['graph']), encoding='ascii')).hexdigest())
            if _LOG_DEBUG:
                _logger.debug(f"RGC checksum {rgc_codes[-1]} for {result[0]['graph']}")
            for fgc in tuple(result[1].values())[1:]:
                fgc_codes.append(md5(bytearray(pformat(fgc['graph']), encoding='ascii')).hexdigest())
                if _LOG_DEBUG:
                    _logger.debug(f"FGC checksum {fgc_codes[-1]} for {fgc['graph']}")
        rgc_counts = Counter(rgc_codes)
        fgc_counts = Counter(fgc_codes)

        for checksum in fgc_counts:
            assert checksum in results['basic_insert_4'].keys()
        for checksum in rgc_counts:
            assert checksum in results['basic_insert_4'].keys()
        if stdev(fgc_counts.values()) > 35.115:
            _logger.debug(f"Standard deviation = {stdev(fgc_counts.values())}")
            unlikely += 1

    if unlikely == 1:
        _logger.info("Suspicious: Random connection probability > 1 in 2149")
    elif unlikely == 2:
        _logger.warn("Very suspicious: Random connection probability > 1 in 4618201")
    elif unlikely == 3:
        _logger.error("Something is wrong: Random connection probability > 1 in 9924513949")
    assert unlikely < 3


def test_basic_insert_5_stats():
    """Check the statistics of case #5 in the basic scenario.

    Case #5 FGC should have the same statistics as case #3.
    Only one RGC result is possible.
    """
    graph = {
        'A': [['I', 1, 2], ['I', 1, 2]],
        'B': [['I', 0, 2], ['I', 1, 2]],
        'O': [['B', 0, 2]],
        'U': [['A', 0, 2]]
    }
    tgc = mGC(igraph=gc_graph(graph), sv=False)
    igc = eGC(inputs=(2, 2), outputs=(2,), sv=False)

    unlikely = 0
    for _ in range(3):
        rgc_codes = []
        fgc_codes = []
        for _ in range(STATS_N):
            result = stablize(None, tgc, igc, 'B')
            #TODO: Make code generation a function and reuse.
            rgc_codes.append(md5(bytearray(pformat(result[0]['graph']), encoding='ascii')).hexdigest())
            if _LOG_DEBUG:
                _logger.debug(f"RGC checksum {rgc_codes[-1]} for {result[0]['graph']}")
            for fgc in tuple(result[1].values())[1:]:
                fgc_codes.append(md5(bytearray(pformat(fgc['graph']), encoding='ascii')).hexdigest())
                if _LOG_DEBUG:
                    _logger.debug(f"FGC checksum {fgc_codes[-1]} for {fgc['graph']}")
        rgc_counts = Counter(rgc_codes)
        fgc_counts = Counter(fgc_codes)

        for checksum in fgc_counts:
            assert checksum in results['basic_insert_5'].keys()
        for checksum in rgc_counts:
            assert checksum in results['basic_insert_5'].keys()
        if stdev(fgc_counts.values()) > 35.115:
            _logger.debug(f"Standard deviation = {stdev(fgc_counts.values())}")
            unlikely += 1

    if unlikely == 1:
        _logger.info("Suspicious: Random connection probability > 1 in 2149")
    elif unlikely == 2:
        _logger.warn("Very suspicious: Random connection probability > 1 in 4618201")
    elif unlikely == 3:
        _logger.error("Something is wrong: Random connection probability > 1 in 9924513949")
    assert unlikely < 3


def test_basic_insert_6_stats():
    """Check the statistics of case #6 in the basic scenario.

    Case #6 FGC should have the same statistics as case #3.
    Only one RGC result is possible.
    """
    graph = {
        'A': [['I', 1, 2], ['I', 1, 2]],
        'B': [['I', 0, 2], ['I', 1, 2]],
        'O': [['B', 0, 2]],
        'U': [['A', 0, 2]]
    }
    tgc = mGC(igraph=gc_graph(graph), sv=False)
    igc = eGC(inputs=(2, 2), outputs=(2,), sv=False)
    def f(): return md5(bytearray(pformat(stablize(None, tgc, igc, 'O')[0]['graph']), encoding='ascii')).hexdigest()
    def func(x): return [f() for _ in range(x)]
    unlikely = 0
    for _ in range(3):
        counts = Counter(func(STATS_N))
        for checksum in counts:
            assert checksum in results['basic_insert_6'].keys()
        if stdev(counts.values()) > 35.115:
            _logger.debug(f"Standard deviation = {stdev(counts.values())}")
            unlikely += 1

    if unlikely == 1:
        _logger.info("Suspicious: Random connection probability > 1 in 2149")
    elif unlikely == 2:
        _logger.warn("Very suspicious: Random connection probability > 1 in 4618201")
    elif unlikely == 3:
        _logger.error("Something is wrong: Random connection probability > 1 in 9924513949")
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
    gc_list = [mGC(eGC(inputs=[2]*randint(1, 8), outputs=[2]*randint(1, 8)), sv=False) for _ in range(10)]
    for _ in range(1000):
        tgc = choice(gc_list)
        igc = choice(gc_list)
        gc_list.append(stablize(None, tgc, igc, choice('ABO'))[0])
        assert gc_list[-1]['igraph'].validate()
