"""Test GC operations.

This test module assumes it has access to a postgresql instance as configured in
data/test_glib_config.json. The user requires database CREATE & DELETE rights.
"""


import pytest
from collections import Counter
from hashlib import md5
from pprint import pformat
from statistics import stdev
from random import randint


@pytest.mark.parametrize("inputs, outputs",
    (
        ((int32, int, int), ('float32',)),
        ((int, str, int), ('int', 'str', 'int')),
        ((int, object, list), ('int', 'dict'))
    )
)
def test_eGC(inputs, outputs):
    """eGC() returns an invalid GC so can only test for basic errors."""
    assert eGC(inputs, outputs)


@pytest.mark.good
def test_basic_insert_2_simple_A(small_gene_pool):
    """Test case #2 of GC insertion.

    Two codons with compatible input and output types are stacked.
    The compatible types mean a steady state exception is avoided
    (hence it is a 'basic' insert). Connectivity within the constraints
    of types and insertion location is random and so several variants
    may be created and one of which is correct.
    """
    _logger.info("Starting test test_basic_insert_2_simple_A")
    small_gene_pool.clear()
    tgc, igc = small_gene_pool[ADDI_SIGNATURE], small_gene_pool[SUBI_SIGNATURE]
    graph = gc_insert(tgc, igc, 'A')[0]['graph']
    code = md5(bytearray(pformat(graph), encoding='ascii')).hexdigest()

    """
    # Debug code
    if not code in results['basic_insert_2'].keys():
        graphs = {}
        for _ in range(100):
            graph = gc_insert(tgc, igc, 'A')[0]['graph']
            code = md5(bytearray(pformat(graph), encoding='ascii')).hexdigest()
            graphs[code] = graph
        print(pformat(graphs))
    """
    assert code in results['basic_insert_2'].keys()


@pytest.mark.good
def test_basic_insert_2_simple_B(small_gene_pool):
    """Test GC data after an insertion is as expected.

    When a GC is created from two others many parameters are updated.
    In this case two codons are stacked.
    """
    _logger.info("Starting test test_basic_insert_2_simple_B")
    small_gene_pool.clear()
    tgc, igc = small_gene_pool[ADDI_SIGNATURE], small_gene_pool[SUBI_SIGNATURE]
    gcs = gc_insert(tgc, igc, 'A')
    rgc_signature = gcs[0]['signature']
    small_gene_pool.add(gcs)
    small_gene_pool.push()
    rgc = small_gene_pool[rgc_signature]
    assert rgc['generation'] == 1
    assert rgc['raw_num_codons'] == 2
    assert rgc['code_depth'] == 2
    assert rgc['num_codes'] == 2


@pytest.mark.good
def test_basic_insert_2_stats(small_gene_pool):
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
    _logger.info("Starting test test_basic_insert_2_stats")
    small_gene_pool.clear()
    tgc, igc = small_gene_pool[ADDI_SIGNATURE], small_gene_pool[SUBI_SIGNATURE]
    func = lambda x: [md5(bytearray(pformat(gc_insert(tgc, igc, 'A')[0]['graph']), encoding='ascii')).hexdigest() for _ in range(x)]
    unlikely = 0
    for iteration in range(3):
        counts = Counter(func(STATS_N))
        for checksum in counts: assert checksum in results['basic_insert_2'].keys()
        if stdev(counts.values()) > 35.115:
            unlikely += 1

    if unlikely == 1: _logger.info("Suspicious: Random connection probability > 1 in 2149")
    elif unlikely == 2: _logger.warn("Very suspicious: Random connection probability > 1 in 4618201")
    elif unlikely == 3: _logger.error("Something is wrong: Random connection probability > 1 in 9924513949")
    assert unlikely < 3


@pytest.mark.good
@pytest.mark.parametrize("a, b", tuple((randint(0, 100), randint(0, 100)) for _ in range(10)))
def test_basic_insert_2_identity(small_gene_pool, a, b):
    """Inserting a graph should not change the function of the original.

    Insertion is a null operation from a functionality perspective.
    Do the test 100 times for a strong statistical likelihood of
    generating all variants and not getting a 'unlucky' correct
    mathematical result.
    """
    _logger.info("Starting test test_basic_insert_2_identity")
    small_gene_pool.clear()
    tgc, igc = small_gene_pool[ADDI_SIGNATURE], small_gene_pool[SUBI_SIGNATURE]
    gcs = gc_insert(tgc, igc, 'A')
    rgc_signature = gcs[0]['signature']
    small_gene_pool.add(gcs)
    rgc = small_gene_pool[rgc_signature]
    assert tgc['_func']((a, b)) == rgc['_func']((a, b))


@pytest.mark.good
def test_basic_insert_3_stats(small_gene_pool):
    """Check the statistics of case #3 in the basic scenario.

    Case #3 has 9 equally probable possiblities.
    """
    _logger.info("Starting test test_basic_insert_3_stats")
    small_gene_pool.clear()
    tgc, igc = small_gene_pool[ADDI_SIGNATURE], small_gene_pool[SUBI_SIGNATURE]
    func = lambda x: [md5(bytearray(pformat(gc_insert(tgc, igc, 'B')[0]['graph']), encoding='ascii')).hexdigest() for _ in range(x)]
    unlikely = 0
    for iteration in range(3):
        counts = Counter(func(int(STATS_N * 2.25)))
        for checksum in counts: assert checksum in results['basic_insert_3'].keys()
        if stdev(counts.values()) > 35.115:
            unlikely += 1

    if unlikely == 1: _logger.info("Suspicious: Random connection probability > 1 in 2149")
    elif unlikely == 2: _logger.warn("Very suspicious: Random connection probability > 1 in 4618201")
    elif unlikely == 3: _logger.error("Something is wrong: Random connection probability > 1 in 9924513949")
    assert unlikely < 3


@pytest.mark.good
@pytest.mark.parametrize("a, b", tuple((randint(0, 100), randint(0, 100)) for _ in range(10)))
def test_basic_insert_3_identity(small_gene_pool, a, b):
    """Inserting a graph should not change the function of the original.

    Insertion is a null operation from a functionality perspective.
    Do the test 10 times to avoid getting an 'unlucky' correct
    mathematical result.
    """
    _logger.info("Starting test test_basic_insert_3_identity")
    small_gene_pool.clear()
    tgc, igc = small_gene_pool[ADDI_SIGNATURE], small_gene_pool[SUBI_SIGNATURE]
    gcs = gc_insert(tgc, igc, 'B')
    rgc_signature = gcs[0]['signature']
    small_gene_pool.add(gcs)
    rgc = small_gene_pool[rgc_signature]
    assert tgc['_func']((a, b)) == rgc['_func']((a, b))


@pytest.mark.good
def test_basic_insert_4_stats(small_gene_pool):
    """Check the statistics of case #4 in the basic scenario.

    Case #4 FGC should have the same statistics as case #2.
    Only one RGC result is possible.
    """
    _logger.info("Starting test test_basic_insert_4_stats")
    small_gene_pool.clear()
    small_gene_pool.add([GEN1_GC])
    tgc = small_gene_pool[GEN1_GC['signature']]
    igc = small_gene_pool[ADDI_SIGNATURE]

    unlikely = 0
    for iteration in range(3):
        rgc_codes = []
        fgc_codes = []
        for _ in range(STATS_N):
            result = gc_insert(tgc, igc, 'A')
            rgc_codes.append(md5(bytearray(pformat(result[0]['graph']), encoding='ascii')).hexdigest())
            fgc_codes.append(md5(bytearray(pformat(result[1]['graph']), encoding='ascii')).hexdigest())
        rgc_counts = Counter(rgc_codes)
        fgc_counts = Counter(fgc_codes)

        for checksum in fgc_counts: assert checksum in results['basic_insert_4'].keys()
        for checksum in rgc_counts: assert checksum in results['basic_insert_4'].keys()
        if stdev(fgc_counts.values()) > 35.115:
            unlikely += 1

    if unlikely == 1: _logger.info("Suspicious: Random connection probability > 1 in 2149")
    elif unlikely == 2: _logger.warn("Very suspicious: Random connection probability > 1 in 4618201")
    elif unlikely == 3: _logger.error("Something is wrong: Random connection probability > 1 in 9924513949")
    assert unlikely < 3


@pytest.mark.good
@pytest.mark.parametrize("a, b", tuple((randint(0, 100), randint(0, 100)) for _ in range(10)))
def test_basic_insert_4_identity(small_gene_pool, a, b):
    """Inserting a graph should not change the function of the original.

    Insertion is a null operation from a functionality perspective.
    Do the test 10 times to avoid getting an 'unlucky' correct
    mathematical result.
    """
    _logger.info("Starting test test_basic_insert_4_identity")
    small_gene_pool.clear()
    small_gene_pool.add([GEN1_GC])
    tgc = small_gene_pool[GEN1_GC['signature']]
    igc = small_gene_pool[ADDI_SIGNATURE]
    gcs = gc_insert(tgc, igc, 'A')
    rgc_signature = gcs[0]['signature']
    small_gene_pool.add(gcs)
    rgc = small_gene_pool[rgc_signature]
    assert tgc['_func']((a, b)) == rgc['_func']((a, b))


@pytest.mark.good
def test_basic_insert_5_stats(small_gene_pool):
    """Check the statistics of case #5 in the basic scenario.

    Case #5 FGC should have the same statistics as case #3.
    Only one RGC result is possible.
    """
    _logger.info("Starting test test_basic_insert_5_stats")
    small_gene_pool.clear()
    small_gene_pool.add([GEN1_GC])
    tgc = small_gene_pool[GEN1_GC['signature']]
    igc = small_gene_pool[ADDI_SIGNATURE]

    unlikely = 0
    for iteration in range(3):
        rgc_codes = []
        fgc_codes = []
        for _ in range(int(STATS_N * 2.25)):
            result = gc_insert(tgc, igc, 'B')
            rgc_codes.append(md5(bytearray(pformat(result[0]['graph']), encoding='ascii')).hexdigest())
            fgc_codes.append(md5(bytearray(pformat(result[1]['graph']), encoding='ascii')).hexdigest())
        rgc_counts = Counter(rgc_codes)
        fgc_counts = Counter(fgc_codes)

        for checksum in fgc_counts: assert checksum in results['basic_insert_5'].keys()
        for checksum in rgc_counts: assert checksum in results['basic_insert_5'].keys()
        if stdev(fgc_counts.values()) > 35.115:
            unlikely += 1

    if unlikely == 1: _logger.info("Suspicious: Random connection probability > 1 in 2149")
    elif unlikely == 2: _logger.warn("Very suspicious: Random connection probability > 1 in 4618201")
    elif unlikely == 3: _logger.error("Something is wrong: Random connection probability > 1 in 9924513949")
    assert unlikely < 3


@pytest.mark.good
@pytest.mark.parametrize("a, b", tuple((randint(0, 100), randint(0, 100)) for _ in range(10)))
def test_basic_insert_5_identity(small_gene_pool, a, b):
    """Inserting a graph should not change the function of the original.

    Insertion is a null operation from a functionality perspective.
    Do the test 10 times to avoid getting an 'unlucky' correct
    mathematical result.
    """
    _logger.info("Starting test test_basic_insert_5_identity")
    small_gene_pool.clear()
    small_gene_pool.add([GEN1_GC])
    tgc = small_gene_pool[GEN1_GC['signature']]
    igc = small_gene_pool[ADDI_SIGNATURE]
    gcs = gc_insert(tgc, igc, 'B')
    rgc_signature = gcs[0]['signature']
    small_gene_pool.add(gcs)
    rgc = small_gene_pool[rgc_signature]
    assert tgc['_func']((a, b)) == rgc['_func']((a, b))


@pytest.mark.good
def test_basic_insert_6_stats(small_gene_pool):
    """Check the statistics of case #5 in the basic scenario.

    Case #6 should have the same statistics as case #2.
    """
    _logger.info("Starting test test_basic_insert_5_stats")
    small_gene_pool.clear()
    small_gene_pool.add([GEN1_GC])
    tgc = small_gene_pool[GEN1_GC['signature']]
    igc = small_gene_pool[ADDI_SIGNATURE]

    func = lambda x: [md5(bytearray(pformat(gc_insert(tgc, igc, 'O')[0]['graph']), encoding='ascii')).hexdigest() for _ in range(x)]
    unlikely = 0
    for iteration in range(3):
        counts = Counter(func(STATS_N))
        for checksum in counts: assert checksum in results['basic_insert_6'].keys()
        if stdev(counts.values()) > 35.115:
            unlikely += 1

    if unlikely == 1: _logger.info("Suspicious: Random connection probability > 1 in 2149")
    elif unlikely == 2: _logger.warn("Very suspicious: Random connection probability > 1 in 4618201")
    elif unlikely == 3: _logger.error("Something is wrong: Random connection probability > 1 in 9924513949")
    assert unlikely < 3


@pytest.mark.good
@pytest.mark.parametrize("a, b", tuple((randint(0, 100), randint(0, 100)) for _ in range(10)))
def test_basic_insert_6_identity(small_gene_pool, a, b):
    """Inserting a graph should not change the function of the original.

    Insertion is a null operation from a functionality perspective.
    Do the test 10 times to avoid getting an 'unlucky' correct
    mathematical result.
    """
    _logger.info("Starting test test_basic_insert_6_identity")
    small_gene_pool.clear()
    small_gene_pool.add([GEN1_GC])
    tgc = small_gene_pool[GEN1_GC['signature']]
    igc = small_gene_pool[ADDI_SIGNATURE]
    gcs = gc_insert(tgc, igc, 'O')
    rgc_signature = gcs[0]['signature']
    small_gene_pool.add(gcs)
    rgc = small_gene_pool[rgc_signature]
    assert tgc['_func']((a, b)) == rgc['_func']((a, b))


@pytest.mark.good
def test_basic_insert_1_stats(small_gene_pool):
    """Check the statistics of case #1 in the basic scenario."""
    _logger.info("Starting test test_basic_insert_1_stats")
    egc = eGC((int, int), (int,))
    small_gene_pool.clear()
    small_gene_pool.add([egc])
    tgc = small_gene_pool[egc['signature']]
    igc = small_gene_pool[ADDI_SIGNATURE]

    func = lambda x: [md5(bytearray(pformat(gc_insert(tgc, igc, 'A')[0]['graph']), encoding='ascii')).hexdigest() for _ in range(x)]
    unlikely = 0
    for iteration in range(3):
        counts = Counter(func(STATS_N * 3))
        for checksum in counts: assert checksum in results['basic_insert_1'].keys()
        if stdev(counts.values()) > 35.115:
            unlikely += 1

    if unlikely == 1: _logger.info("Suspicious: Random connection probability > 1 in 2149")
    elif unlikely == 2: _logger.warn("Very suspicious: Random connection probability > 1 in 4618201")
    elif unlikely == 3: _logger.error("Something is wrong: Random connection probability > 1 in 9924513949")
    assert unlikely < 3


@pytest.mark.good
def test_random_homogeneous_insertion(big_gene_pool):
    """Randomly insert GC's which all have the same type inputs and outputs."""
    _logger.info("Starting test random homogeneous insertion")
    query = {
        "input_types": {
            "operator": "contained by",
            "array": [283] * 8
        },
        "output_types": {
            "operator": "contained by",
            "array": [283] * 8
        },
        "order by": "RANDOM()",
        "limit": 1
    }
    tgc = big_gene_pool.gl([query])[0]
    for _ in range(10):
        igc = big_gene_pool.gl([query])[0]
        gcs = gc_insert(tgc, igc)
        rgc_signature = gcs[0]['signature']
        big_gene_pool.add(gcs)
        tgc = big_gene_pool[rgc_signature]
        big_gene_pool.push()

    # TODO: Need to assert on the inputs & outputs
    #big_gene_pool.save_gc_graph(rgc_signature, join(_log_location, "test_random_homogeneous_insertion_gc_graph"))
    #big_gene_pool.save_codon_graph(rgc_signature, join(_log_location, "test_random_homogeneous_insertion_codon_graph"))
    #git big_gene_pool.gl([query])


@pytest.mark.good
def test_initial_gc(big_gene_pool):
    """Create a GC that meets the input & output criteria."""
    _logger.info("Starting test of initial_gc()")
    register_gene_pool(big_gene_pool)
    signature = initial_gc([128], [128])
    # FIXME: This is not a valid assertion
    assert signature



def test_remove_gc(big_gene_pool):
    """Create a GC that meets the input & output criteria."""
    _logger.info("Starting test of initial_gc()")
    register_gene_pool(big_gene_pool)
    signature = initial_gc([128], [128])
