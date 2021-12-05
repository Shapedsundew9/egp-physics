"""gc_graph verficiation."""


from copy import deepcopy
from json import load
from os.path import dirname, join
from random import randint, random
from numpy.random import choice

import pytest
from egp_physics.ep_type import EP_TYPE_VALUES, INVALID_EP_TYPE_VALUE, asint
from egp_physics.gc_graph import (DESTINATION_ROWS, DST_EP, SOURCE_ROWS,
                                  SRC_EP, conn_idx, const_idx, gc_graph)


_TEST_RESULTS_JSON = 'data/test_gc_graph_results.json'
_VALID_STRUCTURES = (
    ('A', 'O'),
    ('C', 'O'),     # TODO: Is this valid?
    ('I', 'O'),     # TODO: Is this valid?
    ('A', 'O', 'U'),
    ('C', 'O', 'U'),
    ('I', 'O', 'U'),     # TODO: Is this valid?
    ('A', 'C', 'O'),
    ('A', 'I', 'O'),
    ('A', 'B', 'O'),
    ('I', 'C', 'O'),     # TODO: Is this valid?
    ('A', 'C', 'O', 'U'),
    ('A', 'I', 'O', 'U'),
    ('A', 'B', 'O', 'U'),
    ('I', 'C', 'O', 'U'),     # TODO: Is this valid?
    ('A', 'C', 'O', 'B'),
    ('A', 'I', 'O', 'B'),
    ('I', 'O', 'F', 'P'),
    ('A', 'C', 'O', 'B', 'U'),
    ('A', 'I', 'O', 'B', 'U'),
    ('I', 'O', 'F', 'P', 'U'),
    ('A', 'I', 'O', 'F', 'P'),
    ('I', 'C', 'O', 'F', 'P'),
    ('A', 'I', 'O', 'F', 'P', 'U'),
    ('I', 'C', 'O', 'F', 'P', 'U'),
    ('A', 'I', 'O', 'B', 'F', 'P'),
    ('A', 'I', 'O', 'B', 'F', 'P', 'U')
)


# Types are in string format for readability in the results file.
with open(join(dirname(__file__), _TEST_RESULTS_JSON), "r") as results_file:
    results = load(results_file)


def random_type(p=0.0):
    """Choose a random type.

    If a random type is selected the probability of each type is even.
    By default 'int' is returned.

    Args
    ----
    p (float): Probablity that the type is random (otherwise it is an 'int')

    Returns
    -------
    (str) The selected type string.
    """
    if random() < p:
        value = choice(tuple(EP_TYPE_VALUES))
        if value != INVALID_EP_TYPE_VALUE:
            return value
    return asint('builtins_int')


def random_graph(p=0.0, must_be_valid=False):  # noqa: C901
    """Create a random graph.

    The graph is not guaranteed to be valid when p > 0.0. If a destination row requires a type that
    is not present in any valid source row the graph cannot be normalized.

    Args
    ----
    p (float): 0.0 <= p <= 1.0 probability of choosing a random type on each type selection.

    Returns
    -------
    graph
    """
    valid = False
    while not valid:
        graph = gc_graph()
        structure = choice(_VALID_STRUCTURES)
        valid = False
        while not valid:
            destinations = {row: randint(1, 10) for row in structure if row in DESTINATION_ROWS and row not in ('F', 'U', 'P')}
            if 'F' in structure:
                destinations['F'] = 1
            sources = {row: randint(1, 8) for row in structure if row in SOURCE_ROWS and row not in ('U', 'P')}
            destination_types = [random_type(p) for row in destinations.values() for _ in range(row)]
            type_set = set(destination_types)
            valid = sum(sources.values()) >= len(type_set)
        source_types = [random_type(p) for _ in range(sum(sources.values()))]
        indices = choice(sum(sources.values()), len(type_set), replace=False)
        for idx in indices:
            source_types[idx] = type_set.pop()
        for _ in range(len(type_set)):
            source_types[randint(len(source_types))] = type_set

        for row in structure:
            if row not in ('U', 'P'):
                if row in DESTINATION_ROWS and any([src_row in structure for src_row in gc_graph.src_rows[row]]):
                    for i in range(destinations[row]):
                        rtype = destination_types.pop()
                        graph._add_ep([DST_EP, row, i, rtype, []])
                        if row == 'O' and 'P' in structure:
                            graph._add_ep([DST_EP, 'P', i, rtype, []])

                if row in SOURCE_ROWS:
                    for i in range(sources[row]):
                        ep = [SRC_EP, row, i, source_types.pop(), []]
                        if row == 'C':
                            ep.append('int(' + str(randint(-1000, 1000)) + ')')
                            ep[3] = asint('builtins_int')
                        graph._add_ep(ep)

        for _ in range(len(list(filter(graph.dst_filter(), graph.graph.values())))):
            graph.random_add_connection()
        #graph.reindex_row('C')
        graph.normalize()
        valid = graph.validate() or not must_be_valid
    return graph


@pytest.mark.parametrize("i, case", enumerate(results))
def test_graph_validation(i, case):
    """Verification the validate() method correctly functions."""
    gcg = gc_graph(case['graph'])
    assert i == case['i']
    assert case['valid'] == gcg.validate()
    if not case['valid']:
        assert all([e in [g.code for g in gcg.status] for e in case['errors']])


@pytest.mark.parametrize("i, case", enumerate(results))
def test_graph_str(i, case):
    """Verification the __repr__() method is not broken."""
    gcg = gc_graph(case['graph'])
    assert str(gcg)


@pytest.mark.parametrize("i, case", enumerate(results))
def test_graph_draw(i, case):
    """Verification the draw() method is not broken."""
    gcg = gc_graph(case['graph'])
    if case['valid']:
        gcg.draw(join(dirname(__file__), '../logs/gc_graph_' + str(i)))


@pytest.mark.parametrize("i, case", enumerate(results))
def test_graph_internal(i, case):
    """Verification initializing with an internal representation is self consdistent."""
    gcg = gc_graph(case['graph'])
    assert gcg.application_graph() == gc_graph(internal=deepcopy(gcg.save())).application_graph()
    assert gcg.graph == gc_graph(internal=deepcopy(gcg.save())).graph


@pytest.mark.parametrize("i, case", enumerate(results))
def test_graph_conversion(i, case):
    """Verification that converting to internal format and back again is the identity operation."""
    gcg = gc_graph(case['graph'])
    assert i == case['i']
    if case['valid']:
        for k, v in case['graph'].items():
            idx = const_idx.TYPE if k == 'C' else conn_idx.TYPE
            for r in v:
                r[idx] = r[idx]
        assert case['graph'] == gcg.application_graph()


@pytest.mark.parametrize("test", range(100))
def test_remove_connection_simple(test):
    """Verify adding connections makes valid graphs.

    Create a random graph remove some connections & re-normalise.
    To keep it simple all the endpoints have the same type ("int").
    """
    # TODO: These random test cases need to be made static when we are confident in them.
    # Generate them into a JSON file.
    graph = random_graph()
    assert graph.validate()

    # TODO: Split this out into its own test case when the graphs are staticly defined in a JSON file.
    for _ in range(int(len(list(filter(graph.dst_filter(), graph.graph.values()))) / 2)):
        graph.random_remove_connection()
    graph.normalize()
    assert graph.validate()
    # graph.draw(join(_log_location, 'graph_' + str(test)))


@pytest.mark.parametrize("test", range(100))
def test_add_connection(test):
    """Verify adding connections makes valid graphs.

    In this version multiple types endpoint types are used. This can lead to a legitimate invalid
    graph with error codes E01001 or E01004.
    """
    # TODO: These random test cases need to be made static when we are confident in them.
    # Generate them into a JSON file.
    gc = random_graph(0.5)
    if not gc.validate():
        codes = set([t.code for t in gc.status])
        codes.discard('E01001')
        codes.discard('E01004')
        assert not codes


@pytest.mark.parametrize("test", range(100))
def test_stack_simple(test):
    """Verify stacking valid graphs.

    Create two random graphs, gA & gB, and stack them.
    If gB has no inputs it cannot be stacked and the stacking method returns None.
    To keep it simple all the endpoints have the same type ("int"). This
    ensures all validation criteria will be met.
    """
    # TODO: These random test cases need to be made static when we are confident in them.
    # Generate them into a JSON file.
    global none_limit
    if not test:
        none_limit = 5000

    gA = random_graph()
    gB = random_graph()
    gC = gA.stack(gB)

    if gC is None:
        none_limit -= 1
    assert none_limit
    assert gC is None or gC.validate()
    # if not gC is None:
    #    print(gA)
    #    gA.draw('gA')
    #    print(gB)
    #    gB.draw('gB')
    #    print(gC)
    #    gC.draw('gC')
    #    barf()


@pytest.mark.parametrize("test", range(100))
def test_stack(test):
    """Verify stacking valid graphs.

    Create two random graphs, gA & gB, and stack them.
    If gB has no inputs it cannot be stacked and the stacking method returns None.
    In this version multiple types endpoint types are used. This can lead to a legitimate invalid
    stacked graphs which also return as None.
    """
    # TODO: These random test cases need to be made static when we are confident in them.
    # Generate them into a JSON file.
    global none_limit
    if not test:
        none_limit = 500

    gA = random_graph(0.5, True)
    gB = random_graph(0.5, True)
    gC = gA.stack(gB)

    if gC is None:
        none_limit -= 1
    assert none_limit
    assert gC is None or gC.validate()
    # if not gC is None:
    #    print(gA)
    #    gA.draw('gA')
    #    print(gB)
    #    gB.draw('gB')
    #    print(gC)
    #    gC.draw('gC')
    #    barf()
