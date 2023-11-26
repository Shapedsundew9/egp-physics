"""Test GC insertion operations.

See https://docs.google.com/spreadsheets/d/1YQjrM91e5x30VUIRzipNYX3W7yiFlg6fy9wKbMTx1iY/edit?usp=sharing
for insertion cases.

Simple validation is the following:
    1. TGC inputs == RGC inputs
    2. TGC outputs == RGC outputs
    3. IGC appears in the correct location in the RGC/FGC graph
    4. Row C & F are present/not present in the right places.
    5. FGC interfaces are as defined.
    6. gc_graph normalization
    7. gc_graph validation
It is NOT:
    1. Steady state exception validation

Special cases to check:
    a. IGC & TGC are the same instance

The preamble to the test cases codifies the rules in the google spreadsheet insertion case definitions.
The test cases are generated from the rules by means of a Cerberus validator using surebrec.
The validator is used to generate a set of sample graphs for each variant of TGC & IGC structure.
It is easierest to generate the variety needed from the internal graph representation and then convert
to the gc_graph representation.

Repeated generation is avoided by saving the sample graphs to a file and loading them in future runs.
"""
from copy import deepcopy
from itertools import count
from pprint import pformat
from typing import Any, Callable, cast
from random import choice
from itertools import product
import pytest
from logging import DEBUG, INFO, NullHandler, getLogger, Logger
from tqdm import tqdm
from json import load, dump
from os.path import exists, dirname, join
from itertools import combinations

from egp_types.dGC import dGC
from egp_types.ep_type import ep_type_lookup
from egp_types.egp_typing import Row, ROWS, VALID_GRAPH_ROW_COMBINATIONS, ConnectionGraph
from egp_types.gc_graph import gc_graph, random_gc_graph, SRC_EP, DST_EP
from egp_types.internal_graph import random_internal_graph, internal_graph, internal_graph_from_json
from egp_types.graph_validators import base_graph_validator
from egp_types.graph_validators import LIMITED_INTERNAL_GRAPH_SCHEMA, GRAPH_REGISTRY
from egp_physics.egp_typing import NewGCDef
from egp_physics.insertion import _insert_gc
from egp_physics import insertion


# Logging
_logger: Logger = getLogger(__name__)
_logger.addHandler(NullHandler())
_LOG_DEBUG: bool = _logger.isEnabledFor(DEBUG)
getLogger('surebrec').setLevel(INFO)


# Fake reference generator
_REFERENCE: count = count(1)


# Fake GMS with as much as we need to test the insertion.
# Note that a steady state exception should never be raised in this test.
class gms:
    """Fake GMS."""
    next_reference: Callable = lambda self: next(_REFERENCE)
_GMS = gms()


# Codification of the rules in the google spreadsheet insertion case definitions
# {Case: ((TGC must have, TGC must not have), (IGC must have, IGC must not have), above row options)}
# NOTE if F is required I is required.
_CASE_REQS: dict[int, tuple[tuple[str, str], tuple[str, str], str]] = {
    0: (("I", ""), ("O", ""), "I"),
    1: (("", "A"), ("", ""), "ABO"),
    2: (("A", "BF"), ("", ""), "A"),
    3: (("A", "BF"), ("", ""), "BO"),
    4: (("AB", "F"), ("", ""), "A"),
    5: (("AB", "F"), ("", ""), "B"),
    6: (("AB", "F"), ("", ""), "O"),
    7: (("IAF", ""), ("", ""), "A"),
    8: (("IAF", ""), ("", ""), "O"),
    9: (("IABF", ""), ("", ""), "B"),
    10: (("IABF", ""), ("", ""), "P"),
    11: (("O", ""), ("I", ""), "Z")
}


# Generation of all possible TGC & IGC variants for each case.
# {Case: [[TGC variants], [IGC variants], above row options]}
_CASE_VARIANTS: dict[int, list[Any]] = {}
for case, reqs in _CASE_REQS.items():
    for xgc_reqs in reqs[:2]:
        variables: str = "".join(set(cast(str, ROWS)) - set(xgc_reqs[0]) - set(xgc_reqs[1]))
        variants: list[str] = []
        for n in range(len(variables)):
            for s in combinations(variables, n):
                s_str: str = ''.join(sorted(xgc_reqs[0] + "".join(s)))
                # A graph with just an O row is guaranteed to be unstable 
                if s_str in VALID_GRAPH_ROW_COMBINATIONS and s_str != "O":
                    variants.append(s_str)
        _CASE_VARIANTS.setdefault(case, []).append(variants)
    _CASE_VARIANTS[case].append(reqs[2])


# Generate all combinations of TGC & IGC for each case
_CASE_COMBOS: list[tuple[Any, ...]] = [(case, *combo) for case, variants in _CASE_VARIANTS.items() for combo in product(*variants)]
_logger.debug(f"Case combos:\n{pformat(_CASE_COMBOS)}")


# Generate a set of sample graphs for each variant if they have not already been generated and saved to a file.
filename: str = join(dirname(__file__), "data", "test_insertion.json")
if not exists(filename):
    # 1. Generate 100 random src & dst endpoints in the row required
    # 2. Create variants by randomly choosing a subset of required row end points

    # Create NUM_SAMPLES sample graphs for each variant.
    NUM_SAMPLES: int = 10
    _VARIANT_SAMPLES: dict[str, list[gc_graph]] = {}
    for rseed, variant in enumerate(tqdm(VALID_GRAPH_ROW_COMBINATIONS, desc="Generating graphs")):
        gc_graph_list: list[gc_graph] = [
            gc_graph(i_graph=random_internal_graph(variant, (ep_type_lookup["n2v"]["int"],), verify=True, rseed=rseed * NUM_SAMPLES + i, row_stablization=True)) for i in range(NUM_SAMPLES)]
        for gcg in gc_graph_list:
            gcg.normalize()
            assert gcg.validate()
        _VARIANT_SAMPLES[variant] = gc_graph_list

    # Write the samples to a file so we can just load them in future
    _logger.debug(f"Variants created: {sorted(_VARIANT_SAMPLES.keys())}")
    i_graph_json = {k: [v.i_graph.json_obj() for v in l] for k, l in _VARIANT_SAMPLES.items()}
    _logger.debug(f"Variants stored: {sorted(i_graph_json.keys())}")
    with open(filename, "w", encoding="ascii") as f:
        dump(i_graph_json, f, indent=4, sort_keys=True)


with open(filename, "r", encoding="ascii") as f:
    _VARIANT_SAMPLES = {k: [gc_graph(i_graph=internal_graph_from_json(v)) for v in l] for k, l in load(f).items()}

# O row only graphs are guaranteed to be unstable
del _VARIANT_SAMPLES["O"]
_VARIANT_VALID_SUPERSET = tuple(_VARIANT_SAMPLES.keys())
_logger.debug(f"Variants loaded: {sorted(_VARIANT_SAMPLES.keys())}")


#_logger.debug(pformat(_VARIANT_SAMPLES))
# Codification of the RGC and FGC graphs after insertion case.
# {Case: ({RGC Row: (Source xGC, Row) | None}} | None, {FGC Row: (Source xGC, Row) | None}} | None)
# None = Must not exist, Empty tuple = Any
_CASE_RESULTS: dict[int, tuple[dict[Row, tuple[str, ...] | None] | None, ...]] = {
    0: (
        {
            "I": ("IGC", "I"),
            "C": None,
            "F": None,
            "A": ("IGC", "IO"),
            "B": ("TGC", "IO"),
            "O": ("TGC", "O"),
        },
        None
    ),
    1: (
        {
            "I": ("TGC", "I"),
            "C": ("TGC", "C"),
            "F": None,
            "A": ("IGC", "IO"),
            "B": None,
            "O": ("TGC", "O"),
        },
        None
    ),
    2: (
        {
            "I": ("TGC", "I"),
            "C": ("TGC", "C"),
            "F": None,
            "A": ("IGC", "IO"),
            "B": ("TGC", "A"),
            "O": ("TGC", "O"),
        },
        None
    ),
    3: (
        {
            "I": ("TGC", "I"),
            "C": ("TGC", "C"),
            "F": None,
            "A": ("TGC", "A"),
            "B": ("IGC", "IO"),
            "O": ("TGC", "O"),
        },
        None
    ),
    4: (
        {
            "I": ("TGC", "I"),
            "C": ("TGC", "C"),
            "F": None,
            "A": ("FGC", "IO"),
            "B": ("TGC", "B"),
            "O": ("TGC", "O"),
        },
        {
            "I": ("TGC", "I"),
            "C": ("TGC", "C"),
            "F": None,
            "A": ("IGC", "IO"),
            "B": ("TGC", "A"),
            "O": tuple(),
        }
    ),
    5: (
        {
            "I": ("TGC", "I"),
            "C": ("TGC", "C"),
            "F": None,
            "A": ("FGC", "IO"),
            "B": ("TGC", "B"),
            "O": ("TGC", "O"),
        },
        {
            "I": ("TGC", "I"),
            "C": ("TGC", "C"),
            "F": None,
            "A": ("TGC", "A"),
            "B": ("IGC", "IO"),
            "O": tuple(),
        }
    ),
    6: (
        {
            "I": ("TGC", "I"),
            "C": ("TGC", "C"),
            "F": None,
            "A": ("FGC", "IO"),
            "B": ("IGC", "IO"),
            "O": ("TGC", "O"),
        },
        {
            "I": ("TGC", "I"),
            "C": ("TGC", "C"),
            "F": None,
            "A": ("TGC", "A"),
            "B": ("TGC", "B"),
            "O": tuple(),
        }
    ),
    7: (
        {
            "I": ("TGC", "I"),
            "C": ("TGC", "C"),
            "F": ("TGC", "F"),
            "A": ("FGC", "IO"),
            "B": ("TGC", "B"),
            "O": ("TGC", "O"),
        },
        {
            "I": ("TGC", "A"),
            "C": None,
            "F": None,
            "A": ("IGC", "IO"),
            "B": ("TGC", "A"),
            "O": ("TGC", "A"),
        }
    ),
    8: (
        {
            "I": ("TGC", "I"),
            "C": ("TGC", "C"),
            "F": ("TGC", "F"),
            "A": ("FGC", "IO"),
            "B": ("TGC", "B"),
            "O": ("TGC", "O"),
        },
        {
            "I": ("TGC", "A"),
            "C": None,
            "F": None,
            "A": ("TGC", "A"),
            "B": ("IGC", "IO"),
            "O": ("TGC", "A"),
        }
    ),
    9: (
        {
            "I": ("TGC", "I"),
            "C": ("TGC", "C"),
            "F": ("TGC", "F"),
            "A": ("TGC", "A"),
            "B": ("FGC", "IO"),
            "O": ("TGC", "O"),
        },
        {
            "I": ("TGC", "B"),
            "C": None,
            "F": None,
            "A": ("IGC", "IO"),
            "B": ("TGC", "B"),
            "O": ("TGC", "B"),
        }
    ),
    10: (
        {
            "I": ("TGC", "I"),
            "C": ("TGC", "C"),
            "F": ("TGC", "F"),
            "A": ("TGC", "A"),
            "B": ("FGC", "IO"),
            "O": ("TGC", "O"),
        },
        {
            "I": ("TGC", "B"),
            "C": None,
            "F": None,
            "A": ("TGC", "B"),
            "B": ("IGC", "IO"),
            "O": ("TGC", "B"),
        }
    ),
    11: (
        {
            "I": ("TGC", "I"),
            "C": None,
            "F": None,
            "A": ("TGC", "IO"),
            "B": ("IGC", "IO"),
            "O": ("IGC", "O"),
        },
        None
    ),
}


def new_dgc(graph: gc_graph) -> dGC:
    """Return a new dGC with a new reference."""
    return {
        'gc_graph': graph,
        'ref': next(_REFERENCE),
        'gca_ref': next(_REFERENCE),
        'gcb_ref': next(_REFERENCE),
        'ancestor_a_ref': next(_REFERENCE),
        'ancestor_b_ref': next(_REFERENCE)
    }


def new_xgc(xgc_variant: str) -> dGC:
    """Return an t or i gc with a new reference."""
    gca: dGC = new_dgc(choice(_VARIANT_SAMPLES[choice(_VARIANT_VALID_SUPERSET)]))
    gcb: dGC = new_dgc(choice(_VARIANT_SAMPLES[choice(_VARIANT_VALID_SUPERSET)]))
    return {
        'gc_graph': choice(_VARIANT_SAMPLES[xgc_variant]),
        'ref': next(_REFERENCE),
        'gca_ref': gca['ref'],
        'gcb_ref': gcb['ref'],
        'ancestor_a_ref': gca['ref'],
        'ancestor_b_ref': gcb['ref']
    }


# Mock interface proximity select so that any steady state exception
# returns a GC that will not raise further exceptions.
GCG = gc_graph(cast(ConnectionGraph, {"O": [["A", 0, ep_type_lookup["n2v"]["int"]]]}))
def mock_ips(_, __) -> dGC:
    """Mock interface proximity select so that any steady state exception returns a GC that will not raise further exceptions."""
    return new_dgc(GCG)


# Mock steady state exception: Should not happen in this test.
def mock_sse(_, __) -> dGC:
    """Mock steady state exception: Should not happen in this test."""
    assert False, "Steady State Exception should not happen!"


@pytest.mark.parametrize("i_case, tgc_variant, igc_variant, above_row", _CASE_COMBOS)
def test_gc_insert(monkeypatch, i_case, tgc_variant, igc_variant, above_row) -> None:
    """Test all insertion cases for all combinations of IGC & TGC structures."""
    monkeypatch.setattr(insertion, "steady_state_exception", mock_sse)

    tgc: dGC = new_xgc(tgc_variant)
    igc: dGC = new_xgc(igc_variant)

    # The FGC reference must be greater than this else it is not an FGC...
    new_base_ref: int = next(_REFERENCE)

    if _LOG_DEBUG:
        _logger.debug(f"Case: {i_case}, TGC variant: {tgc_variant}, IGC variant: {igc_variant}, Above row: {above_row}")
        _logger.debug(f"TGC:\n {pformat(tgc)}")
        _logger.debug(f"IGC:\n {pformat(igc)}")

    # Insert the IGC into the TGC
    new_gc_definition: NewGCDef = _insert_gc(_GMS, tgc, igc, above_row)  # type: ignore

    # Check the new gc definition is valid
    # 1. TGC inputs == RGC inputs
    # 2. TGC outputs == RGC outputs
    # 3. IGC appears in the correct location in the RGC/FGC graph
    # 4. Row C & F are present/not present in the right places.
    # 5. FGC interfaces are as defined.
    # 6. gc_graph normalization
    # 7. gc_graph validation
    expected_rgc, expected_fgc = _CASE_RESULTS[i_case]
    rgc, fgc_dict = new_gc_definition
    fgc: dGC | None = cast(dGC, tuple(fgc_dict.values())[0]) if fgc_dict else None

    if expected_rgc is not None:
        assert rgc["ref"] > new_base_ref
        for rgc_row, expected_row in expected_rgc.items():
            # If the expected row is None then it must not exist and gc_grpah.has_x is False
            if expected_row is None:
                assert not rgc['gc_graph'].has_row(rgc_row)
            else:
                # The row may exist & if it does it must match the expected row
                if rgc['gc_graph'].has_row(rgc_row):
                    assert fgc is not None
                    match expected_row[0]:
                        case "TGC": xgc: dGC | None = tgc
                        case "IGC": xgc = igc
                        case "FGC": xgc = fgc
                        case _: raise ValueError(f"Unexpected source {expected_row[0]}")
                    if expected_row[1] != "IO":
                        assert rgc["gc_graph"].row_if(rgc_row) == xgc["gc_graph"].row_if(cast(Row, expected_row[1])), f"Expectations: {expected_row}"
                    else:
                        assert rgc["gc_graph"].row_if(rgc_row) == (
                            xgc["gc_graph"].output_if(), xgc["gc_graph"].input_if()), f"Expectations for row {rgc_row}: {expected_row}"

