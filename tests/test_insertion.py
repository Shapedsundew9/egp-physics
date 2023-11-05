"""Test GC insertion operations.

See https://docs.google.com/spreadsheets/d/1YQjrM91e5x30VUIRzipNYX3W7yiFlg6fy9wKbMTx1iY/edit?usp=sharing

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

from egp_types.dGC import dGC
from egp_types.egp_typing import Row, connection_graph_to_json, JSONGraph
from egp_types.gc_graph import gc_graph, random_gc_graph, SRC_EP, DST_EP
from egp_utils.base_validator import base_validator
from egp_types.graph_validators import LIMITED_INTERNAL_GRAPH_SCHEMA, GR
from egp_physics.egp_typing import NewGCDef
from egp_physics.insertion import _insert_gc


# Logging
_logger: Logger = getLogger(__name__)
_logger.addHandler(NullHandler())
_LOG_DEBUG: bool = _logger.isEnabledFor(DEBUG)
getLogger('surebrec').setLevel(INFO)


_REFERENCE: count = count(1)
_ROWS: str = "ICFABO"


# Fake GMS with as much as we need to test the insertion.
# Note that a steady state exception should never be raised in this test.
class gms:
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
        variables: str = "".join(set(_ROWS) - set(xgc_reqs[0]) - set(xgc_reqs[1]))
        variants: list[str] = ["".join(sorted(xgc_reqs[0] + variables[0:v])) for v in range(len(variables) + 1)]
        _CASE_VARIANTS.setdefault(case, []).append(variants)
    _CASE_VARIANTS[case].append(reqs[2])


# Somes cases are invalid
for case, variant_pair in _CASE_VARIANTS.items():
    for variants in variant_pair[:2]:
        for variant in deepcopy(variants):
            if "B" in variant and "A" not in variant:
                variants.remove(variant)
            elif "F" in variant and "O" in variant  and "P" not in variant:
                variants.remove(variant)
            elif "F" in variant and "O" not in variant and "P" in variant:
                variants.remove(variant)
            elif variant in ("O", "P", "Z"):
                variants.remove(variant)


# Generate all combinations of TGC & IGC for each case
_CASE_COMBOS: list[tuple[Any, ...]] = [(case, *combo) for case, variants in _CASE_VARIANTS.items() for combo in product(*variants)]


# Superset of all the possible xGC graph row variants.
_VARIANT_SUPERSET_SET: set[str] = set()
for variants in _CASE_VARIANTS.values():
    for variant in variants[:2]:
        _VARIANT_SUPERSET_SET |= set(variant)
_VARIANT_SUPERSET: tuple[str, ...] = tuple(_VARIANT_SUPERSET_SET)






# Restrict the graph schema to only the integer type and short interface lengths.
_SIMPLE_GRAPH_SCHEMA: dict[str, dict[str, Any]] = deepcopy(LIMITED_INTERNAL_GRAPH_SCHEMA)
for variant in _SIMPLE_GRAPH_SCHEMA["internal_graph"]["valuesrules"]["oneof"]:
    # Don't need to create anything in row U it will be populated, as needed, by normalization.
    for row in variant.values():
        # Max 4 destinations on a row
        if row.get("maxlength", 1) == 256:
            row["maxlength"] = 16
        # All end point types set to 2
        if "ep_type" in row["schema"]["items"]:
            row["schema"]["items"][row["schema"]["items"].index("ep_type")] = {"type": "integer", "allowed": [2]}
validator = base_validator(_SIMPLE_GRAPH_SCHEMA)


# Create 10 sample graphs for each variant.
_DUMMY_LIST: list[int] = [0] * 10
_VARIANT_SAMPLES: dict[str, list[JSONGraph]] = {v: [] for v in _VARIANT_SUPERSET}
samples: int = len(_VARIANT_SUPERSET_SET) * 10
with tqdm(total=samples) as pbar:
    while samples:
        new_graph: gc_graph = random_gc_graph(validator, True)
        rows = set(new_graph.rows[SRC_EP])
        rows.update(new_graph.rows[DST_EP])
        rows.discard("U")
        graph_variant: str = ''.join(sorted(rows))
        if len(_VARIANT_SAMPLES.get(graph_variant, _DUMMY_LIST)) < 10:
            _VARIANT_SAMPLES[graph_variant].append(connection_graph_to_json(new_graph.connection_graph()))
            pbar.update(1)
            samples -= 1
            if samples < 200:
                with open("tests/data/test_insertion.json", "w", encoding="utf-8") as f:
                    dump(_VARIANT_SAMPLES, f, indent=4, sort_keys=True)

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
    gca: dGC = new_dgc(choice(_VARIANT_SAMPLES[choice(_VARIANT_SUPERSET)]))
    gcb: dGC = new_dgc(choice(_VARIANT_SAMPLES[choice(_VARIANT_SUPERSET)]))
    return {
        'gc_graph': choice(_VARIANT_SAMPLES[xgc_variant]),
        'ref': next(_REFERENCE),
        'gca_ref': gca['ref'],
        'gcb_ref': gcb['ref'],
        'ancestor_a_ref': gca['ref'],
        'ancestor_b_ref': gcb['ref']
    }


@pytest.mark.parametrize("i_case, tgc_variant, igc_variant, above_row", _CASE_COMBOS)
def test_gc_insert(i_case, tgc_variant, igc_variant, above_row) -> None:
    """Test all insertion cases for all combinations of IGC & TGC structures."""
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
                assert not getattr(rgc['gc_graph'], 'has_' + rgc_row)
            else:
                # The row may exist & if it does it must match the expected row
                if getattr(rgc['gc_graph'], 'has_' + rgc_row):
                    assert fgc is not None
                    match expected_row[0]:
                        case "TGC": xgc: dGC | None = tgc
                        case "IGC": xgc = igc
                        case "FGC": xgc = fgc
                        case _: raise ValueError(f"Unexpected source {expected_row[0]}")
                    if expected_row[1]!= "IO":
                        assert rgc["gc_graph"].row_if(rgc_row) == xgc["gc_graph"].row_if(cast(Row, expected_row[1]))
                    else:
                        assert rgc["gc_graph"].row_if(rgc_row) == (
                            xgc["gc_graph"].row_if("I"),
                            xgc["gc_graph"].row_if("O")
                        )
