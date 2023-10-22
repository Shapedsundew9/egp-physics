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

from egp_types.dGC import dGC
from egp_types.egp_typing import Row
from egp_types.gc_graph import gc_graph, random_gc_graph
from egp_types.xgc_validator import GRAPH_SCHEMA, base_validator
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
# {Case: [(TGC variants), (IGC variants), above row options]}
_CASE_VARIANTS: dict[int, list[Any]] = {}
for case, reqs in _CASE_REQS.items():
    for xgc_reqs in reqs[:2]:
        variables: str = "".join(set(_ROWS) - set(xgc_reqs[0]) - set(xgc_reqs[1]))
        variants = tuple("".join(sorted(xgc_reqs[0] + variables[0:v])) for v in range(len(variables) + 1))
        _CASE_VARIANTS.setdefault(case, []).append(variants)
    _CASE_VARIANTS[case].append(reqs[2])


# Generate all combinations of TGC & IGC for each case
_CASE_COMBOS: list[tuple[Any, ...]] = [(case, *combo) for case, variants in _CASE_VARIANTS.items() for combo in product(*variants)]


# Superset of all the possible xGC graph row variants.
_VARIANT_SUPERSET_SET: set[str] = set()
for variants in _CASE_VARIANTS.values():
    for variant in variants:
        _VARIANT_SUPERSET_SET |= set(variant)
_VARIANT_SUPERSET: tuple[str, ...] = tuple(_VARIANT_SUPERSET_SET)


# Restrict the graph schema to only the integer type and short interface lengths.
_SIMPLE_GRAPH_SCHEMA: dict[str, dict[str, Any]] = deepcopy(GRAPH_SCHEMA)
for variant in _SIMPLE_GRAPH_SCHEMA["graph"]["oneof_schema"]:
    # Don't need to create anything in row U it will be populated, as needed, by normalization.
    for row in variant.values():
        # Max 4 destinations on a row
        if row.get("maxlength", 1) == 256:
            row["maxlength"] = 4
        # All end point types set to 2
        if "ep_type" in row["schema"]["items"]:
            row["schema"]["items"][row["schema"]["items"].index("ep_type")] = {"type": "integer", "allowed": [2]}


# Creation of a graph schema for each variant and the associated validator.
# {Variant: (Schema, Validator)}
_VARIANT_VALIDATORS: dict[str, tuple[dict[str, Any], Any]] = {}
for variant in _VARIANT_SUPERSET:
    schema: dict[str, Any] = {"graph": {
            "type": "dict",
            "required": True,
            "schema": deepcopy(_SIMPLE_GRAPH_SCHEMA["graph"]["oneof_schema"]["F" in variant])
        }
    }

    # Remove rows that are not in the variant
    for row in _ROWS:
        if row not in variant and row in schema["graph"]["schema"] and row not in "PU":
            del schema["graph"]["schema"][row]
            if row == "O" and "P" in schema["graph"]["schema"]:
                del schema["graph"]["schema"]["P"]
            if row == "A" and "B" in schema["graph"]["schema"]:
                del schema["graph"]["schema"]["B"]

    # Make sure allowed rows do not include any rows that have been removed
    for row, rdef in deepcopy(schema["graph"]["schema"]).items():
        if "allowed" in rdef["schema"]["items"][0]:
            rdef["schema"]["items"][0]["allowed"] = list(set(rdef["schema"]["items"][0]["allowed"]) - (set(_ROWS) - set(variant)))
            # There are no allowed row reference options. So row must be empty.
            if not rdef["schema"]["items"][0]["allowed"]:
                del schema["graph"]["schema"][row]

    # Row U is superflous unless only row I exists.
    if variant != "I" and "U" in schema["graph"]["schema"]:
        del schema["graph"]["schema"]["U"]

    validator = base_validator(schema)
    _VARIANT_VALIDATORS[variant] = (schema, validator)

print(pformat(_VARIANT_VALIDATORS))

# Create 10 sample graphs for each variant.
_VARIANT_SAMPLES: dict[str, tuple[gc_graph, ...]] = {
    variant: tuple(random_gc_graph(validator, True) for _ in range(10))
    for variant, (_, validator) in _VARIANT_VALIDATORS.items()
}


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
