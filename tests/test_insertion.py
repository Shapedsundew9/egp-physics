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
from pprint import pprint
from typing import Any

from egp_types.dGC import dGC
from egp_types.egp_typing import ConnectionGraph
from egp_types.gc_graph import gc_graph, random_gc_graph
from egp_types.xgc_validator import GRAPH_SCHEMA, base_validator


_REFERENCE: count = count(1)
_ROWS: str = "ICFABO"


# Codification of the rules in the google spreadsheet insertion case definitions
# {Case: ((TGC must have, TGC must not have), (IGC must have, IGC must not have))}
_CASE_REQS: dict[int, tuple[tuple[str, str], tuple[str, str]]] = {
    0: (("I", ""), ("O", "")),
    1: (("", "A"), ("", "")),
    2: (("A", "BF"), ("", "")),
    3: (("A", "BF"), ("", "")),
    4: (("AB", "F"), ("", "")),
    5: (("AB", "F"), ("", "")),
    6: (("AB", "F"), ("", "")),
    7: (("AF", ""), ("", "")),
    8: (("AF", ""), ("", "")),
    9: (("ABF", ""), ("", "")),
    10: (("ABF", ""), ("", "")),
    11: (("O", ""), ("I", ""))
}

# Generation of all possible TGC & IGC variants for each case.
# {Case: [(TGC variants), (IGC variants)]}
_CASE_VARIANTS: dict[int, list[tuple[str, ...]]] = {}
for case, reqs in _CASE_REQS.items():
    for xgc_reqs in reqs:
        variables: str = "".join(set(_ROWS) - set(xgc_reqs[0]) - set(xgc_reqs[1]))
        variants = tuple("".join(sorted(xgc_reqs[0] + variables[0:v])) for v in range(len(variables) + 1))
        _CASE_VARIANTS.setdefault(case, []).append(variants)


# Superset of all the possible xGC graph row variants.
_VARIANT_SUPERSET: set[str] = set()
for variants in _CASE_VARIANTS.values():
    for variant in variants:
        _VARIANT_SUPERSET |= set(variant)


# Restrict the graph schema to only the integer type and short interface lengths.
_SIMPLE_GRAPH_SCHEMA: dict[str, dict[str, Any]] = deepcopy(GRAPH_SCHEMA)
for variant in _SIMPLE_GRAPH_SCHEMA["graph"]["oneof_schema"]:
    # Don't need to create anything in row U it will be populated, as needed, by normalization.
    if "U" in variant:
        del variant["U"]
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
    for row in _ROWS:
        if row not in variant and row in schema["graph"] and row != "P":
            del schema["graph"][row]

    # Row I is implicitly defined from the other rows references and so
    # has to be removed as an option if it should not be one (but is required if row F is present)
    if "I" not in variant and "F" not in variant:
        for definition in schema["graph"]["schema"].values():
            if "I" in definition["schema"]["items"][0].get("allowed", {}):
                definition["schema"]["items"][0]["allowed"].remove("I")
    validator = base_validator(schema)
    _VARIANT_VALIDATORS[variant] = (schema, validator)


# Create 10 sample graphs for each variant.
_VARIANT_SAMPLES: dict[str, tuple[gc_graph, ...]] = {
    variant: tuple(random_gc_graph(validator, True) for _ in range(10))
    for variant, (_, validator) in _VARIANT_VALIDATORS.items()
}


def new_dgc() -> dGC:
    """Return a new dGC with a new reference."""
    return {
        'gc_graph': gc_graph(_xGCCG),
        'ref': next(_REFERENCE),
        'gca_ref': next(_REFERENCE),
        'gcb_ref': next(_REFERENCE),
        'ancestor_a_ref': next(_REFERENCE),
        'ancestor_b_ref': next(_REFERENCE)
    }


def new_xgc(variant: str) -> dGC:
    """Return an t or i gc with a new reference."""
    gca: dGC = new_dgc()
    gcb: dGC = new_dgc()
    return {
        'gc_graph': gc_graph(_xGCCG),
        'ref': next(_REFERENCE),
        'gca_ref': gca['ref'],
        'gcb_ref': gcb['ref'],
        'ancestor_a_ref': gca['ref'],
        'ancestor_b_ref': gcb['ref']
    }

_TGC: dGC = new_xgc()
_IGC: dGC = new_xgc()


def test_gc_insert():
    """Test gc_insert."""
