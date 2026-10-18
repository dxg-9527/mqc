"""Interpret the existing MetaCyc table as directions, not numeric capacities.

No bounds are changed here. Rules scores allowed directions; Rules2 tests the
actual offending flux used by MQC's existing iterative closure procedure.
Unresolved reaction orientation must not be guessed from equation order.
"""
from functools import lru_cache
from pathlib import Path
import math
import re

import pandas as pd
from mqc.defaults import ATP_NAME, NADP, NADPH, NAD_NAME, NADH_NAME

RULE = "refer to the Metacyc database to modify the reversibility of this reaction"
# Exclude ambiguous hydroxide/oxonium synonyms from the legacy water list.
ANCHORS = {
    "ATP": frozenset(ATP_NAME), "NADP": frozenset(NADP),
    "NADPH": frozenset(NADPH), "NAD": frozenset(NAD_NAME),
    "NADH": frozenset(NADH_NAME),
    "WATER": frozenset(["h2o", "H2O", "H2O H2O", "Water", "WATER", "water",
                         "H2O-c0", "H2O-e0", "H2O_c0", "H2O_e0", "H2O[cytosol]"]),
}


def directions(bounds):
    """Capacity magnitudes do not determine thermodynamic direction."""
    lower, upper = (float(v) for v in bounds)
    if not math.isfinite(lower) or not math.isfinite(upper) or lower > upper:
        raise ValueError("Invalid reaction bounds")
    return frozenset(s for s, allowed in [(-1, lower < 0), (1, upper > 0)] if allowed)


def base_identifier(value):
    value = str(value).strip()
    if value.startswith("META:"):
        value = value[5:]
    # BiGG-style variants append compartment or explicit compound mappings.
    # RXN-1 must never match RXN-10 simply because it is a substring.
    value = value.split("[", 1)[0]
    if "/" in value:
        match = re.match(r"^(.+?-RXN|RXN-\d+)-", value)
        if match:
            value = match.group(1)
    return value


def compound_tokens(side):
    tokens = set()
    for term in re.split(r"\s+\+\s+", side.strip()):
        term = re.sub(r"^[-+]?\d+(?:\.\d+)?\s+", "", term.strip())
        # Compartment suffix, not a substring search for ATP within dATP.
        term = re.sub(r"_[a-z][a-z0-9]*$", "", term)
        if term:
            tokens.add(term)
    return tokens


def reference_direction(rxn, lower, upper, equation):
    allowed = directions((lower, upper))
    if len(allowed) != 1:
        return allowed, "orientation_independent", ()
    sides = str(equation).split("-->")
    if len(sides) != 2:
        return None, "unparsed_equation", ()
    left, right = map(compound_tokens, sides)
    model_left = set(rxn.get("reactants_mets", []))
    model_right = set(rxn.get("products_mets", []))
    evidence = []
    for compound, aliases in ANCHORS.items():
        model_sides = [s for s, names in [(-1, model_left), (1, model_right)] if aliases & names]
        ref_sides = [s for s, names in [(-1, left), (1, right)] if compound in names]
        if len(model_sides) == len(ref_sides) == 1:
            evidence.append((compound, model_sides[0] * ref_sides[0]))
    orientations = {sign for _, sign in evidence}
    if len(orientations) != 1:
        return None, "conflicting_anchors" if orientations else "no_orientation_anchor", tuple(evidence)
    orientation = next(iter(orientations))
    return frozenset(orientation * sign for sign in allowed), "oriented", tuple(evidence)


@lru_cache(maxsize=4)
def _read_index(path, mtime_ns, size):
    index = {}
    for row in pd.read_excel(path).iloc[:, :4].itertuples(index=False, name=None):
        if pd.isna(row[0]):
            continue
        index.setdefault(base_identifier(row[0]), []).append(tuple(row))
    return {key: tuple(rows) for key, rows in index.items()}


def reference_index(path):
    path = Path(path).resolve()
    stat = path.stat()
    return _read_index(str(path), stat.st_mtime_ns, stat.st_size)


def resolve(rxn, index):
    annotations = rxn.get("annotation", {})
    annotation = annotations.get("biocyc", []) if isinstance(annotations, dict) else []
    values = [annotation] if isinstance(annotation, str) else annotation
    if not isinstance(values, (list, tuple)):
        values = []
    matched = []
    for value in dict.fromkeys(v for v in values if isinstance(v, str)):
        matched.extend(index.get(base_identifier(value), ()))
    if not matched:
        return dict(allowed=None, reason="no_reference", references=[])
    records = []
    for identifier, lower, upper, equation in dict.fromkeys(matched):
        try:
            signs, reason, anchors = reference_direction(rxn, lower, upper, equation)
        except (ValueError, TypeError):
            signs, reason, anchors = None, "invalid_reference_bounds", ()
        records.append(dict(id=identifier, allowed=None if signs is None else sorted(signs),
                            reason=reason, anchors=anchors))
    options = {None if r["allowed"] is None else tuple(r["allowed"]) for r in records}
    if len(options) != 1 or None in options:
        return dict(allowed=None, reason="unresolved_or_disagreeing_reference", references=records)
    return dict(allowed=list(next(iter(options))), reason="resolved", references=records)


def violates(rxn, resolution, flux=None, tolerance=1e-9):
    if resolution["allowed"] is None:
        return False
    allowed = set(resolution["allowed"])
    if flux is None:
        return not directions(rxn["bounds"]).issubset(allowed)
    if not math.isfinite(float(flux)):
        raise ValueError("Non-finite flux is not a valid direction test")
    return abs(flux) > tolerance and (1 if flux > 0 else -1) not in allowed


def apply_rule(model_info, path, reaction_id=None, flux=None):
    index = reference_index(path)
    for rxn in model_info["reactions"]:
        if reaction_id is not None and rxn["id"] != reaction_id:
            continue
        if violates(rxn, resolve(rxn, index), flux):
            key = "rules" if reaction_id is None else "rules2"
            rxn.setdefault(key, {})[RULE] = "true"
            if reaction_id is not None:
                rxn["rule_conflict"] = "true"
