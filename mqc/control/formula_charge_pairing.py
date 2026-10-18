"""Conservative completion of a missing charge using a matched reference pair.

No bounds, stoichiometry, rule weights or gene evidence are changed here.
"""
import math
import re

import cobra
import libsbml
from cobra.core.formula import elements_and_molecular_weights
from cobra.io.sbml import F_REPLACE, F_SPECIE


def finite_number(value):
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def missing_input_charges(file_path, model):
    """Restore explicit legacy charges before classifying true missing values.

    COBRA may return FBC's default zero even when a legacy notes charge is
    present. FBC attributes take precedence; explicit zero is not missing.
    """
    missing = {m.id for m in model.metabolites if not finite_number(m.charge)}
    if str(file_path).lower().endswith(('.json', '.yaml', '.yml', '.mat')):
        return missing  # Preserve explicitly loaded zero in non-SBML formats.
    document = libsbml.readSBMLFromFile(str(file_path))
    sbml_model = document.getModel()
    if sbml_model is None:
        raise ValueError('Cannot establish input charge provenance from SBML')
    for species in sbml_model.getListOfSpecies():
        mid = F_REPLACE[F_SPECIE](species.getId())
        if mid not in model.metabolites:
            continue
        plugin = species.getPlugin('fbc')
        explicit = bool(plugin is not None and plugin.isSetCharge())
        metabolite = model.metabolites.get_by_id(mid)
        if not explicit and species.isSetCharge():
            metabolite.charge = species.getCharge()
            explicit = True
        note_charge = metabolite.notes.get('CHARGE')
        if not explicit and finite_number(note_charge) and float(note_charge).is_integer():
            metabolite.charge = int(float(note_charge))
            explicit = True
        if explicit:
            missing.discard(mid)
        if not explicit:
            missing.add(mid)
    return missing


def complete_reference_pair(metabolite, formula, charge):
    """Called only for missing input charge. Accept exact or H-only formula differences.

    References with unknown/polymeric formulas, fractional charges/counts, or
    non-H atom disagreement are not silently merged with the input identity.
    """
    if not finite_number(charge) or not float(charge).is_integer():
        return False
    formulas = [metabolite.formula, formula]
    elements = []
    for value in formulas:
        if not isinstance(value, str) or not re.fullmatch(r'(?:[A-Z][a-z]?\d*)+', value):
            return False
        counts = cobra.Metabolite(formula=value).elements
        if not counts or any(k not in elements_and_molecular_weights or not finite_number(v)
                             or v <= 0 or not float(v).is_integer() for k, v in counts.items()):
            return False
        elements.append({k: v for k, v in counts.items() if k != 'H'})
    if elements[0] != elements[1]:
        return False
    metabolite.formula = formula
    metabolite.charge = int(charge)
    return True
