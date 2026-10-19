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


def zero_charge_placeholders(model):
    """Recognize an all-zero charge column only with three independent ions.

    Explicit neutral charges in an otherwise charged model remain authoritative.
    Recognition does not assign charges: each repair still requires a matched
    reference formula/charge pair. BiGG, SEED and exact ion names are supported.
    """
    if not model.metabolites or not all(finite_number(m.charge) and float(m.charge) == 0
                                       for m in model.metabolites):
        return set()
    anchors = {
        'H': ({'h', 'cpd00067'}, {'h+', 'proton'}),
        'Na': ({'na1', 'cpd00971'}, {'sodium', 'na+'}),
        'K': ({'k', 'cpd00205'}, {'potassium', 'k+'}),
        'Cl': ({'cl', 'cpd00099'}, {'chloride', 'cl-'}),
        'Mg': ({'mg2', 'cpd00254'}, {'magnesium', 'mg2+'}),
        'Ca': ({'ca2', 'cpd00063'}, {'calcium', 'ca2+'}),
    }
    evidence = set()
    for met in model.metabolites:
        base = re.sub(r'(?:_[a-z]\d*|\[[a-z]\d*\]|@[a-z]\d*)$', '', met.id.lower())
        if met.formula in anchors:
            ids, names = anchors[met.formula]
            if base in ids or (met.name or '').strip().lower() in names:
                evidence.add(met.formula)
    return set(model.metabolites.list_attr('id')) if len(evidence) >= 3 else set()


def missing_input_charges(file_path, model):
    """Restore explicit legacy charges before classifying true missing values.

    COBRA may return FBC's default zero even when a legacy notes charge is
    present. FBC attributes take precedence; explicit zero is not missing.
    """
    missing = {m.id for m in model.metabolites if not finite_number(m.charge)}
    if str(file_path).lower().endswith(('.json', '.yaml', '.yml', '.mat')):
        placeholders = zero_charge_placeholders(model)
        for mid in placeholders:
            model.metabolites.get_by_id(mid).notes['MQC_ZERO_CHARGE_PLACEHOLDER'] = 'all_zero_column_with_three_ion_anchors'
        return missing | placeholders
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
    placeholders = zero_charge_placeholders(model)
    for mid in placeholders:
        model.metabolites.get_by_id(mid).notes['MQC_ZERO_CHARGE_PLACEHOLDER'] = 'all_zero_column_with_three_ion_anchors'
    return missing | placeholders


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
    if metabolite.notes.get('MQC_ZERO_CHARGE_PLACEHOLDER'):
        # A suspect explicit zero alone is insufficient to change protonation.
        if cobra.Metabolite(formula=formulas[0]).elements != cobra.Metabolite(formula=formulas[1]).elements:
            return False
        metabolite.notes['MQC_CHARGE_REPAIR'] = 'all-zero input; exact reference formula; charge 0 -> %s' % int(charge)
    metabolite.formula = formula
    metabolite.charge = int(charge)
    return True
