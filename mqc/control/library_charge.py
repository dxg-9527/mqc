"""Restore exact formula/charge pairs in a known all-zero reference artifact.

This is a loader erratum for the distributed BiGG-style universal library,
whose 24,403 FBC charges were serialized as zero (including H+ and ATP).
Never infer that a user's explicit zero charge is missing. Other libraries,
including the MetaNetX library, are untouched. No formula, stoichiometry,
boundary, GPR, or target-model metabolite is changed by this function.
"""
import gzip
import hashlib
import math
from pathlib import Path
import re

import cobra
import pandas as pd
from cobra.core.formula import elements_and_molecular_weights


KNOWN_ZERO_LIBRARY_SHA256 = '800e0d7afc6f75b47aa7af2907ae540396e0a6a5490259656ce071147758fe81'


def _formula_key(formula):
    if not isinstance(formula, str) or not re.fullmatch(r'(?:[A-Z][a-z]?\d*)+', formula):
        return None
    elements = cobra.Metabolite(formula=formula).elements
    if not elements or any(k not in elements_and_molecular_weights or
                           not math.isfinite(v) or v <= 0 or int(v) != v
                           for k,v in elements.items()):
        return None
    return tuple(sorted(elements.items()))


def repair_library_charges(library, library_path, reference_path):
    """Apply only unique exact-ID and exact-formula pairs, in memory.

The file digest establishes the *known erroneous library*, not a chemical
identity match. Chemical matching separately requires the reference ID and
all element counts (including H) to agree. Conflicting/invalid table rows,
unknown/nonzero existing charges and different protonation states are skipped.
COBRA's subsequent target-model resolution retains existing native chemistry.
"""
    report = dict(status='not_known_zero_charge_library', corrected=0, records=[],
                  library_path=str(library_path), reference_path=str(reference_path),
                  scope='known_reference_artifact_only_not_user_input', skipped={})
    if library_path is None:
        return report
    path = Path(library_path)
    raw = path.read_bytes()
    if raw[:2] == b'\x1f\x8b':
        raw = gzip.decompress(raw)
    digest = hashlib.sha256(raw).hexdigest()
    report['library_sha256_uncompressed'] = digest
    if digest != KNOWN_ZERO_LIBRARY_SHA256:
        return report
    if reference_path is None:
        report['status'] = 'matched_library_but_no_reference_table'
        return report
    reference = Path(reference_path)
    report['reference_sha256'] = hashlib.sha256(reference.read_bytes()).hexdigest()
    table = pd.read_excel(reference)
    if not {'id','formula','charge'} <= set(table.columns):
        raise ValueError('Reference charge table requires id, formula and charge')
    index = {}
    for row in table[['id','formula','charge']].itertuples(index=False, name=None):
        mid, formula, charge = row
        if not isinstance(mid, str):
            continue
        key = _formula_key(formula)
        try:
            number = float(charge)
            valid = math.isfinite(number) and number.is_integer()
        except (TypeError, ValueError):
            valid = False
        pair = (key, int(number)) if key and valid else None
        index.setdefault(mid, set()).add(pair)
    for met in library.metabolites:
        if met.charge != 0:
            reason = 'existing_nonzero_or_unknown_charge'
        elif met.id not in index:
            reason = 'no_exact_reference_id'
        elif len(index[met.id]) != 1 or None in index[met.id]:
            reason = 'invalid_or_conflicting_reference'
        else:
            key, charge = next(iter(index[met.id]))
            if _formula_key(met.formula) != key:
                reason = 'formula_or_protonation_mismatch'
            elif charge == 0:
                reason = 'matched_neutral_species'
            else:
                old_note = met.notes.get('CHARGE')
                met.charge = charge
                met.notes['MQC_LIBRARY_ORIGINAL_CHARGE'] = '0'
                met.notes['MQC_LIBRARY_CHARGE_REFERENCE'] = (
                    met.id+'; '+str(met.formula)+'; '+report['reference_sha256'])
                if old_note is not None:
                    met.notes['MQC_LIBRARY_ORIGINAL_CHARGE_NOTE'] = str(old_note)
                    met.notes['CHARGE'] = str(charge)
                report['records'].append(dict(metabolite=met.id, formula=met.formula,
                    before=0, after=charge, reference_id=met.id, original_charge_note=old_note))
                continue
        report['skipped'][reason] = report['skipped'].get(reason, 0)+1
    report['status'] = 'exact_reference_pairs_restored'
    report['corrected'] = len(report['records'])
    return report
