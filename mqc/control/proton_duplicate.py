"""Reject a NEW H-imbalanced duplicate that closes a pure-proton cycle.

Do not ban all H/charge review items. Require a balanced native counterpart,
identical non-proton species in their original compartments, and compatible
flux directions that cancel that counterpart. No chemistry is auto-corrected.
"""
from mqc.control.gapfill_chemistry import assess_gapfill_chemistry


def nonproton_signature(r, resolved=None):
    coefficients = {}
    for original, c in r.metabolites.items():
        m = (resolved or {}).get(original.id, original)
        if m.elements == {'H': 1} and m.charge == 1:
            continue
        coefficients[m.id] = float(c)
    if len(coefficients) < 2:
        return None
    first = sorted(coefficients)[0]
    scale = coefficients[first]
    return tuple((mid, round(coefficients[mid] / scale, 12)) for mid in sorted(coefficients)), scale


def native_proton_index(model):
    index = {}
    for r in model.reactions:
        if not assess_gapfill_chemistry(r)['full_mass_charge_pass']:
            continue
        signature = nonproton_signature(r)
        if signature is not None:
            key, scale = signature
            index.setdefault(key, []).append((r, scale))
    return index


def proton_duplicate_conflict(candidate, assessment, resolved, index):
    if assessment['status'] != 'hydrogen_charge_or_annotation_review':
        return None
    imbalance = assessment['atom_imbalance']
    # A known net H+ discrepancy, not unknown charge or uncertain annotation.
    if set(imbalance) != {'H'} or assessment['charge_imbalance'] is None:
        return None
    if abs(imbalance['H'] - assessment['charge_imbalance']) > 1e-9:
        return None
    signature = nonproton_signature(candidate, resolved)
    if signature is None: return None
    key, scale = signature
    for native, native_scale in index.get(key, []):
        ratio = scale / native_scale
        for direction in [-1, 1]:
            if direction < 0 and candidate.lower_bound >= 0: continue
            if direction > 0 and candidate.upper_bound <= 0: continue
            opposite = -ratio * direction
            if opposite < 0 and native.lower_bound >= 0: continue
            if opposite > 0 and native.upper_bound <= 0: continue
            return dict(candidate=candidate.id, native=native.id,
                        candidate_direction=direction, native_relative_flux=opposite,
                        net_H_and_charge=imbalance['H'] * direction,
                        reason='known_unbalanced_duplicate_enables_pure_proton_cycle')
    return None
