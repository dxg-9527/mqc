"""Experiment-independent chemistry gate for NEW gap-filling candidates.

Reject known non-hydrogen atom imbalance and unverifiable molecular formulas.
H/charge-only discrepancies remain explicit review items, not full chemical
passes: resolving protonation or compartment coupling is a separate task.
No stoichiometry, formulas, bounds or existing model reactions are changed here.
"""
import math
import re
from mqc.control.pathway_rules import sugar_identity_conflict
from mqc.control.quinone_methylation import incomplete_methylation


ELEMENTS = set("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og".split())
FORMULA = re.compile(r"(?:[A-Z][a-z]?(?:\d+(?:\.\d*)?|\.\d+)?)+")


def incomplete_reference_template(reaction, resolved=None):
    """Quarantine the documented R01734/DBRC lump, not ordinary orphan enzymes.

    iRC1080 Table S2 maps DBRC to KEGG R01734, whose source explicitly says
    'unclear reaction' and 'incomplete reaction'. Atom balancing cannot supply
    the missing biochemical evidence. Match its complete heavy-species
    signature across namespaces, not a strain ID, reaction name, or formula
    alone. Protonation variants and reversed/scaled equations are equivalent.
    Native reactions are not deleted by the new-candidate gate.

    https://rest.kegg.jp/get/R01734
    https://doi.org/10.1038/msb.2011.52 (Supplementary Table S2, DBRC)
    """
    aliases = {
        'dmb': {'dmbzid', 'C03114', 'cpd01997', 'MNXM1531'},
        'rib': {'ribflv', 'C00255', 'cpd00220', 'MNXM270'},
        'nic': {'nicrnt', 'C01185', 'cpd00873', 'MNXM194'},
        'nh4': {'nh4', 'C01342', 'cpd00013', 'MNXM15'},
        'o2': {'o2', 'C00007', 'cpd00007', 'MNXM4'},
        'water': {'h2o', 'C00001', 'cpd00001', 'MNXM2'},
        'formate': {'for', 'C00058', 'cpd00047', 'MNXM39'},
        'pi': {'pi', 'C00009', 'cpd00009', 'MNXM9'},
    }
    heavy = {'dmb': {'C':9,'N':2}, 'rib': {'C':17,'N':4,'O':6},
             'nic': {'C':11,'N':1,'O':9,'P':1}, 'nh4': {'N':1},
             'o2': {'O':2}, 'water': {'O':1}, 'formate': {'C':1,'O':2},
             'pi': {'O':4,'P':1}}
    if len(reaction.metabolites) not in (8, 9):
        return None
    found, compartments, proton_count = {}, set(), 0
    for original, coefficient in reaction.metabolites.items():
        met = (resolved or {}).get(original.id, original)
        if not math.isfinite(coefficient) or not met.compartment:
            return None
        compartment = met.compartment[2:] if met.compartment.startswith('C_') else met.compartment
        compartments.add(compartment)
        if met.elements == {'H': 1} and met.charge == 1:
            proton_count += 1
            continue
        identities = {met.id}
        for comp in {met.compartment, compartment}:
            for suffix in ('_'+comp, '['+comp+']', '@'+comp):
                if met.id.endswith(suffix):
                    identities.add(met.id[:-len(suffix)])
        for key, value in met.annotation.items():
            if any(token in key.lower() for token in ('bigg', 'kegg', 'seed', 'metanetx')):
                values = [value] if isinstance(value, str) else value if isinstance(value, (list,tuple,set)) else []
                identities.update(v.rstrip('/').rsplit('/',1)[-1] for v in values if isinstance(v,str))
        roles = [role for role, names in aliases.items() if identities & names]
        if len(roles) != 1 or roles[0] in found:
            return None
        role = roles[0]
        if {e:n for e,n in met.elements.items() if e != 'H'} != heavy[role]:
            return None
        found[role] = coefficient
    if len(compartments) != 1 or proton_count > 1 or set(found) != set(aliases):
        return None
    scale = -found['dmb']
    if not scale or any(not math.isclose(found[role]/scale, value, abs_tol=1e-9, rel_tol=1e-9)
                       for role, value in {'dmb':-1,'nic':-1,'nh4':-1,'water':-1,
                                           'o2':-3,'formate':3,'pi':1,'rib':1}.items()):
        return None
    return dict(reference='KEGG:R01734 / iRC1080:DBRC',
                reason='documented_unclear_incomplete_conversion_not_validated_enzyme_template',
                source='https://rest.kegg.jp/get/R01734',
                scope='new_gapfill_candidate_only', automatic_replacement=False)


def assess_gapfill_chemistry(reaction, tolerance=1e-9, resolved_metabolites=None):
    """Classify a candidate without converting unknown chemistry into balance.

Single-metabolite boundaries are excluded by this gate; MQC handles boundaries
separately. Fractional, finite stoichiometry is allowed. Unknown charge does not
invalidate known atom conservation, but is explicitly marked as unverified.
An optional ID-to-metabolite mapping predicts COBRA's target-model resolution
without adding unverified candidates to a solver or modifying either model.
"""
    record = {"reaction_id": reaction.id, "allowed": False, "status": "unverifiable_formula",
              "atom_imbalance": {}, "charge_imbalance": None, "formula_issues": {},
              "unknown_charge": [], "full_mass_charge_pass": False}
    if len(reaction.metabolites) <= 1:
        record["status"] = "boundary_or_empty"
        return record
    total, charge = {}, 0.0
    for met, coefficient in reaction.metabolites.items():
        if resolved_metabolites is not None:
            met = resolved_metabolites.get(met.id, met)
        if not math.isfinite(float(coefficient)):
            record["formula_issues"][met.id] = "nonfinite_stoichiometry"
            continue
        try:
            elements = (met.elements if isinstance(met.formula, str) and FORMULA.fullmatch(met.formula)
                        else None)
        except (TypeError, ValueError):
            elements = None
        if (not elements or any(k not in ELEMENTS or not math.isfinite(float(v)) or v < 0
                                for k,v in elements.items())):
            record["formula_issues"][met.id] = "missing_invalid_or_generic_formula"
            continue
        for element, count in elements.items():
            total[element] = total.get(element, 0.0) + coefficient * count
        try:
            numeric_charge = float(met.charge)
        except (TypeError, ValueError):
            numeric_charge = float("nan")
        if not math.isfinite(numeric_charge):
            record["unknown_charge"].append(met.id)
        else:
            charge += coefficient * numeric_charge
    if record["formula_issues"]:
        # Partial atom sums would be misleading; never infer balance from them.
        return record
    reference_issue = incomplete_reference_template(reaction, resolved_metabolites)
    if reference_issue:
        record['status'] = 'documented_incomplete_reference_template'
        record['reference_issue'] = reference_issue
        return record
    record["atom_imbalance"] = {k: float(v) for k,v in total.items() if abs(v) > tolerance}
    if not record["unknown_charge"]:
        record["charge_imbalance"] = float(charge) if abs(charge) > tolerance else 0.0
    methylation_issue = incomplete_methylation(reaction, resolved_metabolites)
    if methylation_issue:
        record['status'] = 'incomplete_quinone_methylation_redox'
        record['redox_issue'] = methylation_issue
        return record
    if any(k != "H" for k in record["atom_imbalance"]):
        record["status"] = "non_hydrogen_atom_imbalance"
        return record
    conflict = sugar_identity_conflict(reaction, resolved_metabolites)
    if conflict:
        record['status'] = 'sugar_hydrolysis_identity_conflict'
        record['identity_conflict'] = conflict
        return record
    record["allowed"] = True
    if record["atom_imbalance"] or record["charge_imbalance"] != 0.0:
        record["status"] = "hydrogen_charge_or_annotation_review"
    else:
        record["status"] = "balanced"
        record["full_mass_charge_pass"] = True
    return record


def record_gapfill_assessment(report, assessment, phase):
    """Keep counts plus all exclusions/review items for source traceability."""
    report["assessments"] += 1
    key = phase + ":" + assessment["status"]
    report["counts"][key] = report["counts"].get(key, 0) + 1
    if assessment["status"] != "balanced":
        report["records"].append(dict(assessment, phase=phase))
