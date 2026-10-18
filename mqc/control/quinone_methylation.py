"""Detect two exact, electron-incomplete legacy quinone methylation lumps.

SAM -> SAH methyl transfer does not supply the missing two electrons when
the quinone ring also changes oxidation state. Require complete, known
chemistry; do not guess a donor, acceptor or replacement substrate. Balanced
quinone->quinone and quinol->quinol lumps and explicitly coupled reactions
are unaffected. No phenotype, reaction ID or strain name is used.

Primary enzymology: https://pubmed.ncbi.nlm.nih.gov/36417754/
The older EC 2.1.1.163 web equation itself mixes oxidation states; that text
is NOT used as a balanced stoichiometric template.
"""
import math

METHYL_REDOX_RULE = 'Quinone methylation changes redox state without an electron partner'
SOURCE = 'https://pubmed.ncbi.nlm.nih.gov/36417754/'
SPECIES = {
    'dmq': ({'2dmmq8', 'cpd15352', 'MNXM2178'},
            {'2-demethylmenaquinone 8', '2-demethylmenaquinone-8'}, 'C50H70O2', 0),
    'dmqh2': ({'2dmmql8', 'cpd15353', 'MNXM558'},
              {'2-demethylmenaquinol 8', '2-demethylmenaquinol-8'}, 'C50H72O2', 0),
    'mq': ({'mqn8', 'cpd15500', 'MNXM509'},
           {'menaquinone 8', 'menaquinone-8'}, 'C51H72O2', 0),
    'mqh2': ({'mql8', 'cpd15499', 'MNXM223'},
             {'menaquinol 8', 'menaquinol-8'}, 'C51H74O2', 0),
    'sam': ({'amet', 'C00019', 'cpd00017', 'MNXM16'},
            {'s-adenosyl-l-methionine'}, 'C15H23N6O5S', 1),
    'sah': ({'ahcys', 'C00021', 'cpd00019', 'MNXM19'},
            {'s-adenosyl-l-homocysteine'}, 'C14H20N6O5S', 0),
    'h': ({'h', 'C00080', 'cpd00067', 'MNXM1'}, {'h+', 'proton'}, 'H', 1),
}


def incomplete_methylation(reaction, resolved=None):
    """Read-only exact signature, including atom and charge residual checks."""
    if len(reaction.metabolites) != 5:
        return None
    roles, compartments, atoms, charge = {}, set(), {}, 0.
    for original, coefficient in reaction.metabolites.items():
        met = (resolved or {}).get(original.id, original)
        if not met.compartment or not math.isfinite(coefficient) or not coefficient:
            return None
        comp = met.compartment
        comp_aliases = [comp] + ([comp[2:]] if comp.startswith('C_') else [])
        base = met.id
        for alias in comp_aliases:
            for suffix in ('_'+alias, '['+alias+']', '@'+alias):
                if base.endswith(suffix):
                    base = base[:-len(suffix)]
                    break
        ids = {base}
        for key, value in met.annotation.items():
            if key.lower() not in {'bigg.metabolite', 'seed.compound', 'kegg.compound',
                                   'metanetx.chemical', 'seed compound', 'kegg compound',
                                   'metanetx (mnx) chemical'}:
                continue
            values = [value] if isinstance(value, str) else value if isinstance(value, (list, tuple, set)) else []
            for entry in values:
                if isinstance(entry, str):
                    ids.update(v.strip().rstrip('/').rsplit('/', 1)[-1] for v in entry.split(','))
        name = (met.name or '').strip().lower()
        matches = [role for role, (aliases, names, _, _) in SPECIES.items()
                   if ids & aliases or name in names]
        if len(matches) != 1 or matches[0] in roles:
            return None
        role = matches[0]
        if (met.formula, met.charge) != SPECIES[role][2:]:
            return None
        roles[role] = coefficient
        compartments.add(comp_aliases[-1])
        for e, n in met.elements.items():
            atoms[e] = atoms.get(e, 0.) + n * coefficient
        charge += met.charge * coefficient
    if len(compartments) != 1 or 'sam' not in roles:
        return None
    scale = -roles['sam']
    for substrate, product, protons, electrons in [('dmq', 'mqh2', -1, -2),
                                                 ('dmqh2', 'mq', 3, 2)]:
        expected = {'sam': -1, 'sah': 1, substrate: -1, product: 1, 'h': protons}
        if set(roles) != set(expected):
            continue
        if any(not math.isclose(roles[k]/scale, v, rel_tol=1e-9, abs_tol=1e-9)
               for k, v in expected.items()):
            continue
        if (any(abs(v/scale) > 1e-8 for v in atoms.values()) or
                not math.isclose(charge/scale, electrons, abs_tol=1e-8)):
            return None
        return dict(reaction=reaction.id, substrate_role=substrate, product_role=product,
                    charge_residual_per_methylation=electrons, atom_residual={},
                    reason='two_electron_change_without_donor_or_acceptor', source=SOURCE,
                    replacement_inferred=False)
    return None


def methylation_bounds(reaction):
    """Quarantine an incomplete equation, retaining the reaction for audit."""
    if incomplete_methylation(reaction) is None:
        return reaction.bounds
    if reaction.lower_bound > 0 or reaction.upper_bound < 0:
        raise ValueError('Forced flux in incomplete quinone methylation: ' + reaction.id)
    return (0., 0.)
