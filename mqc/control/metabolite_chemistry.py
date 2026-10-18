"""Correct verified charge errata, not arbitrary reference disagreements.

CHEBI:4991 explicitly specifies C12H10FeO14, net charge -3.
Some input reconstructions instead encode the same identified species as -12.
Only that exact incorrect combination is changed. Preserve other ionic forms,
explicit zero/unknown charges, formulas, bounds and reaction stoichiometry.

CHEBI:58695 / CHEBI:134283 specify C6H11O9P, charge -2, for D-tagatose
6-phosphate. Repair only the known -3 mismatch for this identified species,
not other hexose phosphates or different protonation states.
"""

SOURCE = 'https://www.ebi.ac.uk/chebi/CHEBI:4991'
ALIASES = {'bigg.metabolite': {'fe3dcit'}, 'seed.compound': {'cpd03725'},
           'kegg.compound': {'C06229'}, 'metanetx.chemical': {'MNXM146388'},
           'chebi': {'CHEBI:4991', '4991'}}
ELEMENTS = {'C':12, 'H':10, 'Fe':1, 'O':14}
TAGATOSE_SOURCE = 'https://www.ebi.ac.uk/chebi/CHEBI:58695'
TAGATOSE_ALIASES = {
    'bigg.metabolite': {'tag6p__D'}, 'seed.compound': {'cpd00805'},
    'kegg.compound': {'C01097'}, 'metanetx.chemical': {'MNXM795'},
    'biocyc': {'META:TAGATOSE-6-PHOSPHATE', 'META:CPD-15826'},
    # Acid/base and ring/open-chain cross-references describe the same
    # stereospecific sugar. The exact molecular formula remains mandatory.
    'chebi': {'CHEBI:13025','CHEBI:17837','CHEBI:21097','CHEBI:4251',
              'CHEBI:58695','CHEBI:47947','CHEBI:134283'},
}


def ferric_dicitrate_charge_error(metabolite):
    if metabolite.charge != -12 or metabolite.elements != ELEMENTS:
        return False
    comp = metabolite.compartment or ''
    base = metabolite.id
    for alias in [comp] + ([comp[2:]] if comp.startswith('C_') else []):
        if not alias:
            continue
        for suffix in ['_' + alias, '[' + alias + ']', '@' + alias]:
            if base.endswith(suffix):
                base = base[:-len(suffix)]
                break
    known = {value for values in ALIASES.values() for value in values}
    matched = base in known
    for key, allowed in ALIASES.items():
        values = metabolite.annotation.get(key, [])
        if isinstance(values, str):
            values = [values]
        if not isinstance(values, (list, tuple, set)):
            return False
        # Do not silently overrule a contradictory explicit identity.
        if any(value not in allowed for value in values):
            return False
        matched = matched or bool(values)
    return matched


def tagatose_6p_charge_error(metabolite):
    """An identity AND formula AND known wrong charge match; never name-only."""
    if metabolite.charge != -3 or metabolite.elements != {'C':6,'H':11,'O':9,'P':1}:
        return False
    comp = metabolite.compartment or ''
    base = metabolite.id
    for alias in [comp] + ([comp[2:]] if comp.startswith('C_') else []):
        if not alias:
            continue
        for suffix in ['_'+alias, '['+alias+']', '@'+alias]:
            if base.endswith(suffix):
                base = base[:-len(suffix)]
                break
    known = {value for values in TAGATOSE_ALIASES.values() for value in values}
    matched = base in known
    for key, allowed in TAGATOSE_ALIASES.items():
        values = metabolite.annotation.get(key, [])
        if isinstance(values, str):
            values = [values]
        if (not isinstance(values, (list, tuple, set)) or
                any(not isinstance(value, str) or value not in allowed for value in values)):
            return False
        matched = matched or bool(values)
    return matched


def repair_known_charges(model):
    """In-memory, idempotent curation with original values and source retained."""
    records = []
    for met in model.metabolites:
        if ferric_dicitrate_charge_error(met):
            before, after, source = -12, -3, SOURCE
        elif tagatose_6p_charge_error(met):
            before, after, source = -3, -2, TAGATOSE_SOURCE
        else:
            continue
        original_note = met.notes.get('CHARGE')
        row = dict(metabolite=met.id, formula=met.formula, before=before, after=after,
                   source=source, original_charge_note=original_note,
                   kind='verified_input_annotation_erratum_not_charge_parser_bug')
        met.charge = after
        # Keep legacy and FBC exports consistent while retaining the old note.
        if original_note is not None:
            met.notes['MQC_ORIGINAL_CHARGE_NOTE'] = str(original_note)
            met.notes['CHARGE'] = str(after)
        met.notes['MQC_CHARGE_ERRATUM'] = '{} -> {}; {}'.format(before, after, source)
        records.append(row)
    return records
