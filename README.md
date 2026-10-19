# MQC

## Optimized cluster validation release

Branch `optimized-20261018` contains the frozen optimized correction implementation
and a portable Slurm batch runner. It does not replace the AWS deployment or the
original `main` branch. It includes the optimized implementation and subsequent
cluster failure fixes; it is no longer byte-identical to the initial three-strain
release. Exact source and reference hashes are in `cluster/release_lock.json`.

For the latest fixes and the 213-model follow-up list, see
[the cluster review release notes](CLUSTER_REVIEW_FIXES_20261019.md).

See [the Chinese cluster guide](cluster/README_zh.md) for model transfer,
environment checks, array submission, failed-job retry, and return files.
The runner invokes correction only, without a detection pre-run, and preserves
stage models, flux tables, final exports, failures and version provenance.
The first cluster run and a 622-model retry used earlier source locks. A full
6,543-model validation of this new release is still pending; do not merge results
from different locks as a same-version cohort.
The original runner tests and e_coli_core correction/export smoke run are
recorded in [the release validation](validation/optimized_cluster_smoke_20261018.json).

MQC detects quantitative defects in genome-scale metabolic models and supports
automated reaction and biomass curation. Detection can be run independently of
correction. The five diagnostic modules assess reducing-equivalent production,
energy production, net-metabolite production, metabolite yields and biomass
synthesis.

This branch contains the optimized MQC runtime derived from the web workflow, its required
reference tables and reaction libraries, and the programs and data for the
manuscript figures. The web interface is available at
<https://mqc.biodesign.ac.cn/>.

## Install

Use a dedicated Python environment. The package metadata supports Python 3.8 or
later, but this batch comparison uses Python 3.8.16 and the pinned dependencies
in `cluster/requirements-validated.txt`. Follow the cluster guide for this study;
the generic installation below does not reproduce its validated environment.

```bash
git clone --branch optimized-20261018 --single-branch https://github.com/dxg-9527/mqc.git mqc-optimized
cd mqc-optimized
python -m pip install .
```

MQC uses COBRApy for metabolic modelling and CPLEX in biomass curation. A CPLEX
installation licensed for the size of the input and merged reaction-library
models is required for a full correction run. The size-limited Community Edition
may be insufficient for these optimization problems. Reference files are
installed with the package; no separate model-library download is needed.

## Run

Run detection only:

```bash
mqc model.xml output/check --mode check
```

Run detection followed by correction:

```bash
mqc model.xml output/correction
```

Equivalent source-checkout entry points are `python main.py` and `python -m mqc`,
with the same arguments. Input formats include SBML, COBRA JSON, YAML and MATLAB.
Use a separate output directory for each run.

The Python interface exposes the same workflows:

```python
from mqc.main import run

report, flags = run("model.xml", "output/check", mode="check")
report, corrected_model, reference_model = run("model.xml", "output/correction")
```

Detection writes `check_result.json` and diagnostic plots. Correction additionally
writes `control_result.json`, `model_info.json`, a corrected model and a reference
model for subsequent curation. Reaction-level intermediate flux tables and
gap-filling models are written when required by the workflow.

The `ModelCheck` and `ModelControl` classes are retained for applications that use
the web workflow's report and reaction-editing interfaces.

Unsupported input structures or unresolved correction failures are retained in
the batch failure report, rather than counted as successful correction.

## Publication figures

The `figures/` directory retains the historical six figures, 41 standalone panels
for Figures 2-6, plotting inputs and their regeneration programs. See
[the figure instructions](figures/README.md) for dependencies and commands.
Regeneration uses the archived numerical inputs; it does not rerun metabolic
model correction. Figure 1 is supplied as the approved framework artwork.

## Repository layout

- `mqc/`: detection, correction, preprocessing and reference data.
- `main.py`: a small command-line launcher.
- `figures/`: figure programs, source tables, full figures and panels.
- `tests/`: installation, entry-point and packaged-resource checks.

The optimized workflow retains model serialization and reloading between
correction modules. Package-relative paths allow execution outside the repository
directory. Cluster tools live in `cluster/` and run correction without a detection pre-run.

## Validation

The AWS equivalence records below describe the original public baseline, not
the optimized branch. Optimized outputs are expected to differ from that
baseline. The cluster runner preserves both successful and unsuccessful cases
for comparison with the manuscript's historical corrected-model cohort.

A paired iCN718 run reproduced the deployed AWS implementation's detection
outputs, correction report, model metadata and final models. The corrected and
reference model files were byte-identical; the correction report differed only
in a translated text fragment. Details are recorded in
[the iCN718 comparison](validation/iCN718_aws_equivalence.json).

Both runs also reproduced an existing SBML export issue: two generated
metabolites have empty compartment identifiers, so the final XML cannot be
reopened by COBRApy. The accompanying `iCN718.json` files load successfully and
have identical model structure and a biomass objective value of 0.0. Use the JSON
output for this case. Those observations apply to the original baseline release.

## License

MQC code is distributed under the MIT License. Bundled fonts retain their own
license in the figure asset directory. Reference data preserve their upstream
identifiers and annotations.
