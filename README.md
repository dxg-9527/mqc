# MQC

MQC detects quantitative defects in genome-scale metabolic models and supports
automated reaction and biomass curation. Detection can be run independently of
correction. The five diagnostic modules assess reducing-equivalent production,
energy production, net-metabolite production, metabolite yields and biomass
synthesis.

This repository contains the MQC runtime used by the web workflow, its required
reference tables and reaction libraries, and the programs and data for the
manuscript figures. The web interface is available at
<https://mqc.biodesign.ac.cn/>.

## Install

Use a dedicated Python environment. The core package supports Python 3.8 or later;
Python 3.10 or later is recommended for new environments.

```bash
git clone https://github.com/dxg-9527/mqc.git
cd mqc
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

The deployed biomass-correction workflow currently expects an ADP term in the
biomass reaction. Models without that term can stop during correction; detection
mode remains available for inspecting their diagnostic results.

## Publication figures

The `figures/` directory contains the current six figures, 41 standalone panels
for Figures 2-6, plotting inputs and their regeneration programs. See
[the figure instructions](figures/README.md) for dependencies and commands.
Regeneration uses the archived numerical inputs; it does not rerun metabolic
model correction. Figure 1 is supplied as the approved framework artwork.

## Repository layout

- `mqc/`: detection, correction, preprocessing and reference data.
- `main.py`: a small command-line launcher.
- `figures/`: figure programs, source tables, full figures and panels.
- `tests/`: installation, entry-point and packaged-resource checks.

The standalone entry point follows the deployed web workflow, including model
serialization and reloading between correction modules. Package-relative paths
allow execution outside the repository directory. Duplicate launchers, local
batch jobs and commented-out code are excluded.

## Validation

A paired iCN718 run reproduced the deployed AWS implementation's detection
outputs, correction report, model metadata and final models. The corrected and
reference model files were byte-identical; the correction report differed only
in a translated text fragment. Details are recorded in
[the iCN718 comparison](validation/iCN718_aws_equivalence.json).

Both runs also reproduced an existing SBML export issue: two generated
metabolites have empty compartment identifiers, so the final XML cannot be
reopened by COBRApy. The accompanying `iCN718.json` files load successfully and
have identical model structure and a biomass objective value of 0.0. Use the JSON
output for this case. The release preserves this deployed behavior.

## License

MQC code is distributed under the MIT License. Bundled fonts retain their own
license in the figure asset directory. Reference data preserve their upstream
identifiers and annotations.
