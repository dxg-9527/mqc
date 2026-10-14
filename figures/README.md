# Manuscript figures

This directory contains the current manuscript artwork and the numerical inputs
used to draw Figures 2-6. It can be copied to another directory and rebuilt there.
No metabolic-model optimization is performed by the figure scripts.

## Contents

| Figure | Content | Build entry point |
| --- | --- | --- |
| 1 | MQC framework | Approved PNG artwork; no numerical plotting program |
| 2 | Model-level errors and co-occurring defects | `code/fig2_panels.py` |
| 3 | Energy, reducing equivalents and biomass outputs | `code/fig3_panels.py` |
| 4 | Reaction-level mechanism in iCN718 | `code/fig4_panels.py` |
| 5 | Literature phenotype comparisons | `code/fig5_panels.py` |
| 6 | Three isolates tested on nine substrates | `code/fig6_panels.py` |

`full_figures/` contains the six complete figures and their legends. `panels/`
contains 41 independent panels for Figures 2-6, supplied in PNG, PDF and SVG;
editable PowerPoint files are also provided for the mechanism panels.
`panel_code_data_manifest.csv` links each panel to its program and input tables.

The `fig4_final/`, `fig5_final/` and `fig6_final/` subdirectories contain the native
drawing programs and inputs used by their respective build entry points. The
figure-specific outputs are copied into `full_figures/` and `panels/` by these
entry points.

## Regenerate

Use a separate Python 3.12 environment for the figure dependencies:

```bash
python -m pip install -r figures/requirements.txt
python figures/code/generate_all.py
python figures/code/validate_release.py
```

LibreOffice (`soffice`) is also required to export the editable Figure 4
PowerPoint artwork to vector PDF. For example, on Ubuntu it can be installed with
`sudo apt-get install libreoffice-impress`. The Liberation Sans fonts used by the
mechanism artwork and their license are bundled under `fig4_final/assets/`.

The complete build regenerates the panels for Figures 2-6 and composes their
full figures. It preserves Figure 1 as the supplied framework artwork.
Individual figure entry points may also be run directly. Numerical and export
checks for Figures 4-6 run as part of their build workflows.

## Data and interpretation

The source tables retain the archived quantitative values used in the
manuscript. English substrate names and provenance labels replace Chinese labels
and machine-specific paths. Rows, numerical values, missing-value indicators and
model identifiers are preserved. A provenance field ending in `(archived source)`
identifies the original input record; it is not a path required for regeneration.

Figure 2 uses 6,543 models. Figure 3 uses the available quantitative outputs for
each metric and stage. Figure 5 retains prediction regressions and archived
zero-valued predictions. Figure 6 includes three isolates, nine substrates and
27 paired comparisons; its separate model survey covers 213 MAG models and 52
substrates. The isolate growth-fit flags remain in the input tables.

Figure 4 uses archived flux values. Its packaged iCN718 model supports a
read-only stoichiometric check; the build does not reproduce the historical
correction run. Historical imbalance annotations are not interchangeable with
the annotations in that packaged model. The final artwork therefore omits exact
imbalance residuals, and the allophanate branch denotes an independently tested
objective.

Hashes in the manifests describe the files in this directory. PDF metadata and
font rendering can vary across software installations, even when plotted values
are identical.
