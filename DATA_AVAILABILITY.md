# Data Availability

The manuscript data release is staged in a Zenodo draft:

- Draft URL: https://zenodo.org/deposit/20537119
- Status: draft, not published. Replace this draft URL with the final DOI after publication.

## Zenodo Data Packages

The Zenodo draft contains paired genome-scale metabolic models before and after MQC quality control:

- `mqc_before_bigg_models.zip` and `mqc_after_bigg_models.zip`
- `mqc_before_vmh_models.zip` and `mqc_after_vmh_models.zip`
- `mqc_before_modelseed_models.zip` and `mqc_after_modelseed_models.zip`
- `mqc_before_carveme_models.zip` and `mqc_after_carveme_models.zip`

Model-level pairing and checksums are provided in:

- `mqc_before_after_model_manifest.csv`
- `mqc_before_after_model_summary.csv`
- `checksums_sha256.tsv`

The draft also contains manuscript source data and reproducibility resources:

- `mqc_manuscript_figure_source_data.zip`
- `mqc_analysis_input_tables.zip`
- `mqc_experimental_validation_raw_data.zip`
- `mqc_reference_resources.zip`
- `mqc_supplemental_source_data_manifest.csv`
- `checksums_sha256_supplemental.tsv`

## GitHub Repository Scope

The GitHub repository should contain source code, lightweight metadata, and reproducibility scripts. Large models, Zenodo release archives, local test data, intermediate analysis directories, generated build artifacts, and access tokens should not be committed to GitHub.

Sensitive files such as `zenodo.txt`, `.pypirc`, and `PyPi_API_token` must remain local only.
