# Zenodo Release Helpers

These helper scripts prepare and upload MQC manuscript data packages to a Zenodo draft deposition.

They do not publish the Zenodo record or release a DOI. Publication should be done manually after reviewing the Zenodo metadata, license, creators, and uploaded file list.

## Scripts

- `zenodo_upload_mqc_models.py`: prepares paired before/after GEM model archives and uploads them to a Zenodo draft.
- `prepare_zenodo_supplemental_data.py`: prepares supplemental source-data archives for manuscript figures, analysis inputs, experimental validation data, and MQC reference resources.

## Paths

By default, both scripts infer the project root from this repository layout. To use a different local checkout, set:

```bash
export MQC_PROJECT_ROOT=/path/to/mqc
```

## Credentials

The scripts expect a local Zenodo token file, typically `zenodo.txt`, at the project root. Token files must never be committed to GitHub.
