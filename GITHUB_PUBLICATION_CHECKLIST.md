# GitHub Publication Checklist

Use this checklist before pushing the manuscript-ready MQC repository to `https://github.com/dengxiao01/mqc.git`.

## Include In GitHub

- Core package code under `mqc/`, especially `mqc/main.py`, `mqc/utils.py`, and `mqc/control/`.
- Packaging and environment files: `setup.py`, `requirements.txt`, `environment.yaml`, `MANIFEST.in`, `LICENSE`, and `README.md`.
- Lightweight package metadata and database-mapping files required for normal operation, if each file is small enough for GitHub.
- Reproducibility documentation such as `DATA_AVAILABILITY.md` and this checklist.
- Manuscript figure-generation scripts after moving or copying them into a tracked directory such as `analysis/manuscript_figures/`.

## Do Not Include In GitHub

- Access tokens or credentials: `zenodo.txt`, `.pypirc`, `PyPi_API_token`, or any `*token*`, `*secret*`, or `*credential*` file.
- Zenodo release artifacts: `zenodo_upload/`, model ZIP archives, checksum files generated for a specific upload, or deposition JSON state files.
- Local model datasets: `mqc/local_test_data/`, `local_test_data/`, and `tmp/`.
- Generated build artifacts: `build/`, `dist/`, `*.egg-info/`, wheel files, and source distributions.
- Large local installers or binaries: CPLEX installers, R source archives, `.dmg`, `.bin`, `.sqlite`, and large `.zip` files.
- Large reference files that exceed GitHub's normal file-size limits, unless they are moved to Git LFS. Prefer Zenodo for these files.

## Current Repository Risks To Resolve

- The current working tree contains many unrelated local changes. Do not run `git add .`.
- `master` tracks the internal `gitlab/master`, while GitHub uses `origin/main`. Prefer preparing a clean branch from `origin/main` or a new `manuscript-release` branch.
- Some large reference XML files in `mqc/summary/` are larger than ordinary GitHub limits. These are included in the Zenodo supplemental package and should not be newly added to GitHub without Git LFS.
- If `.pypirc` or other token files were ever committed with real secrets, rotate those tokens before publication.

## Suggested Commit Groups

- Commit 1: source-code updates only.
- Commit 2: package metadata, README, and data-availability documentation.
- Commit 3: manuscript figure-generation scripts in a clean tracked `analysis/` directory.
- Commit 4, optional: lightweight example data only, if needed for tests and clearly documented.
