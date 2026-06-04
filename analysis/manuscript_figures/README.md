# Manuscript Figure Reproduction Scripts

This directory contains scripts used to organize manuscript source data and redraw manuscript figures.

The corresponding source data are not stored in GitHub. Download them from the Zenodo record described in `DATA_AVAILABILITY.md`, especially:

- `mqc_manuscript_figure_source_data.zip`
- `mqc_analysis_input_tables.zip`
- `mqc_experimental_validation_raw_data.zip`
- `mqc_reference_resources.zip`

## Scripts

- `build_mqc_fig2_mechanistic_atlas.py`: reproduces and refines the Figure 2 mechanistic error clustering analysis.
- `extract_new_features_to_csv.py`: extracts the feature matrix used by the original Figure 2 new-clustering baseline.
- `cluster_new_error_profiles.py`: clusters the original new error profiles for traceability.
- `redraw_mqc_fig3_fig4_fig5_true_data.py`: redraws Figure 3, Figure 4, and Figure 5 panels using true source data.
- `organize_mqc_fig2_5_data.py`: organizes source data for Figures 2-5.
- `organize_mqc_figure_data.py`: earlier source-data organization helper.
- `generate_supplementary_tables.py`: generates supplementary tables.
- `redraw_mqc_figure1_framework.py`: redraws the MQC framework schematic.

Some scripts were originally developed in a local analysis workspace and may contain local default paths. For publication-grade reproducibility, run them after unpacking the Zenodo source-data archives into the expected project-relative locations or adapt the path constants near the top of each script.
