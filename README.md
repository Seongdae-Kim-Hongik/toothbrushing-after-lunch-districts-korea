# Toothbrushing after lunch across districts of South Korea, 2014-2024

Analysis code for the manuscript

> Kim S. *Does toothbrushing after lunch go together with other health behaviours across Korean districts?*

Archived version: Zenodo, https://doi.org/10.5281/zenodo.24359315 (v1.0.0).

This code is released under the MIT licence.

## Contents

This repository contains code only. The processed district-year data (analytic panel, result tables and boundary file) are not included; they are available from the author on reasonable request. All source statistics are public (see Data sources).

| Path | What it is |
|---|---|
| `01_스크립트/` | Analysis scripts (below) |
| `04_분석/BRSH/` | Model settings fixed before fitting (`sets_BRSH_std3.json`, `config_BRSH_std3.json`), the model chain (`run_fh_chain.sh`) and the district adjacency graph (`graph_593.graph`) |
| `04_분석/BRSH/MS/build/` | Scripts that read the result tables and write every number, table and figure in the manuscript and supplementary material |

Folder and some column names are in Korean because the source tables are Korean; each script has a header comment describing its inputs and outputs.

## Data sources

All inputs are public, aggregated statistics. No individual-level data are included.

- Toothbrushing after lunch (age-sex-standardised prevalence, standard error and respondents) by district: Korea Disease Control and Prevention Agency, Korea Community Health Survey, KOSIS table DT_H_OR_BRUSH (2014-2024).
- District characteristics: Korea Community Health Survey (walking, moderate or vigorous physical activity, current smoking, high-risk drinking, rear seat-belt use), public sports facility registers, health-resource statistics (dentists, clinics), public assistance statistics, census (tertiary education), census on establishments (manufacturing employment) and KOSIS local indicators (fixed adjustments and negative controls).

## Pipeline

Set `ST_DATA_ROOT` and `ST_GPT_DATASET` to the local copies of the public source files and the read-only district database, then run from the repository root.

| Step | Script | Output |
|---|---|---|
| 1 | `99BR_brsh_fh_prep.py` | Outcome on the logit scale with sampling variances (Sejong from the published table) |
| 2 | `99c_extra_canonical_pool_brsh.py` | Additional covariates |
| 3 | `80_build_disease_panel_v.py BRSH` (`ST_VER=std3`) | Analytic panel (interpolation, standardisation, variance inflation factors) |
| 4 | `04_분석/BRSH/run_fh_chain.sh` (`81fh_brsh.R`) | Spatio-temporal Fay-Herriot models, random covariate subsets, Knorr-Held interactions, sensitivity analyses |
| 5 | `99BR4_brsh_maps.py`, `04_분석/BRSH/MS/build/*.py` | Maps, numbers, figures, manuscript, supplementary material |

R 4.6.0 with R-INLA (1 thread), Python 3.9 with pandas, geopandas and python-docx.
