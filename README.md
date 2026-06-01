# NULISA-CART-Proteome

Analysis code and supplementary data for:

**Distinct proteomic signatures define CAR T-cell expansion and development of acute toxicities in large B cell lymphoma patients undergoing CD19-directed CAR T-cell therapy**

Blaum EM, Beer JC, Mian A, Pramanik D, Kirkpatrick R, Ariyapala IS, Hao Q, Feng W, Ma X-J, Luo Y, Hill BT, Gupta N

## Overview

This repository contains the R analysis pipeline for a longitudinal proteomic study of 80 relapsed/refractory large B-cell lymphoma (LBCL) patients treated with CD19-directed CAR T-cell therapy. Plasma proteins (n = 204) were measured at six timepoints using the NULISAseq 200-plex Inflammation Panel (Alamar Biosciences) to characterize proteomic signatures associated with cytokine release syndrome (CRS), immune effector cell-associated neurotoxicity syndrome (ICANS), and CAR T-cell expansion.

## Repository structure

```
NULISA-CART-Proteome/
├── scripts/
│   ├── s00_setup.R                          # Configuration, packages, shared utilities
│   ├── s01_cohort_summary.R                 # Table 1, swimmer plots, sample collection
│   ├── s02_LME_models.R                     # Linear mixed-effects models (4 model sets)
│   ├── s03_clustering.R                     # K-means clustering + trajectory plots
│   ├── s04_pathway_enrichment.R             # GO, KEGG, Reactome, Hallmark enrichment
│   ├── s05_heatmaps.R                       # ComplexHeatmap figures with annotations
│   ├── s06_clinical_associations.R          # Fisher's exact, Wilcoxon rank-sum tests
│   ├── s07_cox_time_to_toxicity.R           # Cox PH models, forest plots, cum. incidence
│   ├── s08_survival_PFS_OS.R                # PFS and OS analysis
│   └── s09_decouple_expansion_toxicity.R    # Decoupling analysis (Fig 6)
├── data/
│   ├── private/                             # Patient-level data (.gitignored)
│   └── public/                              # Aggregate outputs, swimmer plot & sample
│       │                                    #   collection data (tracked in repo)
│       ├── swimmer_plot_data.xlsx
│       └── sample_collection.xlsx
├── figures/
│   ├── main/
│   └── extended_data/
├── supplementary_tables/
│   ├── SuppTable1_protein_panel_annotation.xlsx
│   ├── SuppTable2_LME_model_coefficients.xlsx
│   ├── SuppTable3_cluster_assignments.xlsx
│   ├── SuppTable4_pathway_enrichment.xlsx
│   ├── SuppTable5_cox_PH_results.xlsx
│   └── SuppTable6_decoupling_analysis.xlsx
├── renv.lock
├── .gitignore
├── LICENSE
└── README.md
```

## Setup

This project uses [renv](https://rstudio.github.io/renv/) for reproducible package management.

```r
# 1. Clone the repository
# git clone https://github.com/jcbeer/NULISA-CART-Proteome.git

# 2. Open R in the project directory and restore packages
renv::restore()

# 3. Set your project directory path in scripts/s00_setup.R
#    (line ~96: projectDir <- '/your/path/to/NULISA-CART-Proteome')

# 4. Run scripts in order
source("scripts/s00_setup.R")
source("scripts/s01_cohort_summary.R")
# ... etc.
```

**Note on NULISAseqR**: This package is developed by [Alamar Biosciences](https://alamarbio.com/) and is installed from GitHub. The `renv.lock` file records the installation source so `renv::restore()` can install it automatically. For manual installation outside of renv, see instructions in `scripts/s00_setup.R` or the [NULISAseqR documentation](https://nulisaseqr.alamarbio.com/user-guide/index.html#installation).

**Note on ggforestplot**: This package is not available on CRAN and must be installed from GitHub (`NightingaleHealth/ggforestplot`). It is recorded in `renv.lock` and will be installed automatically by `renv::restore()`.

## Data availability

Patient-level proteomic data (NPQ values), clinical metadata, and survival data are available upon request from the corresponding author. The survival data file (`data/private/survival_analysis.xlsx`) contains protected health information (PHI) and cannot be shared publicly.

Aggregate analysis outputs (model coefficients, cluster assignments, pathway enrichment results) are provided in the `data/public/` and `supplementary_tables/` directories. The supplementary tables in this repository correspond to Supplementary Tables 1–6 referenced in the manuscript.

## Scripts

All scripts source `s00_setup.R` for shared configuration. To reproduce the full analysis, run scripts in numerical order. Scripts s03–s09 can be run using the aggregate data outputs in `data/public/` without access to the patient-level data.

| Script | Description | Main outputs |
|--------|-------------|-------------|
| s00 | Setup, packages, constants, helper functions | — |
| s01 | Cohort summary, Table 1, swimmer plots, sample collection | Fig 1, Extended Data Fig 1, Fig 2A |
| s02 | LME models (time, CRS, ICANS, expansion interactions) | Supp Table 2, Extended Data Figs 2a, 3a, 5a, 9a |
| s03 | K-means clustering of LME coefficients + trajectories | Supp Table 3, Figs 2B, 3A, 3D, 5B |
| s04 | Pathway over-representation analysis | Supp Table 4, Figs 2D, 3C, 3F, 5D |
| s05 | Annotated heatmaps (ComplexHeatmap) | Figs 2C, 3B, 3E, 5C, Extended Data Figs 2b, 3b, 5b, 9b |
| s06 | Clinical associations (toxicity, expansion, response) | Extended Data Figs 4, 10, Fig 5E, 5H, 5I |
| s07 | Cox PH: time-to-severe CRS and ICANS | Supp Table 5, Fig 4, Extended Data Figs 6–7 |
| s08 | PFS and OS survival analysis | Fig 5F, 5G |
| s09 | Decoupling expansion from toxicity | Supp Table 6, Fig 6 |

**Note:** Extended Data Fig 8 (flow cytometry gating scheme) is a lab-generated figure and is not produced by this pipeline.

## Key R packages

- **NULISAseqR** (v1.5.0) — LME modeling for NULISAseq data
- **ComplexHeatmap** (v2.24.1) — annotated heatmaps
- **clusterProfiler** (v4.16.0) — pathway enrichment analysis
- **survival** (v3.8-3) — Cox proportional hazards models
- **survminer** (v0.5.1) — survival curve visualization
- **ggsurvfit** — Kaplan-Meier survival curves (PFS/OS)
- **ggforestplot** — forest plots (installed from GitHub: `NightingaleHealth/ggforestplot`)

## Citation

[Manuscript citation to be added upon publication]

## License

MIT License. See [LICENSE](LICENSE) for details.
