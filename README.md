# Resiquimod in Cutaneous T-Cell Lymphoma (CTCL)

Single-cell transcriptomic analysis of resiquimod (TLR7/8 agonist) treatment in cutaneous T-cell lymphoma (CTCL) skin biopsies.

## Overview

This project analyzes two complementary single-cell RNA-seq datasets from CTCL patient skin biopsies:

- **FFPE scRNA-seq**: 10X Chromium Fixed RNA Profiling (Flex) from formalin-fixed paraffin-embedded biopsies of treated lesional skin (Week 0 and Week 8), across 5 patients treated with resiquimod 0.03% or 0.06%
- **Single-nucleus RNA-seq (snRNA-seq)**: 10X Chromium snRNA-seq from fresh-frozen biopsies of untreated lesional and non-lesional skin (Week 0 and Week 24), across 8 patients multiplexed with souporcell for demultiplexing

After quality control and doublet removal, the merged dataset contains ~71,000 cells spanning 12 cell types: keratinocytes (undifferentiated and differentiated), fibroblasts, T cells, myeloid cells, vascular and lymphatic endothelial cells, pericytes, melanocytes, hair follicle cells, mast cells, and dermal adipocytes.

## Repository Structure

```
resiquimod_ctcl/
├── _notebooks/          # Analysis notebooks (numbered in execution order)
│   ├── utils.R          # Shared R helper functions and color palettes
│   ├── 01_Seurat_ffpe.ipynb
│   ├── 01_Seurat_nuc-seq.ipynb
│   ├── 02_run_scrublet.ipynb
│   ├── 03_merge_datasets_label_cell_types.ipynb
│   ├── 04_Figure1_plots.ipynb
│   ├── 05_DE_LesUntxvsNL_Wk0.ipynb
│   ├── 06_HTS_responder_status.ipynb
│   ├── 06_HTS_clonal_tracking.ipynb
│   ├── 07_myeloid_subclusters_TLR.ipynb
│   ├── 08_T_subclusters.ipynb
│   ├── 09_DE_Treated.ipynb
│   └── 10_DE_Untreated.ipynb
├── _data/               # Input data (not committed; see Data Structure below)
├── _figures/            # Output figures (PDFs)
└── _results/            # Output tables and GSEA plots
```

## Analysis Pipeline

All notebooks are in `_notebooks/` and are numbered in order of execution. Each notebook begins with a markdown cell documenting its purpose, inputs, outputs, and dependencies.

| Notebook | Language | Description |
|----------|----------|-------------|
| `01_Seurat_ffpe.ipynb` | R | Seurat preprocessing of FFPE data: QC filtering, Scrublet doublet removal (mixture model approach), normalization, Harmony batch correction, clustering, and cell type annotation (27,448 cells) |
| `01_Seurat_nuc-seq.ipynb` | R | Seurat preprocessing of snRNA-seq data: souporcell-based demultiplexing, QC filtering, Scrublet doublet removal, Harmony integration, clustering, and cell type annotation (43,854 cells) |
| `02_run_scrublet.ipynb` | Python | Scrublet doublet scoring for both FFPE and snRNA-seq count matrices |
| `03_merge_datasets_label_cell_types.ipynb` | R | Integration of FFPE and snRNA-seq datasets via Harmony, joint clustering, and consensus cell type labeling across the merged 71,302-cell dataset |
| `04_Figure1_plots.ipynb` | R | Figure 1 visualizations: UMAPs, cell type proportion analysis, marker gene heatmaps/dot plots, and statistical comparison of lesional vs. non-lesional cell composition |
| `05_DE_LesUntxvsNL_Wk0.ipynb` | R | Differential expression between untreated lesional and non-lesional skin at baseline (Week 0), with GSEA pathway analysis per cell type |
| `06a_HTS_responder_status.ipynb` | R | Treatment responder classification using high-throughput TCR sequencing: quantifies malignant T-cell clone abundance changes from Week 0 to Week 8/24 |
| `06b_HTS_clonal_tracking.ipynb` | R | TCR clonal tracking analysis: benign T-cell recruitment dynamics and correlation with malignant clone clearance |
| `07_myeloid_subclusters_TLR.ipynb` | R | Myeloid cell subclustering and TLR7/8 expression analysis |
| `08_T_subclusters.ipynb` | R | T-cell subclustering with marker gene characterization |
| `09_DE_Treated.ipynb` | R | Differential expression in treated lesions: Week 8 vs. Week 0, responders vs. non-responders comparisons, with GSEA pathway enrichment |
| `10_DE_Untreated.ipynb` | R | Differential expression in untreated lesional skin over time (Week 0 → Week 24), with GSEA pathway enrichment |
| `utils.R` | R | Shared helper functions (see below) |

### Shared utilities (`utils.R`)

`utils.R` is sourced at the top of every R notebook via `source('./utils.R')`. It provides:

- **Color palettes**: `custom.colors` (12 cell types), `color_scheme` (tissue condition)
- **`fig.size(h, w)`**: Set notebook output plot dimensions
- **`load_snRNA_obj()`**: Load the snRNA-seq Seurat object with harmonized cell type labels from the merged dataset
- **`load_ffpe_obj()`**: Load the FFPE Seurat object with harmonized cell type labels
- **`load_gene_sets()`**: Download MSigDB Hallmark, KEGG canonical, and GO:BP gene sets; returns a list with `$hallmark` and `$all`
- **`volcano_plot(markers_ct, tolabel, title, nx, ny)`**: Consistent volcano plot for DE results
- **`run_GSEA_on_list(markers_list, out_dir, label, gene_set, plot_width)`**: Run GSEA across all cell types in a DE results list and save plots
- **`run_gsea_analysis(de_results, cell_type, gene_set)`**: Run GSEA for a single cell type
- **`plot_genes_hex()`** / **`plot_split_hex()`**: Hexbin expression plots on UMAP
- **`plot_gene_violins()`**: Split violin plots by timepoint
- **`find_valley()`**: Gaussian mixture model valley detection (used for scrublet thresholding)

## Environment Setup

### R (v4.4.3)

This project uses [`renv`](https://rstudio.github.io/renv/) for reproducible package management. The exact package versions are recorded in `renv.lock`.

To restore the environment:

```r
install.packages("renv")
renv::restore()  # installs all packages from renv.lock
```

If you prefer to install manually:

```r
install.packages(c("Seurat", "ggplot2", "ggrepel", "data.table", "tidyverse",
                   "Matrix", "viridis", "lme4", "ggpubr", "ggrastr", "cowplot",
                   "mixtools", "rstatix", "RColorBrewer", "patchwork", "scales"))

# Bioconductor packages
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("MAST", "SingleCellExperiment", "limma",
                       "clusterProfiler", "org.Hs.eg.db"))

# GitHub packages
remotes::install_github("immunogenomics/harmony")
remotes::install_github("immunogenomics/presto")

# MSigDB gene sets
install.packages("msigdbr")

#Note
If you run into the error: "inv(): use of LAPACK must be enabled" when running Harmony, please install Harmony with the fix below (per https://github.com/immunogenomics/harmony/issues/284):

devtools::install_github("pati-ni/harmony", ref="lapack-win-fix", force=TRUE)
```

### Python (notebook 02 only)

```bash
pip install scrublet scanpy scipy pandas matplotlib numpy
```

## Data Structure

Raw data is not committed to this repository. Notebooks expect data in the `_data/` directory:

```
_data/
├── ffpe/                          # FFPE Cell Ranger per-sample output directories
│   └── <run_id>_4plex/<run_id>_multi/outs/per_sample_outs/<sample_id>/
├── snRNAseq/                      # snRNA-seq Cell Ranger output + souporcell results
│   ├── <run_id>/outs/filtered_feature_bc_matrix/
│   └── obj_43854cells_labeled_v2.rds   # generated by 01_Seurat_nuc-seq.ipynb
├── HTS/                           # Adaptive Biotechnologies TCR-seq rearrangement TSVs
│   └── CombinedRearrangements_<patient>.tsv
├── obj_27448cells_labeled_v2.rds  # generated by 01_Seurat_ffpe.ipynb
├── obj_merged_71302cells_labeled.rds   # generated by 03_merge_datasets_label_cell_types.ipynb
├── obj_merged_scRNA_labels.csv    # generated by 03_merge_datasets_label_cell_types.ipynb
├── obj_merged_snRNA_labels.csv    # generated by 03_merge_datasets_label_cell_types.ipynb
└── sample_meta_2026.csv           # clinical metadata (patient, timepoint, response scores)
```

Intermediate `.rds` objects are saved between steps so individual notebooks can be re-run without rerunning the full pipeline from scratch.

## Execution Order

```
02_run_scrublet.ipynb  (Python — run first to generate doublet scores)
        ↓
01_Seurat_ffpe.ipynb ──────────────┐
01_Seurat_nuc-seq.ipynb ───────────┤
        ↓                          │
03_merge_datasets_label_cell_types │ (reads both .rds from step above)
        ↓                          │
04_Figure1_plots.ipynb             │
05_DE_LesUntxvsNL_Wk0.ipynb       │ (reads snRNA obj + merged labels)
06a_HTS_responder_status.ipynb     │ (independent — reads only HTS TSVs)
06b_HTS_clonal_tracking.ipynb      │ (independent — reads only HTS TSVs)
07_myeloid_subclusters_TLR.ipynb   │ (reads merged obj)
08_T_subclusters.ipynb             │ (reads merged obj)
09_DE_Treated.ipynb                │ (reads FFPE obj + merged labels)
10_DE_Untreated.ipynb              │ (reads snRNA obj + merged labels)
```

Notebooks 04–10 are largely independent of each other once step 03 is complete.
