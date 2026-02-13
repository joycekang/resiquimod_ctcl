# Resiquimod in Cutaneous T-Cell Lymphoma (CTCL)

Single-cell transcriptomic analysis of resiquimod (TLR7/8 agonist) treatment in cutaneous T-cell lymphoma (CTCL) skin biopsies.

## Overview

This project analyzes two complementary single-cell RNA-seq datasets from CTCL patient skin biopsies:

- **FFPE scRNA-seq**: 10X Chromium Fixed RNA Profiling (Flex) from formalin-fixed paraffin-embedded biopsies of treated lesional skin (Week 0 and Week 8), across 5 patients treated with resiquimod 0.03% or 0.06%
- **Single-nucleus RNA-seq (snRNA-seq)**: 10X Chromium snRNA-seq from fresh-frozen biopsies of untreated lesional and non-lesional skin (Week 0 and Week 24), across 8 patients multiplexed with souporcell for demultiplexing

After quality control and doublet removal, the merged dataset contains ~71,000 cells spanning 12 cell types: keratinocytes (undifferentiated and differentiated), fibroblasts, T cells, myeloid cells, vascular and lymphatic endothelial cells, pericytes, melanocytes, hair follicle cells, mast cells, and B cells.

## Analysis Pipeline

All notebooks are in `_notebooks/` and are numbered in order of execution.

| Notebook | Language | Description |
|----------|----------|-------------|
| `01_Seurat_ffpe.ipynb` | R | Seurat preprocessing of FFPE data: QC filtering, Scrublet doublet removal (mixture model approach), normalization, Harmony batch correction, clustering, and cell type annotation (27,448 cells) |
| `01_Seurat_nuc-seq.ipynb` | R | Seurat preprocessing of snRNA-seq data: souporcell-based demultiplexing, QC filtering, Scrublet doublet removal, Harmony integration, clustering, and cell type annotation (43,854 cells) |
| `02_run_scrublet.ipynb` | Python | Scrublet doublet scoring for both FFPE and snRNA-seq count matrices |
| `03_merge_datasets_label_cell_types.ipynb` | R | Integration of FFPE and snRNA-seq datasets via Harmony, joint clustering, and consensus cell type labeling across the merged 71,302-cell dataset |
| `04_Figure1_plots.ipynb` | R | Figure 1 visualizations: UMAPs, cell type proportion analysis, marker gene heatmaps/dot plots, and statistical comparison of lesional vs. non-lesional cell composition |
| `05_DE_LesUntxvsNL_Wk0.ipynb` | R | Differential expression between untreated lesional and non-lesional skin at baseline (Week 0), with GSEA pathway analysis per cell type |
| `06_HTS_responder_status.ipynb` | R | Treatment responder classification using high-throughput TCR sequencing: quantifies malignant T-cell clone abundance changes from Week 0 to Week 8/24 |
| `06_HTS_clonal_tracking.ipynb` | R | TCR clonal tracking analysis: benign T-cell recruitment dynamics and correlation with malignant clone clearance |
| `07_myeloid_subclusters_TLR.ipynb` | R | Myeloid cell subclustering (7 subtypes including TissueMoDCs, cDCs, Langerhans, pDCs) and TLR7/8 expression analysis |
| `08_T_subclusters.ipynb` | R | T-cell subclustering into 8 subtypes with marker gene characterization |
| `09_DE_Treated.ipynb` | R | Differential expression in treated lesions: Week 8 vs. Week 0, responders vs. non-responders comparisons, with GSEA pathway enrichment |
| `utils.R` | R | Shared helper functions: library loading, color palettes, hex-bin UMAP plotting, violin plots, and GSEA wrapper |

## Key Dependencies

**R (v4.3.1):** Seurat (v5), Harmony, presto, MAST, clusterProfiler, msigdbr, ggplot2, patchwork, cowplot, mixtools, lme4

**Python:** Scrublet, scanpy

## Data Structure

Notebooks expect data in a sibling `_data/` directory with subdirectories for `ffpe/`, `snRNAseq/`, and clinical metadata (`sample_meta_2026.csv`). Raw input is 10X Cell Ranger output (filtered feature-barcode matrices). Intermediate Seurat `.rds` objects are saved between steps.
