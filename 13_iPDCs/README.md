# `devNKI-scripts:13_IPDCS`

## Overview

This folder contains scripts related to **Figure 7 and Supplementary Figure 7** of the associated publication. The primary objective of this analysis is to investigate the **role of MHC class II (MHCII) expression in tumor samples on the efficacy of Immunotherapy (ICB)**, specifically focusing on T cell responses in different experimental setups.

The validation is carried out through two distinct studies:

1.  **Analysis of iPDCs:** Investigating the association between high tumor MHCII expression and CD8 T cell response in iPDC samples (using data from Launonen et al., 2024).
2.  **MHCII Blockade:** Assessing the effect of anti-HLA (MHCII-blocking) treatment on T cell and cytokine responses in tumor samples from OncosysOva cohort.

## File Structure

```
.
├── Launonen2024_iPDCs_flow_cytometry/
│   └── Analysis_of_flow_cytometry.Rmd  # Analyses iPDC T-cell response in flow cytometry
└── MHCII_blockade/
    ├── Cytokines_analysis.R             # LegendPlex cytokine analysis and heatmap
    └── Flow_cytometry_analysis/         # Detailed flow cytometry analysis scripts
        ├── Heatmap.R                    # Generates flow cytometry heatmaps across immmune cells
        ├── Statistics_CD45+.R           # Calculates log2FC statistics for all CD45+ cells
        ├── Statistics_CD3+.R            # Calculates log2FC statistics for all CD3+ T cells
        ├── Statistics_CD3+CD4+.R        # Calculates log2FC statistics for CD4+ T cells
        ├── Statistics_CD3+CD4-.R        # Calculates log2FC statistics for CD4- T cells
        ├── Statistics_CD3+CD8+.R        # Calculates log2FC statistics for CD8+ T cells
        └── Statistics_CD3+CD8-.R        # Calculates log2FC statistics for CD8- T cells
```

## Detailed Analysis Descriptions

-----

### Study 1: MHCII-High Tumors and T-cell Response in iPDCs

  * **Data Source:** Launonen et al., 2024 (10) dataset, focusing on tumor samples categorized as MHCII-high.
  * **Script:** `Launonen2024_iPDCs_flow_cytometry/Analysis_of_flow_cytometry.Rmd`
  * **Input Data:** Flow cytometry gates (singlets/doublets) for cancer, immune, and myeloid populations, and pre-calculated CD8 T cell response statistics (GrzB).
  * **Analysis:** Analyzes the association between tumor MHCII status and the immunotherapy response of CD8 T cells in iPDC samples.
  * **Outputs (Publication Figures):**
      * **Figure 7b:** Heatmap visualizing the associations.
      * **Figure 7c:** Dot plot of T cell response.
      * **Supplementary Figure 7a, b:** Violin plots used to define MHCII-high and low tumor samples.

-----

### Study 2: MHCII Blockade in Primary Tumor Samples

  * **Sample Source:** Frozen tumor-derived samples from **three HGSC patients** with high HLA-DPB1 expression (ONCOSYS-Ova prospective cohort).
  * **Experimental Conditions:**
    1.  **Control** (Untreated)
    2.  **Pembrolizumab** (ICB)
    3.  **Anti-HLA + Pembrolizumab** (MHCII Blockade + ICB)

#### 2A. Flow Cytometry Analysis (T-cell Function)

  * **Scripts:** Files in `MHCII_blockade/Flow_cytometry_analysis/`
  * **Analysis Goal:** Quantify the IFN-$\gamma$ and Granzyme B response (as Median Fluorescent Intensity) in various gated immune populations:
      * CD45+, CD3+, CD3+CD4+, CD3+CD4-, CD3+CD8+, CD3+CD8-
  * **Process:**
    1.  MFI is calculated for IFN-$\gamma$ and Granzyme B for each cell type and condition per patient.
    2.  Three $\log_2$ fold-change ($\log_2\text{FC}$) comparisons are calculated per patient/cell type:
          * Pembrolizumab vs. Control
          * Anti-HLA + Pembrolizumab vs. Pembrolizumab
          * Anti-HLA + Pembrolizumab vs. Control
    3.  The $\log_2\text{FC}$ values are **pooled** by calculating the mean across the three patients (e.g., in `Heatmap.R`).
  * **Outputs (Publication Figures):**
      * **Figure 7d, e:** Heatmaps showing the pooled mean $\log_2\text{FC}$ for IFN-$\gamma$ and Granzyme B.
      * **Supplementary Figure 7c, d:** Violin plots generated from the patient-specific $\log_2\text{FC}$ statistics (e.g., from `Statistics_CD3+CD4+.R`).

#### 2B. Cytokine Profiling (LEGENDplex)

  * **Script:** `MHCII_blockade/Cytokines_analysis.R`
  * **Analysis Goal:** Measure concentrations of various cytokines in the sample media (LEGENDplex).
  * **Process:** Compares cytokine concentrations across the three conditions (Control, Pembrolizumab, Anti-HLA + Pembrolizumab) and calculates the mean $\log_2\text{FC}$ across patients for the same three comparisons as above.
  * **Outputs (Publication Figures):**
      * **Figure 7f:** Heatmap showing the pooled mean $\log_2\text{FC}$ for cytokine concentrations.
