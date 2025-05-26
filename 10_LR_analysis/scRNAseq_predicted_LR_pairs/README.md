# Predicted Ligand-Receptor Pairs in scRNA-seq Cohort

## Overview

This repository provides code for predicting ligand-receptor (LR) interactions from single-cell RNA sequencing (scRNA-seq) data using **MultiNicheNetR**. The analysis pipeline includes data preprocessing, LR interaction prediction, and result visualization comparing samples with high and low cancer MHCII U cell signatures.

### Maintainer
**Gayani Anandagoda**

## Repository Structure

- **`MHCII_Signature_generation.Rmd`**  
  Generates plots to classify samples based on cancer MHCII U cell signature.

- **`1.1_Data_preprocessing.Rmd`**  
  Preprocesses the scRNA-seq data and outputs a filtered Seurat object for downstream analysis.

- **`2.2_LR_Analysis.Rmd`**  
  Predicts ligand-receptor interactions using the filtered Seurat object. Outputs a prioritized table of LR pairs.

- **`3.3_Plot_generation.Rmd`**  
  Visualizes LR predictions. Outputs include circos and bubble plots for selected LR pairs.

- **Wrapper Scripts**  
  - `run_1_LR_Analysis_Script.R`: Automates `2.2_LR_Analysis.Rmd` for specific cell types.  
  - `run_3_Plot_generation_Script.R`: Automates `3.3_Plot_generation.Rmd` for selected cell types.

- **`4_Analysis_of_PrioritisationTbls.Rmd`**  
  Analyzes top ligand-receptor pairs and visualizes differences in median scaled LR product between MHCII-high and MHCII-low samples (e.g., Fig. 5h in the associated publication).

## Usage

### Prerequisites

- **Software**  
  - R & RStudio  
  - Required packages: `Seurat`, `MultiNicheNetR`, `ggplot2`, `circlize`, and others (see Rmd files for specific dependencies)

### Quick Start

1. **Preprocess Data**  
   Run `1.1_Data_preprocessing.Rmd` to generate a filtered Seurat object.

2. **Predict LR Interactions**  
   Use `run_1_LR_Analysis_Script.R` to execute the LR prediction workflow (`2.2_LR_Analysis.Rmd`) for selected cell types.

3. **Generate Plots**  
   Use `run_3_Plot_generation_Script.R` to create visualizations (circos/bubble plots) from the prediction output.

4. **Advanced Analysis**  
   Run `4_Analysis_of_PrioritisationTbls.Rmd` to compare LR interactions between MHCII-high and MHCII-low cancer samples for selected immune cells and tumor pairs
