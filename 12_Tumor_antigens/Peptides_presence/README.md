# MHCII immunopeptidome analysis

## Overview

This folder provides code for the analysis of an immunopeptidomic dataset. The pipeline identifies and visualizes peptides mapped to presented antigen genes and performs pathway enrichment analysis on genes uniquely represented in the MHCII tumor-high peptidome.

### Maintainer
**Aleksandra Shabanova**

## Repository Structure

- **`01_Peptide_presence_analysis.ipynb`**  
  - Identifies peptides mapped to significantly correlated antigen genes across TCGA, scRNA-seq, and GeoMx cohorts tumor fraction.
  - Generates a heatmap showing the number of peptides per gene per patient.
  - Selects genes uniquely mapped to the MHCII peptidome of four tumor-high samples and saves them to a `shared_genes_only_high_based_on_prop.txt` file.

- **`02_Enrichment.Rmd`**  
  - Performs pathway enrichment analysis on the selected genes using the `enrichR` package.
  - Outputs the top 10 overrepresented pathways from chosen databases.

## Usage

### Prerequisites

- **Software**
  - Python
  - R
  - Required R package: `enrichR`

### Quick Start

1. **Peptide Mapping and Visualization**  
   Run `01_Peptide_presence_analysis.ipynb` to generate:
   - A heatmap (e.g., Fig. 6f) showing peptide presence by gene and patient.
   - A `shared_genes_only_high_based_on_prop.txt` file with genes uniquely found in the MHCII tumor-high peptidome.

2. **Pathway Enrichment**  
   Run `02_Enrichment.Rmd` using the `shared_genes_only_high_based_on_prop.txt` file generated in step 1 to obtain enriched pathways.