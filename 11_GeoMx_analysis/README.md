# GeoMX

## Description
Scripts used for analysis of GeoMX data.


## List of scripts

- geomx_for_fer.R: Correlations between abundance of cells and MHCII expression in the tumor cells. Script for making the plots in the manuscript.



### Script: geomx_master.R
This script runs the geomx_processing pipeline and performs pre-processing and analysis steps described below:

Input files:
* **DCC files** - raw expression data for each AOI
* **.pkc file** - provided by GeoMx - information about sequencing probes and library used
* **annotation .xls file** - spreadheed containing all metadata information
* **scRNAseq reference file** in .RDS format - Seurat object with reference scRNAseq dataset, used for deconvolution

Main output files:
* **geomx_qc_norm_batch_eff_rm.RDS** file - GeoMx object containing all the computed expression mtx (raw, normalised, batch corrected)
* **CONDITIONS_expr_mtx_cleaned_vst_harmony/limma_batch_corr.RDS** (name will change based on the parameters used) -
  BayesPrism deconvoluted, normalised and batch corrected cell type specific expression matrices. Contains list of mtx of every cell type
* **bp_res_CONDITIONS_ct_fraction.csv** (names will change based on the parameters used) -
  file with cell type fractions per each AOI computed by BayesPrism
* **ssgsea/_CONDITIONS.csv** (name will change based on the parameters used) - ssGSEA/GSVA scores of provided signatures per each AOI (for full and/or deconvoluted data)

#### pipeline steps:

1. Script: geomx_qc.R
    * filtering AOIs based on basic QC parameters (nr of aligned reads, hight negative control counts)
    * filtering AOIs based on negative probes modelling
    * filtering probes based on geometric mean and grubbs test
    * filtering segments based on LOQ (limit of quantification)
    * calculating gene detection rate  
3. Script: geomx_normalisation.R
    * quantile (Q3) normalisation
    * Deseq2 normalisation
    * vst (variance stabilising transformation) on deseq2 normalised data
    * UMAP and t-SNE projections on all types of normalisation
4. Script: geomx_batch_effect_rmv.R
    * PVCA variance assessment on deseq2 normalised data to find variables responsible for batch effect
    * batch effect correction with harmony
    * PVCA on batch-effect corrected data to assess the correction results
    * UMAP and t-SNE od batch-effect corrected data
5. Script: geomx_deconvolution.R
   * preparing reference scRNAseq dataset (incl removing of low complexity genes)
   * computing cell fractions and ct-specific expression profiles with BayesPrism
   * post-processing of ct-specific expression profiles
       * removing non-reliable predictions
       * vst
       * removing genes with variance = 0 (the same value across all AOIs - artifact of BP + vst)
       * harmony batch effect correction
6. Script: geomx_pathway_analysis.R (on full and/or deconvoluted signal)
   *  removing low complexity genes
   *  calculating GSEA/GSVA for signatures from selected categories of msigdb database
