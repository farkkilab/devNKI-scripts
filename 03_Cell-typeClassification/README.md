# dev-NKIscripts: 3_Cell-typeClassification

## Description
Script used for cell type calling and plotting of results. As input it is used the cell quantification after QCs.

### File list

- caller_NKI_Sep-2023.R: Main script for cell type calling of the TMA cores. It uses internally many of the others scripts in this folder.
- Celltype-markers_logic_NKI_2023.R: List of markers to definy cell types. This script is used by Celltype-markers_logic_NKI_2023.R
- 2b_cell_type_caller.R: Script used to scale data and construct FLOWSOM based celltype classification according to the Celltype_Marker logic. This script is used by caller_NKI_Sep-2023.R
- qc_functions.R: Script used to create heatmaps and UMAPs with the cell type labels produced, for internal inspection of results. This script is used internally by caller_NKI_Sep-2023.R.
- cell_type_caller_functions.R: Scripts with functions used by caller_NKI_Sep-2023.R. It contains functions for data trimming, constructing two gaussian mixture models for detection of cancer cells.
  
- WSI_Tribus_run.ipynb: Scripts of Tribus to perform cell type calling on whole slide validation dataset. 
- WSI_Tribus_visualization.ipynb: Scripts to visualize the result of Tribus analysis. Contains main Fig.3k. 
- logic_table8.xlsx: Logic table required to run Tribus on whole slide validation dataset. 
- WSI_Scimap_gating.ipynb: Scripts to perform cell gating on 2 whole slide images. Those 2 images have special staining patterns and Tribus didn't work well on those. 

- Plot_resultant_cell_type_proportions.Rmd: Main script to plot the results of celltype calling. By this script are produced the **Figures 1 b,c,f,g in the manuscript**.
