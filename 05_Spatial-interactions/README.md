# 5_Spatial-interactions


Scripts used to create Recurrent Cellular neighboorhoods (RNC), explore the results, perform statistical stest of the RCN between samples.

## Description of files

- **RCN_inspection.Rmd**  - Current scripts to inspect final RNC using radious of 46px and creating 30 RCNs. In this script are merged the RCN labels,  stored the labels in Napari format and calculated the proportion of RCN by core

- **RCN_survival_analysis.Rmd** - Script to perform the RCN correlations, survival analysis and differences between molecular profiles.
  
- **RCN_exploration_Figure5ki.ipynb** - Script to calculate fold changes of RCN abundances between IDS and PDS. And to calculate fold changes in immune cell composition between RCN10 and RCN09.

- **Spatial_counts.ipynb** - Python script using modified scimap function to count the number of neighbor cells from each population 

- **WSI_Spatial_counts.ipynb**: Python script using modified scimap function to count the number of neighbor cells from each population on whole slide images.

- **Find_neighboors_of_cancer_cells.ipynb*** - Script to get a list of neighboors of cancer cells, for calculating later expression of relevant functional markers of cancer cells neighboors.

- **Scimap_spatial_analysis_clean.ipynb** - Script to calculate the RCNs used in the manuscript, using different metrics (46px, 92px, 138px). 
