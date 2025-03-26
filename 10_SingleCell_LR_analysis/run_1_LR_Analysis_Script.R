library(rmarkdown)

params_list_1 = list(
  treatment = "all",
  cell_idents = c("Tumor", "CD4_T_cell"), # "B_cell" , "FOXP3_CD4_Treg" , "CD8_T_cell" , "CD4_T_cell","macrophages_CD163_high", "macrophages_ITGAX_high" , "DC"
  batches = "default",
  covariates = "publication_patient_code", 
  input_dir = "/run/user/1356082/gvfs/smb-share:server=group3.ad.helsinki.fi,share=h345/afarkkilab/Projects/NKI/Single_Cell_LR_Annalysis/Input/",
  seurate_obj_name = "launonen_cancercell_allcells_allannot_Filtered_MHCII_HIgh_and_low_paired_only.RDS",
  output_dir = "/run/user/1356082/gvfs/smb-share:server=group3.ad.helsinki.fi,share=h345/afarkkilab/Projects/NKI/Single_Cell_LR_Annalysis/output3/"
)

render("2_LR_Analysis.Rmd", params = params_list_1, output_dir = params_list_1$output_dir)

# treatment : "chemo-naive" , IDS , all(for the whole dataset without filtering)


# Cell types availble in "edited_cell_type" in meta data of "launonen_cancercell_allcells_allannot_Filtered.RDS"
# [1] "FOXP3_CD4_Treg"         "Stroma"                 "CD8_T_cell"             "CD4_T_cell"            
# [5] "NK"                     "Tumor"                  "B_cell"                 "macrophages_CD163_high"
# [9] "DC"                     "Other_immune"           "macrophages_ITGAX_high"


#   batches: batch correction
#   covariates:  add patient Ids if there are paired samples

# At least two samples need to be present per group: Otherwise cannot run the DE analysis 
# For the DE analysis,
# only cell types will be considered if there are at least two samples per group with a sufficient number of cells.
# Check the "abundance_info$abund_plot_sample" plot for this information

