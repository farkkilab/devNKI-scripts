library(rmarkdown)

params_list_1 = list(
  treatment = "all",
  cell_idents = c("Tumor", "B_cell" , "FOXP3_CD4_Treg" , "CD8_T_cell" , "CD4_T_cell","macrophages_CD163_high", "macrophages_ITGAX_high" , "DC"),
  batches = "default",
  covariates = "default"
)

render("2_LR_Analysis.Rmd", params = params_list_1, output_dir = "/run/user/1356082/gvfs/smb-share:server=group3.ad.helsinki.fi,share=h345/afarkkilab/Projects/NKI/Single_Cell_LR_Annalysis/output/")

# treatment : "chemo-naive" , IDS , all(for the whole dataset without filtering)


# Cell types availble in "edited_cell_type" in meta data of "launonen_cancercell_allcells_allannot_Filtered.RDS"
# [1] "FOXP3_CD4_Treg"         "Stroma"                 "CD8_T_cell"             "CD4_T_cell"            
# [5] "NK"                     "Tumor"                  "B_cell"                 "macrophages_CD163_high"
# [9] "DC"                     "Other_immune"           "macrophages_ITGAX_high"


#   batches: batch correction
#   covariates:  add patient Ids if there are paired samples


