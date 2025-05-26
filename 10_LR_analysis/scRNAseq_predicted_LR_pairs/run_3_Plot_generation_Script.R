library(rmarkdown)

params_list = list(
  treatment = "all"
)

render("3_Plot_generation.Rmd", params = params_list, output_dir = "/run/user/1356082/gvfs/smb-share:server=group3.ad.helsinki.fi,share=h345/afarkkilab/Projects/NKI/Single_Cell_LR_Annalysis/output/")

# treatment : "chemo-naive" , IDS , all(for the whole dataset without filtering)