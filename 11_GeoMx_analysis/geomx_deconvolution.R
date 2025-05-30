# README: script for performing deconvolution using SpatialDecon (output: cell_types fractions)
# and BayesPrism (output: cell types fractions and cell-type-specific transcriptional profiles
# which undergo vst normalisation and batch effect correction)

# recommended usage for scRNAseq reference dataset is raw counts
# although not log transofrmation of both sc and bulk is also ok

# define variables --------------------------------------------------------

norm_type <- 'q3_norm' # quantile is best for sd. cannot be in the log scale! (bp works on raw counts)
ct_nr_thr <- 45 # best 45 for batch1 and 2 - rmv cell states lower than thr in scrnaseq

tumor_ct_name <- 'Epithelial cells' # tumor ct label in scrna_anno
adjust_synonym_gene_names <- F # whether or not to adjust synonymical gene names between scRNAsea and GeoMX
# that help rescue typically around 300 genes with synonym names, but sometimes Ensembl not work

# variables to merge the final csv with
meta_names <- c(aoi_id, roi_id, aoi_segment_var, sample_name, main_experimental_condition, 
                main_roi_label, other_vars_bio)

# main experimental conditions for limma batch eff rmv
exp_design <- as.formula(paste('~', aoi_segment_var, '+', main_experimental_condition))

# make dirs and set additional vars ---------------------------------------

dir.create(file.path(output_dir, 'deconvolution'), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'deconvolution', 'spatial_decon', scrna_anno), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'deconvolution', 'bayes_prism', scrna_anno), showWarnings = T, recursive = T)
dir.create(file.path(output_dir, 'deconvolution', 'bayes_prism', scrna_anno, 'hist'), showWarnings = T, recursive = T)

norm_is_log <- ifelse(norm_type %in% c('exprs', 'q3_norm', 'deseq2_norm'), FALSE, TRUE)
if(norm_is_log){stop('norm_type cannot be in the log scale! use "q3_norm" or "deseq2_norm"')}

covname <- ifelse(is.null(cov_design), 'no', gsub(' ', '', as.character(cov_design)[2]))

scrna_ref_cleaned_path <- file.path(output_dir, 'deconvolution', gsub('.RDS', '_cleaned_for_deconv.RDS', basename(scrna_ref_path)))

# prepare scrnaseq reference dataset --------------------------------------

if(!file.exists(scrna_ref_cleaned_path)){
  geomx_obj <- readRDS(geomx_norm_batch_eff_rm_path)
  scrna_ref_obj <- readRDS(scrna_ref_path)
  
  
  if(adjust_synonym_gene_names){
    # repair synonymuous gene names
    length(rownames(geomx_obj@assayData$exprs))
    length(rownames(scrna_ref_obj@assays$RNA@data))
    length(intersect(rownames(geomx_obj@assayData$exprs), rownames(scrna_ref_obj@assays$RNA@data)))
    
    adjusted_genes_rna <- adjust_synonym_genes(rownames(geomx_obj@assayData$exprs), rownames(scrna_ref_obj@assays$RNA@data))
    
    #  make a new assay with renamed genes
    RNA_common_genes <- scrna_ref_obj@assays$RNA
    RNA_common_genes@counts@Dimnames[[1]] <- adjusted_genes_rna
    RNA_common_genes@data@Dimnames[[1]] <- adjusted_genes_rna
    scrna_ref_obj@assays$RNA_common_genes <- RNA_common_genes
    
    length(intersect(rownames(geomx_obj@assayData$exprs), rownames(scrna_ref_obj@assays$RNA@data)))
    length(intersect(rownames(geomx_obj@assayData$exprs), rownames(scrna_ref_obj@assays$RNA_common_genes@data)))
    
    rna_mtx_touse <- 'RNA_common_genes'
  } else{
    rna_mtx_touse <- 'RNA'
  }
  
  # clean cell labels
  scrna_ref_obj@meta.data$cell_type <- ifelse(scrna_ref_obj@meta.data$cell_type == tumor_ct_name, 
                                              'tumor', scrna_ref_obj@meta.data$cell_type)
  
  # cell states - clustering tumor cells by patient
  scrna_ref_obj@meta.data$cell_state <- ifelse(scrna_ref_obj@meta.data$cell_type == 'tumor', 
                                               paste0('tumor_', scrna_ref_obj@meta.data$patient), 
                                               scrna_ref_obj@meta.data$cell_type)
  
  ########################################
  # QC of cell states
  
  # plot.cor.phi (input=t(scrna_ref_obj@assays$RNA@data),
  #               input.labels=scrna_ref_obj@meta.data$cell_state,
  #               title="cell state correlation",
  #               #specify pdf.prefix if need to output to pdf
  #               #pdf.prefix="gbm.cor.cs",
  #               cexRow=0.6, cexCol=0.6,
  #               margins=c(6,6))
  # 
  # dev.off()
  # 
  # plot.cor.phi (input=t(scrna_ref_obj@assays$RNA@data),
  #               input.labels=scrna_ref_obj@meta.data$mid_lvl_ct,
  #               title="cell type correlation",
  #               #specify pdf.prefix if need to output to pdf
  #               #pdf.prefix="gbm.cor.ct",
  #               cexRow=0.5, cexCol=0.5,
  # )
  # 
  # dev.off()
  #################################
  
  # check genes outliers
  scrna_stat <- plot.scRNA.outlier(
    input=t(scrna_ref_obj@assays[[rna_mtx_touse]]@data), #make sure the colnames are gene symbol or ENSMEBL ID
    cell.type.labels=scrna_ref_obj@meta.data$cell_type,
    species="hs", 
    return.raw=TRUE, #return the data used for plotting.
    pdf.prefix= gsub('.RDS', '', scrna_ref_cleaned_path) # specify pdf.prefix if need to output to pdf
  )
  
  
  # filter out outlier genes
  scrna_filt <- cleanup.genes (input=t(scrna_ref_obj@assays[[rna_mtx_touse]]@data),
                               input.type="count.matrix",
                               species="hs", 
                               gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                               exp.cells=5)
  
  dim(t(scrna_ref_obj@assays[[rna_mtx_touse]]@data))
  dim(scrna_filt)
  
  # geomx doesn't have to be filtered since later on they took only intersection of genes
  
  # subset to protein coding genes
  scrna_filt_pc <-  select.gene.type(scrna_filt, gene.type = "protein_coding")
  
  #  make a new assay with filtered genes
  RNA_filt_pc <- scrna_ref_obj@assays[[rna_mtx_touse]]
  RNA_filt_pc@counts <- RNA_filt_pc@counts[rownames(RNA_filt_pc@counts) %in% colnames(scrna_filt_pc),  ]
  RNA_filt_pc@data <- RNA_filt_pc@data[rownames(RNA_filt_pc@data) %in% colnames(scrna_filt_pc),  ]
  scrna_ref_obj@assays[[paste0(rna_mtx_touse, '_filt_pc')]] <- RNA_filt_pc
  
  # save adjusted scRNAseq file
  saveRDS(scrna_ref_obj, file = scrna_ref_cleaned_path)
  
  
} else{
  print(paste0("cleaned reference scRNAseq dataset made from ", scrna_ref_path,
  " already exists under the path: ", scrna_ref_cleaned_path))
}



# deconvolution by bayesprism ---------------------------------------------
# https://github.com/Danko-Lab/BayesPrism/blob/main/tutorial_deconvolution.html

# load cleaned scrna and geomx
geomx_obj <- readRDS(geomx_norm_batch_eff_rm_path)
scrna_ref_obj <- readRDS(scrna_ref_cleaned_path)

# remove cells from cell states with < thr cells in ref scrnaseq
ct_freq <- as.data.frame(table(scrna_ref_obj@meta.data$cell_state))
low_ct_cells <- scrna_ref_obj@meta.data$cell_name[scrna_ref_obj@meta.data$cell_state %in% 
                                                    as.character(ct_freq$Var1[ct_freq$Freq < ct_nr_thr])]

scrna_ref_obj <- scrna_ref_obj[, !colnames(scrna_ref_obj) %in% low_ct_cells]

# make a prism object
prism_obj <- new.prism(
  reference=t(scrna_ref_obj@assays$RNA_filt_pc@data), 
  mixture=t(geomx_obj@assayData$exprs),
  input.type="count.matrix", 
  cell.type.labels = scrna_ref_obj@meta.data[[scrna_anno]], 
  cell.state.labels = scrna_ref_obj@meta.data$cell_state,
  key="tumor",
  outlier.cut=0.01,
  outlier.fraction=0.1,
)

# run bayesprism
bprism_res <- run.prism(prism = prism_obj, n.cores = detectCores() -2)

# save res
saveRDS(bprism_res, file = file.path(output_dir,'deconvolution', 'bayes_prism', 
                                     paste0('bp_res_', scrna_anno, '.RDS')))


# extract and save ct fractions  ------------------------------------------

cell_frac_cv <- as.data.frame(bprism_res@posterior.theta_f@theta.cv)
# mask ct_frac results if cv > 0.2-0.5 (0.1 thr for bulk, 0.5 for Visium, GeoMx should be in the middle)
# histogram from batch 1 suggests 0.2 as thr

ct_frac <- get.fraction (bp=bprism_res,
                         which.theta="final",
                         state.or.type="type")

ct_names <- colnames(ct_frac)

# mask  unreliable results
ct_frac[cell_frac_cv > 0.2] <- NA

ct_frac <- rownames_to_column(as.data.frame(ct_frac), 'dcc_filename')
ct_frac <- left_join(ct_frac, sData(geomx_obj)[, meta_names],
                     by = 'dcc_filename')

fwrite(ct_frac, file.path(output_dir,'deconvolution', 'bayes_prism', 
                          paste0('bp_res_', scrna_anno, '_ct_fraction.csv')))


# normalise deconvolution expr mtx  ---------------------------------------

deconv_ct_list <- lapply(ct_names, function(ct_name){
  print(ct_name)
  cell_frac_cv <- as.data.frame(bprism_res@posterior.theta_f@theta.cv)
  # mask ct_frac results if cv > 0.2-0.5 (0.1 thr for bulk, 0.5 for Visium, GeoMx should be in the middle)
  # histogram from batch 1 suggests 0.2 as thr
  cell_to_rm <- rownames(cell_frac_cv)[cell_frac_cv[[ct_name]] > 0.2]
  
  deconv_ct <- BayesPrism::get.exp(bp=bprism_res,
                                   state.or.type="type",
                                   cell.name=ct_name)
  
  deconv_ct_cleaned <- deconv_ct[!(rownames(deconv_ct) %in% cell_to_rm), ]
  
  plot_expr_distribution(deconv_ct, paste0(ct_name, '_raw'), 
                         file.path(output_dir, 'deconvolution', 'bayes_prism', 
                                   scrna_anno, 'hist', paste0('expr_hist_', ct_name, '_raw.png')), is_log = F)
  
  deconv_ct_cleaned_vst <- tryCatch({
    # do vst normalisation
    deconv_ct_cleaned_vst <- varianceStabilizingTransformation(round(t(deconv_ct_cleaned)))

    # remove genes with 0 variance across whole dataset (artifact from deconv + vst)
    per_gene_variance <- apply(deconv_ct_cleaned_vst, 1, var)
    
    genes_var0 <- names(per_gene_variance)[which(per_gene_variance == 0)]

    deconv_ct_cleaned_vst <- deconv_ct_cleaned_vst[!(rownames(deconv_ct_cleaned_vst) %in% genes_var0),]

    
    plot_expr_distribution(deconv_ct_cleaned_vst, paste0(ct_name, '_vst'), 
                           file.path(output_dir, 'deconvolution', 'bayes_prism', 
                                     scrna_anno, 'hist', paste0('expr_hist_', ct_name, '_vst_0var_rmv.png')), is_log = T)

    return(ct_name = deconv_ct_cleaned_vst)  # Return the result 
  }, error = function(e) {
    print('not enough AOIs with trustable predictions to perform vst. cell type is removed')
    return()
  })
  
  
})

names(deconv_ct_list) <- ct_names

# clean list from ct for which vst was not computed 
deconv_ct_list[sapply(deconv_ct_list, is.null)] <- NULL

saveRDS(deconv_ct_list, file = file.path(output_dir,'deconvolution', 'bayes_prism', 
                                     paste0('bp_res_', scrna_anno, '_expr_mtx_cleaned_vst.RDS')))


# do batch effect correction ----------------------------------------------

# make metadata for batch effect correction
meta_dt <- sData(geomx_obj)[, c(meta_names, primary_batch_var, secondary_batch_var)]

# do batch effect removal with harmony
deconv_batch_rm_harm_list <- lapply(names(deconv_ct_list), function(ct_name){
  print(ct_name)
  deconv_vst <- deconv_ct_list[[ct_name]]
  
  # filter meta if some ROI does not contain given ct
  mata_dt_ct <- meta_dt[meta_dt[[aoi_id]] %in% colnames(deconv_vst), ]
  
  deconv_harmony_res <- t(HarmonyMatrix(deconv_vst, 
                                        meta_data = mata_dt_ct,
                                        vars_use = c(primary_batch_var, secondary_batch_var)))
  
  plot_expr_distribution(deconv_harmony_res, paste0(ct_name, '_harmony_corr'), 
                         file.path(output_dir, 'deconvolution', 'bayes_prism', 
                                   scrna_anno, 'hist', paste0('expr_hist_', ct_name, '_harmony_corr.png')), is_log = T)
  
  return(deconv_harmony_res)
})

names(deconv_batch_rm_harm_list) <- names(deconv_ct_list)


saveRDS(deconv_batch_rm_harm_list, file = file.path(output_dir,'deconvolution', 'bayes_prism', 
                                                    paste0('bp_res_', scrna_anno, '_expr_mtx_cleaned_vst_harmony_batch_corr.RDS')))

# remove batch effect with limma

deconv_batch_rm_limma_list <- lapply(names(deconv_ct_list), function(ct_name){
  print(ct_name)
  deconv_vst <- deconv_ct_list[[ct_name]]
  
  # filter meta if some ROI does not contain given ct
  mata_dt_ct <- meta_dt[meta_dt[[aoi_id]] %in% colnames(deconv_vst), ]
  
  design <- model.matrix(exp_design, data = mata_dt_ct)
  
  prim_batch <-  mata_dt_ct[[primary_batch_var]]
  
  if(!is.null(secondary_batch_var)){
    second_batch <- mata_dt_ct[[secondary_batch_var]]
  }else{second_batch <- NULL} 
  
  if(!is.null(cov_design)){
    cov <- model.matrix(cov_design, data = mata_dt_ct)
  }else{cov <- NULL} 

  
  deconv_limma_res <- limma::removeBatchEffect(deconv_vst, batch = prim_batch, batch2 = second_batch,
                                        covariates = cov, design = design)
  

  
  plot_expr_distribution(deconv_limma_res, paste0(ct_name, '_limma_corr'), 
                         file.path(output_dir, 'deconvolution', 'bayes_prism', 
                                   scrna_anno, 'hist',
                                   paste0('expr_hist_', ct_name, '_limma_batch_corr_', 
                                          primary_batch_var, secondary_batch_var,
                                          '_cov_', covname, '.png')), is_log = T)


  return(deconv_limma_res)
})

names(deconv_batch_rm_limma_list) <- names(deconv_ct_list)

saveRDS(deconv_batch_rm_limma_list, file = file.path(output_dir,'deconvolution', 'bayes_prism', 
                                                     paste0('bp_res_', scrna_anno, '_expr_mtx_cleaned_vst_limma_batch_corr_', 
                                                            primary_batch_var, secondary_batch_var,
                                                            '_cov_', covname, '.RDS')))

# write logs --------------------------------------------------------------

# save logs
writeLines(c('deconvolution logs:',
             '; spatial decon normalisation type : ', norm_type,
             '; minimum cell type number : ', ct_nr_thr,
             '; cell type annotation  : ', scrna_anno,
             '; limma primary batch effect variable : ', primary_batch_var,
             '; limma secondary batch effect variable : ', secondary_batch_var,
             '; limma experimental design : ', as.character(exp_design)[2],
             '; limma covariate : ', as.character(cov_design)[2]), deconv_logs_path)
