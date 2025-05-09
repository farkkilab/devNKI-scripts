# README: script for removing batch effect using limma and/or harmony
# nice explanation
# https://www.biostars.org/p/366403/

# define variables --------------------------------------------------------

# main experimental conditions
exp_design <- as.formula(paste('~', aoi_segment_var, '+', main_experimental_condition))

# all variables to check for variance
batch_vars <- c(primary_batch_var, secondary_batch_var, other_vars_tech, 
                aoi_segment_var, main_roi_label, sample_name, main_experimental_condition, 
                other_vars_bio)

# normalisation used for batch effect correction calculation
norm_type <- 'deseq2_vst' # best to use vst data, eventually deseq2_norm

# PVCA threshold
pct_threshold <- 0.6 

# make dirs and set additional vars ---------------------------------------

dir.create(file.path(output_dir, 'batch_correction'), showWarnings = T, recursive = T)

norm_is_log <- ifelse(norm_type %in% c('exprs', 'q3_norm', 'deseq2_norm'), FALSE, TRUE)
covname <- ifelse(is.null(cov_design), 'no', gsub(' ', '', as.character(cov_design)[2]))

# load geomx object and create expression set -----------------------------

geomx_obj <- readRDS(geomx_norm_path)

# change vars into factors and remove vars with only 1 unique value
batch_vars_filt <- c()
for(colname in batch_vars){
  if(length(unique(pData(geomx_obj)[[colname]])) > 1){
    pData(geomx_obj)[[paste0(colname, '_factor')]] <- as.factor(pData(geomx_obj)[[colname]])
    batch_vars_filt <- c(batch_vars_filt, colname)
  }
}
batch_factors_names <- paste0(batch_vars_filt, '_factor')

# make expression sets for PVCA
phenoData <- new("AnnotatedDataFrame", data=geomx_obj@phenoData@data, 
                 varMetadata=geomx_obj@phenoData@varMetadata)
featureData <- new("AnnotatedDataFrame", data=geomx_obj@featureData@data, 
                 varMetadata=geomx_obj@featureData@varMetadata)

# !! by default deseq2_norm is used to check for initial batch effect by pvca
exprset_deseq2_norm <- ExpressionSet(assayData=geomx_obj@assayData$deseq2_norm, 
                              phenoData = phenoData,
                              featureData = featureData)

if(!norm_is_log){
  # make log2 transformed normalised counts if norm_type not in log scale
  expr_norm_log <- log2(geomx_obj@assayData[[norm_type]] + 1)
} else{
  expr_norm_log <- geomx_obj@assayData[[norm_type]]
}

# check initial batch effect with PVCA ------------------------------------

pvcaObj_ini <- pvcaBatchAssess(exprset_deseq2_norm, batch_factors_names, pct_threshold) 

plot_pvca(pvcaObj_ini, 'before_correction_deseq2_norm', file.path(output_dir, 'batch_correction'))

# remove batch effect with limma ------------------------------------------

design <- model.matrix(exp_design, data = sData(geomx_obj))
prim_batch <-  sData(geomx_obj)[[primary_batch_var]]

if(!is.null(secondary_batch_var)){
  second_batch <- sData(geomx_obj)[[secondary_batch_var]]
}else{second_batch <- NULL} 

if(!is.null(cov_design)){
  cov <- model.matrix(cov_design, data = sData(geomx_obj))
  }else{cov <- NULL} 

limma_res <- limma::removeBatchEffect(expr_norm_log, batch = prim_batch, batch2 = second_batch,
                                      covariates = cov, design = design)

##########################
# https://portals.broadinstitute.org/harmony/articles/quickstart.html
# remove batch effect with harmony
meta_dt <- pData(geomx_obj)[, unique(c(batch_vars_filt, aoi_segment_var, aoi_id))]

harmony_res <- t(HarmonyMatrix(expr_norm_log, 
                               meta_data = meta_dt,
                               vars_use = c(primary_batch_var, secondary_batch_var)))

# check PVCA after batch effect removal -----------------------------------

exprset_after_limma <- ExpressionSet(assayData=limma_res, 
                                         phenoData = phenoData,
                                         featureData = featureData)

exprset_after_harmony <- ExpressionSet(assayData=harmony_res, 
                                     phenoData = phenoData,
                                     featureData = featureData)


pvcaObj_limma <- pvcaBatchAssess(exprset_after_limma, batch_factors_names, pct_threshold) 
pvcaObj_harmony <- pvcaBatchAssess(exprset_after_harmony, batch_factors_names, pct_threshold) 

plot_pvca(pvcaObj_limma, paste0('after_correction_limma_', primary_batch_var, secondary_batch_var,
                                '_cov_', covname), file.path(output_dir, 'batch_correction'))

plot_pvca(pvcaObj_harmony, paste0('after_correction_harmony_', primary_batch_var, secondary_batch_var), 
          file.path(output_dir, 'batch_correction'))

# add limma and harmony res to umap object --------------------------------

geomx_obj@assayData$limma_batch_corr <- limma_res
geomx_obj@assayData$harmony_batch_corr <- harmony_res

# plot counts distribution ------------------------------------------------

plot_expr_distribution(geomx_obj@assayData$limma_batch_corr, 'limma_batch_corr', 
                       file.path(output_dir, 'batch_correction', 
                                 paste0('expr_hist_limma_batch_corr_', 
                                        primary_batch_var, secondary_batch_var,
                                        '_cov_', covname, '.png')), is_log = T)


plot_expr_distribution(geomx_obj@assayData$harmony_batch_corr, 'harmony_batch_corr', 
                       file.path(output_dir, 'batch_correction',
                                 paste0('expr_hist_harmony_batch_corr_', 
                                        primary_batch_var, secondary_batch_var,
                                        '_cov_', covname, '.png')), is_log = T)

# make UMAP and visualise batch-corrected results -------------------------
# TODO run separately for tumor/stroma (incl umap calculation)

geomx_obj <- make_umap_tsne(geomx_obj, 'limma_batch_corr', assay_is_log = T)
geomx_obj <- make_umap_tsne(geomx_obj, 'harmony_batch_corr', assay_is_log = T)

# generate umap and tsne plots and color by variables
for(corr_type in c('limma_batch_corr', 'harmony_batch_corr')){
  for(method in c('UMAP', 'tSNE')){
    for(color_var in batch_vars_filt){
      print(color_var)
      plot_umap_tsne(pData(geomx_obj), method_type = method, 
                     norm_type = corr_type, color_var = color_var,
                     output_name = file.path(output_dir, 'batch_correction', 
                                             paste0(method, '_', corr_type, '_', color_var, '.pdf')))
    }
  }
}

# save geomx_obj with batch eff correction --------------------------------

saveRDS(geomx_obj, file = geomx_norm_batch_eff_rm_path)

# save logs
writeLines(c('batch effect rmv logs:',
             '; normalisation type : ', norm_type,
             '; primary batch effect variable : ', primary_batch_var,
             '; secondary batch effect variable : ', secondary_batch_var,
             '; limma experimental design : ', as.character(exp_design)[2],
             '; limma covariate : ', as.character(cov_design)[2]), 
           file.path(output_dir, 'batch_correction', 'batch_correction_logs.txt'))
