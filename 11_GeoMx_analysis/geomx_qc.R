# README: script for loading GeoMx data and making a basic QC 
# GeoMX vignette
# https://www.bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html

# define variables --------------------------------------------------------

# vals used for sankey, detection rate plots, 
imp_vars <- c(aoi_segment_var, main_roi_label, main_experimental_condition)

# parameters for removing genes based on LOQ
# TODO adjustment may be needed: 10% for batch 1, 5% for batch2
gene_detect_thr <- 0.05 # segment is removed if <5% of genes > LOQ
# TODO adjustments may be needed - thr is very low bcs we expect high biological variability
segment_detect_rate_thr <- 0.01 # genes are removed if its expr > LOQ in less than 1% of segments

# !! 6 samples in batch1 have HighNTC - proibably contaminated wells
# keep_high_NTC = TRUE means that all AOIs belonging to high NTC wells are keeped
# while processing batch1 separately they're keeped 
keep_high_NTC <- ifelse(batch == 'batch1', TRUE, FALSE)

# create dirs -------------------------------------------------------------

dir.create(file.path(output_dir, 'qc'), showWarnings = T, recursive = T)


# load geomx dataset ------------------------------------------------------

geomx_obj <- readNanoStringGeoMxSet(dccFiles = dcc_path, 
                                    pkcFiles = pkc_path, # this goes into fData() - features (probes) annotation
                                    phenoDataFile = anno_path, # this goes into pData() - protocol (samples) annotation
                                    phenoDataSheet = "Sheet1",
                                    phenoDataDccColName = "Sample_ID",
                                    protocolDataColNames = c("Aoi", "Roi"),
                                    experimentDataColNames = c("Panel")) 

pkcs <- annotation(geomx_obj)
modules <- gsub(".pkc", "", pkcs)

# explore
# View(assayData(geomx_obj)$exprs)
# dim(assayData(geomx_obj)$exprs)
# 
# View(pData(geomx_obj))
# View(pData(protocolData(geomx_obj)))
# View(fData(geomx_obj))
# featureType(geomx_obj)
# annotation(geomx_obj)
# 
# View(summary(geomx_obj, MARGIN = 1)) # for probes
# View(summary(geomx_obj, MARGIN = 2)) # for segments (rois)

print(paste('dim of raw dataset is: '))
print(dim(geomx_obj))

# manualfix of NTC --------------------------------------------------------

sdt <- sData(geomx_obj)

# !!! if theres no 'NTC_ID' column added manually, it assumes that each batch 
# have their own NTC
if(!('NTC_ID' %in% colnames(pData(geomx_obj)))){
  sdt$NTC_ID <- apply(sdt, 1, function(x){
    # one NTC/batch
    ntc <- sdt[[aoi_id]][sdt$Slide_Name == 'No Template Control' & 
                           sdt[[batch_var]] == x[[batch_var]] &
                           sdt[[main_batch_var]] == x[[main_batch_var]]]
    return(ntc)
  })
}

sdt$NTC <- apply(sdt, 1, function(x){
  ntc_cnt <- sdt$DeduplicatedReads[sdt[[aoi_id]] == x[['NTC_ID']]]
})

#add to protocolData
identical(rownames(protocolData(geomx_obj)@data), sdt[[aoi_id]])
protocolData(geomx_obj)@data[, c("NTC_ID", "NTC")] <- sdt[, c("NTC_ID", "NTC")]

# remove NTC_ID from pData to prevent duplicated columns
if('NTC_ID' %in% colnames(pData(geomx_obj))){
  pData(geomx_obj) <- pData(geomx_obj)[ , -which(names(pData(geomx_obj)) == 'NTC_ID')]
}

#change 'Area' and 'Nuclei' colnames for correct qc flags
pData(geomx_obj) <- dplyr::rename(pData(geomx_obj), 'area' = 'Area', 'nuclei' = 'Nuclei')

# make overall sankey plot ------------------------------------------------
count_segments <- geomx_obj@phenoData@data[main_roi_label != 'NA' & !is.na(main_roi_label), ]

plot_sankey(count_segments, imp_vars, main_roi_label, 
            file.path(output_dir, 'qc/sankey_slides.png'))


# set and plot basic qc parameters ----------------------------------------
# Shift 0 counts to one -needed for  Q3 norm (but not 100% sure why)
geomx_obj <- shiftCountsOne(geomx_obj, useDALogic = TRUE)

qc_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 75,    # Minimum % of reads aligned (80%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 1,   # Minimum negative control counts (10, 1 in log scale)
       maxNTCCount = 1000,     # Maximum counts observed in NTC well (1000)
       minNuclei = 20,         # Minimum # of nuclei estimated (100) 
       minArea = 1000)         # Minimum segment area (5000)

# set up qc flags for segments
geomx_obj <- setSegmentQCFlags(geomx_obj, qcCutoffs = qc_params)


# rmv NTC segments
geomx_obj <- geomx_obj[, !(geomx_obj$Slide_Name == 'No Template Control')]

qc_results_segment <- protocolData(geomx_obj)[["QCFlags"]]
qc_summary <- qc_summarize(qc_results_segment)

print('qc summary:')
print(qc_summary)

# plot qc histograms
# duplicated cause you have to iterate trough 2 lists of names
QC_histogram(sData(geomx_obj), "Trimmed (%)", aoi_segment_var, qc_params[["percentTrimmed"]], 
             scale_trans = NULL, file.path(output_dir, 'qc/qc_hist_trim.png'))

QC_histogram(sData(geomx_obj), "Stitched (%)", aoi_segment_var, qc_params[["percentStitched"]],
             scale_trans = NULL, file.path(output_dir, 'qc/qc_hist_stich.png'))

QC_histogram(sData(geomx_obj), "Aligned (%)", aoi_segment_var, qc_params[["percentAligned"]],
             scale_trans = NULL, file.path(output_dir, 'qc/qc_hist_align.png'))

QC_histogram(sData(geomx_obj), "Saturated (%)", aoi_segment_var, qc_params[["percentSaturation"]],
             scale_trans = NULL, file.path(output_dir, 'qc/qc_hist_satur.png'))

QC_histogram(sData(geomx_obj), "area", aoi_segment_var, qc_params[["minArea"]], 
             scale_trans = "log10", file.path(output_dir, 'qc/qc_hist_area.png'))

QC_histogram(sData(geomx_obj), "nuclei", aoi_segment_var, qc_params[["minNuclei"]],
             scale_trans = NULL, file.path(output_dir, 'qc/qc_hist_nuclei.png'))


# negative geometric means ------------------------------------------------

# calculate the negative geometric means for each module
negativeGeoMeans <- 
  esBy(negativeControlSubset(geomx_obj), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }) 

protocolData(geomx_obj)[["NegGeoMean"]] <- negativeGeoMeans

# explicitly copy the Negative geoMeans from sData to pData
# TODO this + detaching may be rmv, check if > 1 module
negCols <- paste0("NegGeoMean_", modules)
pData(geomx_obj)[, negCols] <- sData(geomx_obj)[["NegGeoMean"]]

for(ann in paste0("NegGeoMean_", modules)) {
  QC_histogram(sData(geomx_obj), ann, aoi_segment_var, 2, scale_trans = "log10",
               file.path(output_dir, 'qc/qc_hist_neggeomean.png')) #TODO? why exactly thr = 2?
}

# detatch neg_geomean columns ahead of aggregateCounts call
pData(geomx_obj) <- pData(geomx_obj)[, !colnames(pData(geomx_obj)) %in% negCols]

# background modelling based on negative probes ---------------------------

sum(fData(geomx_obj)$Negative) # same as negativeControlSubset(geomx_obj)

# fit poisson distribution
geomx_obj <- fitPoisBG(geomx_obj)

summary(pData(geomx_obj)$sizefact)
summary(fData(geomx_obj)$featfact[fData(geomx_obj)$Negative])

# diagnose Poisson model
set.seed(123)
geomx_diag <- diagPoisBG(geomx_obj)

notes(geomx_diag)$disper

# If the dispersion is >2, one of these factors might be present in the data. 
# We can check for outlier ROIs. People can choose to set outliers to be missing 
# values and rerun the Poisson Background model.

length(which(assayDataElement(geomx_diag, "low_outlier") == 1, arr.ind = TRUE))
length(which(assayDataElement(geomx_diag, "up_outlier") == 1, arr.ind = TRUE))

# Or if a batch effect is assumed, the poisson model can be adjusted to take 
# different groups into account. Here we are grouping the ROIs by slide.

geomx_obj <- fitPoisBG(geomx_obj, groupvar = 'Slide_Name')

set.seed(123)
geomx_diag <- diagPoisBG(geomx_obj, split = TRUE)

notes(geomx_diag)$disper_sp

# remove flagged segments -------------------------------------------------

table(sData(geomx_obj)$NTC)

# 6 samples in batch1 have very high NTC, probably contaminated. remove them or not?
if(keep_high_NTC){
  qc_results_segment <- qc_results_segment[, -which(names(qc_results_segment) == 'HighNTC')]
}

qc_results_segment$qc_status <- apply(qc_results_segment, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})

geomx_obj <- geomx_obj[, qc_results_segment$qc_status == "PASS", ]

# remove segments with neggeomean < 1.5
geomx_obj <- geomx_obj[, sData(geomx_obj)[["NegGeoMean"]] >= 1.5]

print(paste('dim after removing bad quality segments: '))
print(dim(geomx_obj))

# rmv of probes based on geometric mean and grubbs test -------------------

# the geometric mean of that probe’s counts from all segments divided by the geometric mean 
# of all probe counts representing the target from all segments is less than 0.1
# the probe is an outlier according to the Grubb’s test in at least 20% of the segments
# typically don't change this parameters:

geomx_obj <- setBioProbeQCFlags(geomx_obj, 
                                qcCutoffs = list(minProbeRatio = 0.1,
                                                 percentFailGrubbs = 20), 
                                removeLocalOutliers = TRUE)

qc_results_probe <- fData(geomx_obj)[["QCFlags"]]

# summarise probe qc results
qc_probe_df <- data.frame(Passed = sum(rowSums(qc_results_probe[, -1]) == 0),
                          Global = sum(qc_results_probe$GlobalGrubbsOutlier),
                          Local = sum(rowSums(qc_results_probe[, -2:-1]) > 0
                                      & !qc_results_probe$GlobalGrubbsOutlier))

# retain only probes that passed qc (globally)
geomx_obj <- 
  subset(geomx_obj, 
         fData(geomx_obj)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(geomx_obj)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)

print(paste('dim after removing bad quality probes (globally): '))
print(dim(geomx_obj))

# aggregate probes to features --------------------------------------------

# collapse features to targets
geomx_obj <- aggregateCounts(geomx_obj)

print(paste('dim after collapsing features to targets: '))
print(dim(geomx_obj))

# filter based on LOQ per segment and per gene ----------------------------

# The LOQ is calculated based on the distribution of negative control probes 
# and is intended to approximate the quantifiable limit of gene expression per segment. 
# More stable in larger segments. May not be as accurate in segments with low negative probe counts (ex: <2).
# typically use 2 geometric SD above the geometric mean as the LOQ, which is reasonable for most studies. 
# recommend that a minimum LOQ of 2 be used if the LOQ calculated in a segment is below SD threshold.
# typically don't change this parameters

loq_cutoff <- 2
loq_min <- 2

# Calculate LOQ for each segment
LOQ <- data.frame(row.names = colnames(geomx_obj))

for(module in modules){
  LOQ[, module] <-
    pmax(loq_min,
         pData(geomx_obj)[, paste0("NegGeoMean_", module)] * # coming from aggregate_counts
           pData(geomx_obj)[, paste0("NegGeoSD_", module)] ^ loq_cutoff)
}

pData(geomx_obj)$LOQ <- LOQ

# calculate if expr > LOQ per each gene per segment
LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(geomx_obj)$Module == module
  Mat_i <- t(esApply(geomx_obj[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
# ensure ordering 
LOQ_Mat <- LOQ_Mat[fData(geomx_obj)$TargetName, ]

# Save detection rate information to pheno data
# how many genes have been detected in each segment  
pData(geomx_obj)$GenesDetected <- colSums(LOQ_Mat, na.rm = TRUE)
pData(geomx_obj)$GeneDetectionRate <- pData(geomx_obj)$GenesDetected / nrow(geomx_obj)

print(paste("median gene nr is: ", as.character(median(pData(geomx_obj)$GenesDetected))))
print(paste("mean gene nr is: ", as.character(mean(pData(geomx_obj)$GenesDetected))))
print(paste("median gene detection rate is: ", as.character(median(pData(geomx_obj)$GeneDetectionRate))))


sapply(imp_vars, function(vname){
  plot_detection_rate(pData(geomx_obj), vname, 
                      file.path(output_dir, 'qc', paste0('gene_detect_rate_', vname, '.png')))
})

# filter out segments with too low gene detection rate
geomx_obj <- geomx_obj[, pData(geomx_obj)$GeneDetectionRate >= gene_detect_thr]

print(paste('dim after removing segments based on LOQ: '))
print(dim(geomx_obj))

# save to probe data in how many segments the given gene was detected
LOQ_Mat <- LOQ_Mat[, colnames(geomx_obj)]
fData(geomx_obj)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(geomx_obj)$DetectionRate <- fData(geomx_obj)$DetectedSegments / nrow(pData(geomx_obj))
LOQ_Mat <- LOQ_Mat[fData(geomx_obj)$TargetName, ]

# plot detection rate per gene
plot_gene_detection_rate(fData(geomx_obj), file.path(output_dir, 'qc/gene_detection_rate.png'))

# manually include the negative control probe, for downstream use
negativeProbefData <- subset(fData(geomx_obj), CodeClass == "Negative") # 1 bcs already collapsed to targets
neg_probes <- unique(negativeProbefData$TargetName)

# filter out genes detected > LOQ in less then thr nr of segments (1% for now)
geomx_obj <- 
  geomx_obj[fData(geomx_obj)$DetectionRate >= segment_detect_rate_thr |
              fData(geomx_obj)$TargetName %in% neg_probes, ]

print(paste('dim after removing genes based on LOQ: '))
print(dim(geomx_obj))


# additional background modelling for genes -------------------------------

geomx_obj <- BGScoreTest(geomx_obj)

sum(fData(geomx_obj)[["pvalues"]] < 1e-3, na.rm = TRUE)

# save geomx object after QC ----------------------------------------------

print(paste("median gene nr is: ", as.character(median(pData(geomx_obj)$GenesDetected))))
print(paste("mean gene nr is: ", as.character(mean(pData(geomx_obj)$GenesDetected))))
print(paste("median gene detection rate is: ", as.character(median(pData(geomx_obj)$GeneDetectionRate))))

saveRDS(geomx_obj, file = geomx_qc_path)
