
#### Loading libraries ####
setwd("C:/Users/fernpere/")
library(ggplot2)
library(tidyverse)
library(DESeq2)
library(sigQC)
library(survival)
library(survminer)
library(gridExtra)
library(reshape2)



################### Reading counts ###########################

TCGA.bulk.counts <- read.table(file="D:/users/fperez/NKI_TMAs_AF/RNAseq_TCGA/BulkRNA_TCGA_OVA_augmented_star_gene_counts.tsv",
                               sep="\t", header=TRUE)

#Reading de-convoluted expression of EOC
EOC.deconv <- read.table(file="D:/users/fperez/NKI_TMAs_AF/TCGA_PRISM_RNAdeconvolution/EOC_component_clinics.tsv",
                               sep="\t", header=TRUE)

###################### Reading clinical data #######################################

input.folder <- "D:/users/fperez/NKI_TMAs_AF/TCGA_PRISM_RNAdeconvolution/"

clin.dat <- read.csv(file = paste0(input.folder,"TCGA-CDR-SupplementalTableS1_stage-survival.csv"),
                     header = TRUE, row.names = 1)

out.put.folder <- "D:/users/fperez/NKI_TMAs_AF/Analysis_results/05_ML_validations/TCGA_RNA/"
dir.create(out.put.folder)

HRD.tab <- read.table(file=paste0(input.folder,"Categories_Konstantinopoulos_HR_HGSC_germline.csv"),
                      sep=",", header=TRUE)


#Stratify samples according to tumoral molecular profiles
HRD.tab$Molecular.profile <- "HRP"
HRD.tab$Molecular.profile[HRD.tab$categories == "CCNE1 amplification"] <- "CCNE1amp"
HRD.tab$Molecular.profile[HRD.tab$HRDsum >= 54] <- "HRD"
HRD.tab$Molecular.profile[grep("BRCA", HRD.tab$categories)] <- "BRCAloss"
HRD.tab.cat <- HRD.tab[,c("Sample","Molecular.profile")]

# Structuring clinical data
clin.dat <- clin.dat[clin.dat$type == "OV",]
clin.dat$PFI.time <- as.numeric(clin.dat$PFI.time)/30.4
clin.dat$OS.time <- as.numeric(clin.dat$OS.time)/30.4
clin.dat$OS <- as.numeric(clin.dat$OS)
clin.dat$PFI <- as.numeric(clin.dat$PFI)
clin.dat$age <- as.numeric(clin.dat$age_at_initial_pathologic_diagnosis)

clin.dat$age.bin <- NA
clin.dat$age.bin[clin.dat$age < 40] = 0
clin.dat$age.bin[clin.dat$age >= 40 & clin.dat$age < 50] = 1
clin.dat$age.bin[clin.dat$age >= 50 & clin.dat$age < 60] = 2
clin.dat$age.bin[clin.dat$age >= 60 & clin.dat$age < 70] = 3
clin.dat$age.bin[clin.dat$age >= 70 & clin.dat$age < 80] = 4
clin.dat$age.bin[clin.dat$age >= 80] = 5

#Selecting only high patients with high grade tumors 
clin.dat <- clin.dat[clin.dat$histological_grade %in% c("G3","G4"),]

#Ignoring low clinical stages
clin.dat <- clin.dat[!clin.dat$clinical_stage %in% c("[Not Available]", "Stage IA","Stage IB", "Stage IC",
                                                     "Stage IIA","Stage IIB","Stage IIC"),]
clin.dat$bcr_patient_barcode <- gsub("-", ".", clin.dat$bcr_patient_barcode)



############ Normalizing RNA.count.expression using DESeq ############################

#Vector with samples, the variable group is random, will not be really used
samples <- data.frame(Sample = colnames(TCGA.bulk.counts)[-1],
                      group=sample(c("Low","High"),length(colnames(TCGA.bulk.counts)[-1]), replace = TRUE))


#Preparing input for DESeq
m.counts <- as.matrix(TCGA.bulk.counts[,-1])
row.names(m.counts) <- TCGA.bulk.counts[,1]

deseq <- DESeqDataSetFromMatrix(countData = round(m.counts),
                                colData = samples,
                                design = as.formula("~group"))

#Selecting only genes with at least 1 counts per group
deseq <- deseq[rowSums(counts(deseq)) > 1, ]
#Stimating size factor with DESeq
deseq <- estimateSizeFactors(deseq)

#Selecting only rows with at least expression of 1
m.counts <- m.counts[rowSums(m.counts) > 1,]

#Normalizing the expression
norm.counts.per.sample <- as.data.frame(t(m.counts) * deseq$sizeFactor)

rm(TCGA.bulk.counts)
rm(m.counts)
rm(deseq)


########## Merging normalized counts with clinical data ################

norm.counts.clin.prev <- merge(clin.dat, norm.counts.per.sample, by.x="bcr_patient_barcode", by.y="row.names")

HRD.tab.cat$Sample <- gsub("-", ".", HRD.tab.cat$Sample)

norm.counts.clin <- merge(norm.counts.clin.prev, HRD.tab.cat, by.x="bcr_patient_barcode", by.y="Sample")
rm(norm.counts.clin.prev)


################# Preparing dataset for signature dectection #####################

#These are the genes differential expressed using bulk-RNA
tcga.DE.genes <- read.table(file="D:/users/fperez/NKI_TMAs_AF/RNAseq_TCGA/Overexpressed_genes.csv",
                            sep=",", header=TRUE)

#Detecting the first and last gene column in dataset
gene.cols <- which(colnames(norm.counts.clin) %in% c("TSPAN6", "AP006621.6"))

mRNA_expr_matrix = list()
aux <- norm.counts.clin[,-1]
row.names(aux) <- norm.counts.clin[,1]
gene.cols <- gene.cols-1
aux <- aux[,gene.cols[1]:gene.cols[2]]
mRNA_expr_matrix[["TCGA_bulk"]] =  t(aux)

gene_sigs_list = list()

interferon1.genes <- c("IFNA1","IFNA2","IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNA8", "IFNA10", "IFNA13",
                       "IFNA14", "IFNA16", "IFNA17", "IFNA21", "IFNB1", "IFNB1", "IFNW1", "IFNK",
                       "IFNE", "IFNAR1", "IFNAR2", "JAK1", "TYK2", "STAT1", "STAT2", "IRF9")

gene_list_names = c("MHCIIgenes", "MHCIgenes","NanostringDE","Intersected","TCGAbulkDE","Interferon1")
list_of_genes=list(c("HLA-DPA1","HLA-DPA2","HLA-DPA3","HLA-DPB1","HLA-DPB2","HLA-DQA1","HLA-DQA2","HLA-DQB1",
                     "HLA-DQB2","HLA-DQB3","HLA-DRA","HLA-DRB1","HLA-DRB2","HLA_DRB3","HLA_DRB4","HLA-DRB5",
                     "HLA-DRB6","HLA-DRB7","HLA-DRB8","HLA-DRB9","CD74"),
                   c("HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-G","HLA-H","HLA-J","HLA-K","HLA-L",
                     "HLA-N","HLA-P","HLA-S","HLA-T","HLA-U","HLA-V","HLA-W","HLA-X","HLA-Y","HLA-Z"),
                   c("CD74", "CXCL9", "CXCL10", "CXCL11","IDO1", "ADAMDEC1", "SLAMF7", "CCL5", "CD3D"),
                   c("CCL5", "CD74", "CXCL9"),
                   tcga.DE.genes$gene,
                   interferon1.genes)

for (i in 1:length(gene_list_names)){
  genes <- intersect(row.names(mRNA_expr_matrix$TCGA_bulk), list_of_genes[[i]])
  gene_sigs_list[[gene_list_names[i]]] = genes
}

rm(aux)


############ Running program for MHCII gene signature ################
#The result is stored in the out_dir
out.put.dir = 'D:/users/fperez/NKI_TMAs_AF/Analysis_results/05_ML_validations/TCGA_RNA/SignaturesQC_bulk/'

make_all_plots(gene_sigs_list = gene_sigs_list,
               mRNA_expr_matrix = mRNA_expr_matrix,
               doNegativeControl = FALSE,
               out_dir = out.put.dir,
               showResults = FALSE)


#Reading signatures and merging results with dataset
standarization_tables <- list.files(paste0(out.put.dir,"/standardisation_tables"),
                                    pattern = "_TCGA_bulk.txt", full.names = TRUE)

norm.counts.clin2 = norm.counts.clin

for (f in (standarization_tables)){
  signature.name <- strsplit(basename(f), "_")[[1]][3]
  sig <- read.table(file=f, sep="\t", header=TRUE)
  colnames(sig) <- c(paste0(signature.name,".sig"), paste0(signature.name,".sig.scaled"))
  norm.counts.clin2 <- cbind(norm.counts.clin2, sig[norm.counts.clin2$bcr_patient_barcode,])
}



########## Stratifying samples according to MHCII expression and performing survival analysis #############

total.genes.to.test <- c("CD74", "CXCL9", "CXCL10", "CXCL11",
                         "IDO1", "ADAMDEC1", "SLAMF7", "CCL5", "CD3D")

hrp.genes.to.test <- c("CRABP1", "CD2", "SLAMF7", "ADAMDEC1", "CD3E", "CD27", "CTLA4", 
                       "CD3D", "CXCL9", "CXCL10", "CXCL11", "IDO1")

signatures <- paste0(gene_list_names, ".sig")

genes.to.test <- unique(c(hrp.genes.to.test, total.genes.to.test, signatures))

out.put.folder2 <- "D:/users/fperez/NKI_TMAs_AF/Analysis_results/05_ML_validations/TCGA_RNA/Bulk_RNA/"

dir.create(out.put.folder2)


formulas.names <- c("OS","PFI")
formulas.labs <- c("OS probability","PFI probability")
xlim.max.val <- c(150,120)
dat <- norm.counts.clin2

formulas.cox <- c('Surv(OS.time, OS) ~ age.bin + Molecular.profile',
                  'Surv(PFI.time, PFI)~ age.bin + Molecular.profile')

formulas.kapplan <- c('Surv(OS.time, OS) ~  Expression',
                      'Surv(PFI.time, PFI)~ Expression')

for (gene in genes.to.test){
  for (i in 1:length(formulas.kapplan)){
    #For cox model data is stratified by tiles 
    dat$Expression <- ntile(dat[,gene],3)
    colnames(dat)[colnames(dat) == "Expression"] = paste0(gene,"exp")
    
    model <- coxph(as.formula(paste0(formulas.cox[i]," + ", gene,"exp")), data=dat)
    p1 <- ggforest(model, data = dat,  main=formulas.labs[i])
    print(p1)
    ggsave(p1, file=paste0(out.put.folder2, "Bulk_Cox_",gene,"_age_",formulas.names[i],".png"),
           width = 15, height = 10, units = "cm")
    
    #For kapplan data is stratified by median expression 
    dat$Expression <- ifelse(dat[,gene] > median(dat[,gene]), "High","Low")
    
    fit <- survfit(as.formula(formulas.kapplan[i]), data = dat)
    logrank = surv_pvalue(fit, dat)
    pval=logrank$pval.txt
    p <- ggsurvplot(fit, data = dat, risk.table = TRUE, p.val=TRUE, conf.int = TRUE,
                    palette = c("#ff5a36","#4d4dff"),
                    legend.labs=c(paste0(gene,".High"), paste0(gene,".Low")),
                    tables.y.text = FALSE, legend.title="",
                    xlim=c(0,xlim.max.val[i]), break.time.by=30, ylab=formulas.labs[i])
    p$plot <- p$plot + theme(legend.key.width = unit(3, "line"),
                             plot.margin = margin(t = 0,  r = 2, b = 0, l = 2, unit = "cm"),
                             legend.key.height = unit(1.5, "line"))
    p$plot <- p$plot + ggplot2::annotate("text", x = 90, y = 0.85, label =pval, size = 5)
    p$table <- p$table + theme(plot.margin = margin(t = 0,  r = 2, b = 5, l = 2, unit = "cm"))
    plot.to.save <- list(p$plot, p$table)
    ggsave(file=paste0(out.put.folder2,"Bulk_Kapplan_",gene,"_",formulas.names[i],".png"),
           arrangeGrob(grobs = plot.to.save, ncol = 1),
           width = 13, height = 16, units = "cm")
  }
}


##################### Survival using immune expression for HRDs or HRPs ########################
classes <- list(c("BRCAloss","HRD"), 
                c("CCNE1amp","HRP"))

profile.names <- c("HRDs","HRPs")

formulas.names <- c("OS","PFI")
formulas.labs <- c("OS probability","PFI probability")
xlim.max.val <- c(150,120)

formulas.cox <- c('Surv(OS.time, OS) ~ age.bin',
                  'Surv(PFI.time, PFI)~ ')

formulas.kapplan <- c('Surv(OS.time, OS) ~  Expression',
                      'Surv(PFI.time, PFI)~ Expression')

dat2 <- norm.counts.clin2


for (index in 1:length(classes)){
  dat <- dat2[dat2$Molecular.profile %in% classes[[index]],]
  
  for (gene in genes.to.test){
    for (i in 1:length(formulas.kapplan)){
      #For cox model data is stratified by tiles 
      dat$Expression <- ntile(dat[,gene],3)
      colnames(dat)[colnames(dat) == "Expression"] = paste0(gene,"exp")
      
      model <- coxph(as.formula(paste0(formulas.cox[i]," + ", gene,"exp")), data=dat)
      p1 <- ggforest(model, data = dat,  main=formulas.labs[i])
      ggsave(p1, file=paste0(out.put.folder2, "Bulk_Cox_",gene,"_age_",profile.names[index],
                             "_", formulas.names[i],".png"),
             width = 15, height = 10, units = "cm")
      
      #For kapplan data is stratified by median expression 
      dat$Expression <- ifelse(dat[,gene] > median(dat[,gene]), "High","Low")
      
      fit <- survfit(as.formula(formulas.kapplan[i]), data = dat)
      logrank = surv_pvalue(fit, dat)
      pval=logrank$pval.txt
      p <- ggsurvplot(fit, data = dat, risk.table = TRUE, p.val=TRUE, conf.int = TRUE,
                      #palette = c("#ff5a36","#4d4dff"),
                      #legend.labs=c(paste0(gene,".High"), paste0(gene,".Low")),
                      tables.y.text = FALSE, legend.title="",
                      xlim=c(0,xlim.max.val[i]), break.time.by=30, ylab=formulas.labs[i])
      p$plot <- p$plot + theme(legend.key.width = unit(3, "line"),
                               plot.margin = margin(t = 0,  r = 2, b = 0, l = 2, unit = "cm"),
                               legend.key.height = unit(1.5, "line"))
      p$plot <- p$plot + ggplot2::annotate("text", x = 90, y = 0.85, label =pval, size = 5)
      p$table <- p$table + theme(plot.margin = margin(t = 0,  r = 2, b = 5, l = 2, unit = "cm"))
      p$plot <-  p$plot + guides(colour = guide_legend(nrow = 5))
      plot.to.save <- list(p$plot, p$table)
      ggsave(file=paste0(out.put.folder2,"Bulk_Kapplan_",gene,"_", profile.names[index], "_", formulas.names[i],".png"),
             arrangeGrob(grobs = plot.to.save, ncol = 1),
             width = 13, height = 21, units = "cm")
    }
  }
}


################### Exploring the relationship of MHCII expression in cancer cells and immune gen expression
#Remove outlier function
#Merging bulk expression with MHCII deconvoluted expression
EOC.mhcii <- EOC.deconv %>% select(all_of(c("Sample","HLA_DPB1"))) %>% 
                            mutate(MHCII.eoc.strat = ntile(HLA_DPB1, 3))
EOC.mhcii$Sample <- gsub("-",".", EOC.mhcii$Sample)
dat3 <- merge(dat2, EOC.mhcii, by.x="bcr_patient_barcode", by.y="Sample")


#Remove outlier function
out.rem <- function(x){
  iqr.x <- IQR(x)
  quantile3 <- quantile(x, 0.75)
  cut.val <- quantile3 + 1.5 * iqr.x
  x[which(x > cut.val)] <- cut.val
  return(x)
}


total.genes.to.test <- c("CD74", "CXCL9", "CXCL10", "CXCL11",
                         "IDO1", "ADAMDEC1", "SLAMF7", "CCL5", "CD3D")

genes.to.test <- unique(c(total.genes.to.test))

mat <- dat3 %>% select(any_of(c(genes.to.test))) %>% 
  mutate_all(., out.rem) %>% scale() %>% t()


genes.exp.immune <- as.data.frame(t(mat))


df <- data.frame(Sample=dat3$bcr_patient_barcode,  MHCII.eoc.strat=dat3$MHCII.eoc.strat,
                 Molecular.profile=factor(dat3$Molecular.profile, levels=c("BRCAloss","HRD","HRP","CCNE1amp")))
df <- cbind(df, genes.exp.immune)

aux <- df %>% filter(MHCII.eoc.strat != 2) %>%
  mutate(MHCII.eoc = ifelse(MHCII.eoc.strat==1, "Low","High")) %>% 
  select(-MHCII.eoc.strat) %>%
  melt() %>% 
  mutate(Molecular.profile2 = case_when(Molecular.profile %in% c("BRCAloss","HRD") ~ "HRD",
                                        Molecular.profile %in% c("HRP","CCNE1amp") ~ "HRP"))

colnames(aux)[4] <- "Gene"

for (m in unique(aux$Molecular.profile2)){
  print(m)
  for (gen in unique(aux$Gene)){
    d <- aux %>% filter(Gene == gen,
                        Molecular.profile2 == m)  
    w.res  <- wilcox.test(d[d$MHCII.eoc == "High","value"],
                          d[d$MHCII.eoc == "Low","value"])
    print(paste0("Difference between high and low MHCII, for ", gen, " pvalue: ", w.res$p.value))
  }
}

p <- ggplot(aux, aes(x=Gene, y=value, fill=MHCII.eoc)) +
  geom_violin(aes(fill=MHCII.eoc)) +
  geom_boxplot(aes(fill=MHCII.eoc),width=0.15, outlier.shape = NA, position = position_dodge(0.9)) +
  stat_summary(geom = "errorbar", fun = "median", aes(ymax = ..y.., ymin = ..y..), linewidth=0.6, size = 0.8,
               position = position_dodge(0.9)) +
  scale_fill_manual(values = c("yellow3","deepskyblue3")) +
  xlab("") + theme_bw() + ylab("Deconvoluted immune cell expression") +
  facet_wrap(~Molecular.profile2) + ylim(-1.3,2.7) +
  theme(axis.text.x = element_text(angle = 55, hjust=1))
print(p)
ggsave(p, file=paste0(out.put.folder, "Boxplot_gene_expression_MHCII_lowHigh.svg"),
       width = 20, height = 13, units = "cm")


########### TCGA stat activation in EOC #####################

#Merging bulk expression with MHCII deconvoluted expression
EOC.mhcii <- EOC.deconv %>% select(all_of(c("Sample","HLA_DPB1","STAT1","STAT2"))) %>% 
                            mutate(MHCII.eoc.strat = factor(ntile(HLA_DPB1, 3), labels = c("Low","Mid", "High")))
EOC.mhcii$Sample <- gsub("-",".", EOC.mhcii$Sample)
dat3 <- merge(dat2, EOC.mhcii, by.x="bcr_patient_barcode", by.y="Sample")

df <- dat3 %>% select(bcr_patient_barcode, MHCII.eoc.strat, STAT1.x, STAT1.y, Interferon1.sig) 

p <- ggplot(df, aes(x=MHCII.eoc.strat, y=Interferon1.sig)) + geom_boxplot()
print(p)


########## TCGA neo-antigens and immune composition#####################
neoantigens <- read.table(file="D:/users/fperez/NKI_TMAs_AF/TCGA_PRISM_RNAdeconvolution/Neoantigens/TCGA_pMHC_SNV_sampleSummary_MC3_v0.2.8.CONTROLLED_170404.tsv",
           sep="\t", header=TRUE)

cibertsort <- read.table(file="D:/users/fperez/NKI_TMAs_AF/TCGA_PRISM_RNAdeconvolution/Neoantigens/TCGA.Kallisto.fullIDs.cibersort.relative.tsv",
                          sep="\t", header=TRUE)


neoantigens$barcode  <- sapply(neoantigens$barcode, function(x){paste(strsplit(x,"-")[[1]][1:3],collapse = ".")})
cibertsort$SampleID  <- sapply(cibertsort$SampleID, function(x){paste(strsplit(x,"[.]")[[1]][1:3],collapse = ".")})

cibertsort <- cibertsort %>% filter(CancerType=="OV") %>% 
                             select(-P.value, -Correlation, -RMSE, -CancerType) %>% 
                            mutate(T.cells.CD4=(T.cells.CD4.naive + T.cells.CD4.memory.resting + T.cells.CD4.memory.activated),
                                   Dendritic.cells=(Dendritic.cells.resting + Dendritic.cells.activated)) %>% 
                            melt()
colnames(cibertsort)[c(2,3)] <- c("Celltype","Proportion")
cibertsort$Celltype <- as.character(cibertsort$Celltype)


neoantigens.mhcii <- merge(neoantigens, EOC.mhcii, by.x="barcode", by.y="Sample")
cibertsort.mhcii <- merge(cibertsort, EOC.mhcii, by.x="SampleID", by.y="Sample")


cibertsort.mhcii$Celltype[cibertsort.mhcii$Celltype == "T.cells.regulatory..Tregs."] <- "T.cells.regulatory"

cells.interest <- c("T.cells.CD8", "T.cells.CD4", "T.cells.regulatory", 
                     "Macrophages.M1", "Macrophages.M2", "Macrophages.M0")


df <- cibertsort.mhcii %>% filter(Celltype %in% cells.interest)


df$Celltype <- factor(df$Celltype, levels=cells.interest)


my_comparisons <- list( c("Low", "Mid"), c("Mid", "High"), c("Low", "High") )


p <- ggboxplot(df, x = "MHCII.eoc.strat", y = "Proportion") + facet_wrap(~Celltype) +
  theme_bw() + xlab("HLA_DPB1 expression") +
  theme(axis.title = element_text(size=rel(1.2)),
        axis.text = element_text(size=rel(1.2)),
        strip.text = element_text(size=rel(1.15))) + 
  stat_compare_means(comparisons = my_comparisons,
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols=c("****", "***", "**", "*", "ns"))) +
  ylim(0,1) + xlab("HLA_DPB1 expression in EOC component") +
  ylab("CIBERTSORT cell type propotions")
print(p)
ggsave(p,
       file="D:/users/fperez/NKI_TMAs_AF/Analysis_results/05_ML_validations/TCGA_RNA/CIBERSORT_cell_proportions.png",
       width = 13, height = 12, units = "cm")
ggsave(p,
       file="D:/users/fperez/NKI_TMAs_AF/Analysis_results/05_ML_validations/TCGA_RNA/CIBERSORT_cell_proportions.svg",
       width = 13, height = 12, units = "cm")




df <- neoantigens.mhcii %>% select(barcode, MHCII.eoc.strat, numberOfNonSynonymousSNP, numberOfImmunogenicMutation,
                             numberOfBindingPMHC, numberOfBindingExpressedPMHC)
exp.hladpb1 <- neoantigens.mhcii %>% select(barcode, HLA_DPB1)
df <- melt(df)
df <- merge(df, exp.hladpb1)
p<- ggplot(df, aes(x=HLA_DPB1, y=value)) + geom_point() + theme_bw() +
    facet_wrap(~variable) + geom_smooth() + stat_cor()
print(p)
ggsave(p,
       file="D:/users/fperez/NKI_TMAs_AF/Analysis_results/05_ML_validations/TCGA_RNA/Neo_antigens_proportions.svg",
       width = 15, height = 15, units = "cm")


###########Analysis for prediction of survival using CD8.T cells #########

output.folder3 <- "D:/users/fperez/NKI_TMAs_AF/Analysis_results/05_ML_validations/TCGA_RNA/"

cibertsort <- read.table(file="D:/users/fperez/NKI_TMAs_AF/TCGA_PRISM_RNAdeconvolution/Neoantigens/TCGA.Kallisto.fullIDs.cibersort.relative.tsv",
                         sep="\t", header=TRUE)

cibertsort$SampleID  <- sapply(cibertsort$SampleID, function(x){paste(strsplit(x,"[.]")[[1]][1:3],collapse = ".")})


formulas.names <- c("OS","PFI")
formulas.labs <- c("OS probability","PFI probability")
xlim.max.val <- c(150,120)

formulas.cox <- c('Surv(OS.time, OS) ~ age.bin + Molecular.profile',
                  'Surv(PFI.time, PFI)~ age.bin + Molecular.profile')

cibertsort.clin <- merge(cibertsort, clin.dat, by.x="SampleID", by.y="bcr_patient_barcode")
cibertsort.clin2 <- merge(cibertsort.clin, HRD.tab.cat, by.x="SampleID", by.y="Sample")
cibert.mhcii.clin2 <- merge(cibertsort.clin2, EOC.mhcii, by.x="SampleID", by.y="Sample")

dat <- cibert.mhcii.clin2

vars <- c("HLA_DPB1", "T.cells.CD8")

for (variable in vars){
  for (i in 1:length(formulas.cox)){
      dat$Expression <- ntile(dat[,variable],3)
      colnames(dat)[colnames(dat) == "Expression"] = paste0(variable,".strat")
      
      model <- coxph(as.formula(paste0(formulas.cox[i]," + ", paste0(variable,".strat"))), data=dat)
      p1 <- ggforest(model, data = dat,  main=formulas.labs[i])
      print(p1)
      ggsave(file=paste0(output.folder3,"Cox_",variable,"_",formulas.names[i],".png"),
               width = 14, height = 10, units = "cm")
  }
}

dat <- cibert.mhcii.clin2 %>% filter(Molecular.profile == "CCNE1amp" | 
                                       Molecular.profile == "HRP")

formulas.cox <- c('Surv(OS.time, OS) ~ age.bin',
                  'Surv(PFI.time, PFI)~ age.bin')

for (variable in vars){
  for (i in 1:length(formulas.cox)){
    dat$Expression <- ntile(dat[,variable],3)
    colnames(dat)[colnames(dat) == "Expression"] = paste0(variable,".strat")
    
    model <- coxph(as.formula(paste0(formulas.cox[i]," + ", paste0(variable,".strat"))), data=dat)
    p1 <- ggforest(model, data = dat,  main=formulas.labs[i])
    print(p1)
    ggsave(file=paste0(output.folder3,"Cox_",variable,"_",formulas.names[i],"_HRPs.png"),
           width = 14, height = 7, units = "cm")
  }
}

