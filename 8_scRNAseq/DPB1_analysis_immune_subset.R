######## Tasks:
# plot interferan alfa beta gama expression
# separate MHC2 high low
# separate PDS and IDS
# maximize separation of int gama vs alpha or beta
# porportion of M1 and M2 in MHC2 low vs high
# DEGs of EOCs intrapatient of MHC2 high vs low
# Unos comentarios extra para alguno de estos puntos:
# "plot interferan alfa beta gama expression" - Pero quisieramos ver correlacion entre MCHII de cancer cells, 
#                                               con interferon gamma (de immune cells) e Interferon I (alpha/beta, de cancer cells).(and also gene expression)

# "Separate MHCII high low" - Puedes usar la media de tu cohort, o dividir en 3 grupos.
# "Proportion of M1 and M2 in MHC2 low vs high" - Si fuera posible, tambien de CD4+T.cells y CD4+T.regulatory.cells
# 
# -Tengo duda de estos puntos:
#   "DEGs of EOCs intrapatient of MHC2 high vs low" - A que se refiere intra-patient? es factible o correcto tambien hacer inter-patient?
#   "# maximize separation of int gama vs alfa o beta" - No entiendo bien este punto.
# 
# En general todos los analysis si deben ser separados por PDS e IDS. 
# Si vemos el mismo trend, en ambos, se podria hacer un merge de ambos PDS e IDS para aumentar significancia en algunos analisi

## Matías M. Falco
####Finished: 16/03/2024
setwd("/shared/matiasfa/DPB1_fer/")
library(Seurat)
library(ggplot2)
library(UCell)
library(ggpubr)
library(patchwork)
library(ggforce)
library(dplyr)
library(tidyr)
library(data.table)
library(corrplot)
library(RColorBrewer)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(edgeR)
library(limma)
library(svglite)
library(EnvStats)
# mhc2_genes<-c("HLA-DPA1","HLA-DPA2","HLA-DPA3","HLA-DPB1","HLA-DPB2","HLA-DQA1","HLA-DQA2","HLA-DQB1",
#               "HLA-DQB2","HLA-DQB3","HLA-DRA","HLA-DRB1","HLA-DRB2","HLA-DRB3","HLA-DRB4","HLA-DRB5",
#               "HLA-DRB6","HLA-DRB7","HLA-DRB8","HLA-DRB9")


#read list of published data
published<-fread("../../resources/Published_samples_240307.txt",data.table = F)

#read gene signatures
int_gamma<-fread("HALLMARK_INTERFERON_GAMMA_RESPONSE.v2024.1.Hs.grp",skip = 1,data.table = F)[,1]
int_alpha<-fread("HALLMARK_INTERFERON_ALPHA_RESPONSE.v2024.1.Hs.grp",skip = 1,data.table = F)[,1]
mhc2_genes<-c("HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQA2","HLA-DQB1","HLA-DRA","HLA-DRB1","HLA-DRB5","CD74")

#read whole dataset and subset published data (this subset actually has no overlap with the one published for the immune subset)
# exp<-readRDS("../all_samples_preproc/preprocessed_data/norm_filt_data_all_samples_seurat_obj_v2.Rds")
# dim(exp)
# exp<-subset(exp,subset= !((PAX8>1&DCN>1)|(PAX8>1&PTPRC>1)|(PTPRC>1&DCN>1)))#filter cells expressing double markers
# dim(exp)

#read info from cancer cell paper to filter out cells and add metainfo
# # can_cell_cd8<-readRDS("/shared/resources/launonen_cancercell/launonen_cancercell_cd8tcell_ProjecTILs.RDS")
# # write.table(can_cell_cd8@meta.data,"/shared/resources/launonen_cancercell/launonen_cancercell_cd8tcell_ProjecTILs_metadata.tsv",sep = "\t",col.names = T,row.names = F,quote = F)
# # can_cell_macrophages<-readRDS("/shared/resources/launonen_cancercell/launonen_cancercell_macrophage_subtyped.RDS")
# # write.table(can_cell_macrophages@meta.data,"/shared/resources/launonen_cancercell/launonen_cancercell_macrophage_subtyped_metadata.tsv",sep = "\t",col.names = T,row.names = F,quote = F)
# can_cell<-readRDS("/shared/resources/launonen_cancercell/launonen_cancercell_allcells.RDS")
# # write.table(can_cell@meta.data,"/shared/resources/launonen_cancercell/launonen_cancercell_allcells_metadata.tsv",sep = "\t",col.names = T,row.names = F,quote = F)
# 
# can_cell_cd8<-fread("/shared/resources/launonen_cancercell/launonen_cancercell_cd8tcell_ProjecTILs_metadata.tsv",data.table = F)
# rownames(can_cell_cd8)<-can_cell_cd8$cell_name
# can_cell_macrophages<-fread("/shared/resources/launonen_cancercell/launonen_cancercell_macrophage_subtyped_metadata.tsv",data.table = F)
# rownames(can_cell_macrophages)<-can_cell_macrophages$cell_name
# # can_cell<-fread("/shared/resources/launonen_cancercell/launonen_cancercell_allcells_metadata.tsv",data.table = F)
# # rownames(can_cell)<-can_cell$cell_name
# 
# table(can_cell_cd8$cell_name%in%can_cell$cell_name)
# table(can_cell_macrophages$cell_name%in%can_cell$cell_name)
# table(can_cell_cd8$functional.cluster,can_cell$cell_type[rownames(can_cell_cd8)])
# table(can_cell_macrophages$CD163_ITGAX_subtype,can_cell$cell_type[rownames(can_cell_macrophages)])
# #add cd8 info and macrophage info to whole annotation matrix
# can_cell$cell_type_detailed<-can_cell$cell_type
# can_cell$cell_type_detailed[rownames(can_cell_cd8)]<-can_cell_cd8$functional.cluster
# can_cell$cell_type_detailed[rownames(can_cell_macrophages)]<-paste0("macrophages_",can_cell_macrophages$CD163_ITGAX_subtype)
# table(can_cell$cell_type_detailed)
# saveRDS(can_cell,"/shared/resources/launonen_cancercell/launonen_cancercell_allcells_allannot.RDS")
# exp<-readRDS("/shared/resources/launonen_cancercell/launonen_cancercell_allcells_allannot.RDS")

# table(exp$sample,exp$cell_type_detailed)

#annotate signatures
# exp <- AddModuleScore_UCell(exp, features = list(mhc2_sig=mhc2_genes,int_alpha_sig=int_alpha,int_gamma_sig=int_gamma), name = NULL,ncores = 30)
# saveRDS(can_cell,"launonen_cancercell_allcells_allannot_sig_scored_v2.RDS")
exp<-readRDS("launonen_cancercell_allcells_allannot_sig_scored_bulk.RDS")

#plot signature in different cell_types
VlnPlot(exp,features = "mhc2_sig",group.by = "cell_type_detailed", pt.size = 0)+NoLegend()
# ggsave("results/v2/mhc2_signature_per_cellType.png",  # jpg, png, eps, tex, etc.
#        # plot = cm,
#        width = 30, height = 17,
#        units = "cm", # other options c("in", "cm", "mm"),
#        dpi = 600)
VlnPlot(exp,features = "int_alpha_sig",group.by = "cell_type_detailed", pt.size = 0)+NoLegend()
# ggsave("results/v2/int_alpha_sig_signature_per_cellType.png",  # jpg, png, eps, tex, etc.
#        # plot = cm,
#        width = 30, height = 17,
#        units = "cm", # other options c("in", "cm", "mm"),
#        dpi = 600)
VlnPlot(exp,features = "int_gamma_sig",group.by = "cell_type_detailed", pt.size = 0)+NoLegend()
# ggsave("results/v2/int_gamma_sig_signature_per_cellType.png",  # jpg, png, eps, tex, etc.
#        # plot = cm,
#        width = 30, height = 17,
#        units = "cm", # other options c("in", "cm", "mm"),
#        dpi = 600)

#filter cells with < 40 EOCs
eoc_samp<-table(exp$sample,exp$cell_type_detailed)
eoc_samp<-rownames(eoc_samp)[eoc_samp[,"Epithelial cells"]>40]
exp<-subset(exp,subset= sample%in%eoc_samp)

#check high-mid-low MHC2 groups 
exp$MHC2_group<-"mid"
exp$MHC2_group[exp$mhc2_sig<quantile(exp$mhc2_sig,.33)]<-"low"
exp$MHC2_group[exp$mhc2_sig>quantile(exp$mhc2_sig,.66)]<-"high"

#read HRD info
clin<-fread("../../resources/clinical_export_2022-09-26.csv",data.table = F)
#read EOC translator
eoc<-fread("../organoids/preprocessed_data/Patient_publicationID_170402.csv",data.table = F)
eoc<-eoc[eoc$Patient_cohort_code%in%unique(exp$patient),]
table(unique(exp$patient)%in%eoc$Patient_cohort_code)
rownames(eoc)<-eoc$PublicationID
table(eoc$PublicationID%in%clin$`Patient card::Publication code`)

#add metadata to seurat object
clin<-clin[clin$`Patient card::Publication code`%in%eoc$PublicationID,]
clin$HID<-eoc[clin$`Patient card::Publication code`,"Patient_cohort_code"]
rownames(clin)<-clin$HID
clin$OS_class<-"medium"
clin[clin$`Time to death_Days from Dg`>1009&!is.na(clin$`Time to death_Days from Dg`)&clin$`HR signature SBS3 per patient`=="HRP","OS_class"]<-"long"
clin[clin$`Time to death_Days from Dg`<532&!is.na(clin$`Time to death_Days from Dg`)&clin$`HR signature SBS3 per patient`=="HRP","OS_class"]<-"short"
clin[clin$`Time to death_Days from Dg`>1336&!is.na(clin$`Time to death_Days from Dg`)&clin$`HR signature SBS3 per patient`=="HRD","OS_class"]<-"long"
clin[clin$`Time to death_Days from Dg`<777&!is.na(clin$`Time to death_Days from Dg`)&clin$`HR signature SBS3 per patient`=="HRD","OS_class"]<-"short"
clin$OS<-clin$`Time to death_Days from Dg`
clin$BRCA<-clin$`BRCA mutation status`

exp$HRD<-clin[exp$patient,"HR signature SBS3 per patient"]
exp$OS_class<-clin[exp$patient,"OS_class"]
exp$OS<-clin[exp$patient,"Time to death_Days from Dg"]
exp$BRCA<-clin[exp$patient,"BRCA mutation status"]
exp$PDS<-exp$class


#check expression of DPB1 in HRD/HRP; prim/int
x<-subset(exp, subset=cell_type_detailed=="Epithelial cells")
x<-NormalizeData(x)

a<-table(x$sample,x$mhc2_sig>0)
t(apply(a, 1,function(z)z/sum(z)))

#categorize MHC2 at sample level
bulk<-AggregateExpression(object = x,group.by = "sample",slot="counts")
bulk<-bulk$RNA

de.exp<- DGEList(counts=as.matrix(bulk))
keep.exprs <- rowSums(cpm(de.exp)>1)>=2###keep genes that in at least 2 samples  has cpm>1
de.exp <- de.exp[keep.exprs,, keep.lib.sizes=FALSE]
limma<- calcNormFactors(de.exp,method = "TMM")
tmmexp<-cpm(limma,log = T,prior.count = 3)##get the TMM normalized matrix (counts*1000000/(lib.size*norm.factors)),
metabulk<-as.data.frame(x@meta.data%>%select(sample,patient,PDS,HRD)%>%group_by(sample,patient,PDS,HRD)%>%distinct())
rownames(metabulk)<-metabulk$sample

mhc2_scores <- ScoreSignatures_UCell(tmmexp, features=list(mhc2_sig=mhc2_genes))

metabulk$mhc2_group<-"mid"
metabulk[rownames(mhc2_scores)[mhc2_scores>quantile(mhc2_scores,.66)],"mhc2_group"]<-"high" #higher quantile 0.165
metabulk[rownames(mhc2_scores)[mhc2_scores<quantile(mhc2_scores,.33)],"mhc2_group"]<-"low" #lower  quantile 0.06
table(metabulk$PDS,metabulk$mhc2_group)

exp$MHC2_group_bulk<-"mid"
exp$MHC2_group_bulk[exp$sample%in%rownames(metabulk)[metabulk$mhc2_group=="low"]]<-"low"
exp$MHC2_group_bulk[exp$sample%in%rownames(metabulk)[metabulk$mhc2_group=="high"]]<-"high"

# saveRDS(exp,"launonen_cancercell_allcells_allannot_sig_scored_bulk.RDS")

#comparison of int gamma and alpha in low and high groups
plot_list<-list()
exp_list<-list()
gamma_genes<-c("IFNG","IFNGR1","IFNGR2")
# gamma_genes<-c("CD3E","CD3D","CD4","CD8A","KLRF1", "KLRB1", "GNLY", "NKG7")
for(ct in unique(exp$cell_type_detailed)){
  df<-exp@meta.data
  df<-df[df$cell_type_detailed==ct,]
  df<-df%>%filter(MHC2_group_bulk!="mid")
  df$MHC2_group_bulk<-factor(df$MHC2_group_bulk,levels=c("low","high"))
  
  if(length(unique(df$MHC2_group_bulk))==1)next
  p<-ggboxplot(df, x = "MHC2_group_bulk", y = "int_gamma_sig",
               color = "MHC2_group_bulk", palette = "jco",
               add = "none") + stat_compare_means(method = "wilcox.test")+
    ylab("interferon gamma signature score")+ggtitle(ct)+NoLegend()
  plot_list[[ct]]<-p
  
  df<-cbind.data.frame(df,t(exp@assays$SCT@data[gamma_genes,rownames(df)]))
  p<-ggboxplot(df, x = "MHC2_group_bulk", y = gamma_genes,combine = T,
               color = "MHC2_group_bulk", palette = "jco",
               add = "none") + stat_compare_means(method = "wilcox.test")+
    ylab("normalized expression")+ggtitle(ct)+NoLegend()
  exp_list[[ct]]<-p
}
pl<-cowplot::plot_grid(plotlist = plot_list,ncol = 4)
ggsave("results/v2/int_gamma_score_highvslow_bulk_groups.png",plot = pl, width = 10, height = 23)

pl<-cowplot::plot_grid(plotlist = exp_list,ncol = 4)
ggsave("results/v2/int_gamma_gene_exp_highvslow_bulk_groups.png",plot = pl, width = 20, height = 25)



plot_list<-list()
exp_list<-list()
alpha_genes<-c("IFNA1","IFNA2","IFNAR1","IFNAR2","IFNB1")

for(ct in unique(exp$cell_type_detailed)){
  df<-exp@meta.data
  df<-df[df$cell_type_detailed==ct,]
  df<-df%>%filter(MHC2_group_bulk!="mid")
  df$MHC2_group_bulk<-factor(df$MHC2_group_bulk,levels=c("low","high"))
  
  if(length(unique(df$MHC2_group_bulk))==1)next
  p<-ggboxplot(df, x = "MHC2_group_bulk", y = "int_alpha_sig",
               color = "MHC2_group_bulk", palette = "jco",
               add = "none") + stat_compare_means(method = "wilcox.test")+
    ylab("interferon alpha signature score")+ggtitle(ct)+NoLegend()
  plot_list[[ct]]<-p
  
  df<-cbind.data.frame(df,t(exp@assays$SCT@data[alpha_genes,rownames(df)]))
  p<-ggboxplot(df, x = "MHC2_group_bulk", y = alpha_genes,combine = T,
               color = "MHC2_group_bulk", palette = "jco",
               add = "none") + stat_compare_means(method = "wilcox.test")+
    ylab("normalized expression")+ggtitle(ct)+NoLegend()
  exp_list[[ct]]<-p
}
pl<-cowplot::plot_grid(plotlist = plot_list,ncol = 4)
ggsave("results/v2/int_alpha_score_highvslow_bulk_groups.png",plot = pl, width = 10, height = 23)

pl<-cowplot::plot_grid(plotlist = exp_list,ncol = 4)
ggsave("results/v2/int_alpha_gene_exp_highvslow_bulk_groups.png",plot = pl, width = 20, height = 39)

# Proportion of M1 and M2 in MHC2 low vs high

m2<-prop.table(table(exp$MHC2_group_bulk[exp$cell_type_detailed=="macrophages_CD163_high"]))

m1<-table(exp$MHC2_group_bulk[exp$cell_type_detailed=="macrophages_ITGAX_high"])

prop.table(table(exp$cell_type_detailed,exp$MHC2_group_bulk),2)

#proportion boxplot at the sample level

complete_grid <- expand.grid(# Create a grid of all possible combinations of sample, cell_type_detailed, and MHC2_group_bulk
  sample = unique(exp@meta.data$sample),
  cell_type_detailed = unique(exp@meta.data$cell_type_detailed)
)
total_cells <- exp@meta.data %>%# Join the meta data to the complete grid, retaining MHC2_group_bulk based on the original sample
  group_by(sample, cell_type_detailed, MHC2_group_bulk,class) %>%
  summarise(total = n(), .groups = 'drop') %>%
  right_join(complete_grid, by = c("sample", "cell_type_detailed")) %>%
  group_by(sample) %>%
  fill(MHC2_group_bulk, .direction = "downup") %>%  # Fill in MHC2_group_bulk based on the sample
  replace_na(list(total = 0)) %>%
  mutate(proportion = total / sum(total)) %>%
  ungroup()

table(total_cells$sample)


plot_list<-list()
cellTypes<-unique(exp$cell_type_detailed)
cellTypes<-grep("macro|T cell|DC|CD8",cellTypes,value=T)
for(ct in cellTypes){
  # df<-total_cells %>% filter(MHC2_group_bulk!="mid"&cell_type_detailed==ct&class=="primary")
  df<-total_cells %>% filter(MHC2_group_bulk!="mid"&cell_type_detailed==ct)
  
  df$MHC2_group_bulk<-factor(df$MHC2_group_bulk,levels=c("low","high"))
  
  p<-ggboxplot(df, x = "MHC2_group_bulk", y = "proportion",
               color = "MHC2_group_bulk", palette = "jco",
               add = "jitter") + stat_compare_means(method = "wilcox.test")+
    ylab(paste0("proportion"))+ggtitle(ct)+NoLegend()+stat_n_text()
  plot_list[[ct]]<-p
}
pl<-cowplot::plot_grid(plotlist = plot_list,ncol = 4)
# ggsave("results/v2/selected_celltype_proportion_primary_highvslow_bulk_groups.png",plot = pl, width = 10, height = 16)
ggsave("results/v2/selected_celltype_proportion_highvslow_bulk_groups_wilcox.png",plot = pl, width = 10, height = 16)
ggsave("results/v2/selected_celltype_proportion_highvslow_bulk_groups_wilcox.svg",plot = pl, width = 12, height = 16)

##plot correlation between mhc2 and NK and T cell markers in epithelial
nk_t_genes<-c("CD3E","CD4","CD8A","KLRF1", "KLRB1", "GNLY", "NKG7")
plot_list<-list()
for(gene_id in nk_t_genes){
  df<-cbind.data.frame(mhc2_scores,gene=tmmexp[gene_id,rownames(mhc2_scores)])
  
  cor_test <- cor.test(df$mhc2_sig_UCell, df$gene, method = "pearson")
  
  # Create the plot with ggplot
  p <- ggplot(df, aes(x = mhc2_sig_UCell, y = gene)) +
    geom_point(alpha = 0.6) +  # Scatter plot points with some transparency
    geom_smooth(method = "lm", col = "blue", se = TRUE) +  # Add linear regression line
    labs(title = gene_id) +  # Title as the cell type
    stat_cor(aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")), method = "pearson") +  # Add correlation coefficient and p-value
    theme_minimal()
  
  plot_list[[gene_id]]<-p
}

pl<-cowplot::plot_grid(plotlist = plot_list,ncol = 4)
ggsave("results/v2/correlations_int_gamma_vs_mhc2.png",plot = pl, width = 15, height = 19,bg="white")



##plot correlation between mhc2 and interferon alpha and gamma in the different celltypes

plot_list<-list()
for(ct in unique(exp$cell_type_detailed)){
  df<-exp@meta.data %>%
    group_by(cell_type_detailed) %>%filter(cell_type_detailed==ct)
  
  cor_test <- cor.test(df$mhc2_sig, df$int_gamma_sig, method = "pearson")
  
  # Create the plot with ggplot
  p <- ggplot(df, aes(x = mhc2_sig, y = int_gamma_sig)) +
    geom_point(alpha = 0.6) +  # Scatter plot points with some transparency
    geom_smooth(method = "lm", col = "blue", se = TRUE) +  # Add linear regression line
    labs(title = ct) +  # Title as the cell type
    stat_cor(aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")), method = "pearson") +  # Add correlation coefficient and p-value
    theme_minimal()
  
  plot_list[[ct]]<-p
}

pl<-cowplot::plot_grid(plotlist = plot_list,ncol = 4)
ggsave("results/v2/correlations_int_gamma_vs_mhc2.png",plot = pl, width = 15, height = 19,bg="white")



plot_list<-list()
for(ct in unique(exp$cell_type_detailed)){
  df<-exp@meta.data %>%
    group_by(cell_type_detailed) %>%filter(cell_type_detailed==ct)
  
  cor_test <- cor.test(df$mhc2_sig, df$int_alpha_sig, method = "pearson")
  
  # Create the plot with ggplot
  p <- ggplot(df, aes(x = mhc2_sig, y = int_alpha_sig)) +
    geom_point(alpha = 0.6) +  # Scatter plot points with some transparency
    geom_smooth(method = "lm", col = "blue", se = TRUE) +  # Add linear regression line
    labs(title = ct) +  # Title as the cell type
    stat_cor(aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")), method = "pearson") +  # Add correlation coefficient and p-value
    theme_minimal()
  
  plot_list[[ct]]<-p
}

pl<-cowplot::plot_grid(plotlist = plot_list,ncol = 4)
ggsave("results/v2/correlations_int_alpha_vs_mhc2.png",plot = pl, width = 15, height = 19,bg="white")



plot_list<-list()
for(ct in unique(exp$cell_type_detailed)){
  df<-exp@meta.data %>%
    group_by(cell_type_detailed) %>%filter(cell_type_detailed==ct)
  
  cor_test <- cor.test(df$int_gamma_sig, df$int_alpha_sig, method = "pearson")
  
  # Create the plot with ggplot
  p <- ggplot(df, aes(x = int_gamma_sig, y = int_alpha_sig)) +
    geom_point(alpha = 0.6) +  # Scatter plot points with some transparency
    geom_smooth(method = "lm", col = "blue", se = TRUE) +  # Add linear regression line
    labs(title = ct) +  # Title as the cell type
    stat_cor(aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")), method = "pearson") +  # Add correlation coefficient and p-value
    theme_minimal()
  
  plot_list[[ct]]<-p
}

pl<-cowplot::plot_grid(plotlist = plot_list,ncol = 4)
ggsave("results/v2/correlations_int_alpha_vs_gamma.png",plot = pl, width = 15, height = 19,bg="white")

### epithelial mean MHC II against the cell type proportions (maybe the grouping looses a bit too much signal) 
mhc2_scores<-as.data.frame(mhc2_scores)
mhc2_scores$sample<-rownames(mhc2_scores)
df_merged <- merge(total_cells, mhc2_scores, by = "sample", all.x = TRUE)

plot_list<-list()
cellTypes<-unique(exp$cell_type_detailed)
# cellTypes<-grep("macro|T cell|DC|CD8",cellTypes,value=T)
for(ct in cellTypes){
  # df<-total_cells %>% filter(MHC2_group_bulk!="mid"&cell_type_detailed==ct&class=="primary")
  df<-df_merged %>% filter(cell_type_detailed==ct)
  
  cor_test <- cor.test(df$proportion, df$mhc2_sig_UCell, method = "pearson")
  
  # Create the plot with ggplot
  p <- ggplot(df, aes(x = proportion, y = mhc2_sig_UCell)) +
    geom_point(alpha = 0.6) +  # Scatter plot points with some transparency
    geom_smooth(method = "lm", col = "blue", se = TRUE) +  # Add linear regression line
    labs(title = ct) +  # Title as the cell type
    stat_cor(aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")), method = "pearson") +  # Add correlation coefficient and p-value
    theme_minimal()
  
  plot_list[[ct]]<-p
}
pl<-cowplot::plot_grid(plotlist = plot_list,ncol = 4)
# ggsave("results/v2/selected_celltype_proportion_primary_highvslow_bulk_groups.png",plot = pl, width = 10, height = 16)
ggsave("results/v2/celltype_proportion_correlation_pseudobulk.png",plot = pl, width = 10, height = 16,bg="white")

### epithelial mean MHC II (pseudobulk etc) vs mean IFNG across the cell types.

plot_list<-list()
box_list<-list()
cellTypes<-unique(exp$cell_type_detailed)
# cellTypes<-grep("macro|T cell|DC|CD8",cellTypes,value=T)
for(ct in cellTypes){
  #subset ct
  x<-subset(exp, subset=cell_type_detailed==ct)
  tab<-table(x$sample)
  x$samp_cell_count<-as.numeric(tab[x$sample])
  
  if(!any(x$samp_cell_count>20))next
  x<-subset(x, subset = (samp_cell_count >= 20))
  
  if(length(unique(x$sample))<3)next
  
  #categorize MHC2 at sample level
  bulk<-AggregateExpression(object = x,group.by = "sample",slot="counts",)
  bulk<-bulk$RNA
  
  de.exp<- DGEList(counts=as.matrix(bulk))
  keep.exprs <- rowSums(cpm(de.exp)>1)>=2###keep genes that in at least 2 samples  has cpm>1
  de.exp <- de.exp[keep.exprs,, keep.lib.sizes=FALSE]
  limma<- calcNormFactors(de.exp,method = "TMM")
  tmmexp<-cpm(limma,log = T,prior.count = 3)##get the TMM normalized matrix (counts*1000000/(lib.size*norm.factors)),
  
  scores <- ScoreSignatures_UCell(tmmexp, features=list(int_alpha_sig=int_alpha))
  
  scores<-as.data.frame(scores)
  scores$mhc2_sig_UCell<-mhc2_scores[rownames(scores),"mhc2_sig_UCell"]
  scores$MHC2_group_bulk<-metabulk[rownames(scores),"mhc2_group"]
  
  
  cor_test <- cor.test(scores$mhc2_sig_UCell, scores$int_alpha_sig_UCell, method = "pearson")
  
  # Create the plot with ggplot
  p <- ggplot(scores, aes(x = int_alpha_sig_UCell, y = mhc2_sig_UCell)) +
    geom_point(alpha = 0.6) +  # Scatter plot points with some transparency
    geom_smooth(method = "lm", col = "blue", se = TRUE) +  # Add linear regression line
    labs(title = ct) +  # Title as the cell type
    stat_cor(aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")), method = "pearson") +  # Add correlation coefficient and p-value
    theme_minimal()
  scores<-scores%>%filter(MHC2_group_bulk!="mid")
  scores$MHC2_group_bulk<-factor(scores$MHC2_group_bulk,levels=c("low","high"))
  if(length(unique(scores$MHC2_group_bulk))<2)next
  if(sum(scores$MHC2_group_bulk=="low")<2|sum(scores$MHC2_group_bulk=="high")<2)next
  
  b<-ggboxplot(scores, x = "MHC2_group_bulk", y = "int_alpha_sig_UCell",
               color = "MHC2_group_bulk", palette = "jco",
               add = "jitter") + stat_compare_means(method = "wilcox.test")+
    ylab(paste0("interferon alpha\nscore pseudobulk"))+ggtitle(ct)+NoLegend()+stat_n_text()
  
  plot_list[[ct]]<-p
  box_list[[ct]]<-b
}
pl<-cowplot::plot_grid(plotlist = plot_list,ncol = 4)
# ggsave("results/v2/selected_celltype_proportion_primary_highvslow_bulk_groups.png",plot = pl, width = 10, height = 16)
ggsave("results/v2/correlation_int_alpha_vs_mhc2_scores_pseudobulk.png",plot = pl, width = 10, height = 16,bg="white")

pl<-cowplot::plot_grid(plotlist = box_list,ncol = 4)
# ggsave("results/v2/selected_celltype_proportion_primary_highvslow_bulk_groups.png",plot = pl, width = 10, height = 16)
ggsave("results/v2/boxplot_int_alpha_pseudobulk_vs_mhc2_scores_pseudobulk_wilcox.png",plot = pl, width = 10, height = 16,bg="white")

#for gamma
plot_list<-list()
box_list<-list()
cellTypes<-unique(exp$cell_type_detailed)
# cellTypes<-grep("macro|T cell|DC|CD8",cellTypes,value=T)
for(ct in cellTypes){
  #subset ct
  x<-subset(exp, subset=cell_type_detailed==ct)
  tab<-table(x$sample)
  x$samp_cell_count<-as.numeric(tab[x$sample])
  
  if(!any(x$samp_cell_count>20))next
  x<-subset(x, subset = (samp_cell_count >= 20))
  
  if(length(unique(x$sample))<3)next
  
  #categorize MHC2 at sample level
  bulk<-AggregateExpression(object = x,group.by = "sample",slot="counts",)
  bulk<-bulk$RNA
  
  de.exp<- DGEList(counts=as.matrix(bulk))
  keep.exprs <- rowSums(cpm(de.exp)>1)>=2###keep genes that in at least 2 samples  has cpm>1
  de.exp <- de.exp[keep.exprs,, keep.lib.sizes=FALSE]
  limma<- calcNormFactors(de.exp,method = "TMM")
  tmmexp<-cpm(limma,log = T,prior.count = 3)##get the TMM normalized matrix (counts*1000000/(lib.size*norm.factors)),
  
  scores <- ScoreSignatures_UCell(tmmexp, features=list(int_gamma_sig=int_gamma))
  
  scores<-as.data.frame(scores)
  scores$mhc2_sig_UCell<-mhc2_scores[rownames(scores),"mhc2_sig_UCell"]
  scores$MHC2_group_bulk<-metabulk[rownames(scores),"mhc2_group"]
  
  cor_test <- cor.test(scores$mhc2_sig_UCell, scores$int_gamma_sig_UCell, method = "pearson")
  
  # Create the plot with ggplot
  p <- ggplot(scores, aes(x = int_gamma_sig_UCell, y = mhc2_sig_UCell)) +
    geom_point(gamma = 0.6) +  # Scatter plot points with some transparency
    geom_smooth(method = "lm", col = "blue", se = TRUE) +  # Add linear regression line
    labs(title = ct) +  # Title as the cell type
    stat_cor(aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")), method = "pearson") +  # Add correlation coefficient and p-value
    theme_minimal()
  
  scores<-scores%>%filter(MHC2_group_bulk!="mid")
  scores$MHC2_group_bulk<-factor(scores$MHC2_group_bulk,levels=c("low","high"))
  if(length(unique(scores$MHC2_group_bulk))<2)next
  if(sum(scores$MHC2_group_bulk=="low")<2|sum(scores$MHC2_group_bulk=="high")<2)next
  
  b<-ggboxplot(scores, x = "MHC2_group_bulk", y = "int_gamma_sig_UCell",
               color = "MHC2_group_bulk", palette = "jco",
               add = "jitter") + stat_compare_means(method = "wilcox.test")+
    ylab(paste0("interferon gamma\nscore pseudobulk"))+ggtitle(ct)+NoLegend()+stat_n_text()
  
  plot_list[[ct]]<-p
  box_list[[ct]]<-b
}
pl<-cowplot::plot_grid(plotlist = plot_list,ncol = 4)
# ggsave("results/v2/selected_celltype_proportion_primary_highvslow_bulk_groups.png",plot = pl, width = 10, height = 16)
ggsave("results/v2/correlation_int_gamma_vs_mhc2_scores_pseudobulk.png",plot = pl, width = 10, height = 16,bg="white")

pl<-cowplot::plot_grid(plotlist = box_list,ncol = 4)
# ggsave("results/v2/selected_celltype_proportion_primary_highvslow_bulk_groups.png",plot = pl, width = 10, height = 16)
ggsave("results/v2/boxplot_int_gamma_pseudobulk_vs_mhc2_scores_pseudobulk_wilcox.png",plot = pl, width = 10, height = 16,bg="white")


### epithelial mean MHC II (pseudobulk etc) vs mean IFNG across the cell types. At gene expression level
# alpha_genes<-c("IFNA1","IFNA2","IFNAR1","IFNAR2","IFNB1")
# gamma_genes<-c("IFNG","IFNGR1","IFNGR2")
# genes<-c(alpha_genes,gamma_genes)
genes<-c("IFNA1","IFNA2","IFNB1","IFNG")
genes<-c("ADGRG2","GPR64","HE6","TM7LN2","CLSTN3","CS3","KIAA0726","CXADR","CAR","DBN1","D0S117E","DPPA2","PESCRG1","EFEMP2","FBLN4",
         "UNQ200/PRO226","FARSB","FARSLB","FRSB","HSPC173","FAT1","CDHF7","FAT","FMOD","FM","SLRR2E","FNDC1","FNDC2","KIAA1866","MEL4B3","HTRA1",
         "HTRA","PRSS11","ITGB5","LAMA3","LAMNA","LCN2","HNL","NGAL","MSLN","MPF","MUC16","CA125","NEO1","IGDCC2","NGN","PLEC","PLEC1","PRDX2",
         "NKEFB","TDPX1","PTGFRN","CD9P1","EWIF","FPRP","KIAA1436","PTPRF","LAR","PTPRS","RARRES1","PEIG1","TIG1","RCN1","RCN","SQSTM1","ORCA",
         "OSIL","UBA52","UBCEP2","UBB","UBC")
genes<-genes[genes%in%rownames(tmmexp)]

plot_list<-list()
exp_list<-list()
cellTypes<-unique(exp$cell_type_detailed)
# cellTypes<-grep("macro|T cell|DC|CD8",cellTypes,value=T)
for(ct in cellTypes){
  #subset ct
  x<-subset(exp, subset=cell_type_detailed==ct)
  tab<-table(x$sample)
  x$samp_cell_count<-as.numeric(tab[x$sample])
  
  if(!any(x$samp_cell_count>20))next
  x<-subset(x, subset = (samp_cell_count >= 20))
  
  if(length(unique(x$sample))<3)next
  
  #categorize MHC2 at sample level
  bulk<-AggregateExpression(object = x,group.by = "sample",slot="counts",)
  bulk<-bulk$RNA
  
  de.exp<- DGEList(counts=as.matrix(bulk))
  # keep.exprs <- rowSums(cpm(de.exp)>1)>=2###keep genes that in at least 2 samples  has cpm>1
  # de.exp <- de.exp[keep.exprs,, keep.lib.sizes=FALSE]
  limma<- calcNormFactors(de.exp,method = "TMM")
  tmmexp<-cpm(limma,log = T,prior.count = 3)##get the TMM normalized matrix (counts*1000000/(lib.size*norm.factors)),
  
  scores <- cbind(mhc2_scores[colnames(tmmexp),"mhc2_sig_UCell"],melt(t(tmmexp)[,genes]),ct)
  colnames(scores)<-c("mhc2_sig_UCell","sample","gene","gene_norm_exp_pseudobulk","cellType")
  # cor_test <- cor.test(scores$mhc2_sig_UCell, scores$gene_norm_exp_pseudobulk, method = "pearson")
  result <- scores %>%
    group_by(gene) %>%
    summarise(cor_coef = cor.test(mhc2_sig_UCell, gene_norm_exp_pseudobulk)$estimate,
              cor_pval = cor.test(mhc2_sig_UCell, gene_norm_exp_pseudobulk)$p.value)
  
  scores <- scores %>%
    left_join(result, by = "gene")
  
  exp_list[[ct]]<-scores
  
  # Create the plot with ggplot
  p <- ggplot(scores, aes(x = gene_norm_exp_pseudobulk, y = mhc2_sig_UCell)) +
    geom_point(alpha = 0.6) +  # Scatter plot points with some transparency
    geom_smooth(method = "lm", col = "blue", se = TRUE) +  # Add linear regression line
    facet_wrap(~gene) +
    labs(title = ct) +  # Title as the cell type
    stat_cor(aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")), method = "pearson") +  # Add correlation coefficient and p-value
    theme_minimal()
  
  plot_list[[ct]]<-p
  print(ct)
}
pl<-cowplot::plot_grid(plotlist = plot_list,ncol = 4)

# ggsave("results/v2/propelation_interferon_genes_expression_vs_mhc2_scores_pseudobulk.png",plot = pl, width = 30, height = 36,bg="white")
ggsave("results/v2/correlation_interferon_genes_expression_vs_mhc2_scores_pseudobulk_v2.png",plot = pl, width = 30, height = 36,bg="white")
# ggsave("results/v2/correlation_list_15_genes_expression_vs_mhc2_scores_pseudobulk.pdf",plot = p, width = 10, height = 8,bg="white")

write.table(rbindlist(exp_list),"results/v2/correlation_interferon_genes_expression_vs_mhc2_scores_pseudobulk.tsv",col.names = T,quote = F,sep="\t",row.names = F)
write.table(scores,"results/v2/correlation_15_genes_expression_vs_mhc2_scores_pseudobulk_epithelial.tsv",col.names = T,quote = F,sep="\t",row.names = F)

#####DEGs
#pseudobulk analysis
x<-subset(exp, subset=cell_type_detailed=="Epithelial cells")

bulk<-AggregateExpression(object = x,group.by = "sample",slot="counts")
bulk<-bulk$RNA

de.exp<- DGEList(counts=as.matrix(bulk))
keep.exprs <- rowSums(cpm(de.exp)>1)>=2###keep genes that in at least 2 samples  has cpm>1
de.exp <- de.exp[keep.exprs,, keep.lib.sizes=FALSE]
limma<- calcNormFactors(de.exp,method = "TMM")
tmmexp<-cpm(limma,log = T,prior.count = 3)##get the TMM normalized matrix (counts*1000000/(lib.size*norm.factors)),
metabulk<-as.data.frame(x@meta.data%>%select(sample,patient,PDS,HRD)%>%group_by(sample,patient,PDS,HRD)%>%distinct())
rownames(metabulk)<-metabulk$sample

mhc2_scores <- ScoreSignatures_UCell(tmmexp, features=list(mhc2_sig=mhc2_genes))

metabulk$mhc2_group<-"mid"
metabulk[rownames(mhc2_scores)[mhc2_scores>quantile(mhc2_scores,.66)],"mhc2_group"]<-"high" #higher quantile 0.165
metabulk[rownames(mhc2_scores)[mhc2_scores<quantile(mhc2_scores,.33)],"mhc2_group"]<-"low" #lower  quantile 0.06
table(metabulk$PDS,metabulk$mhc2_group)


# metabulk<-x@meta.data%>%select(sample,patient,NACT)%>%group_by(sample,patient,NACT)%>%distinct()
design <- model.matrix(~0+metabulk$mhc2_group+metabulk$patient)
# design <- model.matrix(~0+metabulk$NACT+metabulk$patient)

# colnames(design)<-gsub("metabulk$NACT","",colnames(design),fixed = T)
colnames(design)<-gsub("metabulk$mhc2_group","",colnames(design),fixed = T)
colnames(design)<-gsub("metabulk$patient","",colnames(design),fixed = T)

fit<-voom(limma,design,plot = T)
fit<-lmFit(fit,design)

# cont.wt <- makeContrasts(prim_vs_int=interval-primary ,  levels=design)
cont.wt <- makeContrasts(high_vs_low_mhc2=high-low,  levels=design)
fit2 <- contrasts.fit(fit, cont.wt)

fit2 <- eBayes(fit2)
plotSA(fit2)

de_res<-topTable(fit2,coef = "high_vs_low_mhc2",number = Inf)
table(de_res$adj.P.Val<0.05)
write.table(de_res,file=paste0("results/v2/pseudobulk_high_vs_low_mhc2_DEGs_limma.tsv"),sep = "\t",quote=F )






