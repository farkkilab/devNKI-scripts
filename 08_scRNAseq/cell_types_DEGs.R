# get celltypes signatures (DEGs)
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

exp<-readRDS("launonen_cancercell_allcells_allannot_sig_scored_bulk.RDS")

Idents(exp)<-"cell_type_detailed"
seurat_DEGs<-FindAllMarkers(object = exp, test.use = "LR",logfc.threshold=0.05 )


#####pseudobulk

#get matrix
exp_list<-list()
design_list<-list()
for(ct in unique(exp$cell_type_detailed)){
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
  ct_name<-gsub(" ","_",ct,fixed = T)
  ct_name<-gsub("/","_",ct_name,fixed = T)
  
  colnames(bulk)<-paste0(ct_name,"_",colnames(bulk))
  exp_list[[ct]]<-bulk
  
  design<-x@meta.data%>%select(sample,samp_cell_count)%>%distinct()
  design<-cbind(ct,design,colnames(bulk))
  design_list[[ct]]<-design
  
  print(ct)
}


#calculate DEGs
exp_df<- do.call(cbind, exp_list) 

de.exp<- DGEList(counts=as.matrix(exp_df))
keep.exprs <- rowSums(cpm(de.exp)>1)>=5###keep genes that in at least 2 samples  has cpm>1
de.exp <- de.exp[keep.exprs,, keep.lib.sizes=FALSE]
limma<- calcNormFactors(de.exp,method = "TMM")
# tmmexp<-cpm(limma,log = T,prior.count = 3)##get the TMM normalized matrix (counts*1000000/(lib.size*norm.factors)),

design<-do.call(rbind, design_list) 
colnames(design)<-c("cell_type","sample","cell_count","ID")
rownames(design)<-design$ID


pseudobulk_DEGs<-list()
for(ct in unique(design$cell_type)){
  ds<-design
  ds$comp<-ifelse(ds$cell_type==ct,"cell_type","rest")
  ds <- model.matrix(~0+ds$comp+ds$sample)
  # design <- model.matrix(~0+metabulk$NACT+metabulk$patient)
  
  # colnames(design)<-gsub("metabulk$NACT","",colnames(design),fixed = T)
  colnames(ds)<-gsub("ds$comp","",colnames(ds),fixed = T)
  colnames(ds)<-gsub("ds$sample","",colnames(ds),fixed = T)
  
  fit<-voom(limma,ds,plot = T)
  fit<-lmFit(fit,ds)
  
  # cont.wt <- makeContrasts(prim_vs_int=interval-primary ,  levels=ds)
  cont.wt <- makeContrasts(sig=cell_type - rest,  levels=ds)
  fit2 <- contrasts.fit(fit, cont.wt)
  
  fit2 <- eBayes(fit2)
  plotSA(fit2)
  
  de_res<-topTable(fit2,coef = "sig",number = Inf)
  table(de_res$adj.P.Val<0.05)
  de_res$cell_type<-ct
  de_res<-de_res[de_res$adj.P.Val<0.05,]
  de_res<-cbind(gene_name=rownames(de_res),de_res)
  pseudobulk_DEGs[[ct]]<-de_res

  
  print(ct)
}



res<- do.call(rbind, pseudobulk_DEGs) 
write.table(res,file=paste0("results/pseudobulk_celltype_DEGs_limma.tsv"),sep = "\t",quote=F ,row.names = F,col.names = T)

