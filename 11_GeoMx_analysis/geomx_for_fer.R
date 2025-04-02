# MHC-II in tumor vs T-cell abundance in corresponding stroma

# MHC-II in tumor vs IFNg in corresponding stroma

# MHC-II in tumor vs quirky epitopes in same AOI
library(data.table)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(reshape2)
library(tibble)
library(broom)
library(gridExtra)
##################################################################
input_dir <- "D:/users/fperez/NKI_TMAs_AF/GeoMX/Inputs/"
outp_dir <- 'D:/users/fperez/NKI_TMAs_AF/GeoMX/Outputs/'
dir.create(outp_dir)

#geomx_meta <- fread('~/Documents/phd/st/geomx-processing/results/batch1-1903/geomx_metadata.csv')

deconv_ct_frac <- data.frame(fread(paste0(input_dir,
                                          'bp_res_mid_lvl_ct_ct_fraction.csv')))
#GSEA calculated from full signal
gsea_all <- data.frame(fread(paste0(input_dir,'ssgsea_norm_harmony_batch_corr_all_custom_IFNg_pathways.csv.csv')))

#GSEA calculated from deconvoluted tumor signal
gsea_deconv_tum <- data.frame(fread(paste0(input_dir,'ssgsea_norm_harmony_batch_corr_deconv_tumor_custom_IFNg_pathways.csv.csv')))

#Defining order of cells for plots
ct_of_interest <- c('DCs', 'Macrophages', 'Mast.cells', 'Tcells', 'Bcells', 'NKcells', 'Fibroblasts','Endothelial.cells')

#Selecting other columns of interest
meta_names <- c('Sample', 'Patient', 'NACT_status', 'PFS', 'PFS_months', 'Annotation_cell')

###################################################################################################################
# prepare dfs
gsea_res <- gsea_deconv_tum

gsea_res$Roi <- paste0(gsea_res$Sample, '_', gsea_res$Roi)
gsea_res <- left_join(gsea_res, deconv_ct_frac[, c('dcc_filename', ct_of_interest)]) 

gsea_res_str <- gsea_res[gsea_res$Segment == 'stroma', ]
gsea_res_tumor <- gsea_res[gsea_res$Segment == 'tumor', ]

# TODO instead of this make it wide before and join with paths of interest
tum_mhc <- gsea_res_tumor[gsea_res_tumor$pathway == 'FER_MHC_REVISITED', ]
tum_mhc <- rename(tum_mhc, ssgsea_mhc_tum = 'ssgsea_score') #To rename column name
tum_mhc <- rename_with(tum_mhc, ~paste0(., "_tum"), ct_of_interest) #To rename column name

stroma_ifn <- gsea_res_str[gsea_res_str$pathway == 'HALLMARK_INTERFERON_GAMMA_RESPONSE', ]
stroma_ifn <- rename(stroma_ifn, ssgsea_ifng_str = 'ssgsea_score')
stroma_ifn <- rename_with(stroma_ifn, ~paste0(., "_str"), ct_of_interest)

tum_str_dt <- left_join(tum_mhc[, c(meta_names, 'Roi', 'ssgsea_mhc_tum', "Sample", paste0(ct_of_interest, '_tum'))], 
                        stroma_ifn[, c('Roi', 'ssgsea_ifng_str', paste0(ct_of_interest, '_str'))])


################################
# make scatters

inp_df <- tum_str_dt %>% filter(NACT_status == "post")
xname <- 'ssgsea_mhc_tum'
ynames <- c(paste0(ct_of_interest, '_str'))

output.scatter <- paste0(outp_dir,"/scatter/")
dir.create(output.scatter)

#For Figure 4o
#Scatter plots
list.plots <- list()
list.plots <- lapply(ynames, function(yname){
    #rm NA if needed
    inp_df2 <- inp_df[!is.na(inp_df[[xname]]) & !is.na(inp_df[[yname]]), ]
    #manual_colours <- brewer.pal(length(unique(as.factor(inp_df[[color_colname]]))), 'Paired')
    M1 <- lm(get(yname) ~ exp(get(xname)) + Patient, data = inp_df2)
    pval <- glance(M1)$p.value
    
    pval.text = case_when(pval < 1e-12 ~ "P<1e-10",
                          pval < 1e-6 ~ "P<1e-10",
                          pval < 1e-3 ~ "P<1e-10",
                          pval < 0.05 ~ "P<0.05",
                           pval >= 0.05 ~ "P=N.S",
                          .default = "Non significant")
    my_ylab = strsplit(yname,"_")[[1]][1]
    
    max.y <- max(inp_df2[[yname]]) - sd(inp_df2[[yname]]) * 0.7
    
    p <- ggplot(data = inp_df2, aes(x = get(xname), y = get(yname))) +
          geom_point(alpha = 0.3) + 
          annotate("text",
                   label=paste("lm(y ~ exp(x) + Patient) :", "R^2=", round(summary(M1)$r.squared, 2)),
                   x=0.73, y=max.y, size = unit(3.5, "pt")) +
          xlab("Deconv. CC MHCII signature (Ucell)") +
          ylab(paste0(my_ylab, " inferred proportion")) +
          theme_bw() +
          theme(plot.title = element_text(size=rel(0.9)),
               axis.title.y = element_text(size=rel(1.4)),
               axis.text = element_text(size=rel(1.2)),
               axis.title.x = element_text(size=rel(1.4))) +
          #scale_color_manual(values=manual_colours) +
          #geom_smooth(method='glm', formula= y~x ) +
          stat_poly_line(method = 'lm', formula = y~exp(x)) +
          stat_cor(method = "pearson")
    
    ggsave(p, filename=paste0(output.scatter,'/scatter_', xname, '_vs_', yname, '_tumor_deconvoluted_signal.svg'),
           width = 9, height = 9, units = "cm")
    return(p)
    
})
ggsave(file=paste0(output.scatter, "all_plots.svg"), arrangeGrob(grobs = list.plots, ncol = 3),
       width = 26, height = 24, units = "cm")


###Calculating correlation between inferred abundance of immune in stromal-TSI vs MCHII in tumor-TSI
#Analysis for Figure 4p
corrs <- sapply(ynames, function(yname){
  #rm NA if needed
  inp_df2 <- inp_df[!is.na(inp_df[[xname]]) & !is.na(inp_df[[yname]]), ]
  my_ylab = strsplit(yname,"_")[[1]][1]
  #Correlation pearson
  cor.res <- cor.test(inp_df2[[yname]], inp_df2[[xname]], method = "pearson")
  out <- c(my_ylab, cor.res$p.value, cor.res$estimate)
  names(out) <- c("Celltype","P.val","Cor")
  return(out)
})
corrs.df <- as.data.frame(t(corrs))
corrs.df$Cor <- as.numeric(corrs.df$Cor)
corrs.df$P.val <- as.numeric(corrs.df$P.val)
corrs.df$fdr <- p.adjust(corrs.df$P.val)
corrs.df$class <- "Fold"

#Organizing data for plotting
corrs.df <- corrs.df %>% mutate(P.val.class = 
                              case_when(fdr < 1e-3 ~ 3,
                                        fdr < 0.01 ~ 2,
                                        fdr < 0.05 ~ 1,
                                        fdr >= 0.05 ~ 0))

corrs.df <- corrs.df %>% mutate(Celltype = case_when(Celltype == "Tcells" ~ "T.cells",
                                         Celltype == "Bcells" ~ "B.cells",
                                         Celltype == "NKcells" ~ "NK.cells",
                                         Celltype == "DCs" ~ "Dendritic",
                                         Celltype == "Endothelial.cells" ~ "Endothelial",
                                         .default = Celltype))

corrs.df$Celltype <- factor(corrs.df$Celltype, levels=rev(c("Fibroblasts","Endothelial","Macrophages","Mast.cells","Dendritic","NK.cells","B.cells","T.cells")))

#For Figure 4p
#Generating dot plot
p <- ggplot(corrs.df, aes(x=class, y=Celltype)) + geom_point(aes(col=Cor, size=P.val.class)) + 
      scale_colour_gradient2(low= "blue", mid="grey90", high ="red", midpoint = 0, limits = c(-0.5,0.5)) +
      scale_radius(range = c(2,6), limits = c(1, 3), breaks = c(1, 2, 3),
                   labels = c("<0.05", "<1^-2", "<1^-3")) +
      theme_classic() +
      theme(axis.text.x=element_blank(),
            axis.text.y=element_text(size=rel(1.1)),
            axis.title.x=element_blank(),
            legend.text=element_text(size=rel(1.1)),
            legend.title=element_text(size=rel(1.1)),
            axis.title=element_text(size=rel(1.1)),
            legend.position = "right") +
      ylab("RNA inferred immune cells") +
      labs(col='Corr (Pearson)', size='FDR')
print(p)
ggsave(p, file=paste0(output.scatter, "Corr_immune_proportion_MCHII_in_tumor.svg"), 
       width = 8, height = 7, units = "cm")


##############Correlation of antigen genes with MHCII expression in tumor area ###########################

gene.epitopes <- read.table(file=paste0(input_dir,'deconv_tumor_expr_neoepitope_genes.csv'),
                            sep=",", header=TRUE)

gene.names <- gene.epitopes[,1]
exp.matrix <- gene.epitopes[,-1]
exp.df <- as.data.frame(t(exp.matrix))
colnames(exp.df) <- gene.names

#Merging epitopes with MHCII
tum_mhc$dcc_filename <- gsub("-",".",tum_mhc$dcc_filename)
tum.epitopes.mch <- merge(exp.df, tum_mhc, by.x="row.names", by.y="dcc_filename")

xname <- "ssgsea_mhc_tum"
list.plots <- list()

gene.names
tum.epitopes.mch[[xname]]

list.plots <- lapply(gene.names, function(gene){
  #rm NA if needed
  inp_df2 <- tum.epitopes.mch[!is.na(tum.epitopes.mch[[xname]]) & !is.na(tum.epitopes.mch[[gene]]), ]
  #manual_colours <- brewer.pal(length(unique(as.factor(inp_df[[color_colname]]))), 'Paired')
  M1 <- lm(get(gene) ~ exp(get(xname)) + Patient, data = inp_df2)
  pval <- glance(M1)$p.value
  
  pval.text = case_when(pval < 1e-12 ~ "P<1e-10",
                        pval < 1e-6 ~ "P<1e-10",
                        pval < 1e-3 ~ "P<1e-10",
                        pval < 0.05 ~ "P<0.05",
                        pval >= 0.05 ~ "P=N.S",
                        .default = "Non significant")
  my_ylab = strsplit(gene,"_")[[1]][1]
  
  max.y <- max(inp_df2[[gene]]) - sd(inp_df2[[gene]]) * 0.7
  
  p <- ggplot(data = inp_df2, aes(x = get(xname), y = get(gene))) +
    geom_point(alpha = 0.3) + 
    annotate("text",
             label=paste("lm(y ~ exp(x) + Patient) :", "R^2=", round(summary(M1)$r.squared, 2)),
             x=0.73, y=max.y, size = unit(3.5, "pt")) +
    xlab("Deconv. CC MHCII signature (Ucell)") +
    ylab(paste0(my_ylab, " expression")) +
    theme_bw() +
    theme(plot.title = element_text(size=rel(0.9)),
          axis.title.y = element_text(size=rel(1.4)),
          axis.text = element_text(size=rel(1.2)),
          axis.title.x = element_text(size=rel(1.4))) +
    #scale_color_manual(values=manual_colours) +
    #geom_smooth(method='glm', formula= y~x ) +
    stat_poly_line(method = 'lm', formula = y~exp(x)) +
    stat_cor(method = "pearson")
  
  ggsave(p, filename=paste0(output.scatter,'/scatter_', xname, '_vs_', gene, '_tumor_deconvoluted_signal.svg'),
         width = 9, height = 9, units = "cm")
  return(p)
  
})
ggsave(file=paste0(output.scatter, "Antigen_genes_plots.svg"), arrangeGrob(grobs = list.plots, ncol = 3),
       width = 26, height = 54, units = "cm")


##Calculating correlation between expression of antigen genes vs MCHII in tumor-TSI
corrs <- sapply(gene.names, function(gene){
  #rm NA if needed
  inp_df2 <- tum.epitopes.mch[!is.na(tum.epitopes.mch[[xname]]) & !is.na(tum.epitopes.mch[[gene]]), ]
  #Correlation pearson
  cor.res <- cor.test(inp_df2[[gene]], inp_df2[[xname]], method = "pearson")
  out <- c(gene, cor.res$p.value, cor.res$estimate)
  names(out) <- c("Gene","P.val","Cor")
  return(out)
})

corrs.df <- as.data.frame(t(corrs))
corrs.df$Cor <- as.numeric(corrs.df$Cor)
corrs.df$P.val <- as.numeric(corrs.df$P.val)
corrs.df$fdr <- p.adjust(corrs.df$P.val)
corrs.df$class <- "Fold"

corrs.df$Col <- ifelse(corrs.df$Cor >0, "Pos","Neg")
corrs.df$Col[corrs.df$fdr >= 0.05] <- "NULL"

#For figure 5h
p <- ggplot(corrs.df,aes(x=reorder(Gene, Cor), y=Cor, fill=Col)) + 
        geom_bar(stat = "identity") + coord_flip() +
        scale_fill_manual(values=c("#6599cdff","#d9dee4ff","#d33f49ff"), guide = "none") +
        theme_bw() +
        xlab("Pearson correlation") + ylab("Antigen genes")
ggsave(p, file=paste0(output.scatter, "Antigen_genes_barplot.png"), 
       width = 8, height = 10, units = "cm")