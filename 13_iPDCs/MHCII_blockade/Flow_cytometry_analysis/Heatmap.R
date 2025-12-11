# Author: Kürşat Birgin
# Date: 26./28. November 2025
# Description: Script for Anti-HLA blocking experiment for Fernando

# Install pheatmap if not already installed

if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")

# ===================================================

# Libraries

# ===================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)



# ===================================================

# Folder path

# ===================================================

folder_path <- "E:/F_Revision_Exp/Main_Blocking_Exp/iPDC_Revision_mainExp_day5"

# ===================================================

# List CSV files (all cell types)

# ===================================================

files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)

# ===================================================

# Initialize empty data frame

# ===================================================

marker_data <- data.frame()

# ===================================================

# Patient display mapping

# ===================================================

patient_map <- c("2" = "S134", "3" = "S102", "4" = "S190")

# ===================================================

# Loop over files

# ===================================================

for (file in files) {
  
  df <- read.csv(file, sep = ";", stringsAsFactors = FALSE)
  colnames(df) <- trimws(colnames(df))
  
  # Extract condition, patient, cell type from filename
  filename <- basename(file)
  
  # Condition: everything before "_Patient"
  condition <- sub("_Patient.*", "", filename)
  
  # Patient number: digits after "_Patient"
  patient <- sub(".*_Patient([0-9]+)_.*", "\\1", filename)
  
  # Cell type: between patient ID and ".csv"
  cell_type <- sub(".*_Patient[0-9]+_(.*)\\.csv", "\\1", filename)
  
  
  # Identify marker columns
  
  ifng_col <- colnames(df)[grepl("^Alexa\\.Fluor\\.700", colnames(df))]
  ki67_col <- colnames(df)[grepl("^PE\\.Cy7", colnames(df))]
  grzb_col <- colnames(df)[grepl("^APC\\.A", colnames(df))]
  
  # Build temp data frame
  
  n_rows <- nrow(df)
  
  temp <- data.frame(
    IFNg = if(length(ifng_col) == 1) df[[ifng_col]] else rep(NA, n_rows),
    Ki67 = if(length(ki67_col) == 1) df[[ki67_col]] else rep(NA, n_rows),
    GrzB = if(length(grzb_col) == 1) df[[grzb_col]] else rep(NA, n_rows),
    Condition = rep(condition, n_rows),
    Patient = rep(patient, n_rows),
    Patient_display = rep(patient_map[patient], n_rows),
    CellType = rep(cell_type, n_rows),
    stringsAsFactors = FALSE
  )
  
  
  marker_data <- rbind(marker_data, temp)
}

# ===================================================

# Convert factors and rename conditions

# ===================================================

marker_data$Condition <- factor(marker_data$Condition,
                                levels = c("Untreated", "Pembrolizumab", "AntiHLAPembro"))
levels(marker_data$Condition) <- c("Untreated", "Pembrolizumab", "Anti-HLA + Pembrolizumab")

marker_data$Patient_display <- factor(marker_data$Patient_display)
marker_data$CellType <- factor(marker_data$CellType)

# ===================================================

# Define markers

# ===================================================

markers <- c("IFNg", "Ki67", "GrzB")

# ===================================================

# Compute median per patient x cell type x marker x condition

# ===================================================

median_data <- marker_data %>%
  pivot_longer(cols = all_of(markers), names_to = "Marker", values_to = "Value") %>%
  group_by(Patient_display, CellType, Marker, Condition) %>%
  summarise(Median = median(Value, na.rm = TRUE), .groups = "drop")

# ===================================================

# Compute log2 fold changes

# ===================================================

# Option 1: Keep values as is, but rename columns to match the computation
fold_changes <- median_data %>%
  pivot_wider(names_from = Condition, values_from = Median) %>%
  mutate(
    `Pembrolizumab vs Untreated` = log2(Pembrolizumab / Untreated),
    `Anti-HLA + Pembrolizumab vs Pembrolizumab` = log2(`Anti-HLA + Pembrolizumab` / Pembrolizumab),
    `Anti-HLA + Pembrolizumab vs Untreated` = log2(`Anti-HLA + Pembrolizumab` / Untreated)
  ) %>%
  select(Patient_display, CellType, Marker,
         `Pembrolizumab vs Untreated`, 
         `Anti-HLA + Pembrolizumab vs Pembrolizumab`, 
         `Anti-HLA + Pembrolizumab vs Untreated`)

# Rename fold-change columns for nicer x-axis labels
colnames(fold_changes)[colnames(fold_changes) == "Untreated_vs_Pembro"] <- "Untreated vs Pembrolizumab"
colnames(fold_changes)[colnames(fold_changes) == "Pembro_vs_AntiHLA_Pembro"] <- "Pembrolizumab vs Anti-HLA + Pembrolizumab"
colnames(fold_changes)[colnames(fold_changes) == "Untreated_vs_AntiHLA_Pembro"] <- "Untreated vs Anti-HLA + Pembrolizumab"



# ===================================================

# Prepare heatmap matrix

# ===================================================
library(tibble)

heatmap_matrix <- fold_changes %>%
  unite("Row", Patient_display, CellType, Marker, sep = "_") %>%
  column_to_rownames("Row") %>%
  as.matrix()


# ===================================================

# # Heatmap
# 
# # ===================================================
# 
# pheatmap(
#   heatmap_matrix,
#   color = colorRampPalette(c("blue", "white", "red"))(50),
#   cluster_rows = TRUE,
#   cluster_cols = FALSE,
#   main = "Log2 Fold Changes per Patient, Cell Type, Marker",
#   angle_col = 45
# )

###
library(pheatmap)
library(tibble)

# List of patients

patients <- unique(fold_changes$Patient_display)

# Generate one heatmap per patient

# Generate one heatmap per patient with prettier row names
# Generate one heatmap per patient with prettier row names

for (p in patients) {
  
  heatmap_data <- fold_changes %>%
    filter(Patient_display == p) %>%
    # Make row names prettier
    mutate(
      CellType_pretty = gsub("minus", "-", gsub("plus", "+", CellType)),
      RowName = paste0(Marker, " levels in ", CellType_pretty, " cells")
    ) %>%
    column_to_rownames("RowName") %>%
    select(-Patient_display, -CellType, -Marker, -CellType_pretty) %>% # remove original columns now encoded in rownames
    as.matrix()
  
  pheatmap(
    heatmap_data,
    color = colorRampPalette(c("blue", "white", "red"))(50),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = paste("Log2 Fold Changes - Patient", p),
    angle_col = 45
  )
}

###A heatmap per marker 
library(pheatmap)
library(tibble)

# List of patients and markers

patients <- unique(fold_changes$Patient_display)
markers <- c("IFNg", "Ki67", "GrzB")

# Generate heatmaps per patient AND per functional marker

for (p in patients) {
  for (m in markers) {
    
    heatmap_data <- fold_changes %>%
      filter(Patient_display == p, Marker == m) %>%
      # Make row names prettier (only cell type now, marker is in title)
      mutate(
        CellType_pretty = gsub("minus", "-", gsub("plus", "+", CellType)),
        RowName = paste0(CellType_pretty, " cells")
      ) %>%
      column_to_rownames("RowName") %>%
      select(-Patient_display, -CellType, -Marker, -CellType_pretty) %>%
      as.matrix()
    
    pheatmap(
      heatmap_data,
      color = colorRampPalette(c("blue", "white", "red"))(50),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      main = paste("Log2 Fold Changes - ", m, "levels - Patient", p),
      angle_col = 45
    )
    
  }
}

# ===================================================

# Raw fold-change heatmaps per patient and per marker

# ===================================================

# Compute raw fold changes (no log2)

fold_changes_raw <- median_data %>%
  pivot_wider(names_from = Condition, values_from = Median) %>%
  mutate(
    `Pembrolizumab vs Untreated` = Pembrolizumab / Untreated,
    `Anti-HLA + Pembrolizumab vs Pembrolizumab` = `Anti-HLA + Pembrolizumab` / Pembrolizumab,
    `Anti-HLA + Pembrolizumab vs Untreated` = `Anti-HLA + Pembrolizumab` / Untreated
  ) %>%
  select(Patient_display, CellType, Marker,
         `Pembrolizumab vs Untreated`,
         `Anti-HLA + Pembrolizumab vs Pembrolizumab`,
         `Anti-HLA + Pembrolizumab vs Untreated`)

# Patients and markers

patients <- unique(fold_changes_raw$Patient_display)
markers <- unique(fold_changes_raw$Marker)

# Loop per patient and marker

for (p in patients) {
  for (m in markers) {
    
    heatmap_data <- fold_changes_raw %>%
      filter(Patient_display == p, Marker == m) %>%
      # prettier row names: just the cell type
      mutate(
        CellType_pretty = gsub("minus", "-", gsub("plus", "+", CellType)),
        RowName = paste0(CellType_pretty, " cells")
      ) %>%
      column_to_rownames("RowName") %>%
      select(-Patient_display, -CellType, -Marker, -CellType_pretty) %>%
      as.matrix()
    
    # Optional: define breaks to center on 1 for raw fold changes
    max_val <- max(heatmap_data, na.rm = TRUE)
    min_val <- min(heatmap_data, na.rm = TRUE)
    breaks <- seq(min_val, max_val, length.out = 51)
    
    pheatmap(
      heatmap_data,
      color = colorRampPalette(c("blue", "white", "red"))(50),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      main = paste("Raw Fold Changes -", m, "- Patient", p),
      angle_col = 45,
      breaks = breaks
    )
    
  }
}


# ###combined patients for log2 plot
# 
# library(dplyr)
# library(tidyr)
# library(pheatmap)
# library(tibble)
# 
# Markers to plot

markers_to_plot <- c("IFNg", "GrzB")

# for (m in markers_to_plot) {
# 
#   heatmap_data <- fold_changes %>%
#     filter(Marker == m) %>%
#     # Create prettier row names
#     mutate(
#       CellType_pretty = gsub("minus", "-", gsub("plus", "+", CellType)),
#       # Automated replacements for CD45+ and Myeloids+
#       CellType_pretty = gsub("^CD45\\+$", "CD45", CellType_pretty),
#       CellType_pretty = gsub("^Myeloids\\+$", "Myeloid", CellType_pretty),
#       RowName = paste0(Patient_display, " - ", CellType_pretty, " cells")
#     ) %>%
#     column_to_rownames("RowName") %>%
#     select(-Patient_display, -CellType, -Marker, -CellType_pretty) %>%
#     as.matrix()
#   
#   # Define color palette centered at 0 (log2 fold-change = 0)
#   
#   max_abs <- max(abs(heatmap_data), na.rm = TRUE)
#   breaks <- seq(-max_abs, max_abs, length.out = 101)
#   colors <- colorRampPalette(c("blue", "white", "red"))(100)
#   
#   pheatmap(
#     heatmap_data,
#     color = colors,
#     breaks = breaks,
#     cluster_rows = FALSE,
#     cluster_cols = FALSE,
#     main = paste("Log2 Fold Changes -", m),
#     angle_col = 45
#   )
# }
# 
###ordering
# Desired order of cell types

celltype_order <- c(
  "CD45",
  "CD3+",
  "CD3+ CD8+",
  "CD3+ CD4-",
  "CD3+ CD4+",
  "CD3+ CD8-",
  "Myeloid"
)
# 
# for (m in markers_to_plot) {
#   
#   heatmap_data <- fold_changes %>%
#     filter(Marker == m) %>%
#     # Prettify cell types
#     mutate(
#       CellType_pretty = gsub("minus", "-", gsub("plus", "+", CellType)),
#       CellType_pretty = gsub("([+-])(?=[A-Za-z0-9])", "\\1 ", CellType_pretty, perl = TRUE),
#       CellType_pretty = gsub("^CD45\\+$", "CD45", CellType_pretty),
#       CellType_pretty = gsub("^Myeloids\\+$", "Myeloid", CellType_pretty),
#       # Create row name with patient first
#       RowName = paste0(Patient_display, " - ", CellType_pretty, " cells"),
#       # Factor for ordering by cell type only
#       CellType_order = factor(CellType_pretty, levels = celltype_order)
#     ) %>%
#     # Arrange by cell type order (row names still show patient first)
#     arrange(CellType_order) %>%
#     column_to_rownames("RowName") %>%
#     select(-Patient_display, -CellType, -Marker, -CellType_pretty, -CellType_order) %>%
#     as.matrix()
#   
#   
#   # Define color palette
#   
#   max_abs <- max(abs(heatmap_data), na.rm = TRUE)
#   breaks <- seq(-max_abs, max_abs, length.out = 101)
#   colors <- colorRampPalette(c("blue", "white", "red"))(100)
#   
#   pheatmap(
#     heatmap_data,
#     color = colors,
#     breaks = breaks,
#     cluster_rows = FALSE,
#     cluster_cols = FALSE,
#     main = paste("Log2 Fold Changes -", m),
#     angle_col = 45
#   )
# }
# 

# ###median 
# for (m in markers_to_plot) {
#   
#   heatmap_data <- pooled_fold_changes %>%
#     filter(Marker == m) %>%
#     mutate(
#       # Prettify cell types
#       CellType_pretty = gsub("minus", "-", gsub("plus", "+", CellType)),
#       CellType_pretty = gsub("([+-])(?=[A-Za-z0-9])", "\\1 ", CellType_pretty, perl = TRUE),
#       CellType_pretty = gsub("^CD45\\+$", "CD45", CellType_pretty),
#       CellType_pretty = gsub("^Myeloids\\+$", "Myeloid", CellType_pretty),
#       
#       # Create "RowName" WITH "cells" like before
#       RowName = paste0(CellType_pretty, " cells"),
#       
#       # Ordering factor
#       CellType_order = factor(CellType_pretty, levels = celltype_order)
#     ) %>%
#     arrange(CellType_order) %>%
#     column_to_rownames("RowName") %>%   # use pretty label with "cells"
#     select(-CellType, -Marker, -CellType_pretty, -CellType_order) %>%
#     as.matrix()
#   
#   # Clean any invisible characters
#   rownames(heatmap_data) <- gsub("[^[:print:]]", "", rownames(heatmap_data))
#   
#   # Color configuration
#   max_abs <- max(abs(heatmap_data), na.rm = TRUE)
#   breaks <- seq(-max_abs, max_abs, length.out = 101)
#   colors <- colorRampPalette(c("blue", "white", "red"))(100)
#   
#   pheatmap(
#     heatmap_data,
#     color = colors,
#     breaks = breaks,
#     cluster_rows = FALSE,
#     cluster_cols = FALSE,
#     main = paste("Median log2 fold-change across patients -", m),
#     angle_col = 45
#   )
# }

###mean instead of median

pooled_fold_changes_mean <- fold_changes %>%
  group_by(CellType, Marker) %>%
  summarise(across(
    c(`Pembrolizumab vs Untreated`,
      `Anti-HLA + Pembrolizumab vs Pembrolizumab`,
      `Anti-HLA + Pembrolizumab vs Untreated`),
    ~ mean(.x, na.rm = TRUE)
  ))

for (m in markers_to_plot) {
  
  heatmap_data <- pooled_fold_changes_mean %>%   # <-- now using MEAN
    filter(Marker == m) %>%
    mutate(
      CellType_pretty = gsub("minus", "-", gsub("plus", "+", CellType)),
      CellType_pretty = gsub("([+-])(?=[A-Za-z0-9])", "\\1 ", CellType_pretty, perl = TRUE),
      CellType_pretty = gsub("^CD45\\+$", "CD45", CellType_pretty),
      CellType_pretty = gsub("^Myeloids\\+$", "Myeloid", CellType_pretty),
      RowName = paste0(CellType_pretty, " cells"),
      CellType_order = factor(CellType_pretty, levels = celltype_order)
    ) %>%
    arrange(CellType_order) %>%
    column_to_rownames("RowName") %>%
    select(-CellType, -Marker, -CellType_pretty, -CellType_order) %>%
    as.matrix()
  
  rownames(heatmap_data) <- gsub("[^[:print:]]", "", rownames(heatmap_data))
  
  max_abs <- max(abs(heatmap_data), na.rm = TRUE)
  breaks <- seq(-max_abs, max_abs, length.out = 101)
  colors <- colorRampPalette(c("blue", "white", "red"))(100)
  
  pheatmap(
    heatmap_data,
    color = colors,
    breaks = breaks,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = paste("Mean log2 fold-change across patients -", m),
    angle_col = 45
  )
}

# Compute mean log2 fold-changes across patients
pooled_fold_changes_mean <- fold_changes %>%
  group_by(CellType, Marker) %>%
  summarise(across(
    c(`Pembrolizumab vs Untreated`,
      `Anti-HLA + Pembrolizumab vs Pembrolizumab`,
      `Anti-HLA + Pembrolizumab vs Untreated`),
    ~ mean(.x, na.rm = TRUE)
  ))

# Rename columns for shorter x-axis labels
pooled_fold_changes_mean_short <- pooled_fold_changes_mean %>%
  rename(
    "PE_vs_UT" = `Pembrolizumab vs Untreated`,
    "COMBO_vs_PE" = `Anti-HLA + Pembrolizumab vs Pembrolizumab`,
    "COMBO_vs_UT" = `Anti-HLA + Pembrolizumab vs Untreated`
  )

for (m in markers_to_plot) {
  
  heatmap_data <- pooled_fold_changes_mean_short %>%
    filter(Marker == m) %>%
    mutate(
      CellType_pretty = gsub("minus", "-", gsub("plus", "+", CellType)),
      CellType_pretty = gsub("([+-])(?=[A-Za-z0-9])", "\\1 ", CellType_pretty, perl = TRUE),
      CellType_pretty = gsub("^CD45\\+$", "CD45", CellType_pretty),
      CellType_pretty = gsub("^Myeloids\\+$", "Myeloid", CellType_pretty),
      RowName = paste0(CellType_pretty, " cells"),
      CellType_order = factor(CellType_pretty, levels = celltype_order)
    ) %>%
    arrange(CellType_order) %>%
    column_to_rownames("RowName") %>%
    select(-CellType, -Marker, -CellType_pretty, -CellType_order) %>%
    as.matrix()
  
  # Clean invisible characters
  rownames(heatmap_data) <- gsub("[^[:print:]]", "", rownames(heatmap_data))
  
  # Define color scale with custom midpoint
  max_abs <- max(abs(heatmap_data), na.rm = TRUE)
  breaks <- seq(-max_abs, max_abs, length.out = 101)
  colors <- colorRampPalette(c("#508791", "white", "#D33F49"))(100)
  
  # Plot heatmap
  pheatmap(
    heatmap_data,
    color = colors,
    breaks = breaks,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = paste("Mean log2 fold-change across patients -", m),
    angle_col = 45,
    fontsize = 18,        # title & legend
    fontsize_row = 18,    # y-axis labels
    fontsize_col = 18,    # x-axis labels
    display_numbers = FALSE,
    legend = TRUE,
    border_color = "black"
  )
}


