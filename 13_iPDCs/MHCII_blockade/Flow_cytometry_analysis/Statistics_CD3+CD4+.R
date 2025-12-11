# Author: Kürşat Birgin
# Date: 26./28. November 2025
# Description: Script for Anti-HLA blocking experiment for Fernando

#CD3+CD4+

# ===================================================
# Load necessary library for plotting
# ===================================================
library(ggplot2)
if (!requireNamespace("ggpubr", quietly = TRUE)) install.packages("ggpubr")
library(ggpubr)

# ===================================================
# Set the folder path where CSV files are stored
# ===================================================
folder_path <- "E:/F_Revision_Exp/Main_Blocking_Exp/iPDC_Revision_mainExp_day5"

# ===================================================
# List all CSV files
# ===================================================
files <- list.files(folder_path, pattern = "CD3plusCD4plus\\.csv$", full.names = TRUE)

# ===================================================
# Initialize empty data frame
# ===================================================
marker_data <- data.frame()

# ===================================================
# Loop over each file and extract markers
# ===================================================
for (file in files) {
  
  df <- read.csv(file, sep = ";", stringsAsFactors = FALSE)
  colnames(df) <- trimws(colnames(df))
  
  # Match columns based on fluorescence channel start
  ifng_col <- colnames(df)[grepl("^Alexa\\.Fluor\\.700", colnames(df))]
  ki67_col <- colnames(df)[grepl("^PE\\.Cy7", colnames(df))]
  grzb_col <- colnames(df)[grepl("^APC\\.A", colnames(df))]
  
  # Extract condition and patient
  filename  <- basename(file)
  condition <- sub("_Patient.*", "", filename)
  patient   <- sub(".*_Patient([0-9]+)_.*", "\\1", filename)
  
  # Build temporary data frame with all markers
  temp <- data.frame(
    IFNg      = if(length(ifng_col) == 1) df[[ifng_col]] else NA,
    Ki67      = if(length(ki67_col) == 1) df[[ki67_col]] else NA,
    GrzB      = if(length(grzb_col) == 1) df[[grzb_col]] else NA,
    Condition = condition,
    Patient   = patient
  )
  
  marker_data <- rbind(marker_data, temp)
}

# ===================================================
# Convert factors
# ===================================================
marker_data$Condition <- factor(marker_data$Condition, levels = c("Untreated", "Pembrolizumab", "AntiHLAPembro"))
marker_data$Patient   <- factor(marker_data$Patient)

# ===================================================
# Rename conditions for plotting only (nice display names)
# ===================================================
levels(marker_data$Condition) <- c(
  "Untreated",
  "Pembrolizumab",
  "Anti-HLA + Pembrolizumab"
)


# ===================================================
# Define markers to process
# ===================================================
markers <- c("IFNg", "Ki67", "GrzB")

# ===================================================
# Loop over markers for stats and plots
# ===================================================
for (marker in markers) {
  
  cat("\n====================", marker, "====================\n")
  
  # Subset by condition
  untreated <- marker_data[[marker]][marker_data$Condition == "Untreated"]
  pembrolizumab    <- marker_data[[marker]][marker_data$Condition == "Pembrolizumab"]
  antihlapembro   <- marker_data[[marker]][marker_data$Condition == "Anti-HLA + Pembrolizumab"]
  
  # Skip if non-numeric or all NA
  if (!is.numeric(untreated) || all(is.na(untreated))) next
  
  # Mann-Whitney U tests
  cat("Pembrolizumab vs Untreated:", wilcox.test(pembrolizumab, untreated)$p.value, "\n")
  cat("AntiHLAPembro vs Untreated:", wilcox.test(antihlapembro, untreated)$p.value, "\n")
  cat("Pembrolizumab vs Anti-HLA + Pembrolizumab:", wilcox.test(pembrolizumab, antihlapembro)$p.value, "\n")
  
  # # Histogram
  # p_hist <- ggplot(marker_data, aes_string(x = marker, fill = "Condition")) +
  #   geom_histogram(bins = 30, alpha = 0.6, position = "identity", color = "black") +
  #   labs(title = paste(marker, "histogram (CD3+CD4+ cells, patients pooled)"),
  #        x = marker, y = "Count") +
  #   theme_minimal()
  # print(p_hist)
  
  # Violin + boxplot
  p_violin <- ggplot(marker_data, aes_string(x = "Condition", y = marker, fill = "Condition")) +
    geom_violin(trim = FALSE, scale = "area") +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.4) +
    stat_summary(fun = median, geom = "point", color = "black", size = 2) +
    labs(title = paste(marker, "distribution in CD3+CD4+ cells (patients pooled)"),
         x = "Condition", y = marker) +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2")
  print(p_violin)
  
  # Boxplot with p-values
  p_box <- ggplot(marker_data, aes_string(x = "Condition", y = marker)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    stat_compare_means(
      comparisons = list(
        c("Untreated", "Pembrolizumab"),
        c("Untreated", "Anti-HLA + Pembrolizumab"),
        c("Pembrolizumab", "Anti-HLA + Pembrolizumab")
      ),
      method = "wilcox.test",
      label = "p.signif",   # show significance asterisks
      hide.ns = TRUE        # hide non-significant p-values
    ) +
    labs(title = paste(marker, "levels in CD3+CD4+ cells (Patients pooled)"), y = marker) +
    theme_minimal()
  
  # ---- SAVE PLOT ----
  plot_title <- paste(marker, "levels in CD3+CD4+ cells (Patients pooled)")
  file_name <- paste0(gsub("[^A-Za-z0-9_\\-]", "_", plot_title), ".png")
  file_path <- file.path("E:/F_Revision_Exp/Main_Blocking_Exp/Plots/Statistics Plots pvalue", file_name)
  png(filename = file_path, width = 500, height = 500)
  print(p_box)
  dev.off()
  
  # Also display in R
  print(p_box)
}



# ===================================================
# --- Add-on: Per-patient statistics and plots ---
# ===================================================

patients <- levels(marker_data$Patient)  # get all patient IDs

# Map numeric patient IDs to display names
patient_map <- c(
  "2" = "P1",
  "3" = "P2",
  "4" = "P3"
)

# Add display column to marker_data
marker_data$Patient_display <- patient_map[as.character(marker_data$Patient)]


for (patient_id in patients) {
  
  cat("\n==================== Patient", patient_id, "====================\n")
  
  # Subset data for this patient only
  patient_data <- marker_data[marker_data$Patient == patient_id, ]
  
  for (marker in markers) {
    
    cat("\n---", marker, "---\n")
    
    # Subset by condition for this patient
    untreated <- patient_data[[marker]][patient_data$Condition == "Untreated"]
    pembrolizumab    <- patient_data[[marker]][patient_data$Condition == "Pembrolizumab"]
    antihlapembro   <- patient_data[[marker]][patient_data$Condition == "Anti-HLA + Pembrolizumab"]
    
    # Skip if non-numeric or all NA
    if (!is.numeric(untreated) || all(is.na(untreated))) next
    
    # Mann-Whitney U tests per patient
    cat("Pembrolizumab vs Untreated:", wilcox.test(pembrolizumab, untreated)$p.value, "\n")
    cat("AntiHLAPembro vs Untreated:", wilcox.test(antihlapembro, untreated)$p.value, "\n")
    cat("Pembrolizumab vs AntiHLAPembro:", wilcox.test(pembrolizumab, antihlapembro)$p.value, "\n")
    
    # Violin + boxplot per patient
    #    p_violin <- ggplot(patient_data, aes_string(x = "Condition", y = marker, fill = "Condition")) +
    #      geom_violin(trim = FALSE, scale = "area") +
    #      geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.4) +
    #      stat_summary(fun = median, geom = "point", color = "black", size = 2) +
    #      labs(title = paste(marker, "distribution (Patient", patient_id, ")"),
    #           x = "Condition", y = marker) +
    #      theme_minimal() +
    #      scale_fill_brewer(palette = "Set2")
    #    print(p_violin)
    
    # Optional: boxplot with p-values
    p_box <- ggplot(patient_data, aes_string(x = "Condition", y = marker)) +
      geom_boxplot() +
      geom_jitter(width = 0.2, alpha = 0.5) +
      stat_compare_means(
        comparisons = list(
          c("Untreated", "Pembrolizumab"),
          c("Untreated", "Anti-HLA + Pembrolizumab"),
          c("Pembrolizumab", "Anti-HLA + Pembrolizumab")
        ),
        method = "wilcox.test",
        label = "p.signif",   # show significance asterisks
        hide.ns = TRUE        # hide non-significant p-values
      ) +
      labs(title = paste0(marker, " levels in CD3+CD4+ cells (Patient ", patient_data$Patient_display[1], ")"), y = marker) +
      theme_minimal()
    
    # ---- SAVE PLOT ----
    plot_title <- paste0(marker, " levels in CD3+CD4+ cells (Patient ", patient_data$Patient_display[1], ")")
    file_name <- paste0(gsub("[^A-Za-z0-9_\\-]", "_", plot_title), ".png")
    file_path <- file.path("E:/F_Revision_Exp/Main_Blocking_Exp/Plots/Statistics Plots pvalue", file_name)
    png(filename = file_path, width = 500, height = 500)
    print(p_box)
    dev.off()
    
    # Also display in R
    print(p_box)
  }
}

# # unique(marker_data$Condition)
# # 
# # ###CHECKING Table content 
# # unique(marker_data$Patient)
# # # See all files and what patient is extracted from each
# # for (file in files) {
# #   filename <- basename(file)
# #   patient <- sub(".*_Patient([0-9]+)_.*", "\\1", filename)
# #   print(c(filename, patient))
# }
