###############################################################
# LOAD REQUIRED PACKAGES
###############################################################
library(tidyverse)
library(readxl)

###############################################################
# 1) READ ALL CYTOKINE SHEETS FROM EXCEL INTO ONE DATAFRAME
###############################################################

file <- "samples_undiluted_concentrations_25112025.xlsx"   # <-- change if needed

all_data <- map_df(
  excel_sheets(file),
  ~ read_excel(file, sheet = .x) %>%
    mutate(cytokine = .x)
)

###############################################################
# 2) PARSE SAMPLE NAMES (patient / condition / replicate)
###############################################################

all_data <- all_data %>%
  mutate(
    patient = str_extract(sample, "S\\d+"),
    condition = str_extract(sample, "(UT|PEMBRO|aHLA|COMBO)"),
    replicate_num = as.integer(str_extract(sample, "\\d+$"))
  )

# Correct PEMBRO/UT switch
all_data$condition <- as.character(all_data$condition)
all_data[all_data$condition == "UT", "condition"] <- "PEMBRO_corrected"
all_data[all_data$condition == "PEMBRO", "condition"] <- "UT"
all_data[all_data$condition == "PEMBRO_corrected", "condition"] <- "PEMBRO"
all_data$condition <- as.factor(all_data$condition)

###############################################################
# 3) SET ORDER OF CONDITIONS
###############################################################

condition_order <- c("UT", "PEMBRO", "aHLA", "COMBO")
all_data$condition <- factor(all_data$condition, levels = condition_order)

###############################################################
# 3B — FILTER OPTIONS (CUSTOMIZABLE)
###############################################################

# ---- Choose CYTOKINES to include ----
# Use "ALL" to keep all cytokines, or a vector:
cytokines_to_plot <- c("Perforin","IL-6","IFN-g", "Granzyme-B", "Granzyme-A",
                       "Granolosyn", "FasL", "Fas")
# Example: cytokines_to_plot <- c("IL-2","TNF-alpha","Granzyme-B")

# ---- Choose CONDITIONS to include for RAW + AVG PLOTS ----
conditions_to_plot <- c("UT", "PEMBRO", "COMBO")
# Example: conditions_to_plot <- c("UT","PEMBRO","COMBO")

# ---- Apply filtering ----
if (!identical(cytokines_to_plot, "ALL")) {
  all_data <- all_data %>% filter(cytokine %in% cytokines_to_plot)
}

all_data <- all_data %>% filter(condition %in% conditions_to_plot)
all_data$condition <- factor(all_data$condition, levels = condition_order)

###############################################################
# 4) CREATE AVERAGED DATAFRAME (PER PATIENT × CONDITION)
###############################################################

avg_data <- all_data %>%
  group_by(cytokine, patient, condition) %>%
  summarise(
    concentration_avg = mean(predicted_concentration, na.rm = TRUE),
    MFI_avg = mean(MFI, na.rm = TRUE),
    N = n(),
    .groups = "drop"
  )

avg_data <- avg_data %>%
  filter(condition %in% conditions_to_plot)

avg_data$condition <- factor(avg_data$condition, levels = condition_order)

###############################################################
# 5) PLOTS — RAW DUPLICATES (jittered points)
###############################################################

unique_cytokines <- unique(all_data$cytokine)
plot_list_raw <- list()

for (cyt in unique_cytokines) {
  
  df <- all_data %>% filter(cytokine == cyt)
  
  p <- ggplot(df,
              aes(x = condition,
                  y = predicted_concentration,
                  color = patient)) +
    geom_point(size = 3, alpha = 0.9,
               position = position_jitter(width = 0.15)) +
    theme_bw(base_size = 14) +
    labs(title = paste("RAW replicates –", cyt),
         x = "Condition",
         y = "Predicted Concentration") +
    theme(plot.title = element_text(face = "bold"))
  
  print(p)
  plot_list_raw[[cyt]] <- p
}

###############################################################
# 6) PLOTS — AVERAGED REPLICATES (jittered means)
###############################################################

plot_list_avg <- list()

for (cyt in unique(avg_data$cytokine)) {
  
  df <- avg_data %>% filter(cytokine == cyt)
  
  p <- ggplot(df,
              aes(x = condition,
                  y = concentration_avg,
                  color = patient)) +
    geom_point(size = 4,
               position = position_jitter(width = 0.15)) +
    theme_bw(base_size = 14) +
    labs(title = paste("AVERAGED replicates –", cyt),
         x = "Condition",
         y = "Average Predicted Concentration") +
    theme(plot.title = element_text(face = "bold"))
  
  print(p)
  plot_list_avg[[cyt]] <- p
}

###############################################################
# 7) FOLD-CHANGE CALCULATIONS (log2) — based on avg_data
###############################################################

fc_data <- avg_data %>%
  select(cytokine, patient, condition, concentration_avg) %>%
  pivot_wider(
    names_from = condition,
    values_from = concentration_avg
  ) %>%
  mutate(
    log2FC_PE_vs_UT = log2(PEMBRO / UT),
    log2FC_CO_vs_UT = log2(COMBO / UT),
    log2FC_CO_vs_PE = log2(COMBO / PEMBRO)
  ) %>%
  select(cytokine, patient,
         log2FC_PE_vs_UT, log2FC_CO_vs_UT, log2FC_CO_vs_PE)

# filter cytokines for heatmaps too
if (!identical(cytokines_to_plot, "ALL")) {
  fc_data <- fc_data %>% filter(cytokine %in% cytokines_to_plot)
}

###############################################################
# 8) HEATMAPS — PER PATIENT
###############################################################

patients <- unique(fc_data$patient)
fc_order <- c("log2FC_PE_vs_UT", "log2FC_CO_vs_PE", "log2FC_CO_vs_UT")

plot_list_heatmap <- list()

for (pat in patients) {
  
  df <- fc_data %>%
    filter(patient == pat) %>%
    pivot_longer(
      cols = starts_with("log2FC"),
      names_to = "Comparison",
      values_to = "log2FC"
    ) %>%
    mutate(Comparison = factor(Comparison, levels = fc_order))
  
  p <- ggplot(df, aes(x = Comparison, y = cytokine, fill = log2FC)) +
    geom_tile(color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    
    # >>> paste this block <<<
    scale_x_discrete(labels = c(
      "log2FC_PE_vs_UT" = "PE_vs_UT",
      "log2FC_CO_vs_PE" = "COMBO_vs_PE",
      "log2FC_CO_vs_UT" = "COMBO_vs_UT"
    )) +
    coord_fixed() +
    # <<< end of block <<<
    
    theme_bw(base_size = 16) +
    labs(
      title = paste("Fold-change heatmap –", pat),
      x = "Comparison",
      y = "Cytokine"
    ) +
    theme(
      axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 18),
      plot.title = element_text(face = "bold")
    )
  
  print(p)
  plot_list_heatmap[[pat]] <- p
}

###############################################################
# 9) AVERAGED FOLD-CHANGE HEATMAP (mean across patients)
###############################################################

fc_avg <- fc_data %>%
  group_by(cytokine) %>%
  summarise(
    log2FC_PE_vs_UT = mean(log2FC_PE_vs_UT, na.rm = TRUE),
    log2FC_CO_vs_UT = mean(log2FC_CO_vs_UT, na.rm = TRUE),
    log2FC_CO_vs_PE = mean(log2FC_CO_vs_PE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = starts_with("log2FC"),
    names_to = "Comparison",
    values_to = "log2FC"
  ) %>%
  mutate(Comparison = factor(Comparison, levels = fc_order))


p_avgFC <- ggplot(fc_avg, aes(x = Comparison, y = cytokine, fill = log2FC)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "#508791", mid = "white", high = "#D33F49", midpoint = 0) +
  
  scale_x_discrete(labels = c(
    "log2FC_PE_vs_UT" = "PE_vs_UT",
    "log2FC_CO_vs_PE" = "COMBO_vs_PE",
    "log2FC_CO_vs_UT" = "COMBO_vs_UT"
  )) +
  coord_fixed() +
  
  theme_bw(base_size = 18) +
  labs(
    title = "Average fold-change heatmap (mean across patients)",
    x = "Comparison", 
    y = "Cytokine"
  ) +
  theme(
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 18),
    plot.title = element_text(face = "bold")
  )

print(p_avgFC)

plot_list_heatmap_avg <- list(AvgFC = p_avgFC)

###############################################################
# 10) SAVE PLOTS TO FOLDERS (OPTIONAL)
###############################################################

dir.create("plots_raw_switched_PE_UT_selected_data", showWarnings = FALSE)
dir.create("plots_avg_switched_PE_UT_selected_data", showWarnings = FALSE)
dir.create("plots_heatmaps_switched_PE_UT_selected_data", showWarnings = FALSE)

# RAW
for (name in names(plot_list_raw)) {
  ggsave(
    filename = paste0("plots_raw_switched_PE_UT_selected_data/", name, ".png"),
    plot = plot_list_raw[[name]], width = 6, height = 5
  )
}

# AVERAGED
for (name in names(plot_list_avg)) {
  ggsave(
    filename = paste0("plots_avg_switched_PE_UT_selected_data/", name, ".png"),
    plot = plot_list_avg[[name]], width = 6, height = 5
  )
}

# HEATMAPS PER PATIENT
for (name in names(plot_list_heatmap)) {
  ggsave(
    filename = paste0("plots_heatmaps_switched_PE_UT_selected_data/", name, ".png"),
    plot = plot_list_heatmap[[name]], width = 6, height = 5
  )
}

# AVERAGED HEATMAP
dir.create("plots_heatmaps_avgFC_switched_PE_UT_selected_data", showWarnings = FALSE)

ggsave(
  filename = "plots_heatmaps_avgFC_switched_PE_UT_selected_data/AvgFC_new_color_palette.png",
  plot = plot_list_heatmap_avg[[1]],
  width = 8, height = 5
)

ggsave(
  filename = "plots_heatmaps_avgFC_switched_PE_UT_selected_data/AvgFC_new_color_palette.pdf",
  plot = plot_list_heatmap_avg[[1]],
  width = 8, height = 5,
  device = cairo_pdf               # ensures crisp text & vector output
)

###############################################################
# 11) STATISTICS — Paired Wilcoxon tests per cytokine
#     Comparing fold-changes across 3 patients (paired)
#
#     Three comparisons (all paired):
#       A) PE vs UT        vs   CO vs PE
#       B) PE vs UT        vs   CO vs UT
#       C) CO vs PE        vs   CO vs UT
###############################################################

# Prepare a list of cytokines
cytokines_fc <- unique(fc_data$cytokine)

# Empty dataframes to store results
stats_PEUT_vs_COPE <- data.frame()
stats_PEUT_vs_COUT <- data.frame()
stats_COPE_vs_COUT <- data.frame()

for (cyt in cytokines_fc) {
  
  df <- fc_data %>% filter(cytokine == cyt)
  
  # Extract paired values
  x_PEUT <- df$log2FC_PE_vs_UT
  x_COPE <- df$log2FC_CO_vs_PE
  x_COUT <- df$log2FC_CO_vs_UT
  
  # A) PE vs UT  vs  CO vs PE
  test_A <- tryCatch(
    wilcox.test(x_PEUT, x_COPE, paired = TRUE, exact = FALSE),
    error = function(e) return(NA)
  )
  
  # B) PE vs UT  vs  CO vs UT
  test_B <- tryCatch(
    wilcox.test(x_PEUT, x_COUT, paired = TRUE, exact = FALSE),
    error = function(e) return(NA)
  )
  
  # C) CO vs PE  vs  CO vs UT
  test_C <- tryCatch(
    wilcox.test(x_COPE, x_COUT, paired = TRUE, exact = FALSE),
    error = function(e) return(NA)
  )
  
  # Store results
  stats_PEUT_vs_COPE <- rbind(
    stats_PEUT_vs_COPE,
    data.frame(
      cytokine = cyt,
      p_value = ifelse(is.na(test_A), NA, test_A$p.value)
    )
  )
  
  stats_PEUT_vs_COUT <- rbind(
    stats_PEUT_vs_COUT,
    data.frame(
      cytokine = cyt,
      p_value = ifelse(is.na(test_B), NA, test_B$p.value)
    )
  )
  
  stats_COPE_vs_COUT <- rbind(
    stats_COPE_vs_COUT,
    data.frame(
      cytokine = cyt,
      p_value = ifelse(is.na(test_C), NA, test_C$p.value)
    )
  )
}

# Add FDR-corrected p-values
stats_PEUT_vs_COPE$FDR <- p.adjust(stats_PEUT_vs_COPE$p_value, method = "BH")
stats_PEUT_vs_COUT$FDR <- p.adjust(stats_PEUT_vs_COUT$p_value, method = "BH")
stats_COPE_vs_COUT$FDR <- p.adjust(stats_COPE_vs_COUT$p_value, method = "BH")

###############################################################
# PRINT RESULTS
###############################################################

cat("\n\n### WILCOXON RESULTS — PE_vs_UT vs CO_vs_PE ###\n")
print(stats_PEUT_vs_COPE)

cat("\n\n### WILCOXON RESULTS — PE_vs_UT vs CO_vs_UT ###\n")
print(stats_PEUT_vs_COUT)

cat("\n\n### WILCOXON RESULTS — CO_vs_PE vs CO_vs_UT ###\n")
print(stats_COPE_vs_COUT)

###############################################################
# OPTIONAL — SAVE RESULTS TO CSV
# (uncomment to save)
###############################################################

# dir.create("stats_results", showWarnings = FALSE)
# write.csv(stats_PEUT_vs_COPE, "stats_results/Wilcoxon_PEvsUT_vs_COvsPE.csv", row.names = FALSE)
# write.csv(stats_PEUT_vs_COUT, "stats_results/Wilcoxon_PEvsUT_vs_COvsUT.csv", row.names = FALSE)
# write.csv(stats_COPE_vs_COUT, "stats_results/Wilcoxon_COvsPE_vs_COvsUT.csv", row.names = FALSE)

###############################################################
# 11) STATISTICS — Paired Wilcoxon tests on AVERAGED concentrations
# Comparisons requested:
#   1) UT vs PEMBRO
#   2) UT vs COMBO
#   3) COMBO vs PEMBRO
#
# Uses avg_data (one value per patient × condition → n = 3)
###############################################################

# Reshape avg_data into wide format: one row per patient × cytokine
avg_wide <- avg_data %>%
  select(cytokine, patient, condition, concentration_avg) %>%
  pivot_wider(
    names_from = condition,
    values_from = concentration_avg
  )

# Function to run paired Wilcoxon test safely (returns NA if missing data)
safe_wilcox <- function(x, y) {
  tryCatch(
    {
      if (all(!is.na(x)) & all(!is.na(y))) {
        wilcox.test(x, y, paired = TRUE)$p.value
      } else {
        NA
      }
    },
    error = function(e) NA
  )
}

# Run tests per cytokine
stats_conc <- avg_wide %>%
  group_by(cytokine) %>%
  summarise(
    p_UT_vs_PEMBRO = safe_wilcox(UT, PEMBRO),
    p_UT_vs_COMBO  = safe_wilcox(UT, COMBO),
    p_COMBO_vs_PEMBRO = safe_wilcox(COMBO, PEMBRO),
    .groups = "drop"
  )

# Print the results
print(stats_conc)

# OPTIONAL: save statistics to CSV (uncomment if desired)
# write.csv(stats_conc, "statistics_concentrations_wilcoxon.csv", row.names = FALSE)
###############################################################
# END STATISTICS SECTION
###############################################################

###############################################################
# END OF SCRIPT
###############################################################
