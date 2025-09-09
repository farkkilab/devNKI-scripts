# devNKI-scripts:7\_ML-validations

## Description

Scripts used to validate some of results obtanined using the machine learning pipeline

## File list

- `iPDCs_Flow_cytometry_Fig6kl.Rmd` is an R script used to analyze the association between MHCII-high cancer samples and the immunotherapy response of CD8 T cells in iPDCs samples. The input files consist of gates selected from mixed signals (including both singlets and doublets) of cancer, immune, and myeloid populations for MHCII-high samples, as well as statistics on CD8 T cell responses for GrzB, IFNg, and Ki67 (calculated in Launonen et al., 2024). The outputs include a heatmap for Fig. 6k, a dot plot for Fig. 6l, and Sup. Fig. 6 violin plots.
- MHCII\_RNCs.Rmd: R script to evaluate the expression of MHCII in cancer cells across different Recurrent Cellullar Neighborhoods and per Molecular Profile.
