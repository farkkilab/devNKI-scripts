# devNKI-scripts:2\_Quantification2QCs

## Description

After segmentation. This script quantified by core the median signal intensity per channel per core and other stats. After the quantification, big tables are produced per slide and subsequent QC steps are run.

## File list

- Cell_quantification.m: Script used to quantify signal by cells. As result produce a table per core and a big table per slide, the big table is used in subsequent analysis. 

- lost_cell_detection_cycif.py: Script used to label cells that has lost the DAPI intensity (concordance) over channels. To run this python script first was necessary to install the cycif_analysis_suite package. As input take a big table per slide with all signal intensities, as output add to the table a column with the status (True or False) if the cell was lost.

## Additional notes

- After the cell quantification some cores were removed from the final table "by hand" (linux one-liner), because those were not actual cores.
