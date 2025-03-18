# devNKI-scripts:1\_beforeSegmentation

## Description

Scripts used to get ready the pixel probability maps used for segmentations for cells. 
Also to merge clinical information to stablish the molecular profile by patient.

## File list

- Crop\_core.m: Script used to crop whole slides and extract individual cores per image, selecting channels of interest in medium resolution (uint16). The resultant tif files will be used as input for the `headless_ilastik_running.sh`. This script takes as input the whole slide ometif, and a list of coordinates for croping each core. That list of coordinates was generated using segmentation of cores.

- Crop\_core\_outOfbox.m: Script used to crop whole slides and extract individual cores per image, selecting channels of interest in medium resolution (uint16). The resultant tif files will be used as input for the `headless_ilastik_running.sh`. This script takes as input the whole slide ometif, and a list of coordinates for croping each core. That list of coordinates was generated manually using Fiji for cores overlaping.

- headless\_ilastik\_running.sh: To run ilastik (pixel probabilities maps). This script it is necessary be to run inside Linux, I used the ubuntu-terminal inside windows. As output will be generated an image file (.tiff) with the probabilities of belonging to a class for each pixel in for each of the input images.

- Clinical_data_intersection.Rmd: Script used to merge the clinical data, molecular profule data, ShallowHRD results, core numbers, and then stablish the final Molecular profile per sample.

## Aditional notes

For segmentation it was used CellProfiler v3.1, the corresponding workflow is the file "Dev_project_files/2_segment_ilastik.cpproj" stored in: https://datacloud.helsinki.fi/index.php/s/o7yZPTrYWB6jYNE 

The Ilastik version used was 1.3.2, the Ilastik project used is the file "Dev_project_files/Ilastik_labeling/Cellring_nuclei_background_Selected-Channels_3.ilp" stored in: 
https://datacloud.helsinki.fi/index.php/s/o7yZPTrYWB6jYNE  

R version used: 4.2.3
