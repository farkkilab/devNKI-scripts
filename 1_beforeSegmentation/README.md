# devNKI-scripts:1\_beforeSegmentation

## Description

Scripts used to get ready the pixel probability maps used for segmentations.

## File list

- Crop\_core.m: Script used to crop slides and extract channels of interest in medium resolution (uint16) for each core. The resultant files will be used as input for the `headless_ilastik_running.sh`.


- headless\_ilastik\_running.sh: To run ilastik (pixel probabilities maps). This script it is necessary to run it inside Linux, I used the ubuntu-terminal inside windows. As output will generated and image file (.tiff) with the probabilities of belonging to a class for each pixel in the input images.
