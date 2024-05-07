# RF-based pipeline for binary classification 

### About the project
This repository is composed of the code for Random Forest pipelines and their results' analysis.


Files:

- Code is available in mol_profiles_bin_class_pipeline.py
- Conda environment is in nki_env.yml

### Prerequisites
To run the code:
   Set up a conda environment with the help of nki_env.yml file

   ````
   conda env create -f nki_env.yml
   ````
   ````
   conda activate nki_env
   ````  

### Getting started

1. Prepare a CSV file with desired experiments

  | experiment  | channels | channels_to_transform | channels_to_outliers | channels_to_scale | types_of_cells | classes_column | classes_types | therapies | scaling_type | best_parameters | balanced_acc_train | balanced_acc_test | f1_train | f1_test | most_predictive_features | eliminated_features | permutation_scores | random_seed |
    | :---:  | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |


1. Don't forget to define the datasets that you will use:
   ````
   names = ['EXPERIMENTS_FILE_NAME.csv',...] 
   ````
   ````
   df = pd.read_csv("DATASET_FILENAME.csv")
   ````  
   
2. You are now ready to run the script

***

### Contact
Aleksandra Shabanova aleksandra.shabanova@helsinki.fi