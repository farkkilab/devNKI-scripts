# OLD VERSIONL: RF-based pipeline for binary classification 
For new version see CEFIIRA

Files:

- Code is available in automatised_learning.py
- Conda environment is in nki_env.yml

To run the code:
1. Set up a conda environment with the help of nki_env.yml file

   ````
   conda env create -f nki_env.yml
   ````
   ````
   conda activate nki_env
   ````  
2. Prepare a CSV file with desired experiments

   | experiment  | channels | channels_to_transfrorm | channels_to_outliers | chanenels_to_scale | types_of_cells | classes_column | classes_types | therapies | scaling_type | best_parameters | balanced_acc_train | balanced_acc_test | f1_train | f1_test | most_predictive_features | eliminated_features |
   | :---:  | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |


4. Don't forget to define the datasets that you will use: experiments CSV file in line 502 and full cell dataset - in 504
   ````
   names = ['EXPERIMENTS_FILE_NAME.csv',...] 
   ````
   ````
   df = pd.read_csv("DATASET_FILENAME.csv")
   ````  
   
5. You are now ready to run the script
