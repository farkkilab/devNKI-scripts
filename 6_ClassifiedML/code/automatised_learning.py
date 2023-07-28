# Handling dataframe
import numpy as np
import pandas as pd

# Pre-process data
import scipy.stats as stats
from sklearn.preprocessing import MinMaxScaler, StandardScaler

# Preprocess features
from sklearn.preprocessing import LabelEncoder

# Load models
from sklearn.ensemble import RandomForestClassifier

# Model tuning 
from skopt.space import Integer
from skopt.utils import use_named_args
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_score
from skopt import gp_minimize
from time import time
from skopt.plots import plot_convergence, plot_evaluations, plot_objective

# Model testing
from sklearn.metrics import f1_score, balanced_accuracy_score

# Feature selection 
from sklearn.feature_selection import RFECV

import os
import matplotlib.pyplot as plt
import scikitplot as skplt

# multiprocessing
import multiprocessing as mp
from multiprocessing import freeze_support
import concurrent.futures

class Model:
    def __init__(self, exp_name = None, channels = None, channels_to_transform = None, channels_to_outliers = None, channels_to_scale = None, cells = None, target = None, classes = None,  therapies = None, scaler = None, df = None):

        self.rng = np.random.RandomState(0)

        # Variables for data pre-processing
        self.df = df

        # from experiments.csv:
        self.exp_name = exp_name # Experiment column
        
        # Set the output directory and filename for the figures
        self.output_dir = os.path.dirname(os.path.abspath(__file__))

        # Collect required variables from csv file
        self.features = channels.split(",") # channels column
        self.types_of_cells = cells # types_of_cells column
        self.target = target # classes_types column
        self.classes = classes.split(",") # classes_column column
        self.therapies = therapies # therapies column 

        self.features_to_transform = channels_to_transform.split(",")
        self.features_for_ouliers = channels_to_outliers.split(",")
        self.features_to_scale = channels_to_scale.split(",")

        # Define transformation type. Available types: LOG, LOG2, BOXCOX
        self.transformation = 'LOG2' 
        
        # Define how to remove outliers. Available options: trim_by_slide & remove
        self.operation = 'remove'
        
        # Scale data
        self.scaler_dict = {'MinMax': MinMaxScaler(),'Standard': StandardScaler()}
        # MinMaxScaler scales the data to be within a specific range, usually between 0 and 1
        # StandardScaler scales the data to have a mean of zero and a standard deviation of one
        self.scaler_type = scaler # Scaling type
        self.scale_by = 'slide' # Define how to scale. Available options: whole, patient, slide

        # Variables for machine learning step
        self.categorical_variables = ['Molecular.profile2', 'therapy_sequence']
        self.best_params = None
        self.RF = RandomForestClassifier(random_state = self.rng, class_weight = 'balanced')

        # Define the evaluation procedure
        self.cv = RepeatedStratifiedKFold(n_splits=3, n_repeats=3, random_state=self.rng)
        
        # Save results 
        self.train_results = None
        self.train_time = None
        self.predictions_train = None
        
        self.test_results = {}
        self.predictions_test = None

        self.feature_selection_results = {}

    def choose_features(self):
        self.df = self.df.loc[:,self.features]
    
    def transform_df(self):

        if self.transformation == 'LOG':
            self.df.loc[:, self.features_to_transform] = np.log(self.df.loc[:, self.features_to_transform] + 1)

        elif self.transformation == 'LOG2':
            self.df.loc[:, self.features_to_transform] = np.log2(self.df.loc[:, self.features_to_transform] + 1)
            
        elif self.transformation == 'BOXCOX':
            # transform data & save lambda value
            for feature in self.features_to_transform:
                self.df[feature], _ = stats.boxcox(self.df[feature].values)

        else:
            print("The transformation was incorrectly written")
            
    def remove_outliers(self):
        
        if self.operation == 'trim_by_slide':
            slides = self.df['cycif.slide'].unique()
        
            for slide in slides:
                # Create a mask to filter data belonging to the current slide
                slide_mask = self.df['cycif.slide'] == slide
                
                for feature in self.features_for_ouliers:
                    percentiles = np.percentile(df.loc[slide_mask, feature], [1, 99])
                    
                    # Replace values below the 1st percentile with the 1st percentile value
                    self.df.loc[slide_mask & (self.df[feature] < percentiles[0]), feature] = percentiles[0]
                    # Replace values above the 99th percentile with the 99th percentile value
                    self.df.loc[slide_mask & (self.df[feature] > percentiles[1]), feature] = percentiles[1]

        elif self.operation == 'remove':

            df_sub = self.df.loc[:, self.features_for_ouliers]

            # Identify outliers using the 1st (0.01) and 99th (0.99) percentiles for each feature
            # For each data point, 'lim' will be True if the value is within the range [0.01, 0.99], otherwise False
            lim = np.logical_and(df_sub < df_sub.quantile(0.99, numeric_only=False),
                             df_sub > df_sub.quantile(0.01, numeric_only=False))

            # Data points outside the range [0.01, 0.99] will be replaced with NaN
            self.df.loc[:, self.features_for_ouliers] = df_sub.where(lim, np.nan)
            
            # Drop rows with NaN in numerical columns
            self.df.dropna(subset=self.features_for_ouliers, inplace=True)
        
        else: 
            print("The command for outliers removal was incorrectly written")
    
    def scaler(self, df):
    
        # Make sure that features are numerical!
        
        # Get a scaler from the dictionary of supported scaler types
        scaler = self.scaler_dict.get(self.scaler_type)

        if scaler is None:
            raise ValueError(f"Invalid scaler type: {self.scaler_type}")
        
        # Scale the data according to the specified scaler type
        data = scaler.fit_transform(df)
        
        return data
    
    def scaling(self):
 
        if self.scale_by not in ['whole','patient', 'slide']:
            raise ValueError(f"Invalid order: {self.scale_by}")
        
        if self.scale_by == 'patient':
            patients = self.df['patient'].unique()
            # Iterate through each unique patient ID and scale the specified features for each patient separately
            for patient in patients:
                self.df.loc[self.df['patient'] == patient, self.features_to_scale] = self.scaler(self.df.loc[self.df['patient'] == patient, self.features_to_scale])
        
        elif self.scale_by == 'slide':
            slides = self.df['cycif.slide'].unique()
            # Iterate through each unique slide and scale the specified features for each slide separately
            for slide in slides:
                self.df.loc[self.df['cycif.slide'] == slide, self.features_to_scale] = self.scaler(self.df.loc[self.df['cycif.slide'] == slide, self.features_to_scale])
        
        elif self.scale_by == 'whole':
            # Scale the specified features on the entire DataFrame
            self.df.loc[self.features_to_scale] = self.scaler(self.df.loc[self.features_to_scale])
        
        else:
            print('Incorrectly written scale_by')

    def choose_types_of_cells(self):

        # Remove Others and replace all immune cells as Immune
        self.df = self.df[~self.df["GlobalCellType"].isin(['Others'])]
        self.df["GlobalCellType"] = self.df["GlobalCellType"].replace(to_replace = ['CD8.T.cells', 'CD11c.MY', 'Other.MY', 'CD163.MP', 'B.cells', 'T.regs', 'Other.immune', 'CD4.T.cells', 'CD207.MY', 'CD68.MP', 'CD15.MY'], value = 'Immune')
        
        cell_types = self.df["GlobalCellType"].unique()

        unique_patients = self.df['patient'].unique()
        print("Before thresholding number of patients: ", len(unique_patients))
        
        # Remove patients with too few cells for each chosen cell type
        for cell_type in cell_types:
            # Group df by patient and get the size of each group
            grouped_df = self.df[self.df['GlobalCellType'] == cell_type].groupby('patient').size()

            # Filter patients with 100 or fewer cells and store them
            selected_patients = grouped_df[grouped_df <= 100].index

            print(f"Cell type {cell_type}: {len(selected_patients)} patients with <= 100 cells")

            # Remove patients with 100 and fewer cells
            if len(selected_patients) > 0:
                self.df = self.df[~self.df['patient'].isin(selected_patients)]

        unique_patients = self.df['patient'].unique()
        print("After thresholding number of patients: ", len(unique_patients))

        # Leave only chosen cell types in df
        if self.types_of_cells == "Non-cancer":
            self.df = self.df[~self.df["GlobalCellType"].isin(['Others', 'Cancer'])]

        elif self.types_of_cells in ["Cancer", "Stromal", "Immune"]:
            self.df = self.df[self.df["GlobalCellType"] == self.types_of_cells]

        elif self.types_of_cells not in ["All"]:
            print("Incorrect cell type")
        
        print("Chosen cell types: ", self.df['GlobalCellType'].unique())

        # Remove column GlobalCellType
        self.df = self.df.drop(columns='GlobalCellType')
        
    def choose_classes(self):
        
        self.df = self.df[self.df[self.target].isin(self.classes)]

        print(f"Chosen classes of {self.target} are {self.df[self.target].unique()}")

    def choose_therapy_sequences(self):
        if self.therapies != 'All':
            self.df = self.df[self.df['therapy_sequence'] == self.therapies]
        
        print(f"Chosen therapies are {self.df['therapy_sequence'].unique()}")

    def prepare_categorical_inputs(self):

        # Encode categorical variables 
        self.label_encoder_dict = {}
        for variable in self.categorical_variables:
            
            label_encoder = LabelEncoder()
            
            # Convert the values of the current categorical variable to strings and encode them
            # The encoded values will replace the original values in the DataFrame
            self.df.loc[:, variable] = label_encoder.fit_transform(self.df.loc[:, variable].astype('str').values)
            
            # Save the encoding results for the current variable in 'label_encoder_dict'
            # The dictionary will store the unique class labels and their corresponding encoded values
            self.label_encoder_dict[variable] = {
                'label_codes': label_encoder.classes_.tolist(), 
                'label_values': label_encoder.transform(label_encoder.classes_).tolist()
            }

    def prep_train_test_data(self):
    
        # To receive same number of patients independent of runs
        np.random.seed(33)
        
        # Group df by Molecular profule and therapy
        grouped = self.df.groupby(['Molecular.profile2', 'therapy_sequence'])

        # Assign 80% of patients from each unique group of Molecular profile and therapy to training and rest to test sets
        self.X_train_full = grouped.apply(lambda x: x.loc[x['patient'].isin(np.random.choice(x['patient'].unique(), size=int(0.8*len(x['patient'].unique())), replace=False))])
        self.X_test_full = grouped.apply(lambda x: x.loc[x['patient'].isin(np.setdiff1d(x['patient'].unique(), self.X_train_full['patient']))])
        
        # Assign target columns to y_train and y_test
        self.y_train=self.X_train_full[self.target]
        self.y_test=self.X_test_full[self.target]
        
        # Drop target columns
        self.X_train = self.X_train_full.drop(columns = [self.target, 'patient', 'cycif.slide'])
        self.X_test = self.X_test_full.drop(columns = [self.target, 'patient', 'cycif.slide'])

    def find_RF_by_Baseian(self, X_train, y_train):
        '''
        inputs:
           - X_train: features training set
           - y_train: income training set
        '''

        # Define hyperparameters search space 
        space = [Integer(200, 900, name='min_samples_split'),
                Integer(200, 900, name='min_samples_leaf'),
                Integer(50, 200, name='n_estimators')]
        
        # Define evaluation for each call within hyperparameters' space
        @use_named_args(space)
        def evaluate_model(**params):
            scores = cross_val_score(self.RF, X_train, y_train, scoring="f1_weighted", cv=self.cv, n_jobs = -1)
            print("Standard deviation of cross-validated scores:", np.std(scores))
            return -np.mean(scores)
        
        # Perform Bayesian search on defined hyperparameters search space 
        start = time() # Get start time
        n_calls = 30
        self.train_results = gp_minimize(evaluate_model, space, random_state=self.rng, n_random_starts=5, n_calls = n_calls)
        end = time() # Get end time
        
        # Calculate the training time
        self.train_time = end - start

        print(f"\nBest score: {-self.train_results.fun:.5f} with hyperparameters: {self.train_results.x}")
        print("Training_time: ", self.train_time)
        
        # Store the best hyperparameters
        self.best_params = {'min_samples_split': self.train_results.x[0],
                    'min_samples_leaf': self.train_results.x[1],
                    'n_estimators': self.train_results.x[2]}
        
        # Create RF model with stored best hyperparameters
        self.RF = RandomForestClassifier(**self.best_params, random_state = self.rng, class_weight = 'balanced')

        # Save figures of training 
        name_convergence = self.exp_name + "_convergence_plot.png"
        plot_convergence(self.train_results)
        plt.savefig(os.path.join(self.output_dir+"/"+name_convergence))

        name_evaluations = self.exp_name + "_evaluations_plot.png"
        plot_evaluations(self.train_results)
        plt.savefig(os.path.join(self.output_dir+"/"+name_evaluations))

        name_objective = self.exp_name + "_objective_plot.png"
        plot_objective(self.train_results)
        plt.savefig(os.path.join(self.output_dir+"/"+name_objective))       
        
    def test(self, X_train, y_train, X_test_full, X_test, y_test):
        '''
        inputs:
           - X_train: features training set
           - y_train: income training set
           - X_test_full: features testing set with saved patient and slide
           - X_test: features testing set
           - y_test: income testing set
        '''
        
        # Training step: Fit the learner to the training data
        self.RF.fit(X_train, y_train)
            
        # Get the predictions on the training set(X_train),
        # then get predictions on the test samples(X_test)
        self.predictions_train = self.RF.predict(X_train)
        # Test step:
        self.predictions_test = self.RF.predict(X_test)
        
        # Save the predictions of X_test as df:
        # 1. Convert predictions_test to a pandas series
        predictions_series = pd.Series(self.predictions_test, name='predicted_label')
        # 2. Convert MultiIndex to a flat index in X_test
        X_test_flat = X_test_full.reset_index(drop=True)
        # 3. Concatenate X_test_flat and predictions_series along axis=1 to receive df of test set predictions 
        self.predictions_df = pd.concat([X_test_flat, predictions_series], axis=1)
                
        # Compute balanced accuracy on training samples
        self.test_results['balanced_acc_train'] = balanced_accuracy_score(y_train, self.predictions_train)
            
        # Compute balanced accuracy on the test set
        self.test_results['balanced_acc_test'] = balanced_accuracy_score(y_test, self.predictions_test)
        
        # Compute weighted F1-score on training samples
        self.test_results['f1_train'] = f1_score(y_train, self.predictions_train, average='weighted')
            
        # Compute weighted F1-score on the test set
        self.test_results['f1_test'] = f1_score(y_test, self.predictions_test, average='weighted')
        
        # Compute AUC on the test set if needed
        probas_test = self.RF.predict_proba(X_test) # 2D array with number columns as target classes, where each column contains the predicted probabilities of belonging to the corresponding target class.
        # self.test_results['auc_test'] = roc_auc_score(y_test, probas_test[:,1])
        
        # Plots
        classes_on_plot = self.label_encoder_dict['Molecular.profile2']['label_codes']

        # Precision-recall
        skplt.metrics.plot_precision_recall(y_test, probas_test, classes_on_plot)
        name_prec_recall = self.exp_name + '_precision_recall_curve.png'
        plt.savefig(os.path.join(self.output_dir+"/"+name_prec_recall), dpi=300, bbox_inches='tight')
        
        # ROC curves
        name_roc = self.exp_name + '_roc_curve.png'
        skplt.metrics.plot_roc(y_test, probas_test, classes_on_plot)
        plt.savefig(os.path.join(self.output_dir+"/"+name_roc), dpi=300, bbox_inches='tight')
    
    def calculate_mean_per_patient(self):
        name = self.exp_name + '_patients_predictions.csv'
        patients_predictions = self.predictions_df.groupby('patient')['predicted_label'].mean().reset_index()
        patients_predictions.to_csv(name, index=False)

    def eliminate_features(self, X_train, y_train):
        '''
        inputs:
        - X_train: features training set
        - y_train: income training set
        '''
        
        # Define cross-validated recursive feature elimination 
        # RFECV automatically selects the best number of features based on the provided scoring function of accuracy
        rfecv = RFECV(
            estimator=self.RF,
            min_features_to_select=5,
            step=5,
            n_jobs=-1,
            scoring="accuracy",
            cv=self.cv,
        )
        
        # Run cross-validated recursive feature elimination to select the most predictive features
        rfecv.fit(X_train, y_train)
        
        # Calculate feature importances using mean decrease in impurity (Gini impurity)
        
        # "importances" here will represent the feature importances based on Gini impurity
        importances = rfecv.estimator_.feature_importances_
        
        # Get the support mask indicating which features were selected by RFECV
        support = rfecv.support_

        # Get the list of most predictive and eliminated features
        features = X_train.columns
        eliminated_features = list(features[~support])
        important_features = features[support]

        print("Important features: ", important_features)
        print("Eliminated features: ", eliminated_features)
        
        # Plot
        # Create dictionary with important feature names as keys and importances as values
        feature_importances = dict(zip(important_features, importances))
        importances = dict(sorted(feature_importances.items(), key=lambda x: x[1], reverse=True))

        # Creat the plot for features
        plt.figure()
        plt.title("Weights for Most Predictive Features")
        plt.bar(np.arange(len(important_features)), importances.values(), width = 0.6, align="center", \
            label = "Feature Weight")
        plt.xticks(np.arange(len(important_features)), importances.keys(), rotation ='vertical')
        plt.ylabel("Weight", fontsize = 12)
        plt.xlabel("Feature", fontsize = 12)

        # Save figure
        name_pred_features = self.exp_name + '_predictive_features.png'
        plt.savefig(os.path.join(self.output_dir+"/"+name_pred_features), dpi=300, bbox_inches='tight')

        self.feature_selection_results['importances'] = importances
        self.feature_selection_results['eliminated_features'] = eliminated_features

def prep_df():

    all_cells = pd.read_csv('All_cells_scimap-format.csv')
    all_cells[['cycif.slide','cycif.core.id']] = all_cells['imageid'].apply(lambda x: pd.Series(str(x).split("_core")))

    # Load spatial proportions
    spatial_proportions = pd.read_csv("Spatial-cell-counts_radious46px.csv")

    # Merge cancer cells and their corresponding spatial proportions on CellId
    dataset = pd.merge(all_cells, spatial_proportions, on="CellId")

    # Load the text file into a pandas DataFrame
    slides_explanations = pd.read_csv('All-slides_cycif2samples.txt', sep="	")
    slides_explanations['cycif.core.id']=slides_explanations['cycif.core.id'].str.replace("core", "").astype(dataset['cycif.core.id'].dtype)
    
    # Load patinets' clinical info file into a pandas DataFrame
    patients_data = pd.read_csv('T-CycIF.csv')
    patients_data.rename(columns = {'TnumberAVL':'patient'}, inplace = True)

    # Load molecular profile 2
    moleculer_profiles = pd.read_csv("Molecular_profiles_patients_2023.csv")

    # Merge the two data frames on the slide number and core id columns
    combined_dataset = pd.merge(dataset, slides_explanations[["cycif.slide", "cycif.core.id","patient"]], on=["cycif.slide", "cycif.core.id"])

    final_patients_data = pd.merge(moleculer_profiles, patients_data[["patient","therapy_sequence"]], on=["patient"])

    # Merge the two data frames on patient column
    df = pd.merge(combined_dataset, final_patients_data, on=["patient"])

    # Group all NACT therapy sequences
    df['therapy_sequence'] = df['therapy_sequence'].replace(to_replace='NACT followed by re-debulking', value='NACT', regex=True)
    df['therapy_sequence'] = df['therapy_sequence'].replace(to_replace='PDS followed by NACT', value='NACT', regex=True)
    df['therapy_sequence'] = df['therapy_sequence'].replace(['Primairy debulking', 'Only debulking'], 'PDS')

    return df

def perform_ml_run(args):

    exp_name, channels, channels_to_transform, channels_to_outliers, channels_to_scale, cells, target, classes, therapies, scaler, best_params, acc_train, acc_test, f_train, f_test, pred_features, eliminated_features, name, exps, df = args

    model = Model(exp_name, channels, channels_to_transform, channels_to_outliers, channels_to_scale, cells, target, classes, therapies, scaler, df)

    model.choose_features()
    model.transform_df()
    model.remove_outliers()
    model.scaling()
    model.choose_types_of_cells()
    model.choose_classes()
    model.choose_therapy_sequences()

    model.prepare_categorical_inputs()
    model.prep_train_test_data()
    print(model.X_train.columns)
    model.find_RF_by_Baseian(model.X_train,model.y_train) # .iloc[:1000]
    model.test(model.X_train,model.y_train, model.X_test_full, model.X_test,model.y_test)
    model.calculate_mean_per_patient()
    model.eliminate_features(model.X_train,model.y_train)

    # Collect the results in a dictionary
    results = {
        'exp_name': exp_name,
        'best_parameters': model.best_params,
        'balanced_acc_train': model.test_results['balanced_acc_train'],
        'balanced_acc_test': model.test_results['balanced_acc_test'],
        'f1_train': model.test_results['f1_train'],
        'f1_test': model.test_results['f1_test'],
        'most_predictive_features': model.feature_selection_results['importances'],
        'eliminated_features': model.feature_selection_results['eliminated_features']
    }

    print(exp_name, "is done") 

    return results

if __name__ == '__main__':

    freeze_support()

    names = ['EXPERIMENTS_FILE_NAME.csv', ...] 

    df = pd.read_csv("DATASET_FILENAME.csv")
    df['Molecular.profile2'] = df['Molecular.profile2'].replace('BRCAmut/met', 'BRCAmutmet')
    
    # Get the number of available CPU cores for parallel processing
    max_workers= mp.cpu_count()

    for name in names:

        exps = pd.read_csv(name)

        print(exps.head(1))

        # Use ProcessPoolExecutor for parallel execution of multiple experiments
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:

            print("Start of experiment")
            
            # Perform multiple machine learning runs in parallel using executor.map()
            # Each experiment's details are passed to the perform_ml_run function
            results = list(executor.map(perform_ml_run, [(d['experiment'], d['channels'], d['channels_to_transform'], d['channels_to_outliers'], d['channels_to_scale'], d['types_of_cells'], d['classes_column'], d['classes_types'], d['therapies'], d['scaling_type'], d['best_parameters'], d['balanced_acc_train'], d['balanced_acc_test'], d['f1_train'], d['f1_test'], d['most_predictive_features'], d['eliminated_features'], name, exps, df) for i, d in exps.iterrows()]))

            executor.shutdown()

        # Update the experiment DataFrame with the results obtained from the machine learning runs
        for result in results:
            exp_name = result['exp_name']
            exps.loc[exps['experiment'] == exp_name, 'best_parameters'] = str(result['best_parameters'])
            exps.loc[exps['experiment'] == exp_name, 'balanced_acc_train'] = result['balanced_acc_train']
            exps.loc[exps['experiment'] == exp_name, 'balanced_acc_test'] = result['balanced_acc_test']
            exps.loc[exps['experiment'] == exp_name, 'f1_train'] = result['f1_train']
            exps.loc[exps['experiment'] == exp_name, 'f1_test'] = result['f1_test']
            exps.loc[exps['experiment'] == exp_name, 'most_predictive_features'] = str(result['most_predictive_features'])
            exps.loc[exps['experiment'] == exp_name, 'eliminated_features'] = str(result['eliminated_features'])

        # Save the updated DataFrame to experiments file
        exps.to_csv(name, index=False)

        print(exps.head(5))
