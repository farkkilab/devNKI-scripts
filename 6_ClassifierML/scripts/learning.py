import numpy as np
from numpy import random
import pandas as pd

from sklearn.preprocessing import LabelEncoder

from skopt import gp_minimize
from skopt.space import Real, Integer
from skopt.utils import use_named_args

from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import AdaBoostClassifier, RandomForestClassifier
from sklearn.model_selection import RepeatedStratifiedKFold, GridSearchCV, cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import RFE, RFECV

from sklearn.metrics import make_scorer, precision_score, fbeta_score, accuracy_score, roc_auc_score
from time import time

import visuals as vs
import matplotlib.pyplot as plt

# function to perform stratified random sampling
def prep_train_test_data(df, target, col1, col2, sample_size):
    
    # To receive same number of patients independent of runs
    np.random.seed(0)
    
    # Sample_size - the number of patients to be selected from each sub-group of the data defined by a combination of col1 and col2
    
    # Create an empty dataframe for storing the down-sampled data
    X_train = pd.DataFrame()
    
    # Get the unique values of col1 and col2 in the input dataframe
    all_col1 = df[col1].unique()
    all_col2 = df[col2].unique()
    
    # Loop through the unique values of col1 and col2 to select a sub-group of data for each combination of col1 and col2
    for c1 in all_col1:
        for c2 in all_col2:
            sub_df = df[(df[col1] == c1) & (df[col2] == c2)]
            
            # If the size of the sub-group is larger than or equal to the specified sample size, down-sample the sub-group to the specified sample size
            if len(sub_df.groupby("patient")) >= sample_size:
                X_train = pd.concat([X_train, sub_df.sample(sample_size, replace=True)])
            # If the size of the sub-group is smaller than the specified sample size, keep all the data in the sub-group
            else:
                X_train = pd.concat([X_train, sub_df])
    
    # merge the two data frames on all columns and create a new column 'isnull' to identify the null values
    merged_df = pd.merge(df, X_train, how='outer', indicator=True)

    # select only the rows where the down-sampled version has null values
    X_test = merged_df.loc[merged_df['_merge'] == 'left_only', df.columns]
    
    # Drop if needed the columns
    X_train = X_train.drop(columns = ['patient'])
    X_test = X_test.drop(columns = ['patient'])
    
    # Assign target columns to y_train and y_test
    y_train=X_train[target]
    y_test=X_test[target]
    
    # Drop target columns
    X_train = X_train.drop(columns = [target])
    X_test = X_test.drop(columns = [target])
     

#     # Print the size and number of unique values of each column in the down-sampled dataframe
#     print("The size of down-sampled dataframe: ", X_train.shape)
#     print("Patient number: ", len(X_train["patient"].unique()))
#     print("Molecular profile number: ", len(X_train["Molecular_profile"].unique()))
#     print("Therapy number: ", X_train["therapy_sequence"].unique())
    
    return X_train, y_train, X_test, y_test

def prepare_categorical_inputs(df, variables):
    label_encoder = LabelEncoder()
    for variable in variables:
        df.loc[:, variable] = label_encoder.fit_transform((df.loc[:,variable].astype('str')).values)
    return df

def prepare_targets(y_train, y_test):
    le = LabelEncoder()
    le.fit(y_train)
    y_train_enc = le.transform(y_train)
    y_test_enc = le.transform(y_test)
    return y_train_enc, y_test_enc

def find_DT(X_train, y_train):
   
    # Define the min split based on features in dataset
    x = 2**(X_train.shape[1]-1)
    
    # create a pipeline with a decision tree classifier
    pipeline = Pipeline([
        ('tree', DecisionTreeClassifier())
    ])

    # create a dictionary of hyperparameters to search
    param_grid = {
        'tree__max_depth': [random.randint(5,15), random.randint(5,15), random.randint(5,15)],
        'tree__min_samples_split': [random.randint(10,1000), random.randint(1000,10000), random.randint(10000,20000)],
        'tree__min_samples_leaf': [random.randint(10,1000), random.randint(1000,5000), random.randint(5000,10000)],
    }

    # k is the number of features in the dataset -> min_samples_split=2**(k-1)

    # perform grid search to find best model
    grid = GridSearchCV(pipeline, param_grid, cv=5, n_jobs = -1, return_train_score=True)
    grid.fit(X_train, y_train)

    # create a dataframe of the grid search results
    df_results = pd.DataFrame(grid.cv_results_)[['params', 'rank_test_score', 'mean_test_score', 'std_test_score','mean_train_score', 'std_train_score']]

    # sort the results by rank_test_score
    df_results = df_results.sort_values('rank_test_score')

    df_results.to_csv('results_DT.csv', index=False)


def find_DT_by_Baseian(X_train, y_train):
    
    space = [Integer(10, 15, name='max_depth'),
             Integer(10**2, 10**3, name='min_samples_split'),
             Integer(10**2, 10**3, name='min_samples_leaf')]
#              Real(0, 1, 'uniform', name='max_features')]
    
    # define the evaluation procedure
    cv = RepeatedStratifiedKFold(n_splits=3, n_repeats=3, random_state=1)

    @use_named_args(space)
    def evaluate_model(**params):
        model = DecisionTreeClassifier(**params)
        cv = RepeatedStratifiedKFold(n_splits=3, n_repeats=3, random_state=1)
        scoring = make_scorer(precision_score, average='weighted')
        scores = cross_val_score(model, X_train, y_train, cv=cv, n_jobs=-1, scoring='accuracy')
        print(f"params: {params}, scores: {scores}")
        return -np.mean(scores)
    
    result = gp_minimize(evaluate_model, space, random_state=0, n_random_starts=5)

    print(f"Best score: {-result.fun:.5f} with hyperparameters: {result.x}")
    
    # train the final model using the best hyperparameters
    best_params = {'max_depth': result.x[0],
                   'min_samples_split': result.x[1],
                   'min_samples_leaf': result.x[2]}
    model = DecisionTreeClassifier(**best_params)
    
    return model

    
def find_RF(X_train, y_train):
       
    # Define the min split based on features in dataset
    x = 2**(X_train.shape[1]-1)
    
    # create a pipeline with a random forest classifier
    pipeline = Pipeline([
        ('forest', RandomForestClassifier())
    ])

    # create a dictionary of hyperparameters to search
    param_grid = {
        'forest__max_depth': [5, 10, 15],
        'forest__min_samples_split': [x, x+10, x+50],
        'forest__min_samples_leaf': [2, 4, 8],
    }

    # perform grid search to find best model
    cv = RepeatedStratifiedKFold(n_splits=3, n_repeats=3, random_state=1)
    grid = GridSearchCV(pipeline, param_grid, cv=cv, n_jobs=-1, return_train_score=True)
    grid.fit(X_train, y_train)
    
    # create a dataframe of the grid search results
    df_results = pd.DataFrame(grid.cv_results_)[['params', 'rank_test_score', 'mean_test_score', 'std_test_score', 'mean_train_score', 'std_train_score']]

    # sort the results by rank_test_score
    df_results = df_results.sort_values('rank_test_score')

    df_results.to_csv('results_RF.csv', index=False)
    
def eliminate_features(learner, X_train, y_train, X_test, y_test):
    '''
    inputs:
       - learner: the learning algorithm to be trained and predicted on
       - X_train: features training set
       - y_train: income training set
       - X_test: features testing set
       - y_test: income testing set
    '''
    
    cv = RepeatedStratifiedKFold(n_splits=3, n_repeats=3, random_state=1)
    
    # Init, fit
    rfecv = RFECV(
        estimator=learner,
        min_features_to_select=5,
        step=5,
        n_jobs=-1,
        scoring="accuracy",
        cv=cv,
    )

    results = {}
    
    rfecv.fit(X_train, y_train)
    
    importances = rfecv.estimator_.feature_importances_
    support = rfecv.support_
    features = X_train.columns
    eliminated_features = features[~support]
    important_features = features[support]
    
    # Create dictionary with important feature names as keys and importances as values
    feature_importances = dict(zip(important_features, importances))
    results['importances'] = feature_importances
    
    results['eliminated_features'] = eliminated_features
    
    return results, rfecv.transform(X_train), rfecv.transform(X_test)

def test(learner, X_train, y_train, X_test, y_test):
    '''
    inputs:
       - learner: the learning algorithm to be trained and predicted on
       - X_train: features training set
       - y_train: income training set
       - X_test: features testing set
       - y_test: income testing set
    '''

    results = {}
    
    # Fit the learner to the training data using slicing with 'sample_size'
    learner.fit(X_train, y_train)
        
    # Get the predictions on the test set(X_test),
    #then get predictions on the training samples(X_train) using .predict()
    predictions_test = learner.predict(X_test)
    predictions_train = learner.predict(X_train)
            
    # TODO: Compute accuracy on training samples which is y_train[:300]
    results['acc_train'] = accuracy_score(y_train, predictions_train)
        
    # TODO: Compute accuracy on test set using accuracy_score()
    results['acc_test'] = accuracy_score(y_test, predictions_test)
    
    # TODO: Compute F-score on training samples using fbeta_score()
    results['f_train'] = fbeta_score(y_train, predictions_train, beta = 0.5, average='macro')
        
    # TODO: Compute F-score on the test set which is y_test
    results['f_test'] = fbeta_score(y_test, predictions_test, beta = 0.5, average='macro')
    
    # Compute AUC on test set
    probas_test = learner.predict_proba(X_test) # 2D array with number columns as target classes, where each column contains the predicted probabilities of belonging to the corresponding target class.
    results['auc_test'] = roc_auc_score(y_test, probas_test, multi_class='ovr')
    
    cv = RepeatedStratifiedKFold(n_splits=3, n_repeats=3, random_state=1)

    # Compute cross-validation scores on the full training set
    cv_scores = cross_val_score(learner, X_train, y_train, cv=cv, n_jobs=-1)
    results['cv_scores'] = cv_scores
    
    return results 
    
    
def train_predict(learner, sample_size, X_train, y_train, X_test, y_test):
    '''
    inputs:
       - learner: the learning algorithm to be trained and predicted on
       - X_train: features training set
       - y_train: income training set
       - X_test: features testing set
       - y_test: income testing set
    '''
    # create pipeline
    rfe = RFE(estimator=DecisionTreeClassifier(), n_features_to_select=4)
    pipeline = Pipeline(steps=[('s',rfe),('m',learner)])
    
    # define the evaluation procedure
    cv = RepeatedStratifiedKFold(n_splits=3, n_repeats=3, random_state=1)

    results = {}
    
    # Fit the learner to the training data using slicing with 'sample_size'
    start = time() # Get start time
    pipeline.fit(X_train[:sample_size], y_train[:sample_size])
    end = time() # Get end time
    
    # TODO: Calculate the training time
    results['train_time'] = end - start
        
    # TODO: Get the predictions on the test set(X_test),
    #       then get predictions on the first 300 training samples(X_train) using .predict()
    start = time() # Get start time
    predictions_test = pipeline.predict(X_test)
    predictions_train = pipeline.predict(X_train)
    end = time() # Get end time
    
    # TODO: Calculate the total prediction time
    results['pred_time'] = end - start
            
    # TODO: Compute accuracy on training samples which is y_train[:300]
    results['acc_train'] = accuracy_score(y_train, predictions_train)
        
    # TODO: Compute accuracy on test set using accuracy_score()
    results['acc_test'] = accuracy_score(y_test, predictions_test)
    
    # TODO: Compute F-score on training samples using fbeta_score()
    results['f_train'] = fbeta_score(y_train, predictions_train, beta = 0.5, average='macro')
        
    # TODO: Compute F-score on the test set which is y_test
    results['f_test'] = fbeta_score(y_test, predictions_test, beta = 0.5, average='macro')
    
#     # Compute AUC on test set
#     probas_test = pipeline.predict_proba(X_test) # 2D array with number columns as target classes, where each column contains the predicted probabilities of belonging to the corresponding target class.
#     results['auc_test'] = roc_auc_score(y_test, probas_test, multi_class='ovr')

    # Compute cross-validation scores on the full training set
    cv_scores = cross_val_score(pipeline, X_train, y_train, cv=cv, n_jobs=-1)
    results['cv_scores'] = cv_scores
    
#     # Get feature importances
#     if hasattr(pipeline, 'named_steps'):
#         importances = pipeline.named_steps['m'].feature_importances_
#     else:
#         importances = pipeline.steps[-1][1].feature_importances_
#     results['importances'] = importances
    
#     # Get eliminated features
#     support = pipeline.named_steps['s'].support_
#     eliminated_features = [X_train.columns[i] for i in range(len(support)) if not support[i]]
#     print(eliminated_features)

# Get feature importances and eliminate features that were not selected
    if hasattr(pipeline, 'named_steps'):
        importances = pipeline.named_steps['m'].feature_importances_
        support = pipeline.named_steps['s'].support_
        features = X_train.columns
        eliminated_features = features[~support]
        important_features = features[support]
    else:
        importances = pipeline.steps[-1][1].feature_importances_
        support = pipeline.steps[-2][1].support_
        features = X_train.columns
        eliminated_features = features[~support]
        important_features = features[support]
    
    # Create dictionary with important feature names as keys and importances as values
    feature_importances = dict(zip(important_features, importances))
    results['importances'] = feature_importances
    
    results['eliminated_features'] = eliminated_features
    
    # Success
    print("{} trained on samples.".format(learner.__class__.__name__))
#     print("Params: ", pipeline.get_params())
    
    # Return the results
    return pipeline, results

def plot_results(X_train, y_train, results):
    fig, axs = plt.subplots(len(results), 4, figsize=(30, 5*len(results)), squeeze=False)

    for i, (clf_name, clf_data) in enumerate(results.items()):
        x = np.array(list(clf_data.keys()))
        train_times = np.array([d['train_time'] for d in clf_data.values()])
        pred_times = np.array([d['pred_time'] for d in clf_data.values()])
        acc_train = np.array([d['acc_train'] for d in clf_data.values()])
        acc_test = np.array([d['acc_test'] for d in clf_data.values()])
        f_train = np.array([d['f_train'] for d in clf_data.values()])
        f_test = np.array([d['f_test'] for d in clf_data.values()])
#         auc_test = np.array([d['auc_test'] for d in clf_data.values()])

        axs[i, 0].bar(x, train_times)
        axs[i, 0].set_title(f"{clf_name} - Training Time")
        axs[i, 0].set_xlabel('Iteration')
        axs[i, 0].set_ylabel('Time (s)')

        axs[i, 1].plot(x, acc_train, 'o', label='Training')
        axs[i, 1].plot(x, acc_test, 'o', label='Testing')
        axs[i, 1].set_title(f"{clf_name} - Accuracy")
        axs[i, 1].set_xlabel('Iteration')
        axs[i, 1].set_ylabel('Accuracy')
        axs[i, 1].legend()

        axs[i, 2].plot(x, f_train, 'o', label='Training')
        axs[i, 2].plot(x, f_test, 'o', label='Testing')
        axs[i, 2].set_title(f"{clf_name} - F1 Score")
        axs[i, 2].set_xlabel('Iteration')
        axs[i, 2].set_ylim(0,1)
        axs[i, 2].set_ylabel('F1 Score')
        axs[i, 2].legend()
        
#         axs[i, 3].bar(x, auc_test)
#         axs[i, 3].set_title(f"{clf_name} - Testing AUC")
#         axs[i, 3].set_xlabel('Iteration')
#         axs[i, 3].set_ylim(0,1)
#         axs[i, 3].set_ylabel('Time (s)')
        

        # Display the most important features for last model
        
        importances = np.array([d['importances'] for d in clf_data.values()])[-1]
        importances = dict(sorted(importances.items(), key=lambda x: x[1], reverse=True))
#         keys = sorted(list(importances.keys), reverse = True)
#         values = sorted(list(importances.values), reverse = True)
        
        print(importances.values)

        # Creat the plot for feature
        axs[i, 3].set_title("Weights for Most Predictive Features")
        axs[i, 3].bar(np.arange(4), importances.values(), width = 0.6, align="center", \
              label = "Feature Weight")
        axs[i, 3].bar(np.arange(4) - 0.3, np.cumsum(list(importances.values())), width = 0.2, align = "center", \
              label = "Cumulative Feature Weight")
        axs[i, 3].set_xticks(np.arange(4))
        axs[i, 3].set_xticklabels(importances.keys())
        axs[i, 3].set_ylabel("Weight", fontsize = 12)
        axs[i, 3].set_xlabel("Feature", fontsize = 12)
        axs[i, 3].legend()