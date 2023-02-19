from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.decomposition import PCA
from scipy import stats
from scipy.stats import pearsonr
import numpy as np
import pandas as pd
import plotly.express as px
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


SCALER_DICT = {
    'MinMax': MinMaxScaler(),
    'Standard': StandardScaler()
}

stratifying_range3 = [1, 33.3, 66.7, 99]
outliers_percentiles = [1, 99]


def prep_df(patients_file, cancer_cells_file, features, explanations_file):
    
    # Load patinets' clinical info file into a pandas DataFrame
    patients_data = pd.read_csv(patients_file)
    patients_data.rename(columns = {'TnumberAVL':'patient'}, inplace = True)

    # Load cancer cells file into a pandas DataFrame
    cancer_cells = pd.read_csv(cancer_cells_file)
    cancer_cells = cancer_cells[features]
    cancer_cells.rename(columns = {'Slide':'cycif.slide', 'CoreId':'cycif.core.id'}, inplace = True)
    
    # Load the text file into a pandas DataFrame
    slides_explanations = pd.read_csv(explanations_file, sep="	")
    slides_explanations['cycif.core.id']=slides_explanations['cycif.core.id'].str.replace("core", "").astype(cancer_cells['cycif.core.id'].dtype)
    
    # Merge the two data frames on the slide number and core id columns
    dataset = pd.merge(cancer_cells, slides_explanations[["cycif.slide", "cycif.core.id","patient"]], on=["cycif.slide", "cycif.core.id"])
    
    # Merge the two data frames on patient column
    final_dataset = pd.merge(dataset, patients_data[["patient","Molecular_profile","therapy_sequence","daystoprogression","timelastfu", "finalstatus"]], on=["patient"])
    
    print(final_dataset['cycif.slide'].unique())
    
    # Plot the distribution of 'daystoprogression' variable
    # final_dataset['daystoprogression'].plot(kind='hist', edgecolor='black')
    
#     # Stratify 'daystoprogression' to 3 groups: "short", "medium", "high"
    final_dataset['daystoprogression'] = pd.cut(final_dataset['daystoprogression'], np.percentile(final_dataset['daystoprogression'], stratifying_range3), labels=[0, 1, 2])
    
    
# #     final_dataset['daystoprogression'] = final_dataset['daystoprogression'].apply(np.log2)
#     final_dataset['daystoprogression'], _ = stats.boxcox(final_dataset['daystoprogression'].values)
#     # create strata based on the distribution of 'daystoprogression' variable so that each stratum contains an equal number of observations
#     final_dataset['daystoprogression'] = pd.qcut(final_dataset['daystoprogression'], 3, labels=["short", "medium", "high"])
    
    # Drop all rows (axis index = 0)  with NaN values
    final_dataset=final_dataset.dropna(axis=0)
    
    return final_dataset

# function to perform stratified random sampling
def stratified_sampling(df, col1, col2, sample_size):
    
    # To receive same number of patients independent of runs
    np.random.seed(0)
    
    # Sample_size - the number of patients to be selected from each sub-group of the data defined by a combination of col1 and col2
    
    # Create an empty dataframe for storing the down-sampled data
    down_sampled_df = pd.DataFrame()
    
    # Get the unique values of col1 and col2 in the input dataframe
    all_col1 = df[col1].unique()
    all_col2 = df[col2].unique()
    
    # Loop through the unique values of col1 and col2 to select a sub-group of data for each combination of col1 and col2
    for c1 in all_col1:
        for c2 in all_col2:
            sub_df = df[(df[col1] == c1) & (df[col2] == c2)]
            
            # If the size of the sub-group is larger than or equal to the specified sample size, down-sample the sub-group to the specified sample size
            if len(sub_df.groupby("patient")) >= sample_size:
                down_sampled_df = pd.concat([down_sampled_df, sub_df.sample(sample_size, replace=True)])
            # If the size of the sub-group is smaller than the specified sample size, keep all the data in the sub-group
            else:
                down_sampled_df = pd.concat([down_sampled_df, sub_df])

    # Print the size and number of unique values of each column in the down-sampled dataframe
    print("The size of down-sampled dataframe: ", down_sampled_df.shape)
    print("Patient number: ", len(down_sampled_df["patient"].unique()))
    print("Molecular profile number: ", len(down_sampled_df["Molecular_profile"].unique()))
    print("Therapy number: ", down_sampled_df["therapy_sequence"].unique())
    
    return down_sampled_df


def remove_outliers(df, features):
    
#     if flag == 'trim':
#         print('trim')
        
#     elif flag == 'remove':
#         # create a boolean mask to filter the rows
#         mask = df

#         # extract the subset of rows to filter
#         data = df.loc[features]
        
#         for column in data.columns:
#             # calculate the 1st and 99th percentiles for each column
#             percentiles = dataх[column].apply(lambda x: np.percentile(x, [1, 99]))
#     #         print(percentiles)
#             col_mask = (data[column] >= percentiles[column].iloc[0]) & (data[column] <= percentiles[column].iloc[1])
#             mask &= col_mask
#             num_removed = col_mask.size - col_mask.sum()
#             print(f"Removed {num_removed} rows for column {column}")

#         # remove the rows that contain outliers from the initial DataFrame - trim 
#         df = df.loc[~mask]
    
    
    slides = df['cycif.slide'].unique()
    
    for slide in slides:
        print("SLIDE: ", slide)
        # create a boolean mask to filter the rows
        mask = df['cycif.slide'] == slide
        mask.head(5)

#         # extract the subset of rows to filter
#         data = df.loc[mask, features]

#         # calculate the 1st and 99th percentiles for each column
#         percentiles = data.apply(lambda x: np.percentile(x, [1, 99]))
# #         print(percentiles)

#         for column in data.columns:
#             col_mask = (data[column] >= percentiles[column].iloc[0]) & (data[column] <= percentiles[column].iloc[1])
#             num_removed = col_mask.size - col_mask.sum()
#             print(f"Removed {num_removed} rows for column {column}")
#             mask &= col_mask
            

#         # remove the rows that contain outliers from the initial DataFrame - trim 
#         df = df.loc[~mask]


    
    return df


def transform_df(transformation, df, features):
    
    data = df.copy()
    if transformation == 'LOG':
        # Log-transform the features
        data[features] = data[features].apply(np.log)
    elif transformation == 'LOG2':
        data[features] = data[features].apply(np.log2)
    elif transformation == 'BOXCOX':
        # transform data & save lambda value
        for feature in features:
            
#             # Check that the values are positive
#             positive_cols = features

#             for col in positive_cols:
#                 if (df[col] <= 0).any():
#                     print(f"The column {col} contains non-positive values.")
#                 else:
#                     print(f"All values in column {col} are positive.")
          
            data[feature], _ = stats.boxcox(data[feature].values)
    else:
        print("The transformation was incorrectly written")
    
    return data

def scaler(df,scaler_type):
    
    # Make sure that features are numerical 
    
    scaler = SCALER_DICT.get(scaler_type)
    
    if scaler is None:
        raise ValueError(f"Invalid scaler type: {scaler_type}")
    
    # MinMaxScaler scales the data to be within a specific range, usually between 0 and 1
    # StandardScaler scales the data to have a mean of zero and a standard deviation of one, which is a commonly used scaling method for PCA.
    data = scaler.fit_transform(df)

    return data

def scaling(order, scaler_type, df, features):
    
    if order not in ['whole','patient', 'slide']:
        raise ValueError(f"Invalid order: {order}")
    
#     if order == 'patient' and patient is None:
#         raise ValueError("Missing patient value")
    
#     if order == 'slide' and slide is None:
#         raise ValueError("Missing slide value")
        
    data = df.copy()
    if order == 'patient':
        patients = data['patient'].unique()
        for patient in patients:
            data.loc[data['patient'] == patient, features] = scaler(data.loc[data['patient'] == patient, features], scaler_type)
    elif order == 'slide':
        slides = data['cycif.slide'].unique()
        for slide in slides:
            data.loc[data['cycif.slide'] == slide, features] = scaler(data.loc[data['cycif.slide'] == slide, features], scaler_type)
    elif order == 'whole':
        data[features] = scaler(data[features], scaler_type)
    
    return data


# Can be used only on standardly scaled data
def get_PCA(dataframe, target):
    df = dataframe[['DNA1', 'CD11c', 'CD207', 'GranzymeB', 'CD163', 'CD57', 'CD20', 'CD4', 'CD3d', 'CD8a', 'CD45RO', 'FOXP3', 'PD1', 'pTBK1', 'CD68', 'PDL1_2', 'CD15', 'CD11b', 'yH2AX', 'cPARP1', 'PDL1_488', 'PDL1_555', 'Ki67','Vimentin', 'MHCII', 'CK7', 'MHCI', 'ECadherin', 'aSMA', 'CD31', 'Area', 'Eccentricity', 'Roundness']]
    pca = PCA(33)
    principalComponents = pca.fit_transform(df)
    principal_df = pd.DataFrame(data = principalComponents[:,0:3])
    
    # Create 3D plot 
    # Get the first 3 components
    components = np.transpose(pca.components_[:3, :])

    # Project the data onto the first 3 components
    projected_data = pca.transform(df)[:, :3]

    # Create a Pandas DataFrame for the projected data
    df = pd.DataFrame({'PC1': projected_data[:, 0],
                       'PC2': projected_data[:, 1],
                       'PC3': projected_data[:, 2],
                      'Target': dataframe[target]})
    
#     # Remove the rows with NaN values in the target feature
#     df = df.dropna(subset=['Target'])

    # Create the interactive 3D scatter plot
    fig = px.scatter_3d(df, x='PC1', y='PC2', z='PC3', color = 'Target')
    fig.write_html('figure.html', auto_open=True)

    # Create scree plot
    PC_values = np.arange(pca.n_components_) + 1
    plt.plot(PC_values, pca.explained_variance_ratio_, 'ro-', linewidth=2)
    plt.title('Scree Plot')
    plt.xlabel('Principal Component')
    plt.ylabel('Proportion of Variance Explained')
    plt.show()
    
    # Create residual plot
    sns.residplot(x=principal_df[0], y=principal_df[1], data=principal_df)
    plt.show()
    
def create_corr_matrix(df):
    sns.set_theme(style="white")

    # Compute the correlation matrix
    corr = df.corr()

    # Generate a mask for the upper triangle
    mask = np.triu(np.ones_like(corr, dtype=bool))

    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(11, 9))

    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(230, 20, as_cmap=True)

    # Draw the heatmap with the mask and correct aspect ratio
    sns.heatmap(corr, vmin=-0.7, vmax=0.7, mask=mask, cmap=cmap, center=0,
                square=True, linewidths=.5, cbar_kws={"shrink": .5})
    

def corr_against_target(df,profile_lab):
    
    # calculate the Pearson's correlation coefficient for each feature against target 
    results = []
    for feature in df.columns:
        corr, p_value = pearsonr(df[feature], profile_lab)
        results.append((feature, corr, p_value))

    # sort the results by correlation in descending order
    results.sort(key=lambda x: x[1], reverse=True)

    # print the results
    for result in results:
        print('Feature:', result[0], 'Correlation:', result[1], 'P-value:', result[2])
        