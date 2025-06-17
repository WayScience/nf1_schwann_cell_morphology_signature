#!/usr/bin/env python
# coding: utf-8

# # Train new simple logistic regression model with all four plates
# 
# The original model was trained with more complexity to optimize performance using three plates. A new plate was collected with a derivative of the cell lines used in the first three plates. We are now including that plate to train a simple model where:
# 
# 1. One random well per plate per genotype are excluded as holdout wells (4 plates X 2 genotypes = 8 wells in holdout set)
# 2. 70 training 30 testing splits (stratified by plate and genotype)

# ## Import libraries

# In[1]:


import pathlib
import random

import pandas as pd
import joblib
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.base import clone
from sklearn.metrics import precision_recall_curve
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RandomizedSearchCV, StratifiedKFold
from sklearn.preprocessing import LabelEncoder
from sklearn.utils import parallel_backend
import warnings


# ## Helper functions

# In[2]:


# Helper to get X and y from a dataframe
def get_X_y(df: pd.DataFrame, label_encoder: LabelEncoder) -> tuple:
    """Generate X and y from a dataframe.

    Args:
        df (pd.DataFrame): pandas DataFrame containing the data for the split.
        label_encoder (LabelEncoder): LabelEncoder instance to transform the labels.

    Returns:
        tuple: Returns the X and y data as a tuple.
    """
    meta_cols = df.filter(like="Metadata").columns
    X = df.drop(columns=meta_cols)
    y = label_encoder.transform(df["Metadata_genotype"])
    return X, y


def get_pr_df(
    X: pd.DataFrame,
    y: np.array,
    model: LogisticRegression,
    model_name: str,
    split_name: str,
    plate: str = None,
    institution: str = None,
) -> pd.DataFrame:
    """Return a DataFrame of precision-recall points with attached metadata.

    Args:
        X (pd.DataFrame): X features for the model.
        y (np.array): y labels for the model.
        model (LogisticRegression): Model to use for prediction.
        model_name (str): Model name for the DataFrame (e.g., final or shuffled).
        split_name (str):  Name of the split (e.g., "train", "test", "holdout").
        plate (str, optional): Name of plate. Defaults to None.
        institution (str, optional): Name of institution cell line is derived from. Defaults to None.

    Raises:
        ValueError: If the model does not produce two class probabilities.

    Returns:
       pd.DataFrame: DataFrame containing precision, recall, thresholds, and metadata.
    """
    proba = model.predict_proba(X)
    if proba.shape[1] != 2:
        raise ValueError("Expected binary classification with two class probabilities.")
    pos_class_proba = proba[:, 1]

    precision, recall, thresholds = precision_recall_curve(y, pos_class_proba)
    thresholds = np.append(thresholds, np.nan)

    pr_len = len(precision)
    return pd.DataFrame(
        {
            "model": [model_name] * pr_len,
            "split": [split_name] * pr_len,
            "precision": precision,
            "recall": recall,
            "threshold": thresholds,
            "plate": [plate] * pr_len,
            "institution": [institution] * pr_len,
        }
    )


# ## Find the root of the git repo on the host system

# In[3]:


# Get the current working directory
cwd = pathlib.Path.cwd()

if (cwd / ".git").is_dir():
    root_dir = cwd

else:
    root_dir = None
    for parent in cwd.parents:
        if (parent / ".git").is_dir():
            root_dir = parent
            break

# Check if a Git root directory was found
if root_dir is None:
    raise FileNotFoundError("No Git root directory found.")


# ## Load in all of the feature selected plates (post-filtering/single cell QC) and concat with the common features

# ### Set paths

# In[4]:


# If the data (within the cell painting directory) is stored in a different location, add location here
repo_dir = pathlib.Path(
    root_dir / "/media/18tbdrive/1.Github_Repositories/nf1_schwann_cell_painting_data"
)

# Directory containing the feature selected parquet files (post-QC)
data_dir = (
    repo_dir / "3.processing_features/data/single_cell_profiles/cleaned_sc_profiles"
)


# ## Load in and concat dataframes

# In[5]:


# Define plate names
plate_names = ["Plate_3", "Plate_3_prime", "Plate_5", "Plate_6"]

# Load and filter dataframes
dfs = []
for plate in plate_names:
    file_path = data_dir / f"{plate}_sc_feature_selected.parquet"
    if file_path.exists():
        df = pd.read_parquet(file_path)
        if plate == "Plate_6":
            df = df[df["Metadata_genotype"] != "HET"]  # Drop HET samples from Plate_6
        # Ensure Metadata_Institution exists
        if "Metadata_Institution" not in df.columns:
            df["Metadata_Institution"] = "iNFixion"
        if not df.empty:
            dfs.append(df)
        else:
            print(f"Warning: {file_path} is empty.")
    else:
        print(f"Warning: {file_path} does not exist.")

# Define metadata columns as any column with 'Metadata_' prefix from the first dataframe
if dfs:
    metadata_cols = [col for col in dfs[0].columns if col.startswith("Metadata_")]
    feature_sets = [set(df.columns) - set(metadata_cols) for df in dfs]
    common_features = sorted(set.intersection(*feature_sets))
    selected_cols = metadata_cols + common_features

    combined_df = pd.concat(
        [df.loc[:, df.columns.intersection(selected_cols)] for df in dfs],
        ignore_index=True,
    )

    print(f"Combined dataframe shape: {combined_df.shape}")
else:
    print("No valid dataframes found.")

# Print the number of features in this combined dataframe
print(f"Number of features in combined dataframe: {len(common_features)}")

# Print the combined dataframe
combined_df.head()


# ## Split data into training, testing, and holdout sets

# ### Generate holdout data by selecting one random well per genotype in each plate

# In[6]:


# Set random seed for reproducibility
random.seed(0)

# Select one random well per genotype for each plate, with at least 100 single cells
holdout_indices = []

# Define acceptable range of single cells per well (avoid huge differences in cell counts)
min_cells = 100
max_cells = 200  # adjust as needed

# Select one random well per genotype per plate, within the cell count range
for plate in plate_names:
    plate_df = combined_df[combined_df["Metadata_Plate"] == plate]
    for genotype in plate_df["Metadata_genotype"].unique():
        wells_df = plate_df[
            (plate_df["Metadata_genotype"] == genotype)
            & (
                plate_df["Metadata_number_of_singlecells"]
                .astype(int)
                .between(min_cells, max_cells)
            )
        ]
        wells = wells_df["Metadata_Well"].unique()
        if len(wells) > 0:
            selected_well = random.choice(list(wells))
            print(
                f"Plate: {plate}, Genotype: {genotype}, Selected well: {selected_well}"
            )
            well_indices = wells_df[
                wells_df["Metadata_Well"] == selected_well
            ].index.tolist()
            holdout_indices.extend(well_indices)
        else:
            print(
                f"Plate: {plate}, Genotype: {genotype}, No well with {min_cells}-{max_cells} single cells."
            )

# Create the holdout dataframe using the selected indices
holdout_df = combined_df.loc[holdout_indices].copy()

# Print shape and data of the holdout dataframe
print("Holdout dataframe shape:", holdout_df.shape)
holdout_df.head()


# In[7]:


print(holdout_df["Metadata_number_of_singlecells"].unique())


# ### Generate training and testing split data

# In[8]:


# Drop holdout indices from combined_df and make a copy
train_test_df = combined_df.drop(index=holdout_indices).copy()

# Perform train/test split (70/30) stratified by plate and genotype
train_df, test_df = train_test_split(
    train_test_df,
    test_size=0.3,
    random_state=0,
    stratify=train_test_df[["Metadata_Plate", "Metadata_genotype"]],
)

# Show resulting shapes
print(f"Train shape: {train_df.shape}")
print(f"Test shape: {test_df.shape}")


# In[9]:


print(train_df["Metadata_genotype"].value_counts())


# ## Train simple logistic regression model

# In[10]:


# Extract metadata and feature columns
meta_cols = train_df.filter(like="Metadata").columns
feat_cols = train_df.drop(columns=meta_cols).columns

# Initialize the label encoder
le = LabelEncoder()

# Encode the genotype labels and prepare feature matrices
y_train = le.fit_transform(train_df["Metadata_genotype"])
X_train = train_df.drop(columns=meta_cols)

# Create a shuffled version of the X_train DataFrame
rng = np.random.default_rng(0)
X_train_shuffled = X_train.copy()
for column in X_train_shuffled.columns:
    X_train_shuffled[column] = rng.permutation(X_train_shuffled[column].values)


# In[11]:


# Set folds for k-fold cross validation (default is 5)
straified_k_folds = StratifiedKFold(n_splits=5, shuffle=False)

# Set Logistic Regression model parameters (use default for max_iter)
logreg_params = {
    "penalty": "elasticnet",
    "solver": "saga",
    "max_iter": 1000,
    "n_jobs": -1,
    "random_state": 0,
    "class_weight": "balanced",
}

# Define the hyperparameter search space for RandomizedSearchCV
param_dist = {
    "C": np.logspace(-3, 3, 7),
    "l1_ratio": np.linspace(0, 1, 11),
}

# Set the random search hyperparameterization method parameters
random_search_params = {
    "param_distributions": param_dist,
    "scoring": "f1_weighted",
    "random_state": 0,
    "n_jobs": -1,
    "cv": straified_k_folds,
}


# In[12]:


# Create the output directory if it doesn't exist
model_dir = pathlib.Path("./models")
model_dir.mkdir(parents=True, exist_ok=True)

# Initialize Logistic Regression and RandomizedSearchCV
logreg = LogisticRegression(**logreg_params)
random_search = RandomizedSearchCV(logreg, **random_search_params)

# Prevent convergence warnings
with parallel_backend("multiprocessing"):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=ConvergenceWarning, module="sklearn")

        ########################################################
        # Train the model on non-shuffled (final) training data
        ########################################################
        print("Training model on original (non-shuffled) data...")
        final_random_search = clone(random_search)
        final_random_search.fit(X_train, y_train)
        final_model = final_random_search.best_estimator_
        print("Optimal parameters (final):", final_random_search.best_params_)

        # Save the final model
        joblib.dump(final_model, f"{model_dir}/final_logreg_model.joblib")

        ########################################################
        # Train the model on shuffled training data
        ########################################################
        print("Training model on shuffled data...")
        shuffled_random_search = clone(random_search)
        shuffled_random_search.fit(X_train_shuffled, y_train)
        shuffled_model = shuffled_random_search.best_estimator_
        print("Optimal parameters (shuffled):", shuffled_random_search.best_params_)

        # Save the shuffled model
        joblib.dump(shuffled_model, f"{model_dir}/shuffled_logreg_model.joblib")


# ## Extract PR curve metrics and coefficients per feature

# In[13]:


# Models and data
datasets = [("train", train_df), ("test", test_df), ("holdout", holdout_df)]
models = [("final", final_model), ("shuffled", shuffled_model)]

# 1. PR curves per split across all plates
all_split_results = []

# 2. PR curves per plate
per_plate_results = []

# 3. PR curves per institution for Plate 6 only (test + holdout)
plate6_per_inst_results = []

for model_name, model in models:
    for split_name, df in datasets:
        # 1. PR curve for the entire split
        X_all, y_all = get_X_y(df, le)
        pr_df_all = get_pr_df(
            X_all, y_all, model, model_name, split_name, plate="ALL", institution="ALL"
        )
        all_split_results.append(pr_df_all)

        # 2. PR curves per plate
        for plate_name, plate_df in df.groupby("Metadata_Plate"):
            X_plate, y_plate = get_X_y(plate_df, le)
            pr_df_plate = get_pr_df(
                X_plate,
                y_plate,
                model,
                model_name,
                split_name,
                plate=plate_name,
                institution="ALL",
            )
            per_plate_results.append(pr_df_plate)

            # 3. PR curves per institution *within Plate 6* for test data
            if plate_name == "Plate_6" and split_name in ["test"]:
                for inst_name, inst_df in plate_df.groupby("Metadata_Institution"):
                    X_inst, y_inst = get_X_y(inst_df, le)
                    pr_df_inst = get_pr_df(
                        X_inst,
                        y_inst,
                        model,
                        model_name,
                        split_name,
                        plate=plate_name,
                        institution=inst_name,
                    )
                    plate6_per_inst_results.append(pr_df_inst)

# Combine to separate DataFrames
df_all_split = pd.concat(all_split_results, ignore_index=True)
df_per_plate = pd.concat(per_plate_results, ignore_index=True)
df_plate6_per_institution = pd.concat(plate6_per_inst_results, ignore_index=True)

# Save results to parquet files
metrics_dir = pathlib.Path("./pr_results")
metrics_dir.mkdir(parents=True, exist_ok=True)

df_all_split.to_parquet(metrics_dir / "pr_curve_all_plates.parquet", index=False)
df_per_plate.to_parquet(metrics_dir / "pr_curve_per_plate.parquet", index=False)
df_plate6_per_institution.to_parquet(
    metrics_dir / "pr_curve_plate6_per_institution.parquet", index=False
)


# In[14]:


# Create the output directory if it doesn't exist
coeff_dir = pathlib.Path("./coeff_results")
coeff_dir.mkdir(parents=True, exist_ok=True)

# Extract feature names and coefficients from the final model
feature_names = X_train.columns
coefficients = final_model.coef_.flatten()

# Create a DataFrame with features and their coefficients
coeff_df = pd.DataFrame({"feature": feature_names, "coefficient": coefficients})

# Save to CSV in the coeff_results folder
coeff_df.to_csv(coeff_dir / "final_model_coefficients.csv", index=False)

# Display the first few rows
print(coeff_df.shape)
coeff_df.head()


# ## Load in original model coefficients and outer merge

# In[15]:


# Load in the original model coefficients
original_coeff_df = pd.read_parquet(
    "../../2.evaluate_model/model_evaluation_data/feature_importances_qc.parquet"
)

# Rename columns to 'feature' and 'coefficient'
original_coeff_df = original_coeff_df.rename(
    columns={
        original_coeff_df.columns[0]: "feature",
        original_coeff_df.columns[1]: "coefficient",
    }
)

# Perform an outer merge of the original model and the new model
merged_coefs = pd.merge(
    original_coeff_df,
    coeff_df,
    on="feature",
    how="outer",
    suffixes=("_orig_model", "_new_model"),
)

# Fill NaN values with 0
merged_coefs.fillna(0, inplace=True)

# Save the merged coefficients to a CSV file
merged_coefs.to_csv(
    coeff_dir / "merged_coefficients_original_new_model.csv", index=False
)

# Display the merged dataframe
print(merged_coefs.shape)
merged_coefs.head()

