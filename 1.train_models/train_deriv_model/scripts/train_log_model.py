#!/usr/bin/env python
# coding: utf-8

# # Train logistic regression model with all four plates

# ## Import libraries

# In[ ]:


import pathlib
import random
import warnings
from collections import defaultdict

import joblib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, precision_recall_curve
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import LabelEncoder
from sklearn.utils import parallel_backend


# ## Find the root of the git repo on the host system

# In[2]:


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

# In[3]:


# Set the seed
rng = np.random.default_rng(0)

# If the data (within the cell painting directory) is stored in a different location, add location here
repo_dir = pathlib.Path(
    root_dir / "/media/18tbdrive/1.Github_Repositories/nf1_schwann_cell_painting_data"
)

# Directory containing the feature selected parquet files (post-QC)
data_dir = (
    repo_dir / "3.processing_features/data/single_cell_profiles/cleaned_sc_profiles"
)

# Load in the four plate dataframes from data_dir
plate3df = pd.read_parquet(data_dir / "Plate_3_sc_feature_selected.parquet")
plate3pdf = pd.read_parquet(data_dir / "Plate_3_prime_sc_feature_selected.parquet")
plate5df = pd.read_parquet(data_dir / "Plate_5_sc_feature_selected.parquet")
plate6df = pd.read_parquet(data_dir / "Plate_6_sc_feature_selected.parquet")

# Add Metadata_Institution column to each plate dataframe (except plate 6)
plate3df["Metadata_Institution"] = "iNFixion"
plate3pdf["Metadata_Institution"] = "iNFixion"
plate5df["Metadata_Institution"] = "iNFixion"


# ## Split and processes data

# In[4]:


gene_column = "Metadata_genotype"


def down_sample_by_genotype(_df):
    """
    Return an equal number of cells from each genotype.
    The number of cells in a genotype is the minimum number of cells from all genotypes.

    Parameters
    ----------
    _df: Pandas Dataframe
        The data to be downsampled by the gene_column column.

    Returns
    -------
    The dataframe down-sampled by genotype.
    """

    min_gene = _df[gene_column].value_counts().min()
    return _df.groupby(gene_column, group_keys=False).apply(
        lambda x: x.sample(n=min_gene, random_state=0)
    )


def process_plates(_df):
    """
    Drop rows with nans from the single cell data and remove HET cells.

    Parameters
    ----------
    _df: Pandas Dataframe
        Uncleaned plate data with nans and HET cells to be removed. Contains the column "Metadata_genotype".

    Returns
    -------
    _df: Pandas Dataframe
        Cleaned single cell data by removing nans and HET cells.
    """

    _df.dropna(inplace=True)
    _df = _df.loc[_df[gene_column] != "HET"]
    return _df


def shuffle_data(_X):
    """
    Shuffle the columns of the input dataframe independently.

    Parameters
    ----------
    _X: Pandas Dataframe
        Input feature data for shuffling the columns.
    """

    for column in _X.columns:
        _X[column] = rng.permutation(_X[column])


def create_splits(_wells, _plate):
    """
    Create data splits for model training. The splits are rest (train and validation) and test.

    Parameters
    ----------
    _wells: List(String)
        The well names from which single cells will be used in the test set.

    _plate: Pandas Dataframe
        Single cell data from one of the plate's containing a "Metadata_Well" column.

    Returns
    -------
    Dataframes of the split single cell data.
    """

    return (
        _plate[~_plate["Metadata_Well"].isin(_wells)],
        _plate[_plate["Metadata_Well"].isin(_wells)],
    )


def store_pre_evaluation_data(_X, _y, _metadata, _datasplit):
    """
    Store model data to evaluate performance.

    Parameters
    ----------
    _X: Pandas Dataframe
        Feature dataframe from a given plate and data split.

    _y: Numpy Array
        A numerically-encoded label vector ordered according to _X.

    _metadata: Pandas Dataframe
        Plate name.

    _datasplit: String
        Data split name.
    """

    eval_data[f"probability_{probability_class}"].extend(
        logreg.predict_proba(_X)[:, 1].tolist()
    )
    eval_data["datasplit"].extend([_datasplit] * _X.shape[0])
    eval_data["predicted_genotype"].extend(logreg.predict(_X).tolist())
    eval_data["true_genotype"].extend(_y.tolist())
    for meta_col in _metadata.columns:
        eval_data[meta_col].extend(_metadata[meta_col].tolist())


# In[5]:


plate3df = process_plates(plate3df)
p3_wells = ["C11", "E11", "C3", "F3"]
rest3df, test3df = create_splits(p3_wells, plate3df)
rest3df, test3df = down_sample_by_genotype(rest3df), down_sample_by_genotype(test3df)

plate3pdf = process_plates(plate3pdf)
p3p_wells = ["F11", "G11", "C3", "F3"]
rest3pdf, test3pdf = create_splits(p3p_wells, plate3pdf)
rest3pdf, test3pdf = down_sample_by_genotype(rest3pdf), down_sample_by_genotype(
    test3pdf
)

plate5df = process_plates(plate5df)
p5_wells = ["C9", "E11", "E3", "G3"]
rest5df, test5df = create_splits(p5_wells, plate5df)
rest5df, test5df = down_sample_by_genotype(rest5df), down_sample_by_genotype(test5df)

plate6df = process_plates(plate6df)
p6_wells = ["E9", "B10", "D11", "D4"]
rest6df, test6df = create_splits(p6_wells, plate6df)
rest6df, test6df = down_sample_by_genotype(rest6df), down_sample_by_genotype(test6df)


# In[6]:


# Columns common to all plates
plate_cols = list(
    set(plate5df.columns)
    & set(plate3df.columns)
    & set(plate3pdf.columns)
    & set(plate6df.columns)
)

# Set up combined rest and test dataframes
restdf = pd.concat(
    [
        rest3df[plate_cols],
        rest3pdf[plate_cols],
        rest5df[plate_cols],
        rest6df[plate_cols],
    ],
    ignore_index=True,
).reset_index(drop=True)

testdf = pd.concat(
    [
        test3df[plate_cols],
        test3pdf[plate_cols],
        test5df[plate_cols],
        test6df[plate_cols],
    ],
    ignore_index=True,
).reset_index(drop=True)


# ## Encode genotypes and set up X,y data

# In[7]:


meta_cols = testdf.filter(like="Metadata").columns
feat_cols = testdf.drop(columns=meta_cols).columns


# In[8]:


le = LabelEncoder()

y = le.fit_transform(restdf["Metadata_genotype"])
X = restdf.drop(columns=meta_cols)

y_test = le.fit_transform(testdf["Metadata_genotype"])
X_test = testdf.drop(columns=meta_cols)

# Class for saving probabilities
probability_class = le.inverse_transform([1])[0]


# In[9]:


print(restdf["Metadata_genotype"].value_counts())


# ## Train logistic regression model

# ### Specify parameters for training

# In[10]:


model_dir = pathlib.Path("./models")
model_dir.mkdir(parents=True, exist_ok=True)
final_model_path = pathlib.Path(f"{model_dir}/best_final_logreg_model.joblib")


# In[11]:


logreg_params = {
    "max_iter": 250,
    "random_state": 0,
    "n_jobs": -1,
    "penalty": "l2",
}

# Random sampling range of hyperparameter
param_ranges = {"C": (0, 200)}

# Number of iteration to optimize hyperparameters
rand_iter = 500

# Best accuracy
best_acc = 0

# Initial accuracy
acc = 0

# Number of folds
n_splits = 8

# Generate hyperparameter samples
random_params = {
    i: {key: random.uniform(*param_ranges[key]) for key in param_ranges}
    for i in range(rand_iter)
}


# In[12]:


# Check if best model already exists
model_path = pathlib.Path("./models/best_final_logreg_model.joblib")
if model_path.exists():
    print(f"Model already exists at {model_path}, skipping hyperparameter search.")
else:
    # Store model results for evaluation
    eval_data = defaultdict(list)

    best_acc = 0  # Initialize best accuracy
    best_hp = None  # Track best hyperparameters

    # Iterate through hyperparameters
    for idx, rparams in random_params.items():

        skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=0)

        # Combine parameters in current search with logistic regression parameters
        comb_params = logreg_params | rparams
        acc = 0  # Accuracy accumulator for this parameter set

        # Loop through the folds
        for fold, (train_index, val_index) in enumerate(skf.split(X, y)):

            X_train, X_val = X.iloc[train_index], X.iloc[val_index]
            y_train, y_val = y[train_index], y[val_index]

            X_val_shuf = X_val.copy()
            shuffle_data(X_val_shuf)

            # Prevent the convergence warning in sklearn
            with parallel_backend("multiprocessing"):
                with warnings.catch_warnings():
                    warnings.filterwarnings(
                        "ignore", category=ConvergenceWarning, module="sklearn"
                    )
                    logreg = LogisticRegression(**comb_params)
                    logreg.fit(X_train, y_train)

            # Cumulative accuracy for all folds
            preds = logreg.predict(X_val)
            preds_shuf = logreg.predict(X_val_shuf)
            acc += accuracy_score(y_val, preds)

            store_pre_evaluation_data(
                X_val, y_val, restdf.iloc[val_index][meta_cols], "val"
            )
            store_pre_evaluation_data(
                X_val_shuf, y_val, restdf.iloc[val_index][meta_cols], "shuffled_val"
            )

        # Average accuracy for the folds
        acc = acc / n_splits

        # Store the data with the best performance
        if acc > best_acc:
            best_hparam = eval_data.copy()
            best_acc = acc
            best_hp = rparams

    print(f"Best average validation accuracy = {best_acc}")


# Due to error in the code, the cells were not run in order. We want to make sure they are ran in order, but not rerun the model which takes 7 hours to train. 
# 
# As reference, the best average validation accuracy for the hyperparameter search is **0.9173886180912052**.

# ## Train final optimized model

# In[13]:


# Save path
model_dir = pathlib.Path("./models")
model_dir.mkdir(parents=True, exist_ok=True)
model_path = model_dir / "best_final_logreg_model.joblib"

# Only train if the model doesn't already exist
if model_path.exists():
    print(f"Final model already exists at {model_path}, skipping training.")
else:
    logreg_params = {
        "max_iter": 3000,
        "random_state": 0,
        "n_jobs": -1,
        "penalty": "l2",
    }

    comb_params = logreg_params | best_hp

    logreg = LogisticRegression(**comb_params)
    logreg.fit(X, y)

    joblib.dump(logreg, model_path)


# In[14]:


# Path to saved pre-evaluation results
data_path = pathlib.Path("./data")
data_path.mkdir(parents=True, exist_ok=True)
eval_file = data_path / "nf1_model_pre_evaluation_results.parquet"

if eval_file.exists():
    print(f"Pre-evaluation data already exists at {eval_file}, skipping.")
else:
    X_shuf = X.copy()
    shuffle_data(X_shuf)

    X_test_shuf = X_test.copy()
    shuffle_data(X_test_shuf)

    store_pre_evaluation_data(X, y, restdf[meta_cols], "train")
    store_pre_evaluation_data(X_shuf, y, restdf[meta_cols], "shuffled_train")

    store_pre_evaluation_data(X_test, y_test, testdf[meta_cols], "test")
    store_pre_evaluation_data(X_test_shuf, y_test, testdf[meta_cols], "shuffled_test")

    # Store the evaluation data
    pd.DataFrame(eval_data).to_parquet(eval_file)


# ## Extract PR curve metrics and coefficients per feature

# In[15]:


# Define evaluation metric data
# The "metrics" include precision, recall
eval_mets = {"precision_recall": defaultdict(list)}


def compute_metrics(_df, _plate, _split):
    """
    Parameters
    ----------
    _df: Pandas Dataframe
        Model data to be evaluated.

    _plate: String
        Name of the plate for storing the metrics

    _split: String
        Name of the data split for storing the metric
    """

    y_true = _df[gene_column]
    y_proba = _df["probability_WT"]

    # Store precision and recall data
    precision, recall, _ = precision_recall_curve(y_true, y_proba, pos_label="WT")
    pr_size = precision.shape[0]
    eval_mets["precision_recall"]["precision"].extend(precision.tolist())
    eval_mets["precision_recall"]["recall"].extend(recall.tolist())
    eval_mets["precision_recall"]["plate"].extend([_plate] * pr_size)
    eval_mets["precision_recall"]["datasplit"].extend([_split] * pr_size)


# In[16]:


# Make directory to save evaluation metrics
eval_path = pathlib.Path("./pr_results")
eval_path.mkdir(parents=True, exist_ok=True)

# Path to the single metrics file
metrics_file = eval_path / "precision_recall_final_model.parquet"

# If it exists, skip everything and print a message
if metrics_file.exists():
    print(f"Metrics already exist at {metrics_file}, skipping computation.")
else:
    # Set eval_data as a DataFrame to process
    eval_data = pd.DataFrame(eval_data)

    # Iterate through each data split
    for split in eval_data["datasplit"].unique():

        # Calculate metrics for all plates
        df_temp = eval_data.loc[(eval_data["datasplit"] == split)].copy()
        compute_metrics(df_temp, "all_plates", split)

        # Calculate metrics for each plate
        for plate in eval_data["Metadata_Plate"].unique():
            df_temp = eval_data.loc[
                (eval_data["Metadata_Plate"] == plate)
                & (eval_data["datasplit"] == split)
            ].copy()
            df_temp = down_sample_by_genotype(df_temp)
            compute_metrics(df_temp, plate, split)

        # Calculate metrics for each institution *only* from Plate_6 rows
        eval_data_Plate_6 = eval_data.loc[eval_data["Metadata_Plate"] == "Plate_6"]

        for institution in eval_data_Plate_6["Metadata_Institution"].unique():
            df_temp = eval_data_Plate_6.loc[
                (eval_data_Plate_6["Metadata_Institution"] == institution)
                & (eval_data_Plate_6["datasplit"] == split)
            ].copy()
            df_temp = down_sample_by_genotype(df_temp)

            # Rename institution for saving
            if institution == "iNFixion":
                inst_name = "Plate_6_orig"
            elif institution == "MGH":
                inst_name = "Plate_6_deriv"
            else:
                inst_name = f"Plate_6_{institution}"

            compute_metrics(df_temp, inst_name, split)

    # Save all computed metrics
    for met, met_data in eval_mets.items():
        pd.DataFrame(eval_mets[met]).to_parquet(metrics_file)


# In[17]:


# Load model if not already in memory
if "logreg" not in locals():
    model_path = pathlib.Path("./models/best_final_logreg_model.joblib")
    logreg = joblib.load(model_path)
    print(f"Loaded trained model from {model_path}")

# Create the output directory if it doesn't exist
coeff_dir = pathlib.Path("./coeff_results")
coeff_dir.mkdir(parents=True, exist_ok=True)

# Create a DataFrame with features and their coefficients
coeff_df = pd.DataFrame(
    {"feature": logreg.feature_names_in_, "coefficient": logreg.coef_.reshape(-1)}
)

# Save to CSV in the coeff_results folder
coeff_df.to_csv(coeff_dir / "final_model_coefficients.csv", index=False)

# Display the first few rows
print(coeff_df.shape)
coeff_df.head()


# In[18]:


# Sort coefficients by absolute value, descending
sorted_coeff_df = coeff_df.reindex(
    coeff_df["coefficient"].abs().sort_values(ascending=False).index
)

# Print the sorted coefficients
print(sorted_coeff_df)


# ## Load in original model coefficients and outer merge

# In[19]:


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


# In[ ]:


spearman_corr, p_value = spearmanr(
    merged_coefs["coefficient_orig_model"], merged_coefs["coefficient_new_model"]
)

print(f"Spearman correlation: {spearman_corr:.4f}, p-value: {p_value:.4e}")

