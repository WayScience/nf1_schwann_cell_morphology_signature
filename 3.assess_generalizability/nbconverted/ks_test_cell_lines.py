#!/usr/bin/env python
# coding: utf-8

# # Compute KS-test results per feature between the two different dervivatives of the `ipn02.3 2λ` cell line
# 
# Plate 6 contains two derivatives of the cell line acquired from `iNFixion` and `MGH`.

# ## Import libraries

# In[1]:


import pathlib
import pandas as pd
from scipy import stats
from joblib import load


# ## Set results directory and load in model to get list of the features used

# In[2]:


# Set results directory
results_dir = pathlib.Path("./results")
results_dir.mkdir(exist_ok=True)

# Set data type for the ks test evaluation
data_type = "cleaned"

# Set suffix for data files if using QC or cleaned data
if data_type == "cleaned":
    suffix = "_qc"
else:
    suffix = ""

# Load in model
model = load(
    pathlib.Path(f"../1.train_models/data/trained_nf1_model{suffix}.joblib").resolve(
        strict=True
    )
)
model_features = list(model.feature_names_in_)

len(model_features)


# ## Load in Plate 6 normalized data

# In[3]:


# Set directory to find the plate 6 data from based on data type
directory = (
    "single_cell_profiles/cleaned_sc_profiles"
    if data_type == "cleaned"
    else "single_cell_profiles"
)

# Load in the normalized data
plate_6_norm = pd.read_parquet(
    pathlib.Path(
        f"/media/18tbdrive/1.Github_Repositories/nf1_schwann_cell_painting_data/3.processing_features/data/{directory}/Plate_6_sc_normalized.parquet"
    )
)


# ## Perform KS-test comparing the features between the two cell line derivatives for all genotypes and between genotypes

# In[4]:


# Split data by institution for comparison
institution_1_norm = plate_6_norm[plate_6_norm["Metadata_Institution"] == "iNFixion"]
institution_2_norm = plate_6_norm[plate_6_norm["Metadata_Institution"] == "MGH"]

# Perform KS-test for each feature
all_genotypes_ks_test_results_norm = {}

for column in plate_6_norm.columns:
    if column.startswith("Metadata_"):
        continue
    ks_stat, p_value = stats.kstest(
        institution_1_norm[column], institution_2_norm[column]
    )
    all_genotypes_ks_test_results_norm[column] = {
        "ks_stat": ks_stat,
        "p_value": p_value,
    }

# Convert results to DataFrame for better visualization
all_genotypes_ks_test_results_norm_df = (
    pd.DataFrame(all_genotypes_ks_test_results_norm)
    .T.reset_index()
    .rename(columns={"index": "feature"})
)


# In[5]:


# Split data by institution and WT genotype for comparison
institution_WT_1_norm = plate_6_norm[
    (plate_6_norm["Metadata_Institution"] == "iNFixion")
    & (plate_6_norm["Metadata_genotype"] == "WT")
]
institution_WT_2_norm = plate_6_norm[
    (plate_6_norm["Metadata_Institution"] == "MGH")
    & (plate_6_norm["Metadata_genotype"] == "WT")
]

# Perform KS-test for each feature for the WT genotype
WT_ks_test_results_norm = {}

for column in plate_6_norm.columns:
    if column.startswith("Metadata_"):
        continue
    ks_stat, p_value = stats.kstest(
        institution_WT_1_norm[column], institution_WT_2_norm[column]
    )
    WT_ks_test_results_norm[column] = {"ks_stat": ks_stat, "p_value": p_value}

# Convert results to DataFrame for better visualization
WT_ks_test_results_norm_df = (
    pd.DataFrame(WT_ks_test_results_norm)
    .T.reset_index()
    .rename(columns={"index": "feature"})
)


# In[6]:


# Split data by institution and Null genotype for comparison
institution_Null_1_norm = plate_6_norm[
    (plate_6_norm["Metadata_Institution"] == "iNFixion")
    & (plate_6_norm["Metadata_genotype"] == "Null")
]
institution_Null_2_norm = plate_6_norm[
    (plate_6_norm["Metadata_Institution"] == "MGH")
    & (plate_6_norm["Metadata_genotype"] == "Null")
]

# Perform KS-test for each feature for the Null genotype
Null_ks_test_results_norm = {}

for column in plate_6_norm.columns:
    if column.startswith("Metadata_"):
        continue
    ks_stat, p_value = stats.kstest(
        institution_Null_1_norm[column], institution_Null_2_norm[column]
    )
    Null_ks_test_results_norm[column] = {"ks_stat": ks_stat, "p_value": p_value}

# Convert results to DataFrame for better visualization
Null_ks_test_results_norm_df = (
    pd.DataFrame(Null_ks_test_results_norm)
    .T.reset_index()
    .rename(columns={"index": "feature"})
)


# In[7]:


# Add genotype column to each KS-test results DataFrame
WT_ks_test_results_norm_df["genotype_comparison"] = "WT"
Null_ks_test_results_norm_df["genotype_comparison"] = "Null"
all_genotypes_ks_test_results_norm_df["genotype_comparison"] = "All"

# Combine the two DataFrames
ks_test_results_norm_df = pd.concat(
    [
        WT_ks_test_results_norm_df,
        Null_ks_test_results_norm_df,
        all_genotypes_ks_test_results_norm_df,
    ],
    ignore_index=True,
)

# Print the combined results
print("\nKS-test results for normalized data:")
print(ks_test_results_norm_df.shape)
ks_test_results_norm_df.head()


# ## Add absolute value coefficients per feature from the model to the results (filtering down the data to only the features in the model)

# In[8]:


if data_type == "cleaned":
    # Load in the feature importance data from the QC model
    feat_import_df = pd.read_parquet(
        pathlib.Path(
            "../2.evaluate_model/model_evaluation_data/feature_importances_qc.parquet"
        )
    )
else:
    # Load in the feature importance data from non-QC model
    feat_import_df = pd.read_parquet(
        pathlib.Path(
            "../2.evaluate_model/model_evaluation_data/feature_importances.parquet"
        )
    )
print("Number of features in model:", feat_import_df.shape[0])

# Take the absolute value of the feature importance
feat_import_df["feature_importances"] = feat_import_df["feature_importances"].abs()

# Change the column name from feature_names to feature
feat_import_df = feat_import_df.rename(columns={"feature_names": "feature"})

# Merge the feature importance data with the KS test results
ks_test_results_norm_df = ks_test_results_norm_df.merge(feat_import_df, on="feature")

print(ks_test_results_norm_df.shape)
ks_test_results_norm_df.head()


# ## Split feature names into parts and save results

# In[9]:


# Split the feature column into parts
ks_test_results_norm_df[
    [
        "compartment",
        "feature_group",
        "measurement",
        "channel",
        "parameter1",
        "parameter2",
        "parameter3",
    ]
] = (
    ks_test_results_norm_df["feature"]
    .str.split("_", expand=True)
    .reindex(columns=range(7), fill_value=pd.NA)
)

# Save the results with qc suffix if data is cleaned
if data_type == "cleaned":
    ks_test_results_file = (
        pathlib.Path(results_dir) / "ks_test_derivatives_results_qc.parquet"
    )
else:
    ks_test_results_file = (
        pathlib.Path(results_dir) / "ks_test_derivatives_results.parquet"
    )

ks_test_results_norm_df.to_parquet(ks_test_results_file)

# Display the updated DataFrame
print(ks_test_results_norm_df.shape)
ks_test_results_norm_df.head()


# ## Print rows from the top five feature importances

# In[10]:


ks_test_results_norm_df = ks_test_results_norm_df.sort_values(
    by="feature_importances", ascending=False
)
ks_test_results_norm_df.head()

