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

# Load in model
model = load(pathlib.Path("../1.train_models/data/trained_nf1_model.joblib"))
model_features = list(model.feature_names_in_)

len(model_features)


# ## Load in Plate 6 normalized data

# In[3]:


# Load in the normalized data
plate_6_norm = pd.read_parquet(
    pathlib.Path(
        "/media/18tbdrive/1.Github_Repositories/nf1_schwann_cell_painting_data/3.processing_features/data/single_cell_profiles/Plate_6_sc_normalized.parquet"
    )
)


# ## Perform KS-test comparing the features between the two cell line derivatives

# In[4]:


# Split data by institution for comparison
institution_1_norm = plate_6_norm[plate_6_norm["Metadata_Institution"] == "iNFixion"]
institution_2_norm = plate_6_norm[plate_6_norm["Metadata_Institution"] == "MGH"]

# Perform KS-test for each feature
ks_test_results_norm = {}

for column in plate_6_norm.columns:
    if column.startswith("Metadata_"):
        continue
    ks_stat, p_value = stats.kstest(
        institution_1_norm[column], institution_2_norm[column]
    )
    ks_test_results_norm[column] = {"ks_stat": ks_stat, "p_value": p_value}

# Convert results to DataFrame for better visualization
ks_test_results_norm_df = (
    pd.DataFrame(ks_test_results_norm)
    .T.reset_index()
    .rename(columns={"index": "feature"})
)


# In[5]:


print("\nKS-test results for normalized data:")
ks_test_results_norm_df.head()


# ## Add absolute value coefficients per feature to the results

# In[6]:


feat_import_df = pd.read_parquet(
    pathlib.Path(
        "../2.evaluate_model/model_evaluation_data/feature_importances.parquet"
    )
)

# Take the absolute value of the feature importance
feat_import_df["feature_importances"] = feat_import_df["feature_importances"].abs()

# Change the column name from feature_names to feature
feat_import_df = feat_import_df.rename(columns={"feature_names": "feature"})

# Merge the feature importance data with the KS test results
ks_test_results_norm_df = ks_test_results_norm_df.merge(feat_import_df, on="feature")

ks_test_results_norm_df.head()


# ## Split feature names into parts and save results

# In[7]:


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

# Filter out features not in model_features
ks_test_results_norm_df = ks_test_results_norm_df[
    ks_test_results_norm_df["feature"].isin(model_features)
]

# Save the results
ks_test_results_norm_df.to_parquet(pathlib.Path(f"{results_dir}/ks_test_derivatives_results.parquet"))

# Display the updated DataFrame
print(ks_test_results_norm_df.shape)
ks_test_results_norm_df.head()

