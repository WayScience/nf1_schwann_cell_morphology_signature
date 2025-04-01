#!/usr/bin/env python
# coding: utf-8

# # Evaluate generalizability of the model across holdout plate(s)

# ## Import libraries

# In[1]:


import pandas as pd
import pathlib
from joblib import load
import numpy as np
from sklearn.metrics import (
    accuracy_score,
    precision_recall_curve,
)
import seaborn as sns
import matplotlib.pyplot as plt


# ## Set paths and variables

# In[2]:


# Set data type for the generalizability evaluation
data_cleaned = "cleaned"

# Set suffix for data files if using QC or cleaned data
if data_cleaned == "cleaned":
    suffix = "_qc"
else:
    suffix = ""

# Path to folder holding model and encoder files
model_dir = pathlib.Path("../1.train_models/data")

# Load in the model encoder
le = load(pathlib.Path(f"{model_dir}/trained_nf1_model_label_encoder{suffix}.joblib"))

# Load in the model
model = load(pathlib.Path(f"{model_dir}/trained_nf1_model{suffix}.joblib"))

# Path to results directory
results_dir = pathlib.Path("./results")
results_dir.mkdir(exist_ok=True)

# Set the random seed
rng = np.random.default_rng(0)


# ## Load in plate with two cell lines (Plate 6)

# In[3]:


# Set directory to find the plate 6 data from based on data type
directory = (
    "single_cell_profiles/cleaned_sc_profiles"
    if data_cleaned == "cleaned"
    else "single_cell_profiles"
)

# Read in data from plate 6 with two cell lines
plate6_df = pd.read_parquet(
    pathlib.Path(
        f"/media/18tbdrive/1.Github_Repositories/nf1_schwann_cell_painting_data/3.processing_features/data/{directory}/Plate_6_sc_normalized.parquet"
    )
)

# Count rows before dropping NaNs
initial_count = plate6_df.shape[0]

# Drop rows with NaNs
plate6_df = plate6_df.dropna()

# Count rows after dropping NaNs
final_count = plate6_df.shape[0]

# Print the count of dropped rows
print(f"Dropped rows: {initial_count - final_count}")

# Print shape and head of data
print(plate6_df.shape)
plate6_df.head()


# ## Generate a shuffled dataset from the loaded in plate

# In[4]:


# Shuffle the features randomly, excluding columns that start with "Metadata_"
shuffled_plate6_df = plate6_df.apply(
    lambda x: rng.permutation(x) if not x.name.startswith("Metadata_") else x
)

# Print shape and head of data
print(shuffled_plate6_df.shape)
shuffled_plate6_df.head()


# ## Apply model to final and shuffled versions of the plate data

# In[5]:


# Create list of the metadata columns only
meta_cols = [col for col in plate6_df.columns if "Metadata" in col]

# Define a dictionary to handle both data types
data_dict = {"final": plate6_df, "shuffled": shuffled_plate6_df}

# Initialize a list to store processed dataframes
processed_dfs = []

# Loop through the data dictionary to create probability dataframes
for data_type, data in data_dict.items():
    # Ensure no duplicates in data and reset index
    data = data.drop_duplicates().reset_index(drop=True)

    # Predict probabilities and labels
    probabilities = model.predict_proba(data[model.feature_names_in_])[:, 1]
    predicted_genotype = model.predict(
        data[model.feature_names_in_]
    )  # outputs as binary labels

    # Make a copy of the column to avoid modifying the original dataframe
    true_genotype = data["Metadata_genotype"].copy()

    # Set HET values to 2 explicitly
    true_genotype.loc[true_genotype == "HET"] = 2

    # Use label encoder for the remaining values (excluding HET)
    mask = true_genotype != 2  # Identify rows that are not HET
    true_genotype.loc[mask] = le.transform(true_genotype.loc[mask])

    # Convert dtype to integer
    true_genotype = true_genotype.astype(int)

    # Create a dataframe with probabilities and predictions
    probability_df = pd.DataFrame(
        {
            "probability_WT": probabilities,
            "predicted_genotype": predicted_genotype,
            "true_genotype": true_genotype,
            "data_type": data_type,
        },
        index=data.index,  # Ensure alignment with original data
    )

    # Add metadata columns (reset index to align lengths)
    metadata_df = data[meta_cols].reset_index(drop=True)
    assert len(probability_df) == len(
        metadata_df
    ), "Row count mismatch between probabilities and metadata!"

    full_df = pd.concat([probability_df, metadata_df], axis=1)
    processed_dfs.append(full_df)

# Combine all dataframes
combined_df = pd.concat(processed_dfs, axis=0).reset_index(drop=True)

# Save to Parquet with qc suffix if data is cleaned
if data_cleaned == "cleaned":
    output_file = (
        pathlib.Path(results_dir) / "plate_6_single_cell_probabilities_qc.parquet"
    )
else:
    output_file = (
        pathlib.Path(results_dir) / "plate_6_single_cell_probabilities.parquet"
    )

combined_df.to_parquet(output_file)

# Print shape and head of data
print(combined_df.shape)
combined_df.head()


# ## Split the probability data by Institution

# In[6]:


# Create dictionary with the split dataframes based on Institution
institution_dfs = {
    institution: combined_df[combined_df["Metadata_Institution"] == institution].copy()
    for institution in combined_df["Metadata_Institution"].unique()
}


# ## Generate PR curve results (for pre-visualization)

# In[7]:


precision_recall_data = []

for institution, df in institution_dfs.items():
    for data_type in ["final", "shuffled"]:  # Compute separately for both types
        # Subset for data type and remove the HET cells from evaluation
        subset_df = df[
            (df["data_type"] == data_type) & (df["Metadata_genotype"] != "HET")
        ]

        # Compute precision-recall curve
        precision, recall, _ = precision_recall_curve(
            subset_df["true_genotype"], subset_df["probability_WT"]
        )

        institution_results = pd.DataFrame(
            {
                "Precision": precision[:-1],
                "Recall": recall[:-1],
                "Metadata_Institution": institution,
                "data_type": data_type,
            }
        )

        precision_recall_data.append(institution_results)

# Combine all institution-based PR data
precision_recall_df = pd.concat(precision_recall_data, ignore_index=True)

# Save PR curve data to parquet file with qc suffix if data is cleaned
if data_cleaned == "cleaned":
    pr_curve_file = (
        pathlib.Path(results_dir) / "plate6_precision_recall_final_model_qc.parquet"
    )
else:
    pr_curve_file = (
        pathlib.Path(results_dir) / "plate6_precision_recall_final_model.parquet"
    )

precision_recall_df.to_parquet(pr_curve_file)

print(precision_recall_df.shape)
precision_recall_df.head()


# In[8]:


# Set the style of the plot
sns.set_theme(style="whitegrid")

# Create a figure and axis
plt.figure(figsize=(10, 6))

# Define a color palette based on Set2 (you can adjust n_colors to match your needs)
institution_palette = sns.color_palette("Dark2", n_colors=8)

# Create a mapping dictionary of institutions to specific colors from Set2
institution_color_map = {
    "MGH": institution_palette[2],
    "iNFixion": institution_palette[3],
}

# Plot the data
sns.lineplot(
    data=precision_recall_df,
    x="Recall",
    y="Precision",
    hue="Metadata_Institution",
    style="data_type",
    palette=institution_color_map,
    dashes=True,
)

# Set y-axis limits
plt.ylim(0, 1)

# Add labels and title
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.title("Precision vs Recall for Different Institutions and Data Types")
plt.legend(loc="lower right", bbox_to_anchor=(1, 0))
plt.show()


# ## Generate accuracy scores per institution and data type (final or shuffled)

# In[9]:


# Calculate accuracy per institution and data type (final or shuffled) without the HET cells
accuracy_per_group = (
    combined_df[combined_df["Metadata_genotype"] != "HET"]
    .groupby(["Metadata_Institution", "data_type"])
    .apply(lambda x: accuracy_score(x["true_genotype"], x["predicted_genotype"]))
    .reset_index(name="accuracy")
)

# Save accuracy data to parquet file with qc suffix if data is cleaned
if data_cleaned == "cleaned":
    accuracy_file = pathlib.Path(results_dir) / "plate6_accuracy_final_model_qc.parquet"
else:
    accuracy_file = pathlib.Path(results_dir) / "plate6_accuracy_final_model.parquet"

accuracy_per_group.to_parquet(accuracy_file)

accuracy_per_group


# ## Generate bar plot (for pre-visualization)

# In[10]:


# Set the style of the plot
sns.set_theme(style="whitegrid")

# Create a figure and axis
plt.figure(figsize=(10, 6))

# Create a bar plot
sns.barplot(
    data=accuracy_per_group,
    x="data_type",
    y="accuracy",
    hue="Metadata_Institution",
    palette="Dark2",
    errorbar=None,
)

# Set y-axis limits
plt.ylim(0, 1)

# Add labels and title
plt.xlabel("Genotype")
plt.ylabel("Accuracy")
plt.title("Accuracy per Genotype and Institution")
plt.legend(title="Institution")
plt.show()


# ## Look at how the accuracies break down per genotype

# In[11]:


# Calculate accuracy per genotype, institution and data type (final or shuffled)
accuracy_per_group = (
    combined_df[combined_df["Metadata_genotype"] != "HET"]
    .groupby(["Metadata_genotype", "Metadata_Institution", "data_type"])
    .apply(lambda x: accuracy_score(x["true_genotype"], x["predicted_genotype"]))
    .reset_index(name="accuracy")
)

accuracy_per_group

