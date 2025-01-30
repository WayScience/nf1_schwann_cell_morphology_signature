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
from typing import Tuple
import seaborn as sns
import matplotlib.pyplot as plt


# ## Define helper functions

# In[2]:


def get_X_y_data(
    df: pd.DataFrame, label: str, shuffle: bool = False
) -> Tuple[pd.DataFrame, np.array]:
    """Get X (feature space) and labels (predicting class) from pandas Data frame, including feature names.

    Args:
        df (pd.DataFrame): Data frame containing morphology.
        label (str): Name of the Metadata column being used as the predicting class.
        shuffle (bool, optional): Shuffle the feature columns to get a shuffled dataset. Defaults to False.

    Returns:
        Tuple[pd.DataFrame, np.array]: Returns DataFrame for the feature space (X) with feature names,
                                       and np.array for the predicting class (y).
    """
    # Get the feature columns (excluding 'Metadata' columns and the label column)
    feature_columns = [
        col for col in df.columns if not col.startswith("Metadata") and col != label
    ]

    # Extract feature space (X) as a DataFrame (keep feature names)
    X = df[feature_columns]

    # Extract class label (y) as a NumPy array
    y = df.loc[:, [label]].values
    y = np.ravel(y)  # Flatten y to a 1D array

    # If shuffle is True, shuffle the feature columns independently for training
    if shuffle:
        for column in X.T:
            np.random.shuffle(column)

    return X, y


# ## Set paths and variables

# In[ ]:


# Path to folder holding model and encoder files
model_dir = pathlib.Path("../1.train_models/data")

# Path to results directory
results_dir = pathlib.Path("./results")
results_dir.mkdir(exist_ok=True)

# Load in the model encoder
le = load(pathlib.Path(f"{model_dir}/trained_nf1_model_label_encoder.joblib"))

# Load in the model
model = load(pathlib.Path(f"{model_dir}/trained_nf1_model.joblib"))

# Set the random seed
rng = np.random.default_rng(0)


# ## Load in plate with two cell lines (Plate 6)

# In[4]:


# Read in data from plate 6 with two cell lines
plate6_df = pd.read_parquet(
    pathlib.Path(
        "/media/18tbdrive/1.Github_Repositories/nf1_schwann_cell_painting_data/3.processing_features/data/single_cell_profiles/Plate_6_sc_normalized.parquet"
    )
)

# Remove rows where Metadata_genotype is "HET"
plate6_df = plate6_df[plate6_df["Metadata_genotype"] != "HET"].reset_index(drop=True)

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

# In[5]:


# Shuffle the features randomly, excluding columns that start with "Metadata_"
shuffled_plate6_df = plate6_df.apply(
    lambda x: rng.permutation(x) if not x.name.startswith("Metadata_") else x
)

# Print shape and head of data
print(shuffled_plate6_df.shape)
shuffled_plate6_df.head()


# ## Apply model to final and shuffled versions of the plate data

# In[6]:


# Create list of the metadata columns only
meta_cols = [col for col in plate6_df.columns if "Metadata" in col]

# Define a dictionary to handle both data types
data_dict = {"final": plate6_df, "shuffled": shuffled_plate6_df}

# Initialize a list to store processed dataframes
processed_dfs = []

# Loop through the data dictionary to create probability dataframes
for data_type, data in data_dict.items():
    # Drop rows with Metadata_genotype == "HET"
    data = data[data["Metadata_genotype"] != "HET"]

    # Ensure no duplicates in data and reset index
    data = data.drop_duplicates().reset_index(drop=True)

    # Predict probabilities and labels
    probabilities = model.predict_proba(data[model.feature_names_in_])[:, 1]
    predicted_genotype = model.predict(
        data[model.feature_names_in_]
    )  # outputs as binary labels

    # Convert true labels to binary values (0 or 1) for WT and Null
    true_genotype = le.transform(data["Metadata_genotype"]).tolist()

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

# Save to Parquet (uncomment when needed)
# combined_df.to_parquet(f"{model_dir}/plate_6_single_cell_probabilities.parquet")

# Print shape and head of data
print(combined_df.shape)
combined_df.head()


# ## Split the probability data by Institution

# In[7]:


# Create dictionary with the split dataframes based on Institution
institution_dfs = {
    institution: combined_df[combined_df["Metadata_Institution"] == institution].copy()
    for institution in combined_df["Metadata_Institution"].unique()
}


# ## Generate PR curve results

# In[8]:


precision_recall_data = []

for institution, df in institution_dfs.items():
    for data_type in ["final", "shuffled"]:  # Compute separately for both types
        subset_df = df[df["data_type"] == data_type]

        # Compute precision-recall curve
        precision, recall, _ = precision_recall_curve(
            subset_df["true_genotype"], subset_df["probability_WT"]
        )
        
        institution_results = pd.DataFrame({
            "Precision": precision[:-1],  
            "Recall": recall[:-1],        
            "Metadata_Institution": institution,
            "data_type": data_type,
        })
        
        precision_recall_data.append(institution_results)

# Combine all institution-based PR data
precision_recall_df = pd.concat(precision_recall_data, ignore_index=True)

# Save PR curve data to parquet file
precision_recall_df.to_parquet(f"{results_dir}/plate6_precision_recall_final_model.parquet")

print(precision_recall_df.shape)
precision_recall_df.head()


# In[9]:


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

# Add labels and title
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.title("Precision vs Recall for Different Institutions and Data Types")
plt.legend()
plt.show()


# ## Generate accuracy scores per institution and data type (final or shuffled)

# In[10]:


# Calculate accuracy per institution and data type (final or shuffled)
accuracy_per_group = combined_df.groupby(
    ["Metadata_Institution", "data_type"]
).apply(lambda x: accuracy_score(x["true_genotype"], x["predicted_genotype"])).reset_index(name="accuracy")

# Save accuracy data to parquet file
accuracy_per_group.to_parquet(f"{results_dir}/plate6_accuracy_final_model.parquet")

accuracy_per_group


# In[11]:


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
    palette="Dark2"
)

# Add labels and title
plt.xlabel("Genotype")
plt.ylabel("Accuracy")
plt.title("Accuracy per Genotype and Institution")
plt.legend(title="Institution")
plt.show()


# In[ ]:




