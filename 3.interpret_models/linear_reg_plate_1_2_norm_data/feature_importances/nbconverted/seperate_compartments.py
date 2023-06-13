#!/usr/bin/env python
# coding: utf-8

# # Here the features are seperated according to compartments for later analysis

# In[1]:


#%%--%%| <qGnYViiwRD|SyZ3qa8iz3>
r"""°°°
## Imports
°°°"""


# In[ ]:


import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.preprocessing import LabelEncoder
from sklearn.linear_model import LogisticRegression

from sklearn.metrics import (
    confusion_matrix,
    ConfusionMatrixDisplay,
    precision_score,
    accuracy_score,
)
from sklearn.model_selection import train_test_split

import matplotlib.pyplot as plt
import seaborn as sns
import itertools

from joblib import dump, load


# ## Find the git root Directory

# In[ ]:


# Get the current working directory
cwd = Path.cwd()

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


# ## Create the output path if it doesn't exist

# In[ ]:


output_path = Path("data")

output_path.mkdir(
    parents=True, exist_ok=True
)  # Create the parent directories if they don't exist


# ## Import the model data as a dataframe

# In[ ]:


feature_properties = pd.read_csv(
    root_dir
    / "1.train_models/linear_reg_plate_1_2_norm_data/data"
    / "model_properties.tsv",
    sep="\t",
)


# ## Seperate cell data by channel

# In[ ]:


# A map for comparments to channel
compartment2channel = {"actin": "RFP", "er": "GFP", "nucleus": "DAPI"}

# Create a Dictionary to hold the comparment data as dataframes
compartment_data = {
    compartment: feature_properties[feature_properties["feature"].str.contains(channel)]
    for compartment, channel in compartment2channel.items()
}


# ## Find the features that use more than one compartment

# In[ ]:


# Create a list of all possible comparment pairs
pairs = list(itertools.combinations(compartment_data.keys(), 2))

# Use a list of possible compartments for find compartment agnostic features
pos_compartments = list(compartment2channel.values())

# Get the features that do not belong to any compartment specifically
other_compartment = feature_properties[
    ~feature_properties["feature"].str.contains("|".join(pos_compartments))
]

# Find intersected rows between each pair of compartments
for pair in pairs:

    # Create placeholder dataframes for each comparment in the pair
    df1 = compartment_data[pair[0]]
    df2 = compartment_data[pair[1]]

    # Get duplicate features between the two compartment dataframes
    intersection = pd.concat([df1, df2], axis=0)
    intersection = intersection[intersection.duplicated(subset="feature", keep=False)]

    # Remove any duplicate features already added to the other compartment dataframe
    other_compartment = pd.concat([other_compartment, intersection], axis=0)
    other_compartment = other_compartment.drop_duplicates(subset="feature")


# ## Organize the data according to compartment in one dataframe

# In[ ]:


# Remove features from each compartment dataframe if they are duplicates in other dataframes, or if the features only exist in the other_compartment dataframe
compartment_data = {
    compartment: compartment_data[compartment][
        ~compartment_data[compartment]["feature"].isin(other_compartment["feature"])
    ]
    for compartment, channel in compartment2channel.items()
}

compartment_data["other"] = other_compartment

# Create a compartment column for each comparment dataframe
for compartment, df in compartment_data.items():
    compartment_data[compartment]["compartment"] = len(df) * [compartment]

# Concatenate the rows of DataFrames to create the plot below
concatenated_df = pd.concat(compartment_data.values(), ignore_index=True)


# ## Save the dataframes

# In[ ]:


concatenated_df.to_csv(output_path / "feature_compartments.tsv", sep="\t", index=False)

# Create a dataframe with only significant models
concatenated_df.loc[concatenated_df["p_value"] >= 0.05].to_csv(
    output_path / "significant_feature_compartments.tsv", sep="\t", index=False
)
