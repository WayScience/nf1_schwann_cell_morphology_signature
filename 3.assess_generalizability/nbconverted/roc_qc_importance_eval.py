#!/usr/bin/env python
# coding: utf-8

# ## Perform bootstrapping ROC AUC method to determine if performing QC within our pipelines is important and improves performance 
# 
# In this method, we have trained two models; one with QC'd data and the other without QC'd data (more noise). 
# We apply the models to their respective holdout (`Plate_6`) dataset (e.g., QC'd model applied to QC'd data and no-QC model applied to no QC dataset).
# We filter the Plate 6 dataset to only include the Null and WT cells from the iNFixion cell lines (Null C04 and WT A3), so we are directly comparing the models trained on that specific cell line.
# We use bootstrapping, a method that repeatedly samples the dataset with replacement to create random subsets of the same size, where some cells might be duplicated or excluded, simulating variations in the population.
# We calculate the ROC AUC for each subsample and plot as a histogram.
# 
# Our goal is to evaluate if QC is important enough to perform within our workflows where we see a higher performance in classification than if we performed no QC at all.

# ## Import libraries

# In[1]:


import pathlib
import joblib
from joblib import load
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score
import matplotlib.colors as mcolors
from scipy.stats import ttest_ind

import sys

sys.path.append("../utils")
from roc_eval_utils import bootstrap_roc_auc


# ## Set output directory for ROC AUC figure

# In[2]:


figure_path = pathlib.Path("./figures")
# make directory if it doesn't already exist
figure_path.mkdir(exist_ok=True)


# ## Load in label encoder (no QC and QC encoders have the same mapping so only one is loaded in)

# In[3]:


# load in label encoder
le = load(pathlib.Path("../1.train_models/data/trained_nf1_model_label_encoder.joblib"))

# Print label mapping
label_mapping = {label: le.transform([label])[0] for label in le.classes_}
print(label_mapping)


# ## Extract probabilities from the no QC model applied to the no QC'd holdout plate

# In[ ]:


# Load the trained model
no_QC_model = joblib.load(
    pathlib.Path("../1.train_models/data/trained_nf1_model.joblib")
)

# Load the feature-selected QC plate 4 (ensure it includes both features and labels)
plate_6_no_QC = pd.read_parquet(
    pathlib.Path(
        "/media/18tbdrive/1.Github_Repositories/nf1_schwann_cell_painting_data/3.processing_features/data/single_cell_profiles/Plate_6_sc_normalized.parquet"
    )
)

# Drop rows where Metadata_genotype is HET as the model is not predicting this class
# and would not contribute to the evaluation
plate_6_no_QC = plate_6_no_QC[plate_6_no_QC["Metadata_genotype"] != "HET"]

# Drop rows with any NaNs prior to splitting and getting X and y data
plate_6_no_QC = plate_6_no_QC.dropna()

# Split by institution
plate_MGH_no_QC = plate_6_no_QC[plate_6_no_QC["Metadata_Institution"] == "MGH"]
plate_iNFixion_no_QC = plate_6_no_QC[
    plate_6_no_QC["Metadata_Institution"] == "iNFixion"
]

# Get X and y for MGH
X_MGH_noQC = plate_MGH_no_QC[no_QC_model.feature_names_in_]
y_MGH_noQC = le.transform(plate_MGH_no_QC["Metadata_genotype"])

# Get X and y for iNFixion
X_iNFixion_noQC = plate_iNFixion_no_QC[no_QC_model.feature_names_in_]
y_iNFixion_noQC = le.transform(plate_iNFixion_no_QC["Metadata_genotype"])

# Predict probabilities for the positive class
y_probs_MGH_modelNoQC = no_QC_model.predict_proba(X_MGH_noQC)[:, 1]
y_probs_iNFixion_modelNoQC = no_QC_model.predict_proba(X_iNFixion_noQC)[:, 1]


# ## Extract probabilities from the QC model applied to the QC'd holdout plate

# In[5]:


# Load the trained model
QC_model = joblib.load(
    pathlib.Path("../1.train_models/data/trained_nf1_model_qc.joblib")
)

# Load the normalized QC plate 6(ensure it includes both features and labels)
plate_6_QC = pd.read_parquet(
    pathlib.Path(
        "/media/18tbdrive/1.Github_Repositories/nf1_schwann_cell_painting_data/3.processing_features/data/single_cell_profiles/cleaned_sc_profiles/Plate_6_sc_normalized.parquet"
    )
)

# Drop rows where Metadata_genotype is HET as the model is not predicting this class
# and would not contribute to the evaluation
plate_6_QC = plate_6_QC[plate_6_QC["Metadata_genotype"] != "HET"]

# Drop rows with any NaNs prior to splitting and getting X and y data
plate_6_QC = plate_6_QC.dropna()

# Split by institution
plate_MGH_QC = plate_6_QC[plate_6_QC["Metadata_Institution"] == "MGH"]
plate_iNFixion_QC = plate_6_QC[plate_6_QC["Metadata_Institution"] == "iNFixion"]

# Get X and y for MGH
X_MGH_QC = plate_MGH_QC[QC_model.feature_names_in_]
y_MGH_QC = le.transform(plate_MGH_QC["Metadata_genotype"])

# Get X and y for iNFixion
X_iNFixion_QC = plate_iNFixion_QC[QC_model.feature_names_in_]
y_iNFixion_QC = le.transform(plate_iNFixion_QC["Metadata_genotype"])

# Predict probabilities for the positive class
y_probs_MGH_modelQC = QC_model.predict_proba(X_MGH_QC)[:, 1]
y_probs_iNFixion_modelQC = QC_model.predict_proba(X_iNFixion_QC)[:, 1]


# ## Calculate ROC AUC score from the QC and no QC model and data per cell line

# In[6]:


# Calculate ROC AUC for the iNFixion cell line
aucNoQC_iNFixion = roc_auc_score(y_iNFixion_noQC, y_probs_iNFixion_modelNoQC)
aucQC_iNFixion = roc_auc_score(y_iNFixion_QC, y_probs_iNFixion_modelQC)

print(f"AUC Model 1: {aucNoQC_iNFixion}")
print(f"AUC Model 2: {aucQC_iNFixion}")


# In[7]:


# Calculate ROC AUC for the MGH cell line
aucNoQC_MGH = roc_auc_score(y_MGH_noQC, y_probs_MGH_modelNoQC)
aucQC_MGH = roc_auc_score(y_MGH_QC, y_probs_MGH_modelQC)

print(f"AUC Model 1: {aucNoQC_MGH}")
print(f"AUC Model 2: {aucQC_MGH}")


# ## Perform ROC AUC bootstrapping method for both QC and no QC models and data

# In[ ]:


# No QC model iNFixion cell line
scores_model1_iNFixion = bootstrap_roc_auc(y_iNFixion_noQC, y_probs_iNFixion_modelNoQC)

# QC model iNFixion cell line
scores_model2_iNFixion = bootstrap_roc_auc(y_iNFixion_QC, y_probs_iNFixion_modelQC)

# Compare distributions
t_stat, p_value = ttest_ind(scores_model1_iNFixion, scores_model2_iNFixion)
print(f"T-statistic: {t_stat}, P-value: {p_value}")
print(
    f"Mean ROC AUC for Model No-QC for iNFixion cell line: {np.mean(scores_model1_iNFixion)}"
)
print(
    f"Mean ROC AUC for Model QC for iNFixion cell line: {np.mean(scores_model2_iNFixion)}"
)


# In[9]:


# No QC model for MGH cell line
scores_model1_MGH = bootstrap_roc_auc(y_MGH_noQC, y_probs_MGH_modelNoQC)

# QC model for MGH cell line
scores_model2_MGH = bootstrap_roc_auc(y_MGH_QC, y_probs_MGH_modelQC)

# Compare distributions
t_stat, p_value = ttest_ind(scores_model1_MGH, scores_model2_MGH)
print(f"T-statistic: {t_stat}, P-value: {p_value}")
print(f"Mean ROC AUC for Model No-QC: {np.mean(scores_model1_MGH)}")
print(f"Mean ROC AUC for Model QC: {np.mean(scores_model2_MGH)}")


# ## Generate plot comparing QC and no QC model on the iNFixion cell lines

# In[10]:


# Define darker colors for the mean lines
darker_teal = mcolors.to_rgba("teal", 0.8)
darker_orchid = mcolors.to_rgba("orchid", 0.8)

plt.hist(scores_model1_iNFixion, bins=50, alpha=0.5, label="Model No-QC", color="teal")
plt.hist(scores_model2_iNFixion, bins=50, alpha=0.5, label="Model QC", color="orchid")

# Add vertical lines for the means
plt.axvline(
    np.mean(scores_model1_iNFixion),
    color=darker_teal,
    linestyle="dashed",
    linewidth=2,
)
plt.axvline(
    np.mean(scores_model2_iNFixion),
    color=darker_orchid,
    linestyle="dashed",
    linewidth=2,
)

plt.legend(loc="upper left", fontsize=10)
plt.xlabel("ROC AUC Score", fontsize=12)
plt.ylabel("Frequency", fontsize=12)
plt.title(
    "Bootstrap ROC AUC Distributions For\nQC versus No QC Data",
    fontsize=14,
    fontweight="bold",
)
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

# save figure
plt.savefig(f"{figure_path}/iNFixion_bootstrap_ROC_AUC_QC_versus_no_QC.png", dpi=600)

plt.show()


# ## Generate plot comparing QC and no QC model on the MGH cell lines

# In[11]:


# Define darker colors for the mean lines
darker_teal = mcolors.to_rgba("teal", 0.8)
darker_orchid = mcolors.to_rgba("orchid", 0.8)

plt.hist(scores_model1_MGH, bins=50, alpha=0.5, label="Model No-QC", color="teal")
plt.hist(scores_model2_MGH, bins=50, alpha=0.5, label="Model QC", color="orchid")

# Add vertical lines for the means
plt.axvline(
    np.mean(scores_model1_MGH),
    color=darker_teal,
    linestyle="dashed",
    linewidth=2,
)
plt.axvline(
    np.mean(scores_model2_MGH),
    color=darker_orchid,
    linestyle="dashed",
    linewidth=2,
)

plt.legend(loc="upper left", fontsize=10)
plt.xlabel("ROC AUC Score", fontsize=12)
plt.ylabel("Frequency", fontsize=12)
plt.title(
    "Bootstrap ROC AUC Distributions For\nQC versus No QC Data",
    fontsize=14,
    fontweight="bold",
)
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

# save figure
plt.savefig(f"{figure_path}/MGH_bootstrap_ROC_AUC_QC_versus_no_QC.png", dpi=600)

plt.show()

