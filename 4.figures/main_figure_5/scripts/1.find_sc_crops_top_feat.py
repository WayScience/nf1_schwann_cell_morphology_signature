#!/usr/bin/env python
# coding: utf-8

# # Generate min/max representative single-cell images per top two features with highest coefficients  
# 
# This is based on the absolute values of the coefficients.

# ## Import libraries

# In[1]:


import pathlib
from pprint import pprint

import cv2
import pandas as pd
from typing import List, Dict


# ## Define functions

# In[2]:


# Function for formatting min/max row data frames into dictionaries
def create_sc_dict(dfs: List[pd.DataFrame], names: List[str]) -> dict:
    """Format lists of data frames and names into a dictionary with all relevant metadata to find single-cell images.

    Args:
        dfs (List[pd.DataFrame]): List of data frames each containing a single cell and relevant metadata.
        names (List[str]): List of names corresponding to the data frames.

    Returns:
        dict: Dictionary containing info relevant for finding single-cell crops.
    """
    sc_dict = {}
    for df, name in zip(dfs, names):
        for i, (_, row) in enumerate(df.iterrows()):
            key = f"{name}_{i + 1}"
            sc_dict[key] = {
                "plate": row["Metadata_Plate"],
                "well": row["Metadata_Well"],
                "site": row["Metadata_Site"],
                "location_center_x": row["Metadata_Nuclei_Location_Center_X"],
                "location_center_y": row["Metadata_Nuclei_Location_Center_Y"],
            }
    return sc_dict


# In[3]:


# Function for generating and saving single-cell crops per channel as PNGs
def generate_sc_crops(
    sc_dict: Dict,
    channel_mapping: Dict[int, str],
    images_dir: pathlib.Path,
    output_img_dir: pathlib.Path,
    crop_size: int,
) -> None:
    """Using a dictionary with single-cell metadata info per image set, single-cell crops per channel are generated
    and saved as PNGs in an image set folder.

    Args:
        sc_dict (Dict): Dictionary containing info relevant for finding single-cell crops.
        channel_mapping (Dict[int, str]): Dictionary mapping integer to channel name for generating paths.
        images_dir (pathlib.Path): Directory where illumination corrected images are found.
        output_img_dir (pathlib.Path): Main directory to save each image set single-cell crops
        crop_size (int): Size of the box in pixels (example: setting crop_size as 250 will make a 250x250 pixel crop
        around the single-cell center coordinates)
    """
    for key, info in sc_dict.items():
        # Initialize a list to store file paths for every image set
        file_paths = []

        # Create file paths with well, site, and channel
        for i in range(1, 5):  # Update the range to start from 1
            channel = channel_mapping[i]
            filename = f"{images_dir}/{info['well']}_01_{i}_{info['site']}_{channel}_001_illumcorrect.tiff"
            file_paths.append(filename)

            # Read the image
            channel_image = cv2.imread(filename, cv2.IMREAD_UNCHANGED)

            # Use the location_center_x and location_center_y to create a crop
            center_x = info.get("location_center_x")
            center_y = info.get("location_center_y")

            # Crop dimensions (including crop_size)
            half_crop = crop_size // 2

            # Ensure the center coordinates are valid
            if center_x is not None and center_y is not None:
                # Calculate crop boundaries
                top_left_x = max(int(center_x - half_crop), 0)
                top_left_y = max(int(center_y - half_crop), 0)
                bottom_right_x = min(int(center_x + half_crop), channel_image.shape[1])
                bottom_right_y = min(int(center_y + half_crop), channel_image.shape[0])

                # Perform cropping
                cropped_channel = channel_image[
                    top_left_y:bottom_right_y, top_left_x:bottom_right_x
                ]

                # Ensure the cropped image is of size 250x250
                cropped_channel = cv2.resize(cropped_channel, (crop_size, crop_size))

                # Make directory for the key to keep all channels for an image in one folder
                key_dir = pathlib.Path(f"{output_img_dir}/{key}")
                key_dir.mkdir(exist_ok=True, parents=True)

                # Save the cropped image with single_cell and channel information
                output_filename = pathlib.Path(f"{key_dir}/{key}_d{i}_cropped.png")

                # Check if the file already exists
                if not output_filename.exists():
                    cv2.imwrite(str(output_filename), cropped_channel)
                else:
                    print(f"File {output_filename} already exists. Skipping.")


# ## Set paths and variables

# In[4]:


# Path to cell painting data directory
cell_painting_dir = pathlib.Path(
    "/media/18tbdrive/1.Github_Repositories/nf1_schwann_cell_painting_data"
)

# Images directory for plate 5 (using for finding single-cells)
images_dir = pathlib.Path(
    f"{cell_painting_dir}/1.cellprofiler_ic/Corrected_Images/Corrected_Plate_5"
).resolve(strict=True)

# Output dir for cropped images
output_img_dir = pathlib.Path("./sc_crops")
output_img_dir.mkdir(exist_ok=True)

# Define the size of the cropping box (NxN pixels)
crop_size = 300

# Define a mapping for the suffixes
channel_mapping = {1: "DAPI", 2: "GFP", 3: "CY5", 4: "RFP"}

# Create open list for one row data frames for each top feature per channel per cell type
list_of_dfs = []

# Create open list of names to assign each data frame in a list relating to the feature, channel, and cell type
list_of_names = []


# ## Load in Plate 5 data to generate repesentative images from

# In[5]:


# Load in QC normalized + feature selected data as data frame
plate5_df = pd.read_parquet(
    pathlib.Path(
        f"{cell_painting_dir}/3.processing_features/data/single_cell_profiles/cleaned_sc_profiles/Plate_5_sc_feature_selected.parquet"
    )
)

# Load in QC annotated dataframe to extract neighbors
annot_df = pd.read_parquet(
    pathlib.Path(
        f"{cell_painting_dir}/3.processing_features/data/single_cell_profiles/cleaned_sc_profiles/Plate_5_sc_annotated.parquet"
    ),
    columns=[
        "Metadata_Well",
        "Metadata_Site",
        "Metadata_Nuclei_Number_Object_Number",
        "Cells_Neighbors_NumberOfNeighbors_Adjacent",
    ],
)

plate5_df = plate5_df.merge(
    annot_df,
    on=["Metadata_Well", "Metadata_Site", "Metadata_Nuclei_Number_Object_Number"],
    how="inner",
)

plate5_df.rename(
    columns={
        "Cells_Neighbors_NumberOfNeighbors_Adjacent": "Metadata_Number_of_Cells_Neighbors_Adjacent"
    },
    inplace=True,
)

print(plate5_df.shape)
plate5_df.head()


# ## Load in feature importance data and determine the top two highest coefficients 
# 
# We will be creating image montages for the features with the highest coefficients, which both relate to being important for predicting WT genotype (positive values).
# The third highest value (in terms of absolute value) is a feature that is most important for predicting Null genotype, but we do not montage it here.
# 
# **Note:** Top positive feature means the most important in predicting the WT genotype, most negative is most important in predicting Null genotype.

# In[6]:


# Load in feature importances from QC model
feat_import_df = pd.read_parquet(
    pathlib.Path(
        "../../2.evaluate_model/model_evaluation_data/feature_importances_qc.parquet"
    )
)

# Find the top highest coefficient feature (is positive to related to predicting WT)
top_coeff_feature = feat_import_df.sort_values(
    by="feature_importances", ascending=False
).iloc[0]["feature_names"]

# Find the second highest coefficient feature (is positive to related to predicting WT)
second_top_coeff_feature = feat_import_df.sort_values(
    by="feature_importances", ascending=False
).iloc[1]["feature_names"]

# Find the top negative feature (predicting Null) [NOT INCLUDED AS MONTAGE]
top_Null_feature = feat_import_df.loc[
    feat_import_df["feature_importances"].idxmin(), "feature_names"
]

# Print the features
print(top_coeff_feature)
print(second_top_coeff_feature)
print(top_Null_feature)


# ## Filter plate 5 single-cells to only include isolated cells that are not near the edge of the FOV

# In[7]:


# Filter the DataFrame directly
filtered_plate5_df = plate5_df[
    (plate5_df["Metadata_Number_of_Cells_Neighbors_Adjacent"].isin([0]))
    & (plate5_df["Metadata_Nuclei_Location_Center_X"] > crop_size // 2)
    & (
        plate5_df["Metadata_Nuclei_Location_Center_X"]
        < (plate5_df["Metadata_Nuclei_Location_Center_X"].max() - crop_size // 2)
    )
    & (plate5_df["Metadata_Nuclei_Location_Center_Y"] > crop_size // 2)
    & (
        plate5_df["Metadata_Nuclei_Location_Center_Y"]
        < (plate5_df["Metadata_Nuclei_Location_Center_Y"].max() - crop_size // 2)
    )
]

print(filtered_plate5_df.shape)
filtered_plate5_df.head()


# ### Max single-cells for top highest feature

# In[8]:


# Get data frame with the next top 6 single-cells
max_top_feature = (
    filtered_plate5_df[filtered_plate5_df["Metadata_genotype"] == "WT"]
    .nlargest(6, top_coeff_feature)[
        [
            top_coeff_feature,
            "Metadata_genotype",
            "Metadata_Well",
            "Metadata_Plate",
            "Metadata_Site",
            "Metadata_Number_of_Cells_Neighbors_Adjacent",
            "Metadata_Nuclei_Location_Center_X",
            "Metadata_Nuclei_Location_Center_Y",
        ]
    ]
)

# Append the DataFrame and its name to the lists
list_of_dfs.append(max_top_feature)
list_of_names.append("max_top_feature")

print(max_top_feature.shape)
max_top_feature


# ### Min single-cells for top highest feature

# In[9]:


# Get data frame with the top 3 single-cells from the top WT coefficient
min_top_feature = filtered_plate5_df[
    filtered_plate5_df["Metadata_genotype"] == "Null"
].nsmallest(6, top_coeff_feature)[
    [
        top_coeff_feature,
        "Metadata_genotype",
        "Metadata_Well",
        "Metadata_Plate",
        "Metadata_Site",
        "Metadata_Number_of_Cells_Neighbors_Adjacent",
        "Metadata_Nuclei_Location_Center_X",
        "Metadata_Nuclei_Location_Center_Y",
    ]
]

# Append the DataFrame and its name to the lists
list_of_dfs.append(min_top_feature)
list_of_names.append("min_top_feature")

print(min_top_feature.shape)
min_top_feature


# ### Max single-cells for the second highest feature

# In[10]:


# Get data frame with the top 6 single-cells
max_second_top_feature = filtered_plate5_df[
    filtered_plate5_df["Metadata_genotype"] == "WT"
].nlargest(6, second_top_coeff_feature)[
    [
        second_top_coeff_feature,
        "Metadata_genotype",
        "Metadata_Well",
        "Metadata_Plate",
        "Metadata_Site",
        "Metadata_Number_of_Cells_Neighbors_Adjacent",
        "Metadata_Nuclei_Location_Center_X",
        "Metadata_Nuclei_Location_Center_Y",
    ]
]

# Append the DataFrame and its name to the lists
list_of_dfs.append(max_second_top_feature)
list_of_names.append("max_second_top_feature")

print(max_second_top_feature.shape)
max_second_top_feature


# ### Min single-cells for the second highest feature

# In[11]:


# Get data frame with the top 3 single-cells from the second top Null coefficient
min_second_top_feature = filtered_plate5_df[
    filtered_plate5_df["Metadata_genotype"] == "Null"
].nsmallest(6, second_top_coeff_feature)[
    [
        second_top_coeff_feature,
        "Metadata_genotype",
        "Metadata_Well",
        "Metadata_Plate",
        "Metadata_Site",
        "Metadata_Number_of_Cells_Neighbors_Adjacent",
        "Metadata_Nuclei_Location_Center_X",
        "Metadata_Nuclei_Location_Center_Y",
    ]
]

# Append the DataFrame and its name to the lists
list_of_dfs.append(min_second_top_feature)
list_of_names.append("min_second_top_feature")

print(min_second_top_feature.shape)
min_second_top_feature


# ## Merge feature info into dictionary for processing

# In[12]:


sc_dict = create_sc_dict(dfs=list_of_dfs, names=list_of_names)

# Check the created dictionary for the first two items
pprint(list(sc_dict.items())[:2], indent=4)


# ## Generate single-cell crops 

# In[13]:


generate_sc_crops(
    sc_dict=sc_dict,
    channel_mapping=channel_mapping,
    images_dir=images_dir,
    output_img_dir=output_img_dir,
    crop_size=crop_size,
)

