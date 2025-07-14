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

# Images directory for plate 6 (using for finding single-cells)
images_dir = pathlib.Path(
    f"{cell_painting_dir}/1.cellprofiler_ic/Corrected_Images/Corrected_Plate_6"
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


# ## Load in Plate 6 data to generate repesentative images from both derivatives

# In[5]:


# Load in QC normalized + feature selected data as data frame
plate_df = pd.read_parquet(
    pathlib.Path(
        f"{cell_painting_dir}/3.processing_features/data/single_cell_profiles/cleaned_sc_profiles/Plate_6_sc_feature_selected.parquet"
    )
)

# Load in QC annotated dataframe to extract neighbors
annot_df = pd.read_parquet(
    pathlib.Path(
        f"{cell_painting_dir}/3.processing_features/data/single_cell_profiles/cleaned_sc_profiles/Plate_6_sc_annotated.parquet"
    ),
    columns=[
        "Metadata_Well",
        "Metadata_Site",
        "Metadata_Nuclei_Number_Object_Number",
        "Cells_Neighbors_NumberOfNeighbors_Adjacent",
    ],
)

plate_df = plate_df.merge(
    annot_df,
    on=["Metadata_Well", "Metadata_Site", "Metadata_Nuclei_Number_Object_Number"],
    how="inner",
)

plate_df.rename(
    columns={
        "Cells_Neighbors_NumberOfNeighbors_Adjacent": "Metadata_Number_of_Cells_Neighbors_Adjacent"
    },
    inplace=True,
)

# Drop HET cells
plate_df = plate_df[plate_df["Metadata_genotype"] != "HET"]

print(plate_df.shape)
plate_df.head()


# ## Load in feature importance data and determine the top two highest coefficients 
# 
# We will be creating image montages for the features with the highest coefficients (after absolute value). Below will show the top two ranked features with their sign.
# 
# **Note:** Top positive feature means the most important in predicting the WT genotype, most negative is most important in predicting Null genotype.

# In[6]:


# Load in feature importances from QC model
feat_import_df = pd.read_csv(
    pathlib.Path("../supp_figure_7/coeff_results/final_model_coefficients.csv")
)

# Sort by absolute value of coefficient, descending
sorted_abs_feat_import_df = feat_import_df.reindex(
    feat_import_df["coefficient"].abs().sort_values(ascending=False).index
)

# Find the top two features by absolute coefficient value (keep sign)
top_coeff_feature = sorted_abs_feat_import_df.iloc[0]["feature"]
top_coeff_value = sorted_abs_feat_import_df.iloc[0]["coefficient"]

second_top_coeff_feature = sorted_abs_feat_import_df.iloc[1]["feature"]
second_top_coeff_value = sorted_abs_feat_import_df.iloc[1]["coefficient"]

# Print the features and their signed values
print(f"{top_coeff_feature}: {top_coeff_value}")
print(f"{second_top_coeff_feature}: {second_top_coeff_value}")


# ## Filter single-cells to only include isolated cells that are not near the edge of the FOV

# In[7]:


# Filter the DataFrame directly
filtered_plate_df = plate_df[
    (plate_df["Metadata_Number_of_Cells_Neighbors_Adjacent"].isin([0]))
    & (plate_df["Metadata_Nuclei_Location_Center_X"] > crop_size // 2)
    & (
        plate_df["Metadata_Nuclei_Location_Center_X"]
        < (plate_df["Metadata_Nuclei_Location_Center_X"].max() - crop_size // 2)
    )
    & (plate_df["Metadata_Nuclei_Location_Center_Y"] > crop_size // 2)
    & (
        plate_df["Metadata_Nuclei_Location_Center_Y"]
        < (plate_df["Metadata_Nuclei_Location_Center_Y"].max() - crop_size // 2)
    )
]

print(filtered_plate_df.shape)
filtered_plate_df.head()


# ### Max single-cells for top highest feature (original)

# In[8]:


## Get data frame with the top 3 single-cells for WT genotype from iNFixion institution
max_top_feature_orig = (
    filtered_plate_df[
        (filtered_plate_df["Metadata_genotype"] == "WT")
        & (filtered_plate_df["Metadata_Institution"] == "iNFixion")
    ]
    .sort_values(by=top_coeff_feature, ascending=False)
    .iloc[1:4][
        [
            top_coeff_feature,
            "Metadata_genotype",
            "Metadata_Institution",
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
list_of_dfs.append(max_top_feature_orig)
list_of_names.append("max_top_feature_orig")

print(max_top_feature_orig.shape)
max_top_feature_orig


# ### Max single-cells for top highest feature (derivative)

# In[9]:


## Get data frame with the top 3 single-cells for WT genotype from MGH institution
max_top_feature_deriv = filtered_plate_df[
    (filtered_plate_df["Metadata_genotype"] == "WT")
    & (filtered_plate_df["Metadata_Institution"] == "MGH")
].nlargest(3, top_coeff_feature)[
    [
        top_coeff_feature,
        "Metadata_genotype",
        "Metadata_Institution",
        "Metadata_Well",
        "Metadata_Plate",
        "Metadata_Site",
        "Metadata_Number_of_Cells_Neighbors_Adjacent",
        "Metadata_Nuclei_Location_Center_X",
        "Metadata_Nuclei_Location_Center_Y",
    ]
]

# Append the DataFrame and its name to the lists
list_of_dfs.append(max_top_feature_deriv)
list_of_names.append("max_top_feature_deriv")

print(max_top_feature_deriv.shape)
max_top_feature_deriv


# ### Min single-cells for top highest feature (original)

# In[10]:


## Get data frame with the bottom 3 single-cells for Null genotype from iNFixion institution
min_top_feature_orig = filtered_plate_df[
    (filtered_plate_df["Metadata_genotype"] == "Null")
    & (filtered_plate_df["Metadata_Institution"] == "iNFixion")
].nsmallest(3, top_coeff_feature)[
    [
        top_coeff_feature,
        "Metadata_genotype",
        "Metadata_Institution",
        "Metadata_Well",
        "Metadata_Plate",
        "Metadata_Site",
        "Metadata_Number_of_Cells_Neighbors_Adjacent",
        "Metadata_Nuclei_Location_Center_X",
        "Metadata_Nuclei_Location_Center_Y",
    ]
]

# Append the DataFrame and its name to the lists
list_of_dfs.append(min_top_feature_orig)
list_of_names.append("min_top_feature_orig")

print(min_top_feature_orig.shape)
min_top_feature_orig


# ### Min single-cells for top highest feature (derivative)

# In[11]:


## Get data frame with the bottom 3 single-cells for Null genotype from MGH institution
min_top_feature_deriv = filtered_plate_df[
    (filtered_plate_df["Metadata_genotype"] == "Null")
    & (filtered_plate_df["Metadata_Institution"] == "MGH")
].nsmallest(3, top_coeff_feature)[
    [
        top_coeff_feature,
        "Metadata_genotype",
        "Metadata_Institution",
        "Metadata_Well",
        "Metadata_Plate",
        "Metadata_Site",
        "Metadata_Number_of_Cells_Neighbors_Adjacent",
        "Metadata_Nuclei_Location_Center_X",
        "Metadata_Nuclei_Location_Center_Y",
    ]
]

# Append the DataFrame and its name to the lists
list_of_dfs.append(min_top_feature_deriv)
list_of_names.append("min_top_feature_deriv")

print(min_top_feature_deriv.shape)
min_top_feature_deriv


# ### Max single-cells for the second highest feature (original)

# In[12]:


# Get data frame with the top 3 single-cells for Null genotype from iNFixion institution
max_second_top_feature_orig = filtered_plate_df[
    (filtered_plate_df["Metadata_genotype"] == "Null")
    & (filtered_plate_df["Metadata_Institution"] == "iNFixion")
].nlargest(3, second_top_coeff_feature)[
    [
        second_top_coeff_feature,
        "Metadata_genotype",
        "Metadata_Institution",
        "Metadata_Well",
        "Metadata_Plate",
        "Metadata_Site",
        "Metadata_Number_of_Cells_Neighbors_Adjacent",
        "Metadata_Nuclei_Location_Center_X",
        "Metadata_Nuclei_Location_Center_Y",
    ]
]

# Append the DataFrame and its name to the lists
list_of_dfs.append(max_second_top_feature_orig)
list_of_names.append("max_second_top_feature_orig")

print(max_second_top_feature_orig.shape)
max_second_top_feature_orig


# ### Max single-cells for the second highest feature (derivative)

# In[13]:


# Get data frame with the top 3 single-cells for Null genotype from MGH institution
max_second_top_feature_deriv = filtered_plate_df[
    (filtered_plate_df["Metadata_genotype"] == "Null")
    & (filtered_plate_df["Metadata_Institution"] == "MGH")
].nlargest(3, second_top_coeff_feature)[
    [
        second_top_coeff_feature,
        "Metadata_genotype",
        "Metadata_Institution",
        "Metadata_Well",
        "Metadata_Plate",
        "Metadata_Site",
        "Metadata_Number_of_Cells_Neighbors_Adjacent",
        "Metadata_Nuclei_Location_Center_X",
        "Metadata_Nuclei_Location_Center_Y",
    ]
]

# Append the DataFrame and its name to the lists
list_of_dfs.append(max_second_top_feature_deriv)
list_of_names.append("max_second_top_feature_deriv")

print(max_second_top_feature_deriv.shape)
max_second_top_feature_deriv


# ### Min single-cells for the second highest feature (original)

# In[14]:


# Get data frame with the bottom 3 single-cells for WT genotype from iNFixion institution
min_second_top_feature_orig = filtered_plate_df[
    (filtered_plate_df["Metadata_genotype"] == "WT")
    & (filtered_plate_df["Metadata_Institution"] == "iNFixion")
].nsmallest(3, second_top_coeff_feature)[
    [
        second_top_coeff_feature,
        "Metadata_genotype",
        "Metadata_Institution",
        "Metadata_Well",
        "Metadata_Plate",
        "Metadata_Site",
        "Metadata_Number_of_Cells_Neighbors_Adjacent",
        "Metadata_Nuclei_Location_Center_X",
        "Metadata_Nuclei_Location_Center_Y",
    ]
]

# Append the DataFrame and its name to the lists
list_of_dfs.append(min_second_top_feature_orig)
list_of_names.append("min_second_top_feature_orig")

print(min_second_top_feature_orig.shape)
min_second_top_feature_orig


# ### Min single-cells for the second highest feature (derivative)

# In[15]:


# Get data frame with the bottom 3 single-cells for WT genotype from MGH institution
min_second_top_feature_deriv = filtered_plate_df[
    (filtered_plate_df["Metadata_genotype"] == "WT")
    & (filtered_plate_df["Metadata_Institution"] == "MGH")
].nsmallest(3, second_top_coeff_feature)[
    [
        second_top_coeff_feature,
        "Metadata_genotype",
        "Metadata_Institution",
        "Metadata_Well",
        "Metadata_Plate",
        "Metadata_Site",
        "Metadata_Number_of_Cells_Neighbors_Adjacent",
        "Metadata_Nuclei_Location_Center_X",
        "Metadata_Nuclei_Location_Center_Y",
    ]
]

# Append the DataFrame and its name to the lists
list_of_dfs.append(min_second_top_feature_deriv)
list_of_names.append("min_second_top_feature_deriv")

print(min_second_top_feature_deriv.shape)
min_second_top_feature_deriv


# ## Merge feature info into dictionary for processing

# In[16]:


sc_dict = create_sc_dict(dfs=list_of_dfs, names=list_of_names)

# Check the created dictionary for the first two items
pprint(list(sc_dict.items())[:2], indent=4)


# ## Generate single-cell crops 

# In[17]:


generate_sc_crops(
    sc_dict=sc_dict,
    channel_mapping=channel_mapping,
    images_dir=images_dir,
    output_img_dir=output_img_dir,
    crop_size=crop_size,
)

