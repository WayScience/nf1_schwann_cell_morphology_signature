# Creating supplemental figure 8 - Third top feature image montage

To generate the eighth supplemental figure of the manuscript, there are 4 steps to follow:

1. [1.find_sc_crops_top_feat.ipynb](./1.find_sc_crops_top_feat.ipynb): Find the 6 top maximum and minimum representative single cells for the third highest coefficient feature, which represents the top feature in predicting Null genotype.
2. We manually stack the channels together into one composite image where blue is nuclei, red is actin, green is ER, and magenta is mitochondria. Then, we add 25 uM scales to each crop using 3.1065 uM/pixel in the `Analyze > Set Scale... module` (as identified from the metadata of the raw images). The composite images are saved as PNGs back into the same folder.
3. [2.generate_image_montage.ipynb](./2.generate_image_montage.ipynb): Using the composite single cell crops, we can now merge them together to make image montage figures that combines the single-cell crops per feature.
4. [3.main_figure_8.ipynb](./3.supp_figure_8.ipynb): Patch together the image montages into one main figure.

All steps can be ran with the bash script using the command below:

```bash
source supp_figure_8.sh
```
