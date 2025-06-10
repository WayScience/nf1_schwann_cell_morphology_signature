# Generate manuscript figures

After evaluation results are extracted, we generate figures describing the results of our experiment.
All figure PNGs are found in the [figures](./figures/) folder.

1. [Main figure 1](./main_figure_1/): This figure describes our workflow and displays an image montage of the wildtype and null *NF1* genotype single cells, which are hard to distinguish just by eye.
2. [Main figure 2](./main_figure_2/): This figure shows how subtle the morphological differences are between *NF1* genotypes at both the well-population and single-cell levels, which supports are reasoning to pursue a machine learning methodology.
3. [Main figure 3](./main_figure_3/): This figure shows the results of the model evaluations (precision-recall, accuracy, and confusion matrices) as extracted in the second module of this repository.
4. [Main figure 4](./main_figure_4/): This figure looks at the feature importances of the model when predicting *NF1* genotype.
5. [Main figure 5](./main_figure_5/): This figure shows four image montages that show six examples of single-cells, two for each of the top features for predicting each genotype.
6. [Main figure 6](./main_figure_6/): This figure shows shows the results from applying the model to another cell line (WT and Null) and the performance.
7. [Supplemental figure 3](./supp_figure_3/): This figure is an extension of main figure 2, which facets the plot by plate to show that the subtle differences between *NF1* genotype are consistent.
8. [Supplemental figure 4](./supp_figure_4/): This figure shows the distributions of FOVs across blur (PowerLogLogSlope) and saturation (PercentMaximal) metrics and where the thresholds were assigned to detect poor-quality images.
9. [Supplemental figure 5](./supp_figure_5/): Extension of main figure 6 panel E, where the importance scores are split based on genotype comparison (e.g., WT cells in original derivative versus WT cells in new derivative, etc.).
10. [Supplemental figure 6](./supp_figure_6/): Visualization of single-cell quality control metrics used to filter out poor quality nuclei segmentations, blurry nuclei, and nuclei/cells undergoing mitosis.
11. [Supplemental figure 7](./supp_figure_7/): PR curves with different conditions visualizing the performance of a new logistic regression model trained on model cell lines.
12. [Supplemental figure 8](./supp_figure_8/): Image montage of the third high coefficient/most important feature to the model (original) which also reflects the feature highest in Null cells.
