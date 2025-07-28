# Generate manuscript figures

After evaluation results are extracted, we generate figures describing the results of our experiment.
All figure PNGs are found in the [figures](./figures/) folder.

1. [Main figure 1](./main_figure_1/): This figure describes our workflow and displays an image montage of the wildtype and null *NF1* genotype single cells, which are hard to distinguish just by eye.
2. [Main figure 2](./main_figure_2/): This figure shows how subtle the morphological differences are between *NF1* genotypes at both the well-population and single-cell levels, which supports are reasoning to pursue a machine learning methodology.
3. [Main figure 3](./main_figure_3/): This figure shows the results of the model evaluations (precision-recall, accuracy, and confusion matrices) from the original model and evaluated the performance on a new model trained on both cell line derivatives.
4. [Main figure 4](./main_figure_4/): This figure looks at the feature importances of the model when predicting *NF1* genotype in the new model.
5. [Main figure 5](./main_figure_5/): This figure shows four image montages that show six examples of single-cells, two for each of the top features for predicting each genotype in the new model.
6. Supplemental figure 1: Western blots from each cell line used in the experiment demonstrating expected neurofibromin levels. (not generated with code)
7. [Supplemental figure 2](./supp_figure_2/): This figure shows the platemap layouts of the first three plates in this experiment (Plates A-C).
8. Supplemental figure 3: Visualization of single-cell quality control metrics used to filter out poor quality nuclei segmentations, blurry nuclei, and nuclei/cells undergoing mitosis. (not generated with code)
9. [Supplemental figure 4](./supp_figure_4/): This figure is an extension of main figure 2, which facets the plot by plate to show that the subtle differences between NF1 genotype are consistent.
10. [Supplemental figure 5](./supp_figure_5/): This figure shows shows the results from applying the original model to the other cell line (WT and Null) and the performance.
11. [Supplemental figure 6](./supp_figure_6/): This figgre visualizes of the ks-stat scores that are split based on genotype comparison (e.g., WT cells in original derivative versus WT cells in new derivative, etc.).
12. [Supplemental figure 7](./supp_figure_7/): This figure shows the distributions of FOVs across blur (PowerLogLogSlope) and saturation (PercentMaximal) metrics and where the thresholds were assigned to detect poor-quality images.
