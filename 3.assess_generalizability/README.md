# Assess generalizability of the model
In this module, the model is evaluated on new "holdout" plate morphology data with new cell lines not in the training or validation datasets as well as the same cell lines found in these datasets.
The "holdout" plate dataset includes cell lines with both the wild-type (WT) and Null genotypes present in all datasets.

In the previous modules, we trained two different model types, one with quality controlled profiles and another without.
For the manuscript, we use the model trained on quality control data, but we include the non-QC'd model for one of the assessments below.

Assessments include:

1. Generating precision-recall curves of the final and shuffled models applied to each cell line to evaluate model performance.
2. KS-statistic test to determine how different the single-cell morphology features are between the two cell lines.
3. Calculating area under the curve of the receiver operating characteristic curve using a bootstrapping method, which we apply to both model types (QC and no-QC). This calculation evaluates the importance of QC and generalization of each model type.
