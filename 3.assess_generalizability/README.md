# Assess generalizability of the model
In this module, the model is applied to a new "holdout" plate (another plate of data that the model has never seen), which includes wild-type (WT) and Null genotype cell lines from the original plates the model was trained on and new WT and Null genotype cells lines.

Assessments include:

1. Generating precision-recall curves of the final and shuffled models applied to each cell line to evaluate performance.
2. KS-statistic test to determine how different the features are between the two cell lines.
3. Area under the curve of the receiver operating characteristic curve for a model trained on quality controlled data and not to assess importance and generalization of each model type.
