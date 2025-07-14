# New model training with both cell lines

Given the poor performance of the original model on the new derivative cell lines, we decided to train a new model that includes both cell lines.
We include Plate 6 within the model, and extract PR curve and coefficients to evaluate the performance and what or if any changes occurred in the top features.

All steps can be ran with the bash script using the command below:

```bash
source train_deriv_model.sh
```
