[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13694800.svg)](https://doi.org/10.5281/zenodo.13694800)

# NF1 Schwann cell morphology signature

We aim to characterize cell morphology signatures of neurofibromin in Schwann cells.
We applied a modified Cell Painting assay on two isogenic Schwann cell lines, one of the wildtype genotype (*NF1+/+*) and one of the null genotype (*NF1-/-*), both from the same patient.
The modified assay stains for four organelles: nuclei, endoplasmic reticulum, mitochondria, and actin.
We applied CellProfiler pipelines to perform quality control, illumination correction, segmentation, and feature extraction.

![Figure 1](./4.figures/figures/main_figure_1_workflow.png)
> *Overview of our Cell Painting and data analysis workflow.* 
> (A) Example image montages of the Cell Painting channels and composite image (all channels overlayed) for each NF1 genotype. 
> The scale bar represents 25 μM. 
> (B) Image analysis and data processing pipeline to derive morphology signature of NF1 genotype.

We segmented 20,680 wild-type (WT) and null cells across three plates and utilized 894 significant morphology features representing various organelle shapes and intensity patterns.

We trained a logistic regression binary classifier to predict the NF1 genotype of single cells.
The model shows high performance with with accuracy of 0.85 and 0.79 for the training and testing data splits respectively.

![Figure 3](./4.figures/figures/main_figure_3_model_eval.png)
> *Logistic regression model predicts NF1 genotype in ipn02.3 2λ Schwann cells with high performance.*
> (A) Precision-recall curves comparing the final model applied to shuffled (dashed line) and non-shuffled data (solid line). 
> Applying the model to the shuffled data performed worse than the non-shuffled data, demonstrating a strong biological signal between NF1 genotypes.
> (B) Confusion matrices from the training and testing sets show higher performance across genotypes in non-shuffled data compared to shuffled data.
> (C) Accuracies demonstrate high performance in predicting both genotypes from the training and testing sets compared to shuffled data.

The machine learning model learned from a total of 894 morphology features and assigned weights, or importance scores, per feature.
This combination of morphology features represents a high-dimensional signature of how the NF1 genotype influences cell morphology in otherwise isogenic Schwann cells.

![Figure 4](./4.figures/figures/main_figure_4_feature_importance.png)
> *Significant cell morphology features come from multiple kinds of measurement and across organelles.*
> Top absolute value coefficients demonstrate which measurements/organelles/compartments are important in determining the differences between NF1 genotypes.
> Not one feature is important in this prediction; many features in our models are important for making accurate genotype predictions.
> The red boxes indicate the two most important features based on absolute value. 
> Before absolute value, the model coefficients range from positive to negative, where positive indicates more importance in predicting the NF1 WT A03 (NF1+/+) genotype and vice versa.
> Of the two top features, both are highest in the positive range before absolute value. 
> The highest coefficient is the distribution of the nuclei stain on the edge of nuclei, and the second highest coefficient is the distribution of the mitochondria stain on top of nuclei.

We look to improve upon this preliminary model in the future.
We aim to generate further data which includes the heterozygous genotype (*NF1+/-*).
AS well, we plan to apply an improved model to large-scale drug screens to capture candidate drugs that make NF1 patient cells look healthy.

**NOTE:** All image analysis and image-based profiling pipelines can be found in the [nf1_schwann_cell_painting_data](https://github.com/WayScience/nf1_schwann_cell_painting_data) repository.

## Repository Structure

This analysis is categorized as follows:

| Module                                                    | Purpose                             | Description                                                                                                                                                                            |
| :-------------------------------------------------------- | :---------------------------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [0.data_analysis](./0.data_analysis/)                     | Preliminary analysis of NF1 data    | We find interesting patterns in the NF1 cell morphology data.                                                                                                                          |
| [1.train_models](./1.train_models/)                       | Train NF1 models                    | Optimize the NF1 model by training multiple models with a random search.                                                                                                               |
| [2.evaluate_models](./2.evaluate_models/)                 | Evaluate final NF1 model            | After training the final NF1 model, we extract the model performance and model feature importances.                                                                                    |
| [3.assess_generalizability](./3.assess_generalizability/) | Assess performance on holdout plate | We evaluate the performance of the model when applied to a holdout plate (data the model has never seen before), which include cell lines from the model and newly derived cell lines. |
| [4.figures](./4.figures/)                                 | Generate manuscript figures         | We interpret the results from the model by generating figures to include in the manuscript.                                                                                            |

## Steps to reproduce

### Step 1. Clone the analysis repo

```sh
git clone https://github.com/WayScience/nf1_schwann_cell_morphology_signature.git
```

### Step 2. Change the current path to the repo path

### Step 3. Populate the nf1_data_repo folder from the nf1_schwann_cell_painting_data repo

```sh
git submodule update --init --recursive
```

### Step 4. Create the conda environments

There are two different environments used in this analysis:

1. [environment.yml](./environment.yml): A Python-based environment meant for performing model training, evaluation, and image montage generation.
2. [figure_environment.yml](./figure_environment.yml): An R-based environment meant for generating figures.

```sh
conda env create -f environment.yml
conda env create -f figure_environment.yml
```
