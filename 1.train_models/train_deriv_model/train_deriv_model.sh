#!/bin/bash

# initialize the correct shell for your machine to allow conda to work (see README for note on shell names)
conda init bash
# activate the python based analysis env
conda activate nf1_analysis

# convert all notebooks to script files into the scripts folder
jupyter nbconvert --to script --output-dir=scripts/ *.ipynb

# run the notebook to train a new model and evaluate performance
python scripts/train_log_model.py
