#!/bin/bash
#SBATCH --job-name=gwasbrewer
#SBATCH --mem=16G
#SBATCH --time=01:00:00

# Run the R script
Rscript run_gwasbrewer_input_helpers.R
