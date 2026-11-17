#!/bin/bash
#SBATCH --job-name=gfa_susie
#SBATCH --mem=64G
#SBATCH --time=04:00:00

# Run the R script
Rscript test_gfa_susie.R
