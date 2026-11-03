#!/bin/bash
#SBATCH --job-name=gfa_susie
#SBATCH --mem=32G
#SBATCH --time=04:00:00

# Run the R script
Rscript test2_gfa_susie.R
