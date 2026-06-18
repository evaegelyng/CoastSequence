#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 500G
#SBATCH -c 1
#SBATCH -t 4:00:00

Rscript /HPC/bioinformatics/16s/both_silva/scripts/clean_up_ASV_wise_EES.R