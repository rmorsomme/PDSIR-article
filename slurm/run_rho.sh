#!/bin/bash
#SBATCH -o slurm/out/rho_%a.out
#SBATCH -e slurm/err/rho_%a.err
#SBATCH -p common,statdept
#SBATCH --mem=1G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.morsomme@duke.edu
#SBATCH -a 1-1000
hostname # print hostname
module load R
Export R_LIBS_USER = ~/R/x86_64-pc-linux-gnu-library/4.1
R CMD BATCH --no-save R/analysis-rho.R
