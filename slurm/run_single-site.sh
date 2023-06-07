#!/bin/bash
#SBATCH -o slurm/out/single_site.out
#SBATCH -e slurm/err/single_site.err
#SBATCH -p common,statdept
#SBATCH --mem=4G # GB RAM
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.morsomme@duke.edu
hostname # print hostname
module load R
Export R_LIBS_USER = ~/R/x86_64-pc-linux-gnu-library/4.1
R CMD BATCH --no-save R/analysis-single_site.R