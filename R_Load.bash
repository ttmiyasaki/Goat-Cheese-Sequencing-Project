#!/usr/bin/env bash 
##
#SBATCH -p production
#SBATCH --time=0:5:0
#SBATCH --mem=10G
#SBATCH -e /share/magalab/Tyler/Scripts/R.err

module load R/3.5.0
echo "hello world"
R

source("/share/magalab/Tyler/Scripts/dada2load.r")