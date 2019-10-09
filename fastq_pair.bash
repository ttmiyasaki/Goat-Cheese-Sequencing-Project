#!/usr/bin/env bash 
##
#SBATCH -p production
#SBATCH --time=2-12:0:0
#SBATCH --mem=10G
#SBATCH -o /share/magalab/Tyler/Cheese/Fastq_Paired/fastq_pair.out
#SBATCH -e /share/magalab/Tyler/Cheese/Fastq_Paired/fastq_pair.err

aklog
time fastq_pair -t 13400000 /share/magalab/Tyler/Raw_Data/MiSeqAnalysis2/MiSeqAnalysis/TTM2018_S1_L001_R1_001.fastq /share/magalab/Tyler/Raw_Data/MiSeqAnalysis2/MiSeqAnalysis/TTM2018_S1_L001_R2_001.fastq
