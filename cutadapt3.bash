#!/usr/bin/env bash 
##
#SBATCH -p production
#SBATCH --time=2-10:0:0
#SBATCH --mem=10G
#SBATCH -o /share/magalab/Tyler/Cheese/Demux_CutAdapt3/cutadapt_anchored.out
#SBATCH -e /share/magalab/Tyler/Cheese/Demux_CutAdapt3/cutadapt_anchored.err

source activate cutadapt
cd /share/magalab/Tyler/Cheese/Demux_CutAdapt3
ulimit -n 1200
cutadapt -g file:/share/magalab/Tyler/Cheese/cutadapt_formatted_barcode_file.fa -G file:/share/magalab/Tyler/Cheese/cutadapt_formatted_barcode_file.fa -e 0 --no-indels --untrimmed-output=no_barcodes.R1.fastq --untrimmed-paired-output=no_barcodes.R2.fastq -o {name}_R1.fastq -p {name}_R2.fastq /share/magalab/Tyler/Cheese/Kneaddata2/TTM2018_S1_L001_R1_001.fastq.paired_kneaddata_paired_1.fastq /share/magalab/Tyler/Cheese/Kneaddata2/TTM2018_S1_L001_R1_001.fastq.paired_kneaddata_paired_2.fastq
