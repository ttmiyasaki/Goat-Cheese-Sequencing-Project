#!/bin/bash
##
#SBATCH -p production
#SBATCH --mem=30G
#SBATCH --time=1-18:0:0
#SBATCH -n 15
#SBATCH -o /share/magalab/Tyler/Cheese/Kneaddata2/kneaddata.out
#SBATCH -e /share/magalab/Tyler/Cheese/Kneaddata2/kneaddata.err
module load bowtie2/2.3.4.1
module load fastqc/0.11.7
module load trimmomatic/0.33
module load trf/4.0.9

#note that to run this there are some issues with current code base and Python3.
#The issue is discussed here:  https://bitbucket.org/biobakery/kneaddata/issues/8/error-during-reformatting-of-sequence
#And a fix can be found here to edit the code base: https://bitbucket.org/biobakery/kneaddata/commits/dbc99e91f332
time kneaddata --bypass-trim --input /share/magalab/Tyler/Cheese/Fastq_Paired/TTM2018_S1_L001_R1_001.fastq.paired.fq --input /share/magalab/Tyler/Cheese/Fastq_Paired/TTM2018_S1_L001_R2_001.fastq.paired.fq -t 15 -db /share/magalab/bin/Fastq_Screen_Index/Cow_Genome_ARS-UCD1.2_Index/ -db /share/magalab/bin/Fastq_Screen_Index/Phage_Index/ -db /share/magalab/bin/Fastq_Screen_Index/Cow_Mito_Index/ -db /share/magalab/bin/Fastq_Screen_Index/Human_Index/ -db /share/magalab/Tyler/Cheese/GoatGenome/GoatDB/ --output /share/magalab/Tyler/Cheese/Kneaddata2 --run-fastqc-start --run-fastqc-end
