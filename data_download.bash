#!/usr/bin/env bash 

#SBATCH -J data_download 
#SBATCH -e copy-%j.err.output 
#SBATCH -o copy-%j.output 
#SBATCH -p production 
#SBATCH --mem-per-cpu=20G
#SBATCH -t 01-18:00

wget -r -nH -nc -np -R index.html* "http://slimsdata.genomecenter.ucdavis.edu/Data/vowxhmb0pu/181221_M00384_0143_MS6532755-500V2/MiSeqAnalysis/"
