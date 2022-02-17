#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --account=def-ibradbur
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tonykess@gmail.com
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --kill-on-invalid=yes
#SBATCH --job-name=genomeindexfaidx

#ensure you are in your home directory - this script assumes starting in /home/user
#load variables
source WGSparams_Lump.txt

# Load software modules
module load samtools

#change directory and align all reads
cd $projdir/genome

$bwamem2 index $genome
samtools faidx $genome


