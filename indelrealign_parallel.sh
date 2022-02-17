#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --account=def-ibradbur
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tonykess@gmail.com
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --cpus-per-task=32
#SBATCH --kill-on-invalid=yes
#SBATCH --job-name=indelrealign_parallel.sh

#ensure you are in your home directory - this script assumes starting in /home/user
source WGSparams_Lump.txt

#load packages
module load nixpkgs/16.09
module load gatk/3.7

#export variables to keep parallel happy
export projdir
export genome
export set

##Go to the bams
cd $projdir/align

#get to realigning
cat $projdir/sets/$set | \
   parallel --jobs 32 \
 'java -Xmx4G -jar $EBROOTGATK/GenomeAnalysisTK.jar \
 -T IndelRealigner \
 -R $genome \
 -I {}.deDup.bam  \
 -targetIntervals {}.intervals \
 -o {}.realigned.bam '
