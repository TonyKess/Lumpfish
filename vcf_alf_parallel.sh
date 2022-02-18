#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --account=def-ibradbur
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tonykess@gmail.com
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --cpus-per-task=32
#SBATCH --kill-on-invalid=yes
#SBATCH --job-name=vcf_alf_parallel.sh

#ensure you are in your home directory - this script assumes starting in /home/user
source WGSparams_Lump.txt

##Go to the vcf
cd $projdir/phased

module load vcftools

cat ../sets/Lumpfishpoplist.txt | parallel --jobs 32 '
 vcftools --gzvcf Cyclopterus_lumpus_sizechecked_phased_chromchange.vcf.gz \
 --keep ../sets/Groupfiles/Lumpfish_inds_{}.tsv \
 --freq \
 --out Cyclopterus_lumpus_sizechecked_phased_chromchange_{} '
