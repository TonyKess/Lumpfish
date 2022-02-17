#!/bin/bash
#SBATCH --time=0:30:00
#SBATCH --account=def-ibradbur
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tonykess@gmail.com
#SBATCH --nodes=1
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --kill-on-invalid=yes
#SBATCH --job-name=bcf_to_vcf.sh


#ensure you are in your home directory - this script assumes starting in /home/user
## Load software modules
module load bcftools

source WGSparams_Lump.txt

##Go to the bams
cd $projdir/angsd_out

bcftools view $species\_sizechecked.$chrom.bcf > $species\_sizechecked.$chrom.vcf
