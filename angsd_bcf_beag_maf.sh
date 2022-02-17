#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --account=def-ibradbur
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tonykess@gmail.com
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --cpus-per-task=32
#SBATCH --kill-on-invalid=yes
#SBATCH --job-name=angsd_bcf_beag_maf.sh

#ensure you are in your home directory - this script assumes starting in /home/user
source WGSparams_Lump.txt

#chrom and #bamlist could be specified in the sbatch command so that they can be looped easily

##Go to the bams
cd $projdir/align
 
#use that new angsd
  $angsd \
  -nThreads 32 \
  -doGlf 2 \
  -dobcf 1 \
  -dopost 1 \
  -gl 1 \
  -domajorminor 1 \
  -doMaf 1 \
  -docounts 1 \
  -minMapQ 30 \
  -minQ 20 \
  -minInd 870 \
  -out $projdir/angsd_out/$species\_sizechecked.$chrom \
  -SNP_pval 2e-6 \
  -minMaf 0.01 \
  -dumpCounts 2 \
  -doQsDist 1 \
  -setMinDepth 1000 \
  -uniqueOnly 1 \
  -remove_bads 1 \
  -r $chrom: \
  -only_proper_pairs 1 \
  -bam $bamlist
