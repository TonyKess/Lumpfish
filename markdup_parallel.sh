#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --account=def-ibradbur
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tonykess@gmail.com
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --cpus-per-task=5
#SBATCH --kill-on-invalid=yes
#SBATCH --job-name=markdup_parallel.sh

#ensure you are in your home directory - this script assumes starting in /home/user
source WGSparams_Lump.txt

#export variables to keep parallel happy
export gatk
export projdir

##Go to the bams
cd $projdir/align

#deduplicate in parallel
cat $projdir/sets/$set | \
  parallel --jobs 5 '$gatk   --java-options "-Xmx20G" \
  MarkDuplicates \
  I={}.sorted.bam \
  O={}.deDup.bam M={}_deDupMetrics.txt \
  REMOVE_DUPLICATES=true'

