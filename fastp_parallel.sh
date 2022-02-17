#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --account=def-ibradbur
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tonykess@gmail.com
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --cpus-per-task=32
#SBATCH --kill-on-invalid=yes
#SBATCH --job-name=fastp.sh


source WGSparams_Lump.txt

cd $projdir

#export to make parallel happy
export fastp
export projdir

#run in parallel
cat $projdir/sets/Lumpfish_set12.txt | \
  parallel -j 8 \
  '$fastp -i {}_R1.fastq.gz \
  -I {}_R2.fastq.gz \
  -o $projdir/trim/{}_R1.trimmed.fastq.gz \
  -O $projdir/trim/{}_R2.trimmed.fastq.gz \
  --thread 4 '

