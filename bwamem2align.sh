#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --account=def-ibradbur
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tonykess@gmail.com
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --cpus-per-task=32
#SBATCH --kill-on-invalid=yes
#SBATCH --job-name=bwamem2align

#ensure you are in your home directory - this script assumes starting in /home/user
#specify set in sbatch command using --export=ALL

#load variables
source WGSparams_Lump.txt

## Load software modules
module load samtools

#change directory and align all reads
cd $projdir/trim

while read ind;
  do echo $ind\.bam ;
  RGID=$(echo $ind |  sed 's/i5.*/i5/') ;
  SMID=$(echo $ind | sed 's/NS.*i5.//') ;
  LBID=$(echo $ind | sed 's/.UDP.*//');
  $bwamem2 mem \
  -t 32 \
  -R "@RG\tID:$RGID\tSM:$SMID\tLB:$LBID" \
  $genome \
  $ind\_R1.trimmed.fastq.gz  $ind\_R2.trimmed.fastq.gz\
  | samtools sort -o $projdir/align/$ind\.sorted.bam -T $ind -@ 32 -m 3G ;
  done <  $projdir/sets/$set 
