#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=def-ibradbur
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tony.kess@gmail.com
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --kill-on-invalid=yes
#SBATCH --job-name=samtools_index_afterdedup.sh

## Load software modules
module load samtools

#source 
source WGSparams_Lump.txt

cd $projdir/align
 
while read  ind; 
  do samtools index $ind.deDup.bam ; done < $projdir/sets/$set
