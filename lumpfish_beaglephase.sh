#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=def-ibradbur
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tonykess@gmail.com
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --cpus-per-task=32
#SBATCH --kill-on-invalid=yes
#SBATCH --job-name=lumpfish_beaglephase

cd /home/tkess/scratch/Lumpy/

java -Xmx80g -jar /home/tkess/programs/beagle.r1399.jar \
  gl=/home/tkess/scratch/Lumpy/angsd_out/Cyclopterus_lumpus_sizechecked.$chrom.vcf\
  nthreads=32 out=phased/Cyclopterus_lumpus_sizechecked.$chrom

