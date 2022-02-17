#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --account=def-ibradbur
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tonykess@gmail.com
#SBATCH --nodes=1
#SBATCH --mem=250G
#SBATCH --cpus-per-task=32
#SBATCH --kill-on-invalid=yes
#SBATCH --job-name=PCANGSD_Lumpfish.sh

#ensure you are in your home directory - this script assumes starting in /home/user

## Load python and launch virtual environment
module load python
source ~/python_env/bin/activate

python /home/tkess/programs/pcangsd/pcangsd.py \
  -beagle /home/tkess/scratch/Lumpy/angsd_out/Cyclopterus_lumpus_sizechecked.beagle.gz \
  -o /home/tkess/scratch/Lumpy/PCANGSD_out/Cyclopterus_lumpus_sizechecked \
  -selection \
  -sites_save \
  -snp_weights \
  -admix \
  -tree \
  -tree_samples /home/tkess/scratch/Lumpy/sets/treelist.txt \
  -inbreedSamples \
  -minMaf 0.01 \
  -threads 32 

deactivate
