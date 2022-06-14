#!/usr/bin/bash -l

#PBS -l walltime=100:00:00
#PBS -l mem=32GB
#PBS -l ncpus=4

### Initialisation

cd $PBS_O_WORKDIR

module load rosalind
module load gcc-env/6.4.0
module load samtools/1.12
module load bwa
source ~/pythonenvs/bwameth_py2/bin/activate

#human index
bwameth.py index hg38.fa
# mouse index with lambda
bwameth.py index GRCm39_pluslambda.fa
