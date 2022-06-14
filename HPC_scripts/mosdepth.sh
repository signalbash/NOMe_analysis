#!/usr/bin/bash -l

#PBS -l walltime=100:00:00
#PBS -l mem=16GB
#PBS -l ncpus=2
#PBS -J 1-24

### Initialisation

cd $PBS_O_WORKDIR
# ~/PROJECT/bams/

module load rosalind
module load mosdepth

mkdir -p mosdepth

FILENAME=`sed "${PBS_ARRAY_INDEX}q;d" filenames.txt`

cd mosdepth
mosdepth "${FILENAME}" ../dedup/"${FILENAME}".dedup.bam
