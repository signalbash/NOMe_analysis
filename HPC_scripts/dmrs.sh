#!/usr/bin/bash -l

#PBS -l walltime=100:00:00
#PBS -l mem=32GB
#PBS -l ncpus=4

module purge
module load rosalind
module load R/4.1.0-foss-2021a

cd $PBS_O_WORKDIR

Rscript ./dmrs.R \
  --HCG NOMe.dedup.hcg.bed.gz \
  --sample_info sample_info_mouse_nome.txt \
  --output DMRs \
  --condition heart_chamber.noside
