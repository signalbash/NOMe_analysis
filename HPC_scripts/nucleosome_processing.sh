#!/usr/bin/bash -l

#PBS -l walltime=100:00:00
#PBS -l mem=32GB
#PBS -l ncpus=4
#PBS -J 1-22

module purge
module load rosalind
module load R/4.1.0-foss-2021a

cd $PBS_O_WORKDIR

FILENAME=`sed "${PBS_ARRAY_INDEX}q;d" chroms.txt`

# format BED files into tables
# combines strand data, removes extra information so files are smaller
Rscript ./nucleosome_processing.R \
  --GCH NOME_bychr/NOMe.dedup.gch."${FILENAME}".bed.gz \
  --sample_info sample_info_mouse_nome.txt \
  --output NDRS_"${FILENAME}".Rdata
