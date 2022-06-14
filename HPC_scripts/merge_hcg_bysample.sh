#PBS -l walltime=100:00:00
#PBS -l mem=16GB
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

module load rosalind
module load R/4.1.0-foss-2021a

Rscript ./merge_biscuit_tables.R \
--sample_info sample_info.tsv \
-d . --directory_patterns hcg,formatted,dedup \
-o NOMe.dedup.hcg.bed.gz
