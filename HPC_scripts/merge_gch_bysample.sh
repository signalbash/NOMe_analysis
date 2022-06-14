#PBS -l walltime=100:00:00
#PBS -l mem=16GB
#PBS -l ncpus=1
#PBS -J 1-25

cd $PBS_O_WORKDIR

FILENAME=`sed "${PBS_ARRAY_INDEX}q;d" chroms.txt`
cd NOME_bychr/

module load rosalind
module load R/4.1.0-foss-2021a

Rscript ../merge_biscuit_tables.R \
--sample_info ../sample_info.tsv \
-d . --directory_patterns gch,"${FILENAME}" \
-o NOMe.dedup.gch."${FILENAME}".bed.gz


Rscript ../merge_biscuit_tables.R --sample_info ../sample_info.tsv -d . --directory_patterns gch,chr22,dedup -o NOMe.dedup.gch.chr22.bed.gz
