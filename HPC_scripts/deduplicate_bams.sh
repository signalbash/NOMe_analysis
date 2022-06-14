#PBS -l walltime=100:00:00
#PBS -l mem=42GB
#PBS -l ncpus=16
#PBS -J 1-86

cd $PBS_O_WORKDIR
# ~/PROJECT/bams/sorted/

module load rosalind
module load gcc-env/6.4.0
module load samtools/1.12


FILENAME=`sed "${PBS_ARRAY_INDEX}q;d" filenames.txt`

# deduplicate
samtools rmdup "${FILENAME}".sorted.bam "${FILENAME}".dedup.bam

mkdir -p ../dedup
mv "${FILENAME}".dedup.bam ../dedup/

#flagstat sorted and dedup bams
samtools index ../dedup/"${FILENAME}".dedup.bam
samtools flagstat ../dedup/"${FILENAME}".dedup.bam > ../flagstat/"${FILENAME}".dedup.flagstat
