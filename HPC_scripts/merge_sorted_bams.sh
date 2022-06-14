#!/usr/bin/bash -l

#PBS -l walltime=100:00:00
#PBS -l mem=16GB
#PBS -l ncpus=2
#PBS -J 1-24

cd $PBS_O_WORKDIR
# ~/PROJECT/bams/sorted

module load rosalind
module load gcc-env/6.4.0
module load samtools
module load bedtools

FILENAME=`sed "${PBS_ARRAY_INDEX}q;d" merged.samples.txt`

samtools merge "${FILENAME}".merged.bam "${FILENAME}"*.bam
samtools sort "${FILENAME}".merged.bam -o "${FILENAME}".sortmerge.sorted.bam
samtools index "${FILENAME}".sortmerge.sorted.bam

#remove unsorted BAM
if [ -f "${FILENAME}".sortmerge.sorted.bam.bai ]; then
   rm -f "${FILENAME}".merged.bam
fi

samtools flagstat "${FILENAME}".sortmerge.sorted.bam > ../flagstat/"${FILENAME}".sortmerge.sorted.flagstat


# deduplicate
samtools rmdup "${FILENAME}".sortmerge.sorted.bam "${FILENAME}".merged.dedup.bam

mkdir -p ../dedup
mv "${FILENAME}".merged.dedup.bam ../dedup/

#flagstat sorted and dedup bams
samtools index ../dedup/"${FILENAME}".merged.dedup.bam
samtools flagstat ../dedup/"${FILENAME}".merged.dedup.bam > ../flagstat/"${FILENAME}".merged.dedup.flagstat
