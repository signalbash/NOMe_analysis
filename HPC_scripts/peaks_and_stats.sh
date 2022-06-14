#!/usr/bin/bash -l

#PBS -l walltime=100:00:00
#PBS -l mem=16GB
#PBS -l ncpus=2
#PBS -J 1-27

cd $PBS_O_WORKDIR

FILENAME=`sed "${PBS_ARRAY_INDEX}q;d" ac_filenames.txt`

module load rosalind
module load gcc-env/6.4.0
module load samtools
module load bedtools
module load bcftools/1.12

#samtools rmdup ../../sorted/"${FILENAME}".sorted.bam "${FILENAME}".sam.bam
#samtools index "${FILENAME}".sam.bam

#~/apps/biscuit/biscuit pileup -N -o "${FILENAME}".sam.bam.vcf /u/bsignal/Scratch/indices/biscuit_index/hg38.fa "${FILENAME}".sam.bam
#~/apps/biscuit/biscuit vcf2bed -k 1 -t hcg -e -s ALL "${FILENAME}".sam.bam.vcf > "${FILENAME}".sam.bam.hcg.bed
#~/apps/biscuit/biscuit vcf2bed -k 1 -t gch -e -s ALL "${FILENAME}".sam.bam.vcf > "${FILENAME}".sam.bam.gch.bed

gzip "${FILENAME}".sam.bam.*.bed

#module purge
#module load rosalind
#module load Python/3.7.2-GCCcore-8.2.0
#source ~/pythonenvs/macs2_RPy3/bin/activate

#macs2 callpeak -t "${FILENAME}".sam.bam -c ../200_IgG.dedup.bam ../Orig_IgG.dedup.bam ../RPT_IgG.dedup.bam -f BAM -g hs -p 1e-3 -n "${FILENAME}"

module purge
module load rosalind
module load R/4.1.0-foss-2021a

 Rscript ./in_peak_meth_summary.R --peaks "${FILENAME}"_peaks.narrowPeak \
  --cpg "${FILENAME}".sam.bam.hcg.bed.gz \
  --gpc "${FILENAME}".sam.bam.gch.bed.gz \
  -o "${FILENAME}".sam.bam.txt \
  --blacklist /u/bsignal/Scratch/NOMe/biscuit/dedup/hg38-blacklist.v2.bed.gz
