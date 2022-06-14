#!/usr/bin/bash -l

#PBS -l walltime=100:00:00
#PBS -l mem=16GB
#PBS -l ncpus=2
#PBS -J 1-60

cd $PBS_O_WORKDIR

FILENAME=`sed "${PBS_ARRAY_INDEX}q;d" mingle_filenames.txt`

module load rosalind
module load Python/3.7.2-GCCcore-8.2.0
source ~/pythonenvs/macs2_RPy3/bin/activate

cd broad_2controls

macs2 callpeak -t ../../"${FILENAME}".sam.bam \
-c ../../200_IgG.sam.bam ../../Orig_IgG.sam.bam \
--broad --broad-cutoff 0.05 \
-f BAM -g hs -q 0.05 -n "${FILENAME}"

module purge
module load rosalind
module load R/4.1.0-foss-2021a

Rscript ~/in_peak_meth_summary.R --peaks "${FILENAME}"_peaks.broadPeak \
  --cpg ../../biscuit/black_beds/"${FILENAME}".sam.bam.hcg.bed.gz \
  --gpc ../../biscuit/black_beds/"${FILENAME}".sam.bam.gch.bed.gz \
  -o "${FILENAME}".peakmeth.txt \
  --blacklist /u/bsignal/Scratch/NOMe/biscuit/dedup/hg38-blacklist.v2.bed.gz

module purge
module load rosalind
module load mosdepth
mkdir -p mosdepth
cd mosdepth

mosdepth -b ../"${FILENAME}"_peaks.broadPeak "${FILENAME}" ../../../"${FILENAME}".sam.bam
