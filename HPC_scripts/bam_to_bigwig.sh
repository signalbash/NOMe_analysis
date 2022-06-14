#!/usr/bin/bash -l

#PBS -l walltime=100:00:00
#PBS -l mem=8GB
#PBS -l ncpus=2
#PBS -J 1-60

### Initialisation

cd $PBS_O_WORKDIR

module load rosalind
module load Python/3.9.5-GCCcore-10.3.0

## Using deeptools

FILENAME=`sed "${PBS_ARRAY_INDEX}q;d" sambams.txt`

# get coverage bigwig
 ~/.local/bin/bamCoverage -b "${FILENAME}".sam.bam -of bigwig -bs 10 \
 -bl /u/bsignal/Scratch/NOMe/biscuit/dedup/hg38-blacklist.v2.bed.gz \
 --normalizeUsing RPKM -o "${FILENAME}".sam.cov10.bw

# compute matrix +/- 5kb
 ~/.local/bin/computeMatrix scale-regions -S "${FILENAME}".sam.cov10.bw \
 -R gencode.v39.annotation.gtf.gz \
 -o "${FILENAME}".gencode.matrix.gz \
 -b 5000 -a 5000 \
 -bl /u/bsignal/Scratch/NOMe/biscuit/dedup/hg38-blacklist.v2.bed.gz \
 --smartLabels -p max/2

# plot profile (no controls)
 ~/.local/bin/plotProfile -m "${FILENAME}".gencode.matrix.gz \
 -o "${FILENAME}".profile.pdf \
 --outFileNameData "${FILENAME}".profile.txt \
 --averageType mean \
 --plotType lines
