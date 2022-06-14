#!/usr/bin/bash -l

#PBS -l walltime=100:00:00
#PBS -l mem=16GB
#PBS -l ncpus=2
#PBS -J 1-60

### Initialisation

cd $PBS_O_WORKDIR

module load rosalind
module load Python/3.9.5-GCCcore-10.3.0

FILENAME=`sed "${PBS_ARRAY_INDEX}q;d" sambams.txt`

## 200/200/RPT

~/.local/bin/bamCompare -b1 "${FILENAME}".sam.bam \
-b2 200_IgG.sam.bam \
-o "${FILENAME}".200_IgG.bw \
-of bigwig -bs 10 \
-bl /u/bsignal/Scratch/NOMe/biscuit/dedup/hg38-blacklist.v2.bed.gz

~/.local/bin/computeMatrix scale-regions \
-S "${FILENAME}".200_IgG.bw \
-R gencode.v39.annotation.gtf.gz \
-o "${FILENAME}".200_IgG.gencode.matrix.gz \
-b 1000 -a 1000 \
-bl /u/bsignal/Scratch/NOMe/biscuit/dedup/hg38-blacklist.v2.bed.gz \
--smartLabels -p 2

~/.local/bin/plotProfile \
-m "${FILENAME}".200_IgG.gencode.matrix.gz \
-o "${FILENAME}".200_IgG.profile.pdf \
--outFileNameData "${FILENAME}".200_IgG.profile.txt \
--averageType mean \
--plotType lines
