#!/usr/bin/bash -l

#PBS -l walltime=100:00:00
#PBS -l mem=32GB
#PBS -l ncpus=4
#PBS -J 1-24

cd $PBS_O_WORKDIR

#cd ~/PROJECT/biscuit

FILENAME=`sed "${PBS_ARRAY_INDEX}q;d" filenames.txt`

module purge
module load rosalind
module load gcc-env/6.4.0
module load samtools
module load bedtools
module load bcftools/1.12

#remove unsorted BAM
# pileup, NOME options
~/apps/biscuit/biscuit pileup -q 4 -u -N -m 10 -o "${FILENAME}".dedup.vcf \
  /u/bsignal/Scratch/indices/biscuit_index/GRCm39_pluslambda.fa \
  ../bams/dedup/"${FILENAME}".dedup.bam

# compress
bgzip "${FILENAME}".dedup.vcf
tabix -p vcf "${FILENAME}".dedup.vcf.gz

# make CpG (HCG) BED file
# -k 1 = one read minimum (default=3)
~/apps/biscuit/biscuit vcf2bed -k 1 -t hcg -e "${FILENAME}".dedup.vcf.gz  > "${FILENAME}".dedup.hcg.bed
# remove blacklisted regions
bedtools intersect -v -a "${FILENAME}".dedup.hcg.bed \
-b ~/Scratch/mouse_heart_nome/mm10-blacklist.v2.Liftover.mm39.bed.txt \
  -wa > "${FILENAME}".dedup.hcg.black.bed

mv "${FILENAME}".dedup.hcg.black.bed "${FILENAME}".dedup.hcg.bed
gzip -f "${FILENAME}".dedup.hcg.bed


# make GpC (GCH) BED file
~/apps/biscuit/biscuit vcf2bed -k 1 -t gch -e "${FILENAME}".dedup.vcf.gz  > "${FILENAME}".dedup.gch.bed
bedtools intersect -v -a "${FILENAME}".dedup.gch.bed \
-b ~/Scratch/mouse_heart_nome/mm10-blacklist.v2.Liftover.mm39.bed.txt \
  > "${FILENAME}".dedup.gch.black.bed

mv "${FILENAME}".dedup.gch.black.bed "${FILENAME}".dedup.gch.bed
gzip -f "${FILENAME}".dedup.gch.bed

module purge
module load rosalind
module load R/4.1.0-foss-2021a

# format BED files into tables
# combines strand data, removes extra information so files are smaller
Rscript ~/format_biscuit_singles.R "${FILENAME}".dedup.gch.bed.gz --header "${FILENAME}"
Rscript ~/format_biscuit_singles.R "${FILENAME}".dedup.hcg.bed.gz --header "${FILENAME}"
