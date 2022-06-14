#PBS -l walltime=100:00:00
#PBS -l mem=16GB
#PBS -l ncpus=4
#PBS -J 1-24

cd $PBS_O_WORKDIR
#~PROJECT/fastq/

module load rosalind
module load Python/3.7.4-GCCcore-8.3.0
source ~/pythonenvs/trim3/bin/activate
#mkdir -p trimmed

FILENAME=`sed "${PBS_ARRAY_INDEX}q;d" filenames.txt`

module load fastqc
fastqc -t 4 "${FILENAME}"_r1.fastq.gz
fastqc -t 4 "${FILENAME}"_r2.fastq.gz

mkdir -p fastqc_reports
mv *fastqc.* fastqc_reports/

# trim fastq files
~/apps/TrimGalore-0.6.6/trim_galore --paired -j 4 \
  "${FILENAME}"_r1.fastq.gz \
  "${FILENAME}"_r2.fastq.gz

# moved to trimmed folder
mv "${FILENAME}"_r1_val_1.fq.gz trimmed/
mv "${FILENAME}"_r2_val_2.fq.gz trimmed/

mkdir -p trimming_reports
mv *trimming_report.txt trimming_reports/

# rename files to get rid of _r2_val_2
rename '_r1_val_1' '_1' trimmed/*
rename '_r2_val_2' '_2' trimmed/*

# run fastqc
fastqc -t 4 trimmed/"${FILENAME}"_1.fq.gz
fastqc -t 4 trimmed/"${FILENAME}"_2.fq.gz

mkdir -p trimmed/fastqc_reports
mv trimmed/*fastqc.* trimmed/fastqc_reports/
