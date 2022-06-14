#PBS -l walltime=100:00:00
#PBS -l mem=42GB
#PBS -l ncpus=16
#PBS -J 1-86

cd $PBS_O_WORKDIR
# ~/PROJECT/bams

module load rosalind
module load gcc-env/6.4.0
module load samtools/1.12
module load bwa
source ~/pythonenvs/bwameth_py2/bin/activate

FILENAME=`sed "${PBS_ARRAY_INDEX}q;d" filenames.txt`

bwameth.py --threads 16 \
 --reference ../indices/bwameth_index/hg38.fa \
 ../fastq/trimmed/"${FILENAME}"_1.fq.gz \
 ../fastq/trimmed/"${FILENAME}"_2.fq.gz | samtools view -Sb -F 4 - > "${FILENAME}".bam

samtools sort "${FILENAME}".bam -o "${FILENAME}".sorted.bam

mkdir -p sorted
mkdir -p flagstat

mv "${FILENAME}".sorted.bam sorted/

samtools index sorted/"${FILENAME}".sorted.bam
samtools flagstat sorted/"${FILENAME}".sorted.bam > flagstat/"${FILENAME}".sorted.flagstat

#remove unsorted BAM
if [ -f sorted/"${FILENAME}".sorted.bam.bai ]; then
    rm -f "${FILENAME}".bam
fi
