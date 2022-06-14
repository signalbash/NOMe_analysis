#PBS -l walltime=100:00:00
#PBS -l mem=8GB
#PBS -l ncpus=1

cd $PBS_O_WORKDIR
#~PROJECT/fastq/

for x in *1.fastq.gz; do echo $x >> fastq_read_numbers.txt && zcat "$x" | wc -l >> fastq_read_numbers.txt ; done
for x in *1.fq.gz; do echo $x >> fastq_read_numbers.txt && zcat "$x" | wc -l >> fastq_read_numbers.txt ; done

echo 'sample,fastq_lines' > fastq_file_lines.csv
paste -s -d ',\n' fastq_read_numbers.txt >> fastq_file_lines.csv
