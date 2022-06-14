#PBS -l walltime=100:00:00
#PBS -l mem=8GB
#PBS -l ncpus=1
#PBS -J 1-24

cd $PBS_O_WORKDIR

mkdir -p NOME_bychr/
FILENAME=`sed "${PBS_ARRAY_INDEX}q;d" filenames.txt`

## contains list of chr1-->chr22,chrX,chrY,chrM
cat chroms.txt | while read chrom
do
  ./split_bed_bychrom.sh "${FILENAME}".dedup.gch.formatted.bed.gz $chrom
done

gzip "${FILENAME}"*chr*.bed
mv "${FILENAME}"*chr*.bed.gz NOME_bychr/
