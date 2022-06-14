bed=$1
file="${bed##*/}"
bedout="${file%%.bed.gz}"

chrom=$2

gunzip -c $1 | head -n 1 > $bedout.$2.bed

if [ $2 = "chr1" ]; then
  gunzip -c $1| awk '$1 == "chr1" {print $0}'  >> $bedout.chr1.bed
elif [ $2 = "chr2" ]; then
 gunzip -c $1 | awk '$1 == "chr2" {print $0}'  >> $bedout.chr2.bed
else
 gunzip -c $1 | grep $2 >> $bedout.$2.bed
fi

gzip $bedout.chr2.bed
