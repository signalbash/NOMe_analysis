# Scripts for processing NOMe-Seq data

Most scripts are formatted to work on a HPC cluster running the PBS job scheduler.
Some scripts use array jobs to submit a single job per sample.

## Making indices for bwameth and biscuit

`make_bwameth_index.sh`
`make_biscuit_index.sh`

## Trimming fastq files & QC checks
TrimGalore-0.6.6
Python 3.7.4

`trim_and_qc.sh`

Trims fastq files with default parameters, and runs fastQC on trimmed files.

**Makes:**

fastq/trimmed/$SAMPLE_1.fq.gz
fastq/trimmed/$SAMPLE_2.fq.gz

fastqc reports

`count_reads.sh`
Counts all reads in untrimmed samples.
Makes:
fastq_file_lines.csv

## Alignment to genome
We can now align reads to a reference genome. For this step we're using bwameth (https://github.com/brentp/bwa-meth), which builds off bwa-mem (standard aligner). bwameth converts the genome into two seperate genomes - one with C->T conversion, and the other with G->A conversion. The fastq reads also undergo this conversion, so instead of using the typical 3-letter nucleotide alphabet to align, we are only using 3 (for fwd/rev reads).

`bwameth2dedup.sh`
**Requires:**

indices/bwameth_index/GRCm39_pluslambda.fa(bwameth index)
../fastq/trimmed/$SAMPLE_1.fq.gz (trimmed fastq)
../fastq/trimmed/$SAMPLE_2.fq.gz (trimmed fastq)

**Makes:**

bams/sorted/"${FILENAME}".sorted.bam
bams/flagstat/"${FILENAME}".sorted.flagstat
bams/dedup/"${FILENAME}".dedup.bam
bams/flagstat/"${FILENAME}".dedup.flagstat


## Read depth profiling with mosdepth

`mosdepth_array.sh`

**Makes:**

mosdepth/"${FILENAME}".mosdepth.summary.txt
**Requires:**

dedup/"${FILENAME}".dedup.bam

## Methylation calling
We use biscuit to pileup reads and call methylation status (i.e. is the C converted to a T or not).
https://huishenlab.github.io/biscuit/docs/pileup.html
We can either call all samples at once into a single output file (note that this is LARGE when uncompressed... >1TB but compression makes it a lot smaller. This can also take over a day to run due to sample size), or individual files. We have opted for individual in this case.

We use the -N flag as this is NOMe data (helps with downstream splitting into HCG/GCH)

We then convert the VCF to bed files with either CpG methylation (typical methylation; HCG) or GpC methylation (NOMe methylation; HCG), and remove any regions that overlap blacklisted regions.

converting to bed can require some extra IMPORTANT flags

* -k 1 : minimum coverage of 1 (these are some low coverage samples so we don't really want to throw out data)
* -t [hcg/gch] : extract all hCG or GCh nts (so both NOMe GpC and CpG sites are covered in seperate files)
* -e : show context (so we can tell if this is a GpC/CpG/other)
* -s ALL : extract information for **ALL** samples (important if running on a merged VCF)
* -m 10 : minimum quality score of 10 required to count reads (default=40)

We also remove blacklisted regions and convert biscuit bed files to tables with C counts and coverage counts (instead of methylation percent). This is to help with downstream processing.

`pileup_to_bed.sh`

**Requires:**

indices/biscuit_index/GRCm39_pluslambda.fa
bams/dedup/"${FILENAME}".dedup.bam

**Makes:**

"${FILENAME}".merged.dedup.vcf.gz
"${FILENAME}".merged.dedup.hcg.bed.gz
"${FILENAME}".merged.dedup.gch.bed.gz
"${FILENAME}".merged.dedup.hcg.formatted.bed.gz
"${FILENAME}".merged.dedup.gch.formatted.bed.gz

**merging biscuit tables**

## Merging methylation calls for all samples
**splitting by chromosome**

HCG files *can* be combined into a larger table of all samples, as number of CpG sites is lower than GpC sites. However, given high enough coverage, this may not be the case and files may need to be split up by chromosome to process in the available RAM.
Using grep for unambiguous chromosome names (e.g. chr21) and awk for ambiguous chromsomes (e.g. chr2). grep is much faster than awk.

`./split_gch_bychrom.sh`
[calls] `split_bed_bychrom.sh`

Requires:
"${FILENAME}".dedup.hcg.formatted.bed.gz
Makes:
NOME_bychr/"${FILENAME}".dedup.hcg.formatted.chr[1-22].bed.gz

**merging**

`merge_gch_bysample.sh` (runs on single chromosomes)
`merge_hcg_bysample.sh`


**Requires:**

`merge_biscuit_tables.R`

**Makes:**

NOMe.dedup.hcg.bed.gz
NOMe.dedup.gch.bed.gz


## Differential methylation (DMRs and diff. CpG sites)
`dmrs.sh`
calls
`dmrs.R`

Calls DMRs (Bsseq)



## Nucleosome Occupancy

`nucleosome_processing.sh`
calls
`nucleosome_processing.R`
and runs on a **per-chromosome basis**, otherwise the data is too big.

Calls NDRs per-sample and per-condition.
Calls based on aaRon findNDRs(), but with updated data returned and window definition.
window width=140, rolling width=10

#nucleosome_scripts.R ??
`nucleosome_functions.R`

## Differential Nucleosome Occupancy (DiNO)
`diff_nucleosome_functions.R`
Contains altered functions for calling regions of differential nucleosome occupancy, based off methylsig.
Currently takes NDR regions as input to pre-filter for regions of interest.
Differential occupancy is called as the most 'significant' region within a 10kb window. This prevents over-representation of semi-overlapping windows.




## Annotate Regions
