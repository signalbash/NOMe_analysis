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

Makes:
fastq/trimmed/$SAMPLE_1.fq.gz
fastq/trimmed/$SAMPLE_2.fq.gz

fastqc reports

`count_reads.sh`
Counts all reads in untrimmed samples.
Makes:
fastq_file_lines.csv

## Alignment to genome
We can now align reads to a reference genome. For this step we're using bwameth (https://github.com/brentp/bwa-meth), which builds off bwa-mem (standard aligner). bwameth converts the genome into two seperate genomes - one with C->T conversion, and the other with G->A conversion. The fastq reads also undergo this conversion, so instead of using the typical 3-letter nucleotide alphabet to align, we are only using 3 (for fwd/rev reads).

`bwameth_array.sh`
Requires:
indices/bwameth_index/hg38.fa (bwameth index)
../fastq/trimmed/$SAMPLE_1.fq.gz (trimmed fastq)
../fastq/trimmed/$SAMPLE_2.fq.gz (trimmed fastq)

Makes:
bams/sorted/"${FILENAME}".sorted.bam
flagstat/"${FILENAME}".sorted.flagstat

This should align, then sort and index, and run flagstat on the resulting BAM file.

**Merging of samples sequenced across multiple lanes**
We merge after alignment as it means the alignment step is performed much quicker, and on par with times for smaller files.
Merging is done through samtools merge. These files need to be re-sorted following merging. Bams are then deduplicated (PCR duplicates have issues when calling methylation ratios), indexed, and stats gathered using samtools flagstat.

`merge_sorted_bams.sh`

This takes all bam files with a common prefix (in merged.samples.txt) and merges them.

Requires:
bams/sorted/"${FILENAME}"*.bam

Makes:
bams/sorted/"${FILENAME}".sortmerge.sorted.bam (merged bam, only sorted and not deduplicated)
bams/flagstat/"${FILENAME}".sortmerge.sorted.flagstat

bams/dedup/"${FILENAME}".merged.dedup.bam (merged and deduplicated bam)
bams/flagstat/"${FILENAME}".merged.dedup.flagstat

**Deduplication of bams if merging is not required**
`deduplicate_bams.sh`

Makes:
bams/dedup/"${FILENAME}".dedup.bam
bams/flagstat/"${FILENAME}".dedup.flagstat


## Read depth profiling with mosdepth

`mosdepth_array.sh`

Makes:

Requires:

### Alignment to (lambda) genome
Same as above, but aligned to the lambda genome as internal controls.

## Methylation calling
We use biscuit to pileup reads and call methylation status (i.e. is the C converted to a T or not).
https://huishenlab.github.io/biscuit/docs/pileup.html
We can either call all samples at once into a single output file (note that this is LARGE when uncompressed... >1TB but compression makes it a lot smaller. This can also take over a day to run due to sample size), or individual files. We have opted for individual in this case.

We use the -N flag as this is NOMe data (helps with downstream splitting into HCG/GCH)

We then convert the VCF to bed files with either CpG methylation (typical methylation; HCG) or GpC methylation (NOMe methylation; HCG), and remove any regions that overlap blacklisted regions.

converting to bed can require some extra IMPORTANT flags
-k 1 : minimum coverage of 1 (these are some low coverage samples so we don't really want to throw out data)
-t [hcg/gch] : extract all hCG or GCh nts (so both NOMe GpC and CpG sites are covered in seperate files)
-e : show context (so we can tell if this is a GpC/CpG/other)
-s ALL : extract information for **ALL** samples (important if running on a merged VCF)

We also remove blacklisted regions and convert biscuit bed files to tables with C counts and coverage counts (instead of methylation percent). This is to help with downstream processing.

`pileup_to_bed.sh`

Requires:
indices/biscuit_index/hg38.fa
bams/dedup/"${FILENAME}".merged.dedup.bam

Makes:
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

`./split_GCH_bed_bychrom.sh`

Requires:
"${FILENAME}".merged.dedup.hcg.formatted.bed.gz
Makes:
NOME_bychr/"${FILENAME}".merged.dedup.hcg.formatted.chr[1-22].bed.gz

**merging**
merge_gch_bysample.sh
merge_hcg_bysample.sh
Requires:
`merge_biscuit_tables.R`


# MINGLE-specific analyses

## Coverage
`bam_to_bigwig.sh`

Uses deeptools bam Coverage in 10nt blocks. Normalised using RPKM.
Make profiles around gencode transcripts and plot to PDF.
Requires:
"${FILENAME}".sam.bam
Makes:
"${FILENAME}".sam.cov10.bw
"${FILENAME}".gencode.matrix.gz
"${FILENAME}".profile.pdf

After this, manually check controls etc. if they look typical for peaks.

`bam_to_bigwig_control.sh`

Same as above, but uses an IgG control to control for read coverage.

## Call Peaks
`peaks_and_stats_narrow.sh`
`peaks_and_stats_broad.sh`

Calls peaks with MACS2, runs summary R script (`in_peak_meth_summary.R`) to extract CpG and GpC methylation levels from within peaks.
Requires:
"${FILENAME}".sam.bam
Makes:
"${FILENAME}"_peaks.broadPeak
OR
"${FILENAME}"_peaks.narrowPeak
"${FILENAME}".sam.bam.txt

### Mosdepth coverage in peaks


### Check methylation levels in peaks





# Scripts for analysing NOMe-Seq data

## QC-checks


## Plotting functions


## More QC: PCAs, etc.



## Differential methylation (DMRs and diff. CpG sites)
"${FILENAME}".merged.dedup.hcg.formatted.bed.gz
"${FILENAME}".merged.dedup.gch.formatted.bed.gz

## Nucleosome Occupancy



## Differential Nucleosome Occupancy (DiNO)


## Annotate Regions
