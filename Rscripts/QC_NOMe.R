library(data.table)
library(tidyverse)
library(halpme)
setwd("~/Projects/UTAS/NOMe/merged_data")

source("~/Projects/UTAS/NOMe_analysis/Rscripts/QC_functions.R")

multiqc_bin_location = "/Users/bsignal/opt/miniconda3/bin/multiqc"
sample_info = fread("sample_info_nome.tsv", data.table = F)
seqid2sampleid = fread("../sample_info_all.tsv", data.table = F)


seq_reports = read_in_seqreports(report_files = c("../seq_round2/SequencingProjectReport.html",
                                   "../SequencingProjectReport.html",
                                   "../SequencingProjectReport_dec1.html" ), 
                                 sample_info = sample_info,
                                 seqid2sampleid = seqid2sampleid)

total_seqs_raw = read_in_fastq_lines(summary_file = "../fastq_file_lines.csv", 
                                     sample_info = sample_info,
                                     seqid2sampleid = seqid2sampleid)

total_seqs_trimmed = read_in_fastqc_stats(directory = "fastqc", 
                                          sample_info = sample_info, 
                                          multiqc_bin_location="/Users/bsignal/opt/miniconda3/bin/multiqc")

### AUTO SEARCH FOR FLAGSTAT FILES
samtools_flagstat = read_in_flagstat(multiqc_filetype = "multiqc_samtools_flagstat.txt",
                                     multiqc_bin_location = "/Users/bsignal/opt/miniconda3/bin/multiqc")

table(samtools_flagstat$directory)
paste0(samtools_flagstat$directory,"/", samtools_flagstat$sample)

mosdepth_cov = read_in_mosdedpth()
check_sample_hasdata(colnames(mosdepth_cov), sample_info$sample_id)
## mitochondrial / total reads
mito_bias = mosdepth_cov[mosdepth_cov$chrom == "chrM",-c(1,2)] / mosdepth_cov[mosdepth_cov$chrom == "total",-c(1,2)]
autosome_coverage = colSums(mosdepth_cov[mosdepth_cov$chrom %in% c(paste0("chr", c(1:22))), -c(1,2)]) / 
  sum(mosdepth_cov$length[mosdepth_cov$chrom %in% c(paste0("chr", c(1:22)))])


#####






#untrimmed_stats = filter(untrimmed_stats, Seq %in% c("NOMe"))

aggregate(m_seqs~sample_number, untrimmed_stats, function(x) sum(as.numeric(x)))

untrimmed_stats %>% mutate(pd = as.numeric(gsub("%", "", percent_dups))) %>% 
  mutate(seqs=as.numeric(m_seqs)*1e7) %>% 
ggplot( aes(x=pd, y=seqs, col=sample_number, shape=seq_round)) + geom_point()

untrimmed_stats %>% mutate(pd = as.numeric(gsub("%", "", percent_dups))) %>%
   mutate(seqs=as.numeric(m_seqs)*1e7) %>%
  ggplot(aes(x=seq_round, y=pd)) + geom_boxplot()

untrimmed_stats %>% mutate(pd = as.numeric(gsub("%", "", percent_dups))) %>%
  mutate(seqs=as.numeric(m_seqs)*1e7) %>%
  ggplot( aes(x=pd, y=seqs, fill=sample_number)) + geom_point(aes(shape=seq_round)) + scale_y_log10()+geom_mark_ellipse(alpha=0.1)

 untrimmed_stats %>% mutate(pd = as.numeric(gsub("%", "", percent_dups))) %>%
  mutate(seqs=as.numeric(m_seqs)*1e7) %>% 
   #filter(seq_round !="mingle") %>% 
   filter(!(seq_round=="round2" & Seq=="MINGLE")) %>% 
   filter(sample_number %in% untrimmed_stats$sample_number[untrimmed_stats$seq_round == "mingle"]) %>% 
  ggplot( aes(y=pd, x=sample_number, fill=sample_number)) + geom_boxplot() + geom_point(aes(shape=seq_round)) + ggeasy::easy_rotate_x_labels()


#### BAM Alignment numbers

#bam_stats.sorted = fread("flagstat/sorted/multiqc_data/multiqc_samtools_flagstat.txt", data.table = F)
bam_stats.dedup = fread("mosdepth/multiqc_data/multiqc_samtools_flagstat.txt", data.table = F)

#untrimmed_stats$bam_reads = bam_stats.sorted$flagstat_total[match(untrimmed_stats$sample_id, bam_stats.sorted$Sample)] /2
untrimmed_stats$bam_dedup_reads = bam_stats.dedup$flagstat_total[match(untrimmed_stats$sample_id, bam_stats.dedup$Sample)] /2
#untrimmed_stats$percent_duplicated_bam = 1 - (untrimmed_stats$bam_dedup_reads  / untrimmed_stats$bam_reads)

ggplot(untrimmed_stats, aes(x=as.numeric(m_seqs)*1e6, y=bam_dedup_reads)) + geom_point() + geom_abline(slope=1)
#ggplot(untrimmed_stats, aes(x=percent_duplicated_bam, y=bam_reads, col=chip)) + geom_point()

## LM: duplication tends to increase with increased depth.
untrimmed_stats %>% filter(bam_reads < 100e6) %>%
  ggplot( aes(x=percent_duplicated_bam, y=bam_reads)) + geom_point(aes(col=chip)) + geom_smooth(method="lm")

filter(untrimmed_stats, ad_pathology %in% c("Low", "Intermediate", "High")) %>% aggregate(percent_duplicated_bam~chip, ., function(x) summary(as.numeric(x)))
filter(untrimmed_stats, ad_pathology %in% c("Low", "Intermediate", "High")) %>% aggregate(bam_reads~chip, ., function(x) summary(as.numeric(x)))
filter(untrimmed_stats, ad_pathology %in% c("Low", "Intermediate", "High")) %>% aggregate(bam_dedup_reads, ., function(x) summary(as.numeric(x)))


##### MOSDEPTH
## chrom cov

bam_coverage = read_mosdepth_files(list.files("mosdepth/dedupbam/", "summary", full.names = T))
bam_coverage[bam_coverage$chrom == "chrM",-c(1,2)] / bam_coverage[bam_coverage$chrom == "total",-c(1,2)]

##### BISULPHITE CONVERSION EFFICIENCY
## and number aligned to lambda genome // percent total reads
##### BISULPHITE CONVERSION EFFICIENCY
## and number aligned to lambda genome // percent total reads

lambda_stats.dedup = fread("lambda/flagstat/multiqc_data/multiqc_samtools_flagstat.txt", data.table = F)

meth_percent_lambda_hcg = check_total_methylation_levels("lambda/NOMe_lambda_merged.hcg.black.merged.bed.gz")
meth_percent_lambda_gch = check_total_methylation_levels("lambda/NOMe_lambda_merged.gch.black.merged.bed.gz")

lambda_conversion = full_join(meth_percent_lambda_hcg,meth_percent_lambda_gch, by = "sample", suffix = c(".hcg", ".gch"))



