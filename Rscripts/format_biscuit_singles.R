## reformat single sample VCF->BED biscuit files
## Merge + and - strand data from GCH Biscuit BED files and format with colnames/types
suppressPackageStartupMessages({
library(data.table)
library(tidyverse)
library(halpme)
})

library("argparse")
parser <- ArgumentParser()
parser$add_argument('bed_input', type="character",
                    help='BED file from biscuit of either GCH/HCG methylation')

parser$add_argument('--header', type="character",
                    help = 'File names from the VCF header')

args <- parser$parse_args()


read_lines = 1e7
buffer = 100

bed_file = args$bed_input
#bed_file = "biscuit/NM_RA_m1to5_rep1.dedup.hcg.bed.gz"
bed_out = gsub(".bed", ".formatted.bed", bed_file)

if(file.exists(args$header)){
  VCF_sample_header = read.table(args$header, quote="\"", comment.char="")
  VCF_sample_header = VCF_sample_header[1,]
}else{
  VCF_sample_header = args$header
  #VCF_sample_header = "NM_RA_m1to5_rep1"
}
#VCF_sample_header=read.table("../seq_round2/testing_biscuit/NOMe_combined.N.vcf.samplenames.txt", quote="\"", comment.char="")
counter = 1
bed_ended = F
while(!bed_ended){

  bed_read_start = ((counter-1)*read_lines)
  bed = fread(bed_file, data.table = F, header = F, skip = bed_read_start, nrows = read_lines+buffer)
  message(bed[,1] %>% unique() %>% paste(collapse = " "))

  if(nrow(bed) < read_lines){
    bed_ended = T
  }
  #bed_og = fread(bed_file, data.table = F, header = F)

  ### test col 7 == 5mer
  bed_7 = bed[,7][1:min(10, nrow(bed))]
  if(all(nchar(bed_7) == 5) & all(unlist(strsplit(bed_7, "")) %in% c("A", "C", "G", "T"))){
    ## Bed format with 7 columns
    colnames(bed)[1:7] = c("seqnames", "start", "end", "strand", "c", "c_context", "c_5mer")
    last_col = 7
  }else{
    ## Bed format with 3 columns
    colnames(bed)[1:3] = c("seqnames", "start", "end")
    message("BED file does not contain C context, cannot merge strands")
    last_col = 3
  }

  sample_names = as.character(VCF_sample_header)
  if(length(sample_names)*2 == length(colnames(bed)[-c(1:last_col)])){
    premerged = FALSE
    colnames(bed)[-c(1:last_col)] = as.character(VCF_sample_header) %>% halpme::strv_split(., "[.]", 1) %>% rep(., each=2) %>% paste0(., c(".meth_percent", ".cov"))
  }else if(length(sample_names)*3 == length(colnames(bed)[-c(1:last_col)])){
    premerged = TRUE
    colnames(bed)[-c(1:last_col)] = as.character(VCF_sample_header) %>% halpme::strv_split(., "[.]", 1) %>% rep(., each=3) %>% paste0(., c(".meth_percent", ".cov", ".cg"))
    bed = bed[,str_sub(colnames(bed), -3,-1)!=".cg"]
  }else{
    stop("sample names from VCF header file are not a multiple of BED file methylation columns. Please check that these match.")
  }


  #### convert to C_counts
  sample_names = as.character(colnames(bed)[-c(1:last_col)]) %>%
    str_replace_all(".cov", "") %>% str_replace_all(".meth_percent", "") %>%
    unique()

  for(sample in sample_names){
    if(counter == 1){
      message(sample)
      message(paste0("converting methylation percentages to counts"))
    }
    if(any(colnames(bed) == paste0(sample, ".meth_percent"))){
      index = c(which(colnames(bed) == paste0(sample, ".meth_percent")),
                which(colnames(bed) == paste0(sample, ".cov")))

      bed[,index[1]]
      has_reads = which(bed[,index[2]] > 0)

      C_count = round(as.numeric(bed[,index[1]][has_reads])*bed[,index[2]][has_reads])
      bed[,index[1]][has_reads] = C_count
      bed[,index[1]][-has_reads] = 0

      colnames(bed)[index][1] = paste0(sample, ".C")
    }

    if(any(colnames(bed) == paste0(sample, ".C"))){
      bed[,which(colnames(bed) == paste0(sample, ".C"))] = as.numeric(bed[,which(colnames(bed) == paste0(sample, ".C"))])
    }
  }


  # chrom_order = factor(c(paste0("chr", c(1:99)), paste0("chr", c("X", "Y", "M"))),
  #                      levels = c(paste0("chr", c(1:99)), paste0("chr", c("X", "Y", "M"))))
  #
  # # add in any undefined chromosomes
  # if(any(!(unique(bed$seqnames) %in% chrom_order))){
  #   chrom_order_add = unique(bed$seqnames)[!(unique(bed$seqnames) %in% chrom_order)]
  #   chrom_order_keep = unique(bed$seqnames)[(unique(bed$seqnames) %in% chrom_order)]
  #
  #   m = match(chrom_order_add, unique(bed$seqnames))
  #   if(m[1] == 1){
  #     chrom_order = factor(c(chrom_order_add, paste0("chr", c(1:99)), paste0("chr", c("X", "Y", "M"))),
  #                          levels = c(chrom_order_add, paste0("chr", c(1:99)), paste0("chr", c("X", "Y", "M"))))
  #   }else if(any(m == length(unique(bed$seqnames)))){
  #     chrom_order = factor(c(paste0("chr", c(1:99)), paste0("chr", c("X", "Y", "M")), chrom_order_add),
  #                          levels = c(paste0("chr", c(1:99)), paste0("chr", c("X", "Y", "M")), chrom_order_add))
  #   }else{
  #     ## REDUNDANT WITH IF().... KEEP FOR TESTING
  #     chrom_order = factor(c(chrom_order_add, paste0("chr", c(1:99)), paste0("chr", c("X", "Y", "M"))),
  #                          levels = c(chrom_order_add, paste0("chr", c(1:99)), paste0("chr", c("X", "Y", "M"))))
  #
  #   }
  #
  # }


  #### Combine C/G strands
  if(last_col == 7 & premerged == FALSE){
    bed.fwd = bed[bed$strand == "C",]
    bed.rev = bed[bed$strand == "G",]
    if(nrow(bed.rev) > 0){
      message("Merging C and G strand methylation")
      chrom_order = unique(bed.fwd[,1])
      chrom_order = factor(chrom_order, levels = chrom_order)

      bed.rev$start = bed.rev$start -1
      bed.fwd$end = bed.fwd$end + 1

      m = match(bed.rev$start, bed.fwd$start)

      cov_cols = which(str_sub(colnames(bed), -4,-1) == ".cov")
      c_cols = which(str_sub(colnames(bed), -2,-1) == ".C")

      bed.fwd[m[!is.na(m)], cov_cols] = bed.fwd[m[!is.na(m)], cov_cols] + bed.rev[!is.na(m), cov_cols]
      bed.fwd[m[!is.na(m)], c_cols] = bed.fwd[m[!is.na(m)], c_cols] + bed.rev[!is.na(m), c_cols]
      bed.fwd = rbind(bed.fwd, bed.rev[which(is.na(m)),])
      rm(bed.rev)

      bed.fwd$seqnames= factor(bed.fwd$seqnames, levels = chrom_order)
      bed.fwd = arrange(bed.fwd, seqnames, start)

    }
    bed.fwd = bed.fwd[,-which(colnames(bed.fwd) %in% c("c", "c_context", "c_5mer"))]
    bed.fwd$strand = "+"

    bed.fwd$seqnames = factor(bed.fwd$seqnames, levels = chrom_order)

    if(counter > 1){
      ## remove first few rows that have already been written
      new_start = which(as.character(bed.fwd$seqnames) == as.character(last_row$seqnames) & bed.fwd$start == last_row$start & bed.fwd$end == last_row$end)
      if(length(new_start) > 0){
      bed.fwd = bed.fwd[(new_start+1):nrow(bed.fwd),]
      }else{
        warning("first row of bed chunk not found in previous chunk buffer... check buffer size is OK")
      }
    }
    last_row = bed.fwd[nrow(bed.fwd)-2,]

    if(counter == 1 & bed_ended == T){
      bed.fwd = arrange(bed.fwd, seqnames, start)
    }

    if(counter == 1){
      write_delim(bed.fwd[1:(nrow(bed.fwd)-2),], bed_out,  delim = "\t")
    }else{
      write_delim(bed.fwd[1:(nrow(bed.fwd)-2),], bed_out,  delim = "\t", append = T, col_names = F)
    }
    counter = counter+1



  }else{

    if(counter > 1){
      ## remove first few rows that have already been written
      new_start = which(bed$seqnames == last_row$seqnames & bed$start == last_row$start & bed$end == last_row$end)
      if(length(new_start) > 0){
        bed.fwd = bed.fwd[(new_start+1):nrow(bed.fwd),]
      }else{
        warning("first row of bed chunk not found in previous chunk buffer... check buffer size is OK")
      }
    }
    last_row = bed[nrow(bed)-2,]
    if(counter == 1){
      write_delim(bed, bed_out,  delim = "\t")
    }else{
      write_delim(bed, bed_out,  delim = "\t", append = T, col_names = F)
    }
  }

}


