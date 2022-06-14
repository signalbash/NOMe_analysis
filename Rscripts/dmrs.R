library(bsseq)
library(dmrseq)
library(data.table)
library(tidyverse)
library(halpme)
library(BiocParallel)
library(GenomicRanges)


if(interactive()){
  source("~/Projects/UTAS/NOMe_analysis/Rscripts/nucleosome_functions.R")
  source("~/Projects/UTAS/NOMe/scripts/nucleosome_scripts.R")

  sample_to_treat = read.delim("sample_info_mouse_nome.txt")


  cpg_file = paste0("biscuit/NOMe.dedup.hcg.bed.gz")
  outfile = "DMRs.csv"

  condition = "heart_chamber.noside"
  condition_1 = "V"
  condition_2 = "A"
  n_cores = 1

}else{
  source("~/nome_scripts/nucleosome_functions.R")
  source("~/nome_scripts/nucleosome_scripts.R")

  library("argparse")
  parser <- ArgumentParser()
  parser$add_argument('--sample_info', type="character",
                      help='CSV/TSV file containing sample_ids and categories')
  parser$add_argument('--HCG', type="character", default = "./",
                      help='file with HCG (CpG) biscuit formatted bed')
  parser$add_argument('--condition', type="character",
                      help='which condition from --sample_info to compare')
  parser$add_argument('--condition_1', type="character", default=NULL,
                      help='specify which factor to compare (i.e) control')
  parser$add_argument('--condition_2', type="character", default=NULL,
                      help='specify which factor to compare (i.e) case')
  parser$add_argument('-t', '--threads', type="integer", default=1,
                      help='Number of threads for multicore processing [default %(default)s]')
  parser$add_argument('--output', '-o', type="character", default="DMRs",
                      help='output prefix [default %(default)s]')

  args <- parser$parse_args()
  sample_to_treat = read.delim(args$sample_info)
  cpg_file = args$HCG
  outfile = args$output

  condition = args$condition
  condition_1 = args$condition_1
  condition_2 = args$condition_2
  n_cores = args$threads


}



if(is.null(condition_1) | is.null(condition_2)){
  message("no specific conditions given... making all possible comparisons...")
  combinations = make_combinations(sample_to_treat, condition, levels_to_compare=unique(sample_to_treat[,condition]))
}else{
  combinations = make_combinations(sample_to_treat, condition, levels_to_compare=c(condition_1, condition_2))
}


cpg_table = fread(cpg_file, data.table = F)
cpg_table = cpg_table[!is.na(cpg_table[,1]),]

all_chroms=  as.character(unique(cpg_table[,1]))
all_chroms = grep("chr", all_chroms, value=T)
cpg_table = cpg_table[which(cpg_table[,1] %in% all_chroms),]

for(i in seq_along(combinations[,1])){

  dmrs = get_dmrs(condition_1 = combinations$Var1[i],
                  condition_2 = combinations$Var2[i],
                  sample_to_treat = sample_to_treat,
                  cpg_table = cpg_table,
                  n_cores = n_cores,
                  condition = condition)

  colnames(dmrs) = gsub("group1", combinations$Var1[i], colnames(dmrs))
  colnames(dmrs) = gsub("group2", combinations$Var2[i], colnames(dmrs))


  write_csv(dmrs, paste0(outfile, "_", combinations$ordered[i], ".csv"))

}
