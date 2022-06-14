library(data.table)
library(stringr)
library(tidyverse)
library(Repitools)
library(GenomicRanges)
library(halpme)
library(aaRon)
library(Gviz)
library(ggbio)


if(interactive()){
  source("~/Projects/UTAS/NOMe/scripts/findNDRs2.R")
  source("~/Projects/UTAS/NOMe_analysis/Rscripts/diff_nucleosome_functions.R")
  source("~/Projects/UTAS/NOMe_analysis/Rscripts/nucleosome_functions.R")
  sample_to_treat = read.delim("lowqual/sample_info_nome.tsv")
  
  gch_file = paste0("lowqual/NOMe.dedup.gch.chr21.bed.gz")
  chrom = "chr21"
  outfile = "NDRS_chr21.rdata"
  
  condition = "ad_pathology"
  condition_1 = NULL
  condition_2 = NULL
  n_cores = 1

}else{

  source("~/nome_scripts/findNDRs2.R")
  source("~/nome_scripts/diff_nucleosome_functions.R")
  source("~/nome_scripts/nucleosome_functions.R")


  library("argparse")
  parser <- ArgumentParser()
  parser$add_argument('--sample_info', type="character",
                      help='CSV/TSV file containing sample_ids and categories')
  parser$add_argument('--GCH', type="character", default = "./",
                      help='file with GCH biscuit formatted bed')
  parser$add_argument('--chrom', "-c", type="character", default="all",
                      help='chromosome to filter BED file for. [default %(default)s]')
  parser$add_argument('--output', '-o', type="character", default="NDRs.rdata",
                      help='RDATA output file name [default %(default)s]')
  parser$add_argument('--condition', type="character",
                      help='which condition from --sample_info to compare')
  parser$add_argument('--condition_1', type="character", default=NULL,
                      help='specify which factor to compare (i.e) control')
  parser$add_argument('--condition_2', type="character", default=NULL,
                      help='specify which factor to compare (i.e) case')
  

  args <- parser$parse_args()


  sample_to_treat = read.delim(args$sample_info)
  gch_file = args$GCH
  chrom = args$chrom
  outfile = args$output
  
  condition = args$condition
  condition_1 = args$condition_1
  condition_2 = args$condition_2

}

if(tolower(tools::file_ext(outfile)) != "rdata"){
  outfile = paste0(outfile, ".Rdata")
}


gch_table = fread(gch_file, data.table = F)
if(chrom %in% gch_table[,1] & length(unique(gch_table[,1]) > 1)){
  gch_table = gch_table[gch_table[,1]==chrom,]
}else if(chrom != "all"){
  message("not filtering GCH for ", chrom)
  if(!(chrom %in% gch_table[,1])) message(chrom, " not in first column")
  if( length(unique(gch_table[,1]) == 1)) message("only 1 unique value in first column: ", gch_table[1,1])
}

### NDRS in single samples
ndrs = ndrs_from_gchtable(gch_table, p.cutoff = 10)

if(is.null(condition_1) | is.null(condition_2)){
  message("no specific conditions given... making all possible comparisons...")
  combinations = make_combinations(sample_to_treat, condition, levels_to_compare=unique(sample_to_treat[,condition]))
}else{
  combinations = make_combinations(sample_to_treat, condition, levels_to_compare=c(condition_1, condition_2))
}

ndrs_merged = ndrs_by_condition(gch_table, sample_to_treat, condition=condition, levels_to_compare = unique(sample_to_treat[,condition]), p.cutoff = 10, windowBy = 10, windowWidth=140)

diff_nucl.merged = list()
diff_nucl.unmerged = list()

for(i in seq_along(combinations[,1])){
  
  diff_nucl.merged[[i]] = diff_no(gch_table, regions=ndrs_merged, sample_to_treat, condition=condition, levels_to_compare = c(combinations$Var1[i], combinations$Var2[i]), peak_span = 1000)
  diff_nucl.unmerged[[i]] = diff_no(gch_table, regions=ndrs, sample_to_treat, condition=condition, levels_to_compare = c(combinations$Var1[i], combinations$Var2[i]), peak_span = 1000)

}

diff_nucl.merged = do.call("c", unlist(diff_nucl.merged))
diff_nucl.unmerged = do.call("c", unlist(diff_nucl.unmerged))

save(ndrs, ndrs_merged, diff_nucl.merged, diff_nucl.unmerged, file=outfile)
