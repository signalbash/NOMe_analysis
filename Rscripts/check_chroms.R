library(stringr)

args = commandArgs(trailingOnly=TRUE)

file_chrs = read.delim(args[1], header = F)[,1]
file_chrs = unlist(stringr::str_split(file_chrs, "chr"))
file_chrs = file_chrs[file_chrs!=""]
file_chrs = paste0("chr", file_chrs)


all_chrs = paste0("chr", c(1:22, "X","Y"))

if(any((all_chrs %in% file_chrs) == FALSE)){

  missing = all_chrs[!(all_chrs %in% file_chrs)]
  message("missing:")
  for(i in seq_along(missing)){
    message(missing[i])
  }

}else{
  message("all chroms found in file")
}
