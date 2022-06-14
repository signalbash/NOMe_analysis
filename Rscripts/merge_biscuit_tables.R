suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(halpme)
})

library("argparse")
parser <- ArgumentParser()
parser$add_argument('--sample_info', type="character",
                    help='CSV/TSV file containing sample_ids and categories')
parser$add_argument('--directory', '-d', type="character", default = "./",
                    help='directory with single sample BED files')
parser$add_argument('--directory_patterns', type="character", default="",
                    help='list of terms in BED file names, seperated by , with no spaces')
parser$add_argument('--combine_samples',action="store_true", default=FALSE,
                    help='Combine C/T counts for all samples into one')
parser$add_argument('--file_headers',action="store_true", default=FALSE,
                    help='Do files have a header line?')
parser$add_argument('--split_size', default=0,
                    help='Process files in chunks based of number of rows from the first sample. Any value greater than 0 will turn this option on [default %(default)s]. Experimental.')
parser$add_argument('--output', '-o', type="character",
                    help='Output file name')


args <- parser$parse_args()





sample_info = args$sample_info
directory = args$directory
directory_patterns = args$directory_patterns
output = args$output
combine_samples = args$combine_samples
split_size = args$split_size
read_headers = args$file_headers



# split patterns into a vector
if(directory_patterns != ""){
  directory_patterns = unlist(str_split(directory_patterns, ","))
}
message(paste(directory_patterns, collapse = "; "))
tryCatch(
  {setwd(directory)},
  error = function(e){
    message(paste("cannot set working directory to ", directory))
    message(e)
  }
)

cpg_files = list.files(pattern=".bed")

for(p in directory_patterns){
  cpg_files = cpg_files[grep(p, cpg_files)]
}

message(paste("merging data from ", paste(cpg_files, collapse = ", ")))
cpg_samples = strv_split(cpg_files, "[.]", 1)


check_bedheads = function(bed_file){
  x = fread(bed_file, nrows = 2, data.table = F, header = F)
  if(any(x[,1] %in% c("seqnames", "chr", "chromosome", "start", "end", "strand"))){
    bedhead = T
  } else{
    bedhead = F
  }
  return(bedhead)
}

if(split_size == 0){
  for(i in seq_along(cpg_files)){
    
    if(i==1){
      
      bed_check = check_bedheads(cpg_files[i])
      if(read_headers != bed_check){
        message("It appears you have selected the wrong file_headers argument")
        if(bed_check) message(cpg_files[i], " has the first row as column names")
        if(!bed_check) message(cpg_files[i], " has no column names")
        read_headers = bed_check
        message("Now switching to correct file_headers. Please check all files are in the same format.")
      }
      
      nome_cpg = fread(cpg_files[i], data.table = F, header = F, skip = ifelse(read_headers, 1, 0))
    }else{
      #message(i)
      #nome_cpg = full_join(nome_cpg, fread(cpg_files[i], data.table = F, header = F), by=c("seqnames", "start", "end", "strand"))
      if(combine_samples == FALSE){
        nome_cpg = full_join(nome_cpg, fread(cpg_files[i], data.table = F, header = F, skip = ifelse(read_headers, 1, 0)), by=c("V1", "V2", "V3", "V4"))
      }else{
        nome_cpg = full_join(nome_cpg, fread(cpg_files[i], data.table = F, header = F, skip = ifelse(read_headers, 1, 0)), by=c('V1', "V2", "V3", "V4"))
        
        nome_cpg[,5] = replace_na(nome_cpg[,5],0) + replace_na(nome_cpg[,7],0)
        nome_cpg[,6] = replace_na(nome_cpg[,6],0) + replace_na(nome_cpg[,8],0)
        
        nome_cpg = nome_cpg[,c(1:6)]
        
        ### reorder by chrom/start
        nome_cpg[,1] = factor(nome_cpg[,1], levels = chrom_order)
        
      }
      nome_cpg = arrange(nome_cpg, V1, V2)
    }
    
    for(j in 5:ncol(nome_cpg)){
      nome_cpg[,j] = replace_na(nome_cpg[,j], 0)
    }
    
  }
  
  #for(j in 5:ncol(nome_cpg)){
  #  nome_cpg[,j] = replace_na(nome_cpg[,j], 0)
  #}
  colnames(nome_cpg) = c("seqnames", "start", "end", "strand", paste0(rep(cpg_samples, each=2), ".", c("C", "cov")))
  write_delim(nome_cpg, output,  delim = "\t", col_names = T)
}else{
  
  ##### NEED TO FIX HEADER READING AS PER VERSION OF CODE ABOVE
  
  # keep a vector with chromosome order as read-in (may not be chr 1,2,3,--> etc.)
  chrom_order = vector()
  read_lines = split_size
  ### initialise start_index_list
  start_index_list = (rep(list(1), length(cpg_files)))
  
  counter = 1
  bed_ended = F
  while(!bed_ended){

    bed_read_start = ((counter-1)*read_lines)+1
    chrom_order = vector()
    for(i in seq_along(cpg_files)){
      #message(i)
      
      ## read in n-lines of the first file...
      ## this is always constant, all other files are reduced/expanded to cover the same regions as the first file n-lines
      if(i == 1){
        nome_cpg = fread(cpg_files[i], data.table = F, header = F, skip = bed_read_start, nrows = read_lines)
        start_index_list[[i]] = bed_read_start+read_lines+1
        chrom_order = unique(c(chrom_order, unique(nome_cpg[,1])))
        final_line = tail(nome_cpg,1)
        if(nrow(nome_cpg) < read_lines){
          bed_ended = T
        }
        
        ## add table is longer
        ## remove any extra chromosomes
        max_per_chrom = aggregate(V2 ~ V1, nome_cpg, max)
        # hard core max length to 1000M for all full chroms
        max_per_chrom$V2[max_per_chrom$V1 %in% chrom_order[-length(chrom_order)]] = 1e9
        
      }else{
        #message(start_index_list[[i]])
        #message(cpg_files[i])
        nome_add = fread(cpg_files[i], data.table = F, header = F, skip = start_index_list[[i]], nrows = read_lines)
        ## chrom order as read-in
        chrom_order = unique(c(chrom_order, unique(nome_add[,1])))
        
        if(nrow(nome_add) < read_lines){
          second_bed_ended = T
        }else{
          second_bed_ended = F
        }
        
        ## read more until catches up
        extra_counter = 2
        while(tail(nome_add[,2], 1) < final_line[,2] & match(tail(nome_add[,1],1), chrom_order) <= match(final_line[,1], chrom_order) & second_bed_ended==F){
          
          #bed_read_start_extra = ((extra_counter-1)*read_lines)+1
          bed_read_start_extra = start_index_list[[i]]+((extra_counter-1)*read_lines)+1
          
          nome_add = rbind(nome_add,
                           fread(cpg_files[i], data.table = F, header = F, skip = bed_read_start_extra, nrows = read_lines))
          
          chrom_order = unique(c(chrom_order, unique(nome_add[,1])))
          
          extra_counter = extra_counter+1
        }
        
        if(second_bed_ended == F){
          keep = which(nome_add[,2] <= max_per_chrom[,2][match(nome_add[,1], max_per_chrom[,1])])
          start_index_list[[i]] = start_index_list[[i]] + max(keep)
          nome_add = nome_add[keep,]
        }
        
        if(combine_samples == FALSE){
          nome_cpg = full_join(nome_cpg, nome_add, by=c("V1", "V2", "V3", "V4"))
        }else{
          nome_cpg = full_join(nome_cpg, nome_add, by=c('V1', "V2", "V3", "V4"))
          
          nome_cpg[,5] = replace_na(nome_cpg[,5],0) + replace_na(nome_cpg[,7],0)
          nome_cpg[,6] = replace_na(nome_cpg[,6],0) + replace_na(nome_cpg[,8],0)
          
          nome_cpg = nome_cpg[,c(1:6)]
          
          ### reorder by chrom/start
          nome_cpg[,1] = factor(nome_cpg[,1], levels = chrom_order)
          nome_cpg = arrange(nome_cpg, V1, V2)
        }
        
      }
      
    }
    
    
    if(counter == 1){
      write_delim(nome_cpg, output,  delim = "\t", col_names = F)
    }else{
      write_delim(nome_cpg, output,  delim = "\t", append = T, col_names = F)
    }
    counter = counter+1
  }
}




