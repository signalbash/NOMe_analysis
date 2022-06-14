library(stringr)
library(tidyverse)
library(GenomicRanges)
library(halpme)
library(data.table)
library(zoo)
library(ggpmisc)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


gene_annotation = rtracklayer::import("~/Downloads/gencode.v39.annotation.gtf.gz")
TSS = promoters(gene_annotation[gene_annotation$type=="transcript" & gene_annotation$transcript_type == "protein_coding" & gene_annotation$transcript_support_level %in% c("1","2")], 0, 1)

gpc_files = list.files("biscuit/NOME_bychr/", full.names = T)

for(i in seq_along(gpc_files)){
  if(i == 1) first_window = T
  nome_gpc = fread(gpc_files[i], data.table = F)
  
  try({
    gpc_windows = convert_biscuit2methwindows(nome_gpc, window_width = 200, window_spacing = 20, min_window_cov = 1)
    
    if(first_window == T){
      gpc_windows_all = gpc_windows
      first_window = F
    }else{
      suppressWarnings({
        gpc_windows_all = c(gpc_windows_all, gpc_windows)
      })
    }
  })
}


#hcg = fread("merged_data/NOMe_merged.hcg.black.merged.bed", data.table = F)
#hcg.grange = makeGRangesFromDataFrame(hcg)

buffer = 1000
centers = TSS


gch_file = gpc_files[1]

merge_meth_sums = function(meth.x, meth.y){
  col_match = match(colnames(meth.x), colnames(meth.y))
  col_match = col_match[colnames(meth.x) != "centered_dist"]
  meth.x[,col_match] = meth.x[,col_match] + meth.y[,col_match]
  return(meth.x)
}

merge_meth_list = function(meth_list.x, meth_list.y){
  for(i in 1:length(meth_list.x)){
    meth_list.x[[i]] = merge_meth_sums(meth_list.x[[i]], meth_list.y[[i]])
  }
  return(meth_list.x)
}

gch_summary_centered = function(gch_file, centers=TSS, buffer=1000, chunk_size=1e6, verbose_rows = T){
  start(centers) = start(centers)-buffer
  end(centers) = end(centers)+(buffer)
  file_ended = FALSE
  counter = 1
  while(!file_ended){
    
    # read in file in chunks
    file_read_start = ((counter-1)*chunk_size)
    if(counter == 1){
      C_sums = NULL
      cov_sums = NULL
      meth = fread(gch_file, data.table = F, header = T, skip = file_read_start, nrows = chunk_size-1, showProgress = T)
      meth_head = colnames(meth)
    }else{
      meth = fread(gch_file, data.table = F, header = F, skip = file_read_start, nrows = chunk_size-1)
      colnames(meth) = meth_head
    }
    if(nrow(meth) < (chunk_size-1)){
      file_ended = T
    }
    
    # convert to grange
    # find distance to TSS
    meth.grange = makeGRangesFromDataFrame(meth)
    ol = findOverlaps(centers, meth.grange, ignore.strand=T)
    meth = meth[(ol@to),]
    meth$distance_to_center = distance(TSS[ol@from], meth.grange[ol@to])
    meth$TSS_center = start(TSS[ol@from])
    # flip for stranded data
    meth$start_is_greater = meth$start+1 > start(TSS[ol@from])
    meth$centered_dist = meth$distance_to_center
    flip = which((meth$strand == "+" & meth$start_is_greater == FALSE) | 
                   (meth$strand == "-" & meth$start_is_greater == TRUE) )
    meth$centered_dist[flip] =  meth$centered_dist[flip]*-1
    # filter out anything too far away
    meth = meth[meth$distance_to_center <= 1000,]
    
    meth= arrange(meth, centered_dist)
    
    # C sums for each sample normalised to TSS center
    C_sums = meth %>% dplyr::select(ends_with(c(".C", "centered_dist"))) %>% 
      plyr::rbind.fill(C_sums) %>% 
      aggregate(.~centered_dist, ., sum)
    # cov sums for each sample normalised to TSS center
    cov_sums = meth %>% dplyr::select(ends_with(c(".cov", "centered_dist"))) %>% 
      plyr::rbind.fill(cov_sums) %>% 
      aggregate(.~centered_dist, ., sum)
    
    if(verbose_rows) message(paste0("read rows ", as.integer(file_read_start), " to ", as.integer(file_read_start + chunk_size)))
    
    if(counter == 1){
      C_sums_all = C_sums
      cov_sums_all = cov_sums
    }else{
      C_sums_all = merge_meth_sums(C_sums_all, C_sums)
      cov_sums_all = merge_meth_sums(cov_sums_all, cov_sums)
    }
    counter = counter+1
  }
  return(list(C_sums_all, cov_sums_all))
}

for(i in seq_along(gpc_files)){
  
  if(i== 1){
    meth_list = gch_summary_centered(gpc_files[i], centers=TSS, buffer=1000, verbose_rows = F)
  }else{
    meth_list = merge_meth_list(meth_list, gch_summary_centered(gpc_files[i], centers=TSS, buffer=1000, verbose_rows = F))
  }
  message(gpc_files[i])
}



  M = as.matrix(meth_list[[1]][,-which(colnames(meth_list[[1]]) =="centered_dist")])
  Cov = as.matrix(meth_list[[2]][,-which(colnames(meth_list[[2]]) =="centered_dist")])
  
  sampleNames = gsub("[.]C","", colnames(M))
  methyl_percents = M/Cov
  colnames(methyl_percents) = sampleNames
  methyl_percents = as.data.frame(methyl_percents)
  methyl_percents$distance_to_center = C_sums$centered_dist
  
  long_meth = pivot_longer(methyl_percents, !distance_to_center, names_to="sample", values_to="meth_percent")
  long_meth.df = as.data.frame(long_meth)
  long_meth.df$distance_to_center = as.numeric(as.character(long_meth.df$distance_to_center))
  
  arrange(long_meth.df, distance_to_center, sample) %>% 
    filter(sample=="01-31_NOMe") %>% 
    mutate(ma2 =rollapply(meth_percent, 100, mean, align="right", fill=NA)) %>% 
    ggplot(aes(x=distance_to_center, y=ma2, col=sample, group=sample)) + geom_line()
  #+ stat_valleys(span=40)
  
  
  arrange(long_meth.df, distance_to_center, sample) %>% 
    filter(sample=="01-31_NOMe") -> test

 #ggpmisc::find_peaks()
 
 long_meth.df[splus2R::peaks(long_meth.df$meth_percent, span=40),]
  
}





write_delim(C_sums, "Csums_gch_TSS_centered.txt")
write_delim(cov_sums, "covsums_gch_TSS_centered.txt")