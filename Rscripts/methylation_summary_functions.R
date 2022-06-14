convert_biscuit2methylkit = function(nome_cpg, min_cov=3, out_folder="./methylkit"){

  system(command = paste0("mkdir -p ", out_folder))

  samps <- colnames(nome_cpg) %>% grep("[.]cov",., value =T) %>% gsub("[.]cov", "", .)

  for(i in seq_along(samps)){
    keep_cols = c(grep(samps[i], colnames(nome_cpg)))

    nome_cpg.p = nome_cpg[,c(1,2,4, keep_cols)]
    nome_cpg.p = nome_cpg.p[nome_cpg.p[,5] >= min_cov,]

    colnames(nome_cpg.p) = c("chr", "base", "strand_sign",  "c_count","coverage")
    nome_cpg.p$strand = "F"
    nome_cpg.p$strand[nome_cpg.p$strand_sign == "-"] = "R"
    nome_cpg.p$freqC = (nome_cpg.p$c_count/nome_cpg.p$coverage)*100
    nome_cpg.p$freqT = 100 - nome_cpg.p$freqC
    nome_cpg.p$freqC = format(round(nome_cpg.p$freqC, 2), nsmall=2)
    nome_cpg.p$freqT = format(round(nome_cpg.p$freqT, 2), nsmall=2)
    nome_cpg.p$strand_sign = nome_cpg.p$strand
    nome_cpg.p$strand = NULL
    colnames(nome_cpg.p)[3] = "strand"
    nome_cpg.p$c_count= NULL
    nome_cpg.p = cbind(chrBase = paste(nome_cpg.p$chr, nome_cpg.p$base, sep="."), nome_cpg.p)
    write_delim(nome_cpg.p, paste0(out_folder, "/", samps[i], "_methylkit.txt.gz"), delim='\t')
  }
}

total_meth_inoverlap = function(meth_table, overlap_range, progress=T, min_window_cov=10){

  colnames(meth_table)[1:4] = c("seqnames", "start", "end", "strand")
  x = makeGRangesFromDataFrame(meth_table, keep.extra.columns = T)


  #cpg.grange = makeGRangesFromDataFrame(meth_table,keep.extra.columns = T)
  #ol.prom = findOverlaps(prom_locs, cpg.grange)
  cov_cols = colnames(meth_table)[which(str_sub(colnames(meth_table), -4,-1) == ".cov")]
  c_cols = gsub(".cov", ".C", cov_cols)

  pb = progress_bar$new(format = "Getting window coverage for samples [:bar] :percent eta: :eta",
                        clear = FALSE, total = length(cov_cols), width = 60)

  for(i in seq_along(cov_cols)){
    if(progress==T) pb$tick(1)
    #message("processing ", cov_cols[i])
    overlap_range$cov = overlapSums(x, overlap_range, values(x)[[cov_cols[i]]])
    overlap_range = overlap_range[which(overlap_range$cov >= min_window_cov)]
    colnames(elementMetadata(overlap_range))[length(colnames(elementMetadata(overlap_range)))] = cov_cols[i]
  }

  if(length(overlap_range)==0){
    warning("No overlap_range with sufficient coverage, returning NULL")
    return(NULL)
  }else{
    ## calculate meth_percent in remaining covered overlap_range
    pb = progress_bar$new(format = "Getting window methylation for samples [:bar] :percent eta: :eta",
                          clear = FALSE, total = length(cov_cols), width = 60)

    for(i in seq_along(cov_cols)){
      if(progress==T) pb$tick(1)
      cov = overlapSums(x, overlap_range, values(x)[[cov_cols[i]]])
      C = overlapSums(x, overlap_range, values(x)[[c_cols[i]]])
      elementMetadata(overlap_range)[which(colnames(elementMetadata(overlap_range)) == cov_cols[i])] = C/cov
      colnames(elementMetadata(overlap_range))[which(colnames(elementMetadata(overlap_range)) == cov_cols[i])] = gsub(".cov", ".meth_percent", cov_cols[i])
    }

    return(overlap_range)
  }




}



library(progress)


## Convert tables to reduced windows with high variance.
## Used for PCAs on CpG sites and GpC sites

convert_biscuit2methwindows = function(meth_table, window_width=300, window_spacing=50, min_window_cov=1, progress=T){

  # original order of samples
  cov_cols_original = colnames(meth_table)[which(str_sub(colnames(meth_table), -4,-1) == ".cov")]

  # find sample with lowest coverage and start from there for making windows.
  lowest_cov = apply(meth_table[1:min(10e6, nrow(meth_table)),cov_cols_original], 2, sum)
  order_by_cov = stack(lowest_cov) %>% arrange(values)

  ### Make windows covering genome
  colnames(meth_table)[1:4] = c("seqnames", "start", "end", "strand")


  x = makeGRangesFromDataFrame(meth_table, keep.extra.columns = T)
  # calculate 'seqlengths' based on max end values for each chrom processed
  seqlengths(x) =
    aggregate(end~seqnames, meth_table, max)%>%
      mutate(seqnames=factor(seqnames, levels = names(seqlengths(x)))) %>%
      arrange(seqnames) %>% pull(end, seqnames)
  # x must have seqlengths defined
  stopifnot(all(!is.na(seqlengths(x))))
  windows <- genomeBlocks(seqlengths(x), width=window_width, spacing=window_spacing)
  # Only keep windows that overlap x
  windows <- windows[windows %over% x]

  cov_cols = as.character(order_by_cov$ind)
  c_cols = gsub(".cov", ".C", cov_cols)

  ## Calculate coverage of each window, by sample (lowest cov first)
  ## and filter for windows with >= minimum coverage
  ## so all windows should have coverage in all samples
  pb = progress_bar$new(format = "Getting window coverage for samples [:bar] :percent eta: :eta",
                        clear = FALSE, total = length(cov_cols), width = 60)

  for(i in seq_along(cov_cols)){
    if(progress==T) pb$tick(1)
    #message("processing ", cov_cols[i])
    windows$cov = overlapSums(x, windows, values(x)[[cov_cols[i]]])
    windows = windows[windows$cov >= min_window_cov]
    colnames(elementMetadata(windows))[length(colnames(elementMetadata(windows)))] = cov_cols[i]
  }

  if(length(windows)==0){
    warning("No windows with sufficient coverage, returning NULL")
    return(NULL)
  }else{
    ## calculate meth_percent in remaining covered windows
    pb = progress_bar$new(format = "Getting window methylation for samples [:bar] :percent eta: :eta",
                          clear = FALSE, total = length(cov_cols), width = 60)

    for(i in seq_along(cov_cols)){
      if(progress==T) pb$tick(1)
      cov = overlapSums(x, windows, values(x)[[cov_cols[i]]])
      C = overlapSums(x, windows, values(x)[[c_cols[i]]])
      elementMetadata(windows)[which(colnames(elementMetadata(windows)) == cov_cols[i])] = C/cov
      colnames(elementMetadata(windows))[which(colnames(elementMetadata(windows)) == cov_cols[i])] = gsub(".cov", ".meth_percent", cov_cols[i])
    }

    ## reduce overlapping windows by finding window in range with highest variance
    windows_reduced = reduce_windows_by_maxvar(windows)
    elementMetadata(windows_reduced) = elementMetadata(windows_reduced)[,match(cov_cols_original, cov_cols)]

    return(windows_reduced)
  }

}

reduce_windows_by_maxvar = function(windows){

  ## reduce overlapping windows by finding window in range with highest variance
  meth_percent = as.data.frame(elementMetadata(windows))

  smol_windows = reduce(windows)
  smol_index = as.data.frame(findOverlaps(windows, smol_windows))

  meth_percent$var = apply(meth_percent, 1 ,var)
  meth_percent$smol_index = smol_index$subjectHits
  meth_percent$index = smol_index$queryHits
  max_var = aggregate(var ~ smol_index, meth_percent, max)

  # needs distinct() to get rid of multi-matching smol_index and var.
  max_var_join = left_join(max_var, meth_percent, by=c("smol_index", "var")) %>%
    distinct(smol_index, var, .keep_all = TRUE)

  windows_reduced = windows[max_var_join$index]

  return(windows_reduced)

}

convert_biscuit2methwindows_chunk = function(bed_file, chunk_size=1e7, buffer=500, verbose_rows=F,
                                             window_width = 200, window_spacing = 20, min_window_cov = 1){
  file_ended = FALSE
  counter = 1
  while(!file_ended){
    # read in file in chunks
    file_read_start = ((counter-1)*chunk_size)
    if(counter == 1){
      nome_gpc = fread(bed_file, data.table = F, header = T, skip = file_read_start, nrows = (chunk_size-1)+buffer)
      meth_head = colnames(nome_gpc)
      # convert to non-overlapping windows (picking the window with the highest variance)
      gpc_windows = convert_biscuit2methwindows(nome_gpc, window_width = window_width,
                                                window_spacing = window_spacing, min_window_cov = min_window_cov)

    }else{
      nome_gpc = fread(gch_file, data.table = F, header = F, skip = file_read_start, nrows = (chunk_size-1)+buffer)
      colnames(nome_gpc) = meth_head
      # convert to non-overlapping windows (picking the window with the highest variance)
      gpc_windows2 = convert_biscuit2methwindows(nome_gpc, window_width = window_width,
                                                 window_spacing = window_spacing, min_window_cov = min_window_cov)
      # combine with previous chunks, same merging technique as used internally above.
      gpc_windows = reduce_windows_by_maxvar(c(gpc_windows, gpc_windows2))
    }
    # if less rows read than the chunk size+buffer, the file has ended.
    if(nrow(nome_gpc) < ((chunk_size-1)+buffer)){
      file_ended = T
    }
    counter = counter+1
    if(verbose_rows) message(paste0("read rows ", as.integer(file_read_start), " to ", as.integer(file_read_start + chunk_size)))

  }
  return(gpc_windows)
}


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

gch_summary_centered = function(gch_file, centers=TSS, centers_buffer=1000, chunk_size=1e6, verbose_rows = F){
  start(centers) = start(centers)-centers_buffer
  end(centers) = end(centers)+(centers_buffer)
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
