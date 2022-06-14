add_max_seqlengths = function(x){
  seqlengths(x) =
    aggregate(end~seqnames, as.data.frame(x)[,c(1:5)], max)%>%
     mutate(seqnames=factor(seqnames, levels = names(seqlengths(x)))) %>%
    arrange(seqnames) %>% pull(end, seqnames)
  return(x)
}

findNDRs <- function(x, samples, p.cutoff=5, windowWidth=100, windowBy=20, minSize=140, verbose=T) {
  stopifnot(all(c("Sample", "C", "cov") %in% colnames(samples)))
  fast.chisq <- compiler::cmpfun(function(x, p) {
    n <- rowSums(x)
    E <- cbind(n * p[1], n * p[2])
    STATISTIC <- rowSums((x - E)^2/E)
    PARAMETER <- 1
    PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
    return(PVAL)
  })

  # x must have seqlengths defined
  stopifnot(all(!is.na(seqlengths(x))))
  # Rounding for nicer ranges ending on a multiple of windowBy (i.e. every 10/20, instead of ranges ending in 3??)
  plyr::round_any(seqlengths(x), windowBy, f=ceiling)
  # Define windows for pvalue calculation
  windows <- genomeBlocks(seqlengths(x), width=windowWidth, spacing=windowBy)
  # Only keep windows that overlap x
  windows <- windows[windows %over% x]
  windows = windows[width(windows) == windowWidth]
  shift_up = plyr::round_any(end(windows)[1], f=ceiling, windowBy) - end(windows)[1]
  windows = GenomicRanges::shift(windows, shift_up)

  # Find NDRs
  NDRs <- GRangesList(lapply(1:nrow(samples), function(i) {


    message(paste0("Processing ", samples$Sample[i]))
    # Counts of Cs and Ts in each window
    if(verbose) message(" - Counting Cs")
    windows.counts <- cbind("C"=overlapSums(x, windows, values(x)[[samples$C[i]]]))
    if(verbose) message(" - Counting Ts")
    windows.counts <- cbind(windows.counts, "T"=overlapSums(x, windows, values(x)[[samples$cov[i]]])-windows.counts[,"C"])

    # Genomic background ratio to compare to
    Csum <- sum(values(x)[[samples$C[i]]])
    Tsum <- sum(values(x)[[samples$cov[i]]])-sum(values(x)[[samples$C[i]]])
    C.T <- c(Csum,Tsum)/(Csum+Tsum)

    # calculate chisq pvalue
    if(verbose) message(" - Calculating chisq pvalues")
    windows.pvals <- -log10(fast.chisq(windows.counts, C.T))
    windows.pvals[windows.pvals==Inf] <- max(windows.pvals[!windows.pvals==Inf])

    # find regions that meet the cutoff on pvalue and size
    if(verbose) message(" - Finding significant regions")

    windows.p = windows[which(windows.pvals>p.cutoff)]
    windows.p$pvals = windows.pvals[which(windows.pvals>p.cutoff)]
    windows.cutoff <- GenomicRanges::reduce(windows.p)

    # Counts of Cs and Ts in each window
    windows.counts <- cbind("C"=overlapSums(x, windows.cutoff, values(x)[[samples$C[i]]]))
    windows.counts <- cbind(windows.counts, "T"=overlapSums(x, windows.cutoff, values(x)[[samples$cov[i]]])-windows.counts[,"C"])
    windows.cutoff$C = windows.counts[,1]
    windows.cutoff$T = windows.counts[,2]
    windows.cutoff$C.T.diff = C.T[1] - (windows.cutoff$C / (windows.cutoff$C+windows.cutoff$T))
    windows.cutoff$direction = ifelse(windows.cutoff$C.T.diff > 0 , "nucleosome_depleted", "nucleosome_enriched")
    windows.cutoff$pval_log10<- -log10(fast.chisq(windows.counts, C.T))
    windows.cutoff$pval_log10[windows.cutoff$pval_log10==Inf] <- max(windows.cutoff$pval_log10[!windows.cutoff$pval_log10==Inf])
    windows.cutoff <- windows.cutoff[width(windows.cutoff)>=minSize]


    windows.cutoff$sample = samples$Sample[i]
    return(windows.cutoff)

  }))

  names(NDRs) <- samples$Sample
  seqlevels(NDRs) <- seqlevels(x)
  seqlengths(NDRs) <- seqlengths(x)
  NDRs
}

ndrs_from_gchtable = function(gch_table, p.cutoff = 10, windowBy = 10, windowWidth=140 ){

  samps <- colnames(gch_table) %>% grep("[.]cov",., value =T) %>% gsub("[.]cov", "", .)
  samples <- data.frame(Sample=samps,
                        C=paste0(samps, ".C"),
                        cov=paste0(samps, ".cov"), stringsAsFactors=FALSE, row.names=samps)

  gch_grange = makeGRangesFromDataFrame(gch_table, keep.extra.columns = TRUE)
  gch_grange = add_max_seqlengths(gch_grange)

  NDRs <- findNDRs(gch_grange, samples, p.cutoff = p.cutoff,
                           windowBy = windowBy, windowWidth=windowWidth, verbose = F)


  #### All NDRs by SAMPLE (not condition specific)
  NDRs = unlist(NDRs)
  return(NDRs)

}

make_combinations = function(sample_to_treat, condition, levels_to_compare){
  sample_to_treat = make_condition_levels(sample_to_treat, condition)

  combinations = expand.grid(levels_to_compare, levels_to_compare)
  combinations = combinations[combinations[,1]!=combinations[,2],]
  combinations$ordered = apply(combinations, 1, function(x) (paste0(sort(x), collapse = "_vs_")))
  combinations = combinations[combinations$Var1 == gsub("(_vs_).*", "", combinations$ordered),]
  combinations = combinations[!duplicated(combinations$ordered),]
  
  return(combinations)
}


ndrs_by_condition = function(gch_table, sample_to_treat, condition, levels_to_compare,
                             p.cutoff = 10, windowBy = 10, windowWidth=140){


  samps <- colnames(gch_table) %>% grep("[.]cov",., value =T) %>% gsub("[.]cov", "", .)
  samples <- data.frame(Sample=samps,
                        C=paste0(samps, ".C"),
                        cov=paste0(samps, ".cov"), stringsAsFactors=FALSE, row.names=samps)
  gch_grange = makeGRangesFromDataFrame(gch_table, keep.extra.columns = TRUE)
  gch_grange = add_max_seqlengths(gch_grange)

  sample_to_treat = make_condition_levels(sample_to_treat, condition)
  combinations = make_combinations(sample_to_treat, condition, levels_to_compare)

  ## make NDRs for each cat.
  ndr_factors = unique(c(combinations$Var1, combinations$Var2))
  merged_table = gch_table[,c(1:4)]
  for(i in seq_along(ndr_factors)){
    merged_samples = sample_to_treat$sample_id[sample_to_treat[,paste0(condition, ".factor")] == ndr_factors[i]]
    merged_samples = merged_samples[which(merged_samples %in% samples$Sample)]
    merged_table$total.C = apply(gch_table[,which(colnames(gch_table) %in% samples$C[samples$Sample %in% merged_samples])], 1, sum)
    merged_table$total.cov = apply(gch_table[,which(colnames(gch_table) %in% samples$cov[samples$Sample %in% merged_samples])], 1, sum)
    colnames(merged_table) = str_replace_all(colnames(merged_table), "total", as.character(ndr_factors[i]))
  }
  merged_ndrs = ndrs_from_gchtable(merged_table, p.cutoff = p.cutoff, windowBy = windowBy, windowWidth = windowWidth)
  return(merged_ndrs)

}


diff_no = function(gch_table, regions, sample_to_treat, condition, levels_to_compare, peak_span=10000){


  ## setup & formatting
  samps <- colnames(gch_table) %>% grep("[.]cov",., value =T) %>% gsub("[.]cov", "", .)
  samples <- data.frame(Sample=samps,
                        C=paste0(samps, ".C"),
                        cov=paste0(samps, ".cov"), stringsAsFactors=FALSE, row.names=samps)
  gch_grange = makeGRangesFromDataFrame(gch_table, keep.extra.columns = TRUE)
  gch_grange = add_max_seqlengths(gch_grange)

  sample_to_treat = make_condition_levels(sample_to_treat, condition)
  combinations = make_combinations(sample_to_treat, condition, levels_to_compare)

  ## Get significant regions
  signif_dnd = lapply(1:nrow(combinations), function(i){

    case_samples = sample_to_treat$sample_id[sample_to_treat[,paste0(condition, ".factor")] == combinations$Var1[i]]
    case_samples = case_samples[which(case_samples %in% samples$Sample)]

    control_samples = sample_to_treat$sample_id[sample_to_treat[,paste0(condition, ".factor")] == combinations$Var2[i]]
    control_samples = control_samples[which(control_samples %in% samples$Sample)]

    windows = range2window(gch_grange, window_width = 140, window_spacing = 20,
                           case = case_samples, control=control_samples,
                           signif.range = regions, extend = 100)

    dnd = methylsig_dnd(windows, peak_span = peak_span, case_samples = case_samples, control_samples = control_samples)
    colnames(elementMetadata(dnd))[which(colnames(elementMetadata(dnd)) == "meth_case")] = paste0("meth_", combinations$Var1[i])
    colnames(elementMetadata(dnd))[which(colnames(elementMetadata(dnd)) == "meth_control")] = paste0("meth_", combinations$Var2[i])
    dnd$comp = combinations$ordered[i]
    return(dnd)
  })

  return(signif_dnd)

}

get_dmrs = function(condition_1=0, 
                    condition_2=2, 
                    sample_to_treat,
                    cpg_table,
                    n_cores = 1,
                    condition = "ad_pathology",
                    min_coverage = 3,
                    min_samples_per_condition = 3,
                    qcutoff = c(0.01, 0.99),
                    meth_diff = 0.1){
  
  
  ## Convert to Factor levels (0/1/2)
  sample_to_treat = sample_to_treat %>%
    make_condition_levels(condition)
  treat_col = paste0(condition, ".factor")
  
  # filter for samples in each comparison group
  cpg_table.filtered = cpg_table[,tolower(colnames(cpg_table)) %in% tolower(c("seqnames", "start", "strand",
                                                                              paste0(sample_to_treat$sample_id[tolower(sample_to_treat[,treat_col])==tolower(condition_1)], ".c"),
                                                                              paste0(sample_to_treat$sample_id[tolower(sample_to_treat[,treat_col])==tolower(condition_1)], ".cov"),
                                                                              paste0(sample_to_treat$sample_id[tolower(sample_to_treat[,treat_col])==tolower(condition_2)], ".c"),
                                                                              paste0(sample_to_treat$sample_id[tolower(sample_to_treat[,treat_col])==tolower(condition_2)], ".cov")
  ))]
  
  
  cpg_table.filtered$seqnames = factor(cpg_table.filtered$seqnames, levels = paste0("chr", c(1:22, "X", "Y", "M")))
  cpg_table.filtered = arrange(cpg_table.filtered, seqnames, start)
  
  M = as.matrix(cpg_table.filtered[,grep("[.]C", colnames(cpg_table.filtered))])
  Cov = as.matrix(cpg_table.filtered[,grep("[.]cov", colnames(cpg_table.filtered))])
  
  sampleNames = gsub("[.]C","", colnames(M))
  sampleConditions = sample_to_treat[,treat_col][match(sampleNames, sample_to_treat$sample_id)]
  
  colnames(M) = NULL
  colnames(Cov) = NULL
  
  ## Create BSseq object
  if(any(cpg_table.filtered$strand == "-")){
    ## merge strand data if any exists for the - strand
    BS.cpg.f = BSseq(chr = cpg_table.filtered$seqnames[cpg_table.filtered$strand=="+"],
                     pos = cpg_table.filtered$start[cpg_table.filtered$strand=="+"],
                     M = M[cpg_table.filtered$strand=="+",], Cov=Cov[cpg_table.filtered$strand=="+",],
                     sampleNames = sampleNames)
    
    BS.cpg.r = BSseq(chr = cpg_table.filtered$seqnames[cpg_table.filtered$strand=="-"],
                     pos = cpg_table.filtered$start[cpg_table.filtered$strand=="-"] - 1L,
                     M = M[cpg_table.filtered$strand=="-",], Cov=Cov[cpg_table.filtered$strand=="-",],
                     sampleNames = sampleNames)
    
    BS.cpg <- combine(BS.cpg.f, BS.cpg.r)
    BS.cpg <- collapseBSseq(BS.cpg, group = c("a", "a"), type = "integer")
  }else{
    BS.cpg = BSseq(chr = cpg_table.filtered$seqnames,
                   pos = cpg_table.filtered$start,
                   M = M, Cov=Cov,
                   sampleNames = sampleNames)
  }
  # remove objects covered within BS.cpg... these can be very large tables...
  rm(M, Cov, cpg_table.filtered)
  gc()
  
  
  pData(BS.cpg)$Condition = sampleConditions
  
  ## Create BSmooth object
  ## TIME WARNING... This takes a very long time, which is why we save the data
  if(!file.exists(paste0(condition, ".", condition_1, "_", condition_2, "bsmooth.rdata"))){
    BS.cpg.fit <- BSmooth(
      BSseq = BS.cpg,
      BPPARAM = MulticoreParam(workers = n_cores),
      verbose = TRUE)
    # save the data in case we want to change something downstream later.. 
    save(BS.cpg.fit, file=paste0(condition, ".", condition_1, "_", condition_2, "bsmooth.rdata"))
  }else{
    load(paste0(condition, ".", condition_1, "_", condition_2, "bsmooth.rdata"))
  }
  
  
  sam_names=names(round(colMeans(getCoverage(BS.cpg)), 1))
  
  # filter for min coverage
  BS.cov <- getCoverage(BS.cpg.fit)
  keepLoci.ex <- which(rowSums(BS.cov[, BS.cpg.fit$Condition == condition_1] >= min_coverage) >= min_samples_per_condition &
                         rowSums(BS.cov[, BS.cpg.fit$Condition == condition_2] >= min_coverage) >= min_samples_per_condition)
  BS.cpg.fit.filter <- BS.cpg.fit[keepLoci.ex,]
  
  # make groups for stat calculation
  g1 = (sample_to_treat$sample_id)[sample_to_treat[,treat_col] == condition_1]
  g1 = g1[g1 %in% sampleNames(BS.cpg.fit.filter)]
  g2 = (sample_to_treat$sample_id)[sample_to_treat[,treat_col] == condition_2]
  g2 = g2[g2 %in% sampleNames(BS.cpg.fit.filter)]
  
  BS.cpg.tstat <- BSmooth.tstat(BS.cpg.fit.filter,
                                group1 = g1,
                                group2 = g2,
                                estimate.var = "group2",
                                local.correct = TRUE,
                                verbose = TRUE)
  
  # get significant regions
  ps <- getStats(BS.cpg.tstat)[,"tstat.corrected"]
  dmrs3 = dmrFinder(BS.cpg.tstat, qcutoff = qcutoff)
  dmrs <- subset(dmrs3, n >= 3 & abs(meanDiff) >= meth_diff)
  
  # sort into numerical chrom order (not by so 1-2-3 not 1-10-11 etc.)
  dmrs$chr_n = factor(gsub("chr", "", dmrs$chr), levels=c(1:22, "X", "Y", "M"))
  dmrs = arrange(dmrs, chr_n, start)
  dmrs$chr_n = NULL
  
  return(dmrs)
}

make_condition_levels = function(sample_to_treat, condition){
  
  if(is.numeric(sample_to_treat[,condition])){
    
    sample_to_treat = sample_to_treat[!duplicated(sample_to_treat$sample_number),]
    sample_to_treat = sample_to_treat[!(is.na(sample_to_treat[,condition])),]
    sample_to_treat[,condition] = as.numeric(as.character(sample_to_treat[,condition]))
    sample_to_treat = sample_to_treat[order(sample_to_treat[,condition]),]
    order(sample_to_treat[,condition])
    
    group_sizes = ceiling(nrow(sample_to_treat)/3)
    
    sample_to_treat$new_levels = 1
    sample_to_treat$new_levels[1:group_sizes] = 0
    sample_to_treat$new_levels[(nrow(sample_to_treat)-group_sizes+1):nrow(sample_to_treat) ] = 2
    
    
  }else{
    sample_to_treat$new_levels = sample_to_treat[,condition]
  }
  colnames(sample_to_treat)[which(colnames(sample_to_treat) == "new_levels")] = paste0(condition, ".factor")
  
  return(sample_to_treat)
}
