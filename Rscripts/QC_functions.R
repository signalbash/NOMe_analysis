### Functions for reading in FASTQC files


sub_seqids = function(names, sample_info){
  for(i in 1:nrow(sample_info)){
    names = gsub(sample_info$sequencing_id[i], sample_info$sample_id[i], names)
  }
  return(names)
}

check_sample_hasdata = function(file_vector, sample_list = sample_info$sample_id){
  sample_vector = vector()
  for(i in seq_along(file_vector)){
    w = which(unlist(lapply(sample_list, function(x) grepl(x, file_vector[i]))))
    if(length(w) > 0) sample_vector[i] =sample_list[w]
  }

  missing_samples = sample_list[which(!sample_list %in% sample_vector)]
  if(length(missing_samples) >0) message(paste0("missing data for: ", paste(missing_samples, collapse = ", ")))
}


read_in_fastqc_stats = function(directory, sample_info, multiqc_bin_location){

  # read in multiqc fileand/or create if not already done
  fastqc_files = list.files(directory)
  if(all(grepl("multiqc", fastqc_files) == FALSE)){
    system(paste0(multiqc_bin_location, " ", directory, "/ -f -o ", directory, "/"))
  }

  fastqc_stats = fread(paste0(directory, "/multiqc_data/multiqc_fastqc.txt"), data.table = F)
  fastqc_stats$Sample = gsub("_R1_val_1", "_1", fastqc_stats$Sample) %>% gsub("_R2_val_2", "_2", .)
  fastqc_stats$Filename = gsub("_R1_val_1", "_1", fastqc_stats$Filename) %>% gsub("_R2_val_2", "_2", .)

  fastqc_stats = fastqc_stats %>% janitor::clean_names()

  fastqc_stats$filename = sub_seqids(fastqc_stats$filename, sample_info = sample_info)
  fastqc_stats$sample = sub_seqids(fastqc_stats$sample, sample_info = sample_info)

  # add in sample_ids (from sample_info table)
  # and get lane/split number
  fastqc_stats$sample_id = NA
  fastqc_stats$split_number = NA

  for(i in 1:nrow(fastqc_stats)){
    w = which(unlist(lapply(sample_info$sample_id, function(x) grepl(x, fastqc_stats$sample[i]))))
    fastqc_stats$sample_id[i] = sample_info$sample_id[w]
    split_and_read = gsub(sample_info$sample_id[w], "", fastqc_stats$sample[i])
    fastqc_stats$split_number[i] = gsub("_[1,2]$", "", split_and_read) %>% gsub("^_", "", .)
  }

  fastqc_stats$read_direction = str_sub(fastqc_stats$sample, -1,-1)

  # check if both read directions
  missing_samples = sample_info$sample_id[which(!sample_info$sample_id %in% fastqc_stats$sample_id)]
  if(length(missing_samples) >0) message(paste0("missing all fastqc data for: ", paste(missing_samples, collapse = ", ")))

  fwd_samples = fastqc_stats$sample[fastqc_stats$read_direction == 1] %>% gsub("[_][1,2]$", "", .) %>% sort()
  rev_samples = fastqc_stats$sample[fastqc_stats$read_direction == 2] %>% gsub("[_][1,2]$", "", .) %>% sort()
  if(!assertthat::are_equal(fwd_samples, rev_samples)){
    missing_rev = fwd_samples[which(!(fwd_samples %in% rev_samples))]
    missing_fwd = rev_samples[which(!(rev_samples %in% fwd_samples))]
    if(length(missing_fwd) >0) message(paste0("missing read_1 fastqc data for: ", paste(missing_fwd, collapse = ", ")))
    if(length(missing_rev) >0) message(paste0("missing read_2 fastqc data for: ", paste(missing_rev, collapse = ", ")))
  }

  ### ADD IN ANY EXTRA DATA YOU WANT HERE
  total_seqs = aggregate(total_sequences ~ sample_id+read_direction, fastqc_stats, sum) %>% arrange(sample_id, read_direction)

  return(total_seqs)
}


### fastq fine lines
# >> made with ???

read_in_fastq_lines = function(summary_file, sample_info, seqid2sampleid){

  total_reads = fread(summary_file, data.table = F)
  total_reads$sample = gsub("_R[1,2].fastq.gz$", "", total_reads$sample)
  total_reads$sample = sub_seqids(total_reads$sample, seqid2sampleid)

  total_reads$sample_id = NA
  total_reads$split_number = NA

  for(i in 1:nrow(total_reads)){
    w = which(unlist(lapply(sample_info$sample_id, function(x) grepl(x, total_reads$sample[i]))))
    if(length(w)>0){
      total_reads$sample_id[i] = sample_info$sample_id[w]
      split_and_read = gsub(sample_info$sample_id[w], "", total_reads$sample[i])
      total_reads$split_number[i] = gsub("^_", "", split_and_read)
    }
  }

  # check if both read directions
  missing_samples = sample_info$sample_id[which(!sample_info$sample_id %in% total_reads$sample_id)]
  if(length(missing_samples) >0) message(paste0("missing all fastqc data for: ", paste(missing_samples, collapse = ", ")))

  total_reads = total_reads[total_reads$sample_id %in% sample_info$sample_id,]
  total_reads$n_reads = total_reads$fastq_lines / 4

  total_reads_sum = aggregate(n_reads ~ sample_id, total_reads, sum)

  return(total_reads_sum)

}

library(rvest)
read_in_seqreports = function(report_files, sample_info, seqid2sampleid, add_date=TRUE){

  read_reports = NULL
  for(i in seq_along(report_files)){
    if(add_date){
      report_html_text = xml2::read_html(report_files[i]) %>%
        html_text() %>%
        str_trim() %>%
        unlist()
      date_locs = str_locate_all(report_html_text, "generated on")[[1]][1,]
      date = str_sub(report_html_text, date_locs[2]+2, date_locs[2]+20) %>% strv_split(",", 1)
    }

    read_reports = plyr::rbind.fill(read_reports, cbind(date, htmltab::htmltab(report_files[i], which = 1) %>% janitor::clean_names()))
  }
  read_reports$sequencing_id = strv_split(read_reports$sample_name, "_", 4)
  read_reports$lane = strv_split(read_reports$sample_name, "_", 2)
  read_reports$sample_id = seqid2sampleid$sample_id[match(read_reports$sequencing_id, seqid2sampleid$sequencing_id)]

  # check if both read directions
  missing_samples = sample_info$sample_id[which(!sample_info$sample_id %in% read_reports$sample_id)]
  if(length(missing_samples) >0) {
    message(paste0("missing sequencing report data for: ", paste(missing_samples, collapse = ", ")))
    message("please consider removing these samples from `sample_info` to prevent further warnings")
    message(paste0("`sample_info = sample_info[-which(sample_info$sample_id %in% c(", '"', paste(missing_samples, collapse = '", "'),'"',  ")),]`"))
  }

  read_reports = read_reports[read_reports$sample_id %in% sample_info$sample_id,]

  rm_columns = apply(read_reports, 2, function(x) length(which(is.na(x))) /length(x))
  read_reports = read_reports[,rm_columns < 0.25]

  if("percent_dups" %in% colnames(read_reports)){
    read_reports$percent_dups = gsub("%", "", read_reports$percent_dups) %>% as.numeric()
    read_reports$percent_dups = read_reports$percent_dups/100
  }
  if("percent_gc" %in% colnames(read_reports)){
    read_reports$percent_gc = gsub("%", "", read_reports$percent_gc) %>% as.numeric()
    read_reports$percent_gc = read_reports$percent_gc/100
  }
  if("length" %in% colnames(read_reports)){
    read_reports$length = gsub(" bp", "", read_reports$length) %>% as.numeric()
  }
  if("percent_failed" %in% colnames(read_reports)){
    read_reports$percent_failed = gsub("%", "", read_reports$percent_failed) %>% as.numeric()
    read_reports$percent_failed = read_reports$percent_failed/100
  }
  if("m_seqs" %in% colnames(read_reports)){
    read_reports$m_seqs = as.character(read_reports$m_seqs) %>% as.numeric()
    read_reports$n_seqs = read_reports$m_seqs * 1000000
  }
  return(read_reports)

}

check_subdirs_for_multiqc_data = function(dir = ".",
                                          dir_search_string,
                                          filetype_suffix,
                                          multiqc_filetype = "multiqc_samtools_flagstat.txt",
                                          multiqc_bin_location = "/Users/bsignal/opt/miniconda3/bin/multiqc",
                                          multiqc=TRUE){


  matched_dirs = list.dirs(path = dir, full.names = TRUE, recursive = TRUE) %>% grep(dir_search_string, ., value = T)
  #all_flags = NULL
  dirs_with_data = vector()
  for(i in seq_along(matched_dirs)){

    w = which(str_sub(list.files(matched_dirs[i]), -1*nchar(filetype_suffix), -1) == filetype_suffix)

    if(length(w) > 0){
      if(multiqc == TRUE){
        if(any(list.files(matched_dirs[i]) == "multiqc_data")){
          if(any(list.files(paste0(matched_dirs[i], "/multiqc_data/")) == multiqc_filetype)){
            dirs_with_data = c(dirs_with_data, paste0(matched_dirs[i], "/multiqc_data/"))
          }else{
            message("no ",  gsub(".txt", "", gsub("multiqc_" , "", multiqc_filetype)), " multiqc data for: ", matched_dirs[i])
          }
        }else{
          message("Trying to run multiqc in ", matched_dirs[i])
          system(paste0(multiqc_bin_location, " ", matched_dirs[i],  "/ -f -o ", matched_dirs[i], "/"))
          if(any(list.files(paste0(matched_dirs[i], "/multiqc_data/")) == multiqc_filetype)){
            dirs_with_data = c(dirs_with_data, paste0(matched_dirs[i], "/multiqc_data/"))
          }else{
            message("no ",  gsub(".txt", "", gsub("multiqc_" , "", multiqc_filetype)), " multiqc data for: ", matched_dirs[i])
          }
        }
      }else{
          dirs_with_data = c(dirs_with_data, list.files(matched_dirs[i], full.names = T)[w])
      }

    }


  }

  return(dirs_with_data)

}


read_in_flagstat = function(dir = ".",
                            multiqc_filetype = "multiqc_samtools_flagstat.txt",
                            multiqc_bin_location = "/Users/bsignal/opt/miniconda3/bin/multiqc"){
  flagstat_dirs = check_subdirs_for_multiqc_data(dir, "flagstat", ".flagstat", multiqc_filetype, multiqc_bin_location)
  all_flags = NULL
  for(i in seq_along(flagstat_dirs)){
    flagstat = fread(paste0(flagstat_dirs[i], multiqc_filetype)) %>% janitor::clean_names()
    flagstat$directory = flagstat_dirs[i]
    all_flags = plyr::rbind.fill(all_flags, flagstat)
  }
  return(all_flags)
}


read_in_peakmeth = function(dir = ".", sample_info){
  peakmeth_dirs = check_subdirs_for_multiqc_data(dir,
                                                 "peak",
                                                 ".peakmeth.txt", multiqc = F)
  all_peakmeth = NULL
  for(i in seq_along(peakmeth_dirs)){
    peakmeth = fread(peakmeth_dirs[i], data.table = F)
    peakmeth$directory = peakmeth_dirs[i]
    all_peakmeth = plyr::rbind.fill(all_peakmeth, peakmeth)
  }


  all_peakmeth$sample_id = NA
  for(i in seq_along(all_peakmeth$sample)){
    w = which(unlist(lapply(sample_info$sample_id, function(x) grepl(x, all_peakmeth$sample[i]))))
    if(length(w) > 0) all_peakmeth$sample_id[i] = sample_info$sample_id[w][1]
  }
   # check if both read directions
  missing_samples = sample_info$sample_id[which(!sample_info$sample_id %in% all_peakmeth$sample_id)]
  if(length(missing_samples) >0) {
    message(paste0("missing sequencing report data for: ", paste(missing_samples, collapse = ", ")))
  }


  return(all_peakmeth)
}


read_in_mosdepth = function(dir = ".", return_regions = FALSE, return_stats = c("length", "bases", "cov"), autosome_only = TRUE){

  mosdepth_files = check_subdirs_for_multiqc_data(dir,dir_search_string =  "mosdepth",
                                                  filetype_suffix = ".mosdepth.summary.txt",
                                                  multiqc=FALSE)

  for(i in seq_along(mosdepth_files)){
    mosdepth_summary = fread(mosdepth_files[i], data.table = F)
    if(autosome_only){

      mosdepth_summary = mosdepth_summary[mosdepth_summary$chrom %in% c(paste0("chr", c(1:99), "_region"), paste0("chr", c(1:99))),]
      mosdepth_summary = rbind(mosdepth_summary,
                                 data.frame(chrom = "total",
                                            length = sum(mosdepth_summary$length[!grepl("region", mosdepth_summary$chrom)]),
                                            bases = sum(mosdepth_summary$bases[!grepl("region", mosdepth_summary$chrom)]),
                                            mean = NA, min=NA, max=NA),
                                 data.frame(chrom = "total_region",
                                            length = sum(mosdepth_summary$length[grepl("region", mosdepth_summary$chrom)]),
                                            bases = sum(mosdepth_summary$bases[grepl("region", mosdepth_summary$chrom)]),
                                            mean = NA, min=NA, max=NA)
                                 )
    }


    if(return_regions){
      mosdepth_summary = mosdepth_summary[grepl("_region", mosdepth_summary$chrom),]
      mosdepth_summary$cov = as.numeric(mosdepth_summary$bases) / as.numeric(mosdepth_summary$length)
    }else{
      mosdepth_summary = mosdepth_summary[!grepl("_region", mosdepth_summary$chrom),]
    }

    if(i == 1){
      multi_mosdepth = mosdepth_summary[,c(1,2, which(colnames(mosdepth_summary) %in% c("bases", "cov")))]
    }else{
      if(return_regions){
        multi_mosdepth = left_join(multi_mosdepth, mosdepth_summary[,c(1,2, which(colnames(mosdepth_summary) %in% c("bases", "cov")))], by=c("chrom"))
      }else{
        multi_mosdepth = left_join(multi_mosdepth, mosdepth_summary[,c(1,2, which(colnames(mosdepth_summary) %in% c("bases", "cov")))], by=c("chrom", "length"))
      }
    }
  }
  if(return_regions){
    colnames(multi_mosdepth)[-1] = paste0(rep(gsub(paste0(".mosdepth.summary.txt", "$"), "", mosdepth_files), each = 3), c(".length", ".bases", ".cov"))
    multi_mosdepth = multi_mosdepth[,which(apply(multi_mosdepth, 2, function(x) length(which(is.na(x))) / length(x)) != 1)]
    multi_mosdepth = multi_mosdepth[,c(1, which(tools::file_ext(colnames(multi_mosdepth)) %in% return_stats))]
  }else{
    colnames(multi_mosdepth)[-c(1:2)] = gsub(paste0(".mosdepth.summary.txt", "$"), "", mosdepth_files)
  }

  for(j in 2:ncol(multi_mosdepth)){
    multi_mosdepth[,j] = as.numeric(multi_mosdepth[,j])
  }

  return(multi_mosdepth)
}

read_mosdepth_frip = function(dir = ".", autosome_only = TRUE){

  mosdepth_files = check_subdirs_for_multiqc_data(dir,dir_search_string =  "mosdepth",
                                                  filetype_suffix = ".mosdepth.summary.txt",
                                                  multiqc=FALSE)

  frip = data.frame(mosdepth_files, frip=NA)
  for(i in seq_along(mosdepth_files)){

    mosdepth_summary = fread(mosdepth_files[i], data.table = F)
    region_data = mosdepth_summary[grepl("_region", mosdepth_summary$chrom), ] %>% mutate(chrom = gsub("_region", "", chrom))
    genome_data = mosdepth_summary[!grepl("_region", mosdepth_summary$chrom), ]
    if(nrow(region_data) > 0){
      mosdepth_summary = full_join(genome_data, region_data, by="chrom", suffix = c(".genome", ".region"))

      if(autosome_only){
        mosdepth_summary = mosdepth_summary[mosdepth_summary$chrom %in% paste0("chr", c(1:99)),]
      }else if(!by_chrom){
        mosdepth_summary = mosdepth_summary[mosdepth_summary$chrom %in% "chr_total",]
      }
      frip$frip[i] = sum(mosdepth_summary$bases.region) / sum(mosdepth_summary$bases.genome)
    }
  }
  frip = frip[!is.na(frip$frip),]
  frip$sample_id = basename(frip$mosdepth_files) %>% gsub(".mosdepth.summary.txt", "", .)

  return(frip)
}


library(progress)

check_total_methylation_levels = function(bed_file, read_lines=1e6){

  counter = 1
  bed_ended = F
  read_header = T

  total_rows = fpeek::peek_count_lines(bed_file)

  reps_needed = ceiling(total_rows / read_lines)
  if(reps_needed > 3){
    pb = progress_bar$new(format = "Reading in bed file [:bar] :percent eta: :eta",
                          clear = FALSE, total = reps_needed, width = 60)
  }


  bed_head = fread(bed_file, data.table = F,
                  header = T, nrows = 1)

  while(!bed_ended){

    bed_read_start = ((counter-1)*read_lines)+1
    meth_percent = fread(bed_file, data.table = F,
                     header = F, skip = bed_read_start,
                     nrows = read_lines)
    colnames(meth_percent) = colnames(bed_head)

    samples = colnames(meth_percent) %>% grep(".C", .,value = T) %>% gsub(".C", "", .)

    if(counter == 1){
      C_totals = colSums(meth_percent[,paste0(samples, ".C")])
      cov_totals = colSums(meth_percent[,paste0(samples, ".cov")])
    }else{
      C_totals = C_totals + colSums(meth_percent[,paste0(samples, ".C")])
      cov_totals = cov_totals + colSums(meth_percent[,paste0(samples, ".cov")])
    }

    if(nrow(meth_percent) < read_lines) bed_ended = T
    counter = counter +1
    if(reps_needed > 3) pb$tick(1)
    #C_totals/cov_totals
  }

  c_percent = C_totals/cov_totals
  c_percent = stack(c_percent)[,c(2,1)]
  colnames(c_percent) = c("sample", "meth_percent")
  c_percent$sample = gsub(".C$", "", c_percent$sample)

  return(c_percent)
}

