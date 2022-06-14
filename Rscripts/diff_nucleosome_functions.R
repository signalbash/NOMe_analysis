weight_function <- function(u) (1-u^2)^3

.derivative_phi <- function(phi, local_c, local_t, mu, weight) {
  derivative = 0
  indicator_c = local_c > 0
  indicator_t = local_t > 0
  indicator_ct = local_c + local_t > 0
  
  if(nrow(local_c) == 1) {
    derivative =
      sum( indicator_c * ( mu * (digamma((mu * phi) + local_c + 1e-100) - digamma(mu * phi + 1e-100)) ) ) +
      sum( indicator_t * ((1 - mu) * (digamma( ((1 - mu) * phi) + local_t + 1e-100) - digamma( ((1-mu) * phi) + 1e-100))) ) -
      sum( indicator_ct * (digamma(phi + local_c + local_t + 1e-100) - digamma(phi)) )
  } else {
    for(g in seq(ncol(local_c))) {
      derivative = derivative +
        sum( indicator_c[,g] * (weight * mu[,g] * (digamma(mu[,g] * phi + local_c[,g] + 1e-100) - digamma(mu[,g] * phi + 1e-100))) ) +
        sum( indicator_t[,g] * (weight * (1 - mu[,g]) * (digamma((1 - mu[,g]) * phi + local_t[,g] + 1e-100) - digamma((1 - mu[,g]) * phi + 1e-100))) ) -
        sum( indicator_ct[,g] * (weight * (digamma(phi + local_c[,g] + local_t[,g] + 1e-100) - digamma(phi))) )
    }
  }
  
  derivative
}

.derivative_mu <- function(mu, local_c, local_t, phi, weight) {
  derivative = 0
  indicator_c = local_c > 0
  indicator_t = local_t > 0
  
  if(nrow(local_c) == 1) {
    derivative =
      sum( indicator_c * (digamma(mu * phi + local_c + 1e-100) - digamma(mu * phi + 1e-100)) ) -
      sum( indicator_t * (digamma((1 - mu) * phi + local_t + 1e-100) - digamma((1 - mu) * phi + 1e-100)) )
  } else {
    for(g in seq(ncol(local_c))) {
      derivative = derivative +
        sum( indicator_c[,g] * (weight * (digamma(mu * phi + local_c[,g]+ 1e-100) - digamma(mu * phi + 1e-100))) ) -
        sum( indicator_t[,g] * (weight * (digamma((1 - mu) * phi + local_t[,g] + 1e-100) - digamma((1 - mu) * phi + 1e-100))) )
    }
  }
  
  derivative
}

.log_likelihood  <- function(mu, phi, local_c, local_t, weight) {
  llik = 0
  indicator_c = local_c > 0
  indicator_t = local_t > 0
  
  if(nrow(local_c) == 1) {
    llik = llik +
      sum( indicator_c * (lgamma(mu * phi + local_c + 1e-100) - lgamma(mu * phi + 1e-100)) ) +
      sum( indicator_t * (lgamma((1 - mu) * phi + local_t + 1e-100) - lgamma((1 - mu) * phi + 1e-100)) )
  } else {
    for(g in seq(ncol(local_c))) {
      llik = llik +
        sum( indicator_c[,g] * (weight * (lgamma(mu * phi + local_c[,g] + 1e-100) - lgamma(mu * phi + 1e-100))) ) +
        sum( indicator_t[,g] * (weight * (lgamma((1 - mu) * phi + local_t[,g] + 1e-100) - lgamma((1 - mu) + 1e-100))) )
    }
  }
  
  2*llik
}


result_map = function(locus_idx, cov_mat, meth_mat){
  
  # Do not use local information when the local_window_size is 0
  local_loci_idx = locus_idx
  local_weights = 1
  
  # Collect Cov and M matrices for all the loci in the window
  # Rows are loci and columns are samples
  local_cov = matrix(cov_mat[local_loci_idx, ], nrow = 1)
  local_meth = matrix(meth_mat[local_loci_idx, ], nrow = 1)
  local_unmeth = local_cov - local_meth
  
  # Collect the correct rows of meth_est
  local_meth_est = matrix(meth_est[local_loci_idx, ], nrow = 1)
  
  #####################################
  
  ### Compute the degrees of freedom for the locus
  if(all(disp_groups)) {
    df_subtract = 2
  } else {
    df_subtract = 1
  }
  df = pmax(rowSums(local_cov[, disp_groups_idx, drop = FALSE] > 0) - df_subtract, 0)
  # Compute the degrees of freedom to be used in the test for differential methylation
  df = sum(df * local_weights)
  
  #####################################
  
  if(df > 1) {
    ### Common disp_groups calculation
    # This returns a singleton numeric
    if(derivative_phi(
      phi = max_inverse_disp,
      local_c = local_meth[, disp_groups_idx, drop = FALSE],
      local_t = local_unmeth[, disp_groups_idx, drop = FALSE],
      mu = local_meth_est[, disp_groups_idx, drop = FALSE],
      weight = local_weights) >= 0) {
      
      disp_est = max_inverse_disp
    } else if(derivative_phi(
      phi = min_inverse_disp,
      local_c = local_meth[, disp_groups_idx, drop = FALSE],
      local_t = local_unmeth[, disp_groups_idx, drop = FALSE],
      mu = local_meth_est[, disp_groups_idx, drop = FALSE],
      weight = local_weights) <= 0){
      
      disp_est = min_inverse_disp
    } else {
      disp_est = stats::uniroot(
        f = derivative_phi,
        interval = c(min_inverse_disp, max_inverse_disp),
        local_meth[, disp_groups_idx, drop = FALSE],
        local_unmeth[, disp_groups_idx, drop = FALSE],
        local_meth_est[, disp_groups_idx, drop = FALSE],
        local_weights)$root
    }
    
    #####################################
    
    ### Common group means calculation
    # This returns a numeric vector (control, case, control + case) with the mu est
    group_meth_est_list = list(control_idx, case_idx, c(control_idx, case_idx))
    group_meth_est = rep(0, length(group_meth_est_list))
    for(group_idx in seq_along(group_meth_est_list)) {
      if(sum(local_meth[, group_meth_est_list[[group_idx]], drop = FALSE]) == 0) {
        # If there are no local C reads, methylation is 0
        group_meth_est[group_idx] = 0
      } else if (sum(local_unmeth[, group_meth_est_list[[group_idx]], drop = FALSE]) == 0) {
        # If there are no local T reads, methylation is 1
        group_meth_est[group_idx] = 1
      } else {
        # Otherwise, do something fancier
        group_meth_est[group_idx] = stats::uniroot(
          f = derivative_mu,
          interval = c(min_meth, max_meth),
          local_meth[, group_meth_est_list[[group_idx]], drop = FALSE],
          local_unmeth[, group_meth_est_list[[group_idx]], drop = FALSE],
          disp_est,
          local_weights)$root
      }
    }
    
    #####################################
    
    ### log Likelihood ratio calculation
    log_lik_ratio =
      log_likelihood(
        mu = group_meth_est[1],
        phi = disp_est,
        local_c = local_meth[, control_idx, drop = FALSE],
        local_t = local_unmeth[, control_idx, drop = FALSE],
        weight = local_weights) +
      log_likelihood(
        mu = group_meth_est[2],
        phi = disp_est,
        local_c = local_meth[, case_idx, drop = FALSE],
        local_t = local_unmeth[, case_idx, drop = FALSE],
        weight = local_weights) -
      log_likelihood(
        mu = group_meth_est[3],
        phi = disp_est,
        local_c = local_meth[, c(control_idx, case_idx), drop = FALSE],
        local_t = local_unmeth[, c(control_idx, case_idx), drop = FALSE],
        weight = local_weights)
    
    #####################################
    
    locus_data = c(
      disp_est = disp_est,
      log_lik_ratio = log_lik_ratio,
      meth_control = group_meth_est[1]*100,
      meth_case = group_meth_est[2]*100,
      meth_all = group_meth_est[3]*100,
      df = df + 2)
  } else {
    # Not enough degrees of freedom, return NAs, these will be removed
    # with a message to the user with how many
    locus_data = c(
      disp_est = NA,
      log_lik_ratio = NA,
      meth_control = NA,
      meth_case = NA,
      meth_all = NA,
      df = df)
  }
  
  return(locus_data)
}



methylsig_dnd = function(
  windows,
  peak_span = 200,
  case_samples,
  control_samples,
  disp_groups = c('case' = TRUE, 'control' = TRUE),
  t_approx = TRUE,
  n_cores = 1,
  min_cov = 10,
  min_sample = 8){
  
  # Constants
  min_disp = 1e-6
  min_inverse_disp = 0.001
  max_inverse_disp = max(1/max(min_disp, 1e-6), min_inverse_disp)
  min_meth = 0
  max_meth = 1
  local_window_size=0

  
  local_weight_function = weight_function

  #####################################
  
  #####################################
  
  #case = grep("low", samples$Sample, value = T)
  #control = grep("high", samples$Sample, value = T)
  case = case_samples
  control = control_samples

  
  cov_mat = as.matrix(mcols(windows)[,grep(".cov$", colnames(mcols(windows)))])
  meth_mat = as.matrix(mcols(windows)[,grep(".C$", colnames(mcols(windows)))])
  
  control_idx = match(paste0(control, ".C"), colnames(meth_mat))
  case_idx = match(paste0(case, ".C"), colnames(meth_mat))
  
  filtered_rows = which(rowSums(cov_mat[,c(case_idx, control_idx)] >= min_cov) >= min_sample)
  windows = windows[filtered_rows]
  cov_mat = cov_mat[filtered_rows, ]
  meth_mat = meth_mat[filtered_rows, ]
  
  gr = granges(windows)
  num_loci = length(gr)
  
  
  if(all(disp_groups)) {
    disp_groups_idx = c(case_idx, control_idx)
  } else if (disp_groups['case'] & !disp_groups['control']) {
    disp_groups_idx = case_idx
  } else if (!disp_groups['case'] & disp_groups['control']) {
    disp_groups_idx = control_idx
  }
  
  #####################################
  
  # Estimate meth per locus within each group. The same value is used for all samples within the same group.
  # Note, the approach is to sum reads over all samples per group per locus
  meth_est = matrix(0, ncol = ncol(cov_mat), nrow = num_loci)
  meth_est[, case_idx] = base::rowSums(meth_mat[, case_idx, drop=F]) / (base::rowSums(cov_mat[, case_idx, drop=F]) + 1e-100)
  meth_est[, control_idx] = base::rowSums(meth_mat[, control_idx, drop=F]) / (base::rowSums(cov_mat[, control_idx, drop=F]) + 1e-100)

  #####################################
  
  result = do.call(rbind, parallel::mclapply(seq_along(gr), function(locus_idx){
    
    ### Deal with local information (or not)
    if(local_window_size != 0) {
      # NOTE: It is much faster to work with subsets of the result of start()
      # than it is to work with subsets of GRanges.
      
      # Get the indices which are within the local_window_size, but also limit to 5 CpGs on either side
      # NOTE, local information is only used with cytosine/CpG resolution data so start() is valid.
      # If regions were allowed, we would have to pay attention to which side we're on and use start()/end()
      local_loci_idx = intersect(
        which(abs(start(gr)[locus_idx] - start(gr)) < local_window_size),
        max(1, locus_idx - 5):min(num_loci, locus_idx + 5))
      
      if(length(local_loci_idx) == 1) {
        # Do not use local information when there is only one local locus
        local_loci_idx = locus_idx
        local_weights = 1
        
        # Collect Cov and M matrices for all the loci in the window
        # Rows are loci and columns are samples
        local_cov = matrix(cov_mat[local_loci_idx, ], nrow = 1)
        local_meth = matrix(meth_mat[local_loci_idx, ], nrow = 1)
        local_unmeth = local_cov - local_meth
        
        # Collect the correct rows of meth_est
        local_meth_est = matrix(meth_est[local_loci_idx, ], nrow = 1)
      } else {
        # We need to scale the loci in the window onto the interval [-1, 1] because
        # that is the domain of the local_weight_function.
        # This is a vector of the distances of the local loci to the loci of interest (domain)
        local_loci_norm = (start(gr)[local_loci_idx] - start(gr)[locus_idx]) / (local_window_size + 1)
        
        # Calculate the weights
        # Each is a vector of values of the weight function (range)
        local_weights = local_weight_function(local_loci_norm)
        
        # Collect Cov and M matrices for all the loci in the window
        # Rows are loci and columns are samples
        local_cov = cov_mat[local_loci_idx, ]
        local_meth = meth_mat[local_loci_idx, ]
        local_unmeth = local_cov - local_meth
        
        # Collect the correct rows of meth_est
        local_meth_est = meth_est[local_loci_idx, ]
      }
    } else {
      # Do not use local information when the local_window_size is 0
      local_loci_idx = locus_idx
      local_weights = 1
      
      # Collect Cov and M matrices for all the loci in the window
      # Rows are loci and columns are samples
      local_cov = matrix(cov_mat[local_loci_idx, ], nrow = 1)
      local_meth = matrix(meth_mat[local_loci_idx, ], nrow = 1)
      local_unmeth = local_cov - local_meth
      
      # Collect the correct rows of meth_est
      local_meth_est = matrix(meth_est[local_loci_idx, ], nrow = 1)
    }
    
    #####################################
    
    ### Compute the degrees of freedom for the locus
    if(all(disp_groups)) {
      df_subtract = 2
    } else {
      df_subtract = 1
    }
    df = pmax(rowSums(local_cov[, disp_groups_idx, drop = FALSE] > 0) - df_subtract, 0)
    # Compute the degrees of freedom to be used in the test for differential methylation
    df = sum(df * local_weights)
    
    #####################################
    
    if(df > 1) {
      ### Common disp_groups calculation
      # This returns a singleton numeric
      if(.derivative_phi(
        phi = max_inverse_disp,
        local_c = local_meth[, disp_groups_idx, drop = FALSE],
        local_t = local_unmeth[, disp_groups_idx, drop = FALSE],
        mu = local_meth_est[, disp_groups_idx, drop = FALSE],
        weight = local_weights) >= 0) {
        
        disp_est = max_inverse_disp
      } else if(.derivative_phi(
        phi = min_inverse_disp,
        local_c = local_meth[, disp_groups_idx, drop = FALSE],
        local_t = local_unmeth[, disp_groups_idx, drop = FALSE],
        mu = local_meth_est[, disp_groups_idx, drop = FALSE],
        weight = local_weights) <= 0){
        
        disp_est = min_inverse_disp
      } else {
        disp_est = stats::uniroot(
          f = .derivative_phi,
          interval = c(min_inverse_disp, max_inverse_disp),
          local_meth[, disp_groups_idx, drop = FALSE],
          local_unmeth[, disp_groups_idx, drop = FALSE],
          local_meth_est[, disp_groups_idx, drop = FALSE],
          local_weights)$root
      }
      
      #####################################
      
      ### Common group means calculation
      # This returns a numeric vector (control, case, control + case) with the mu est
      group_meth_est_list = list(control_idx, case_idx, c(control_idx, case_idx))
      group_meth_est = rep(0, length(group_meth_est_list))
      for(group_idx in seq_along(group_meth_est_list)) {
        if(sum(local_meth[, group_meth_est_list[[group_idx]], drop = FALSE]) == 0) {
          # If there are no local C reads, methylation is 0
          group_meth_est[group_idx] = 0
        } else if (sum(local_unmeth[, group_meth_est_list[[group_idx]], drop = FALSE]) == 0) {
          # If there are no local T reads, methylation is 1
          group_meth_est[group_idx] = 1
        } else {
          # Otherwise, do something fancier
          group_meth_est[group_idx] = stats::uniroot(
            f = .derivative_mu,
            interval = c(min_meth, max_meth),
            local_meth[, group_meth_est_list[[group_idx]], drop = FALSE],
            local_unmeth[, group_meth_est_list[[group_idx]], drop = FALSE],
            disp_est,
            local_weights)$root
        }
      }
      
      #####################################
      
      ### log Likelihood ratio calculation
      log_lik_ratio =
        .log_likelihood(
          mu = group_meth_est[1],
          phi = disp_est,
          local_c = local_meth[, control_idx, drop = FALSE],
          local_t = local_unmeth[, control_idx, drop = FALSE],
          weight = local_weights) +
        .log_likelihood(
          mu = group_meth_est[2],
          phi = disp_est,
          local_c = local_meth[, case_idx, drop = FALSE],
          local_t = local_unmeth[, case_idx, drop = FALSE],
          weight = local_weights) -
        .log_likelihood(
          mu = group_meth_est[3],
          phi = disp_est,
          local_c = local_meth[, c(control_idx, case_idx), drop = FALSE],
          local_t = local_unmeth[, c(control_idx, case_idx), drop = FALSE],
          weight = local_weights)
      
      #####################################
      
      locus_data = c(
        disp_est = disp_est,
        log_lik_ratio = log_lik_ratio,
        meth_control = group_meth_est[1]*100,
        meth_case = group_meth_est[2]*100,
        meth_all = group_meth_est[3]*100,
        df = df + 2)
    } else {
      # Not enough degrees of freedom, return NAs, these will be removed
      # with a message to the user with how many
      locus_data = c(
        disp_est = NA,
        log_lik_ratio = NA,
        meth_control = NA,
        meth_case = NA,
        meth_all = NA,
        df = df)
    }
    
    return(locus_data)
  }, mc.cores = n_cores))  
  #####################################
  
  # Build GRanges version of result
  result_gr = gr
  mcols(result_gr) = result
  
  ###TEMP FILTERING
  #result_gr = result_gr[start(result_gr) < 3051720]
  
  # Calculate pvalue
  if(t_approx) {
    result_gr$pvalue = stats::pt(-sqrt(pmax(result_gr$log_lik_ratio, 0)), result_gr$df) * 2
  } else {
    result_gr$pvalue = stats::pchisq(pmax(result_gr$log_lik_ratio, 0), 1, lower.tail = FALSE)
  }
  
  # Calculate meth_diff and set very small differences to 0
  result_gr$meth_diff = (result_gr$meth_case - result_gr$meth_control)
  result_gr$meth_diff[abs(result_gr$meth_diff) < 0.01] = 0
  
  dist_to_next = lead(start(result_gr)) - start(result_gr)
  
  #peak_span = 200
  min_step = min(dist_to_next, na.rm = T)
  peak_span_step = peak_span / min_step

  buffer_starts = seq(min(start(result_gr)), max(start(result_gr)), by=min_step)
  
  rm = which(buffer_starts %in% start(result_gr))
  buffer_starts = buffer_starts[-rm]
  buffer_ends = buffer_starts + (median(width(result_gr))-1)
  chrom = as.character(seqnames(result_gr)[1])
  
  buffer_ranges = makeGRangesFromDataFrame(data.frame(seqnames = chrom, start=buffer_starts, end=buffer_ends))
  
  ## add in buffer for first peakspan/2 elements? endbehaviour=1 should fix
  
  result_gr.buff = c(result_gr, buffer_ranges)
  result_gr.buff = sort(result_gr.buff)
  peak_sites = splus2R::peaks(-log(result_gr.buff$pvalue), span = peak_span_step+1, endbehavior = 1)
  
  result_gr.peak = result_gr.buff[peak_sites]
  
  result_gr.peak$fdr = stats::p.adjust(result_gr.peak$pvalue, method = 'BH')
  # Correct for multiple testing
  #result_gr$fdr = stats::p.adjust(result_gr$pvalue, method = 'BH')
  
  # Assign direction of hypermethylation (NOTE, this is not "significance")
  result_gr.peak$direction = ifelse(result_gr.peak$meth_diff >= 0, case[1], control[1])
  
  #####################################
  
  # Order output columns and attach to GRanges
  col_order = c(
    'meth_case',
    'meth_control',
    'meth_diff',
    'direction',
    'pvalue',
    'fdr',
    'disp_est',
    'log_lik_ratio',
    'df'
  )
  mcols(result_gr.peak) = mcols(result_gr.peak)[, col_order]
  
  #####################################
  
  # Check for NA results and indicate how many loci were dropped because of
  # a lack of available degrees of freedom
  insufficient_df = result_gr.peak$df == 1
  
  if(sum(insufficient_df) > 0) {
    result_gr.peak = result_gr.peak[!insufficient_df]
    message(sprintf('%s loci were dropped due to insufficient degrees of freedom (df = 1).', sum(insufficient_df)))
  }
  
  #result_gr = result_gr[which.max(result_gr$log_lik_ratio)]
  
  return(result_gr.peak)
}


filter_to_peak = function(result_range, peak_value, peak_span, peak_direction="high"){
  
  dist_to_next = lead(start(result_range)) - start(result_range)
  min_step = min(dist_to_next, na.rm = T)
  
  peak_span_step = peak_span / min_step
  
  buffer_starts = seq(min(start(result_range)), max(start(result_range)), by=min_step)
  
  rm = which(buffer_starts %in% start(result_range))
  buffer_starts = buffer_starts[-rm]
  buffer_ends = buffer_starts + (median(width(result_range))-1)
  
  ### FIX CHROMOSOMES
  chrom = as.character(seqnames(result_range)[1])
  
  buffer_ranges = makeGRangesFromDataFrame(data.frame(seqnames = chrom, start=buffer_starts, end=buffer_ends))
  
  result_range.buff = c(result_range, buffer_ranges)
  result_range.buff = sort(result_range.buff)
  vals_for_peak = result_range.buff[,peak_value]
  if(peak_direction %in% c("-", "low", "negative", "neg")){
    vals_for_peak = vals_for_peak*-1
  }
  peak_sites = splus2R::peaks(vals_for_peak, span = peak_span_step+1, endbehavior = 1)
  
  result_range.peak = result_range.buff[peak_sites]

}

range2window=function(nome.grange, 
                      signif.range = NULL,
                      window_width=100, window_spacing=20,
                      case=NULL, control=NULL,
                      min_cov=10, min_sample=8, extend=F){
  
  if(extend & !is.null(signif.range)){
    signif.range.ext = signif.range
    start(signif.range.ext) = start(signif.range.ext) -window_width
    end(signif.range.ext) = end(signif.range.ext) +window_width
    ol = findOverlaps(signif.range.ext, nome.grange)
    nome.grange.ol = nome.grange[ol@to]
    windows <- genomeBlocks(seqlengths(nome.grange.ol), width=window_width, spacing=window_spacing)
    #windows = windows[findOverlaps(nome.grange.ol, windows)@to]
    windows <-windows[windows %over% nome.grange.ol]
  }else{
    nome.grange.ol = nome.grange
    windows <- genomeBlocks(seqlengths(nome.grange), width=window_width, spacing=window_spacing)
    windows <-windows[windows %over% nome.grange]
  }
  
  
  for(i in 1:length(c(case, control))){
    #message(i)
    message(paste0("Processing ", c(case, control)[i]))
    # Counts of Cs and Ts in each window
    #message(" - Counting Cs")
    windows.counts <- cbind("C"=overlapSums(nome.grange.ol, windows, values(nome.grange.ol)[[paste0(c(case, control)[i], ".C")]]))
    #message(" - Counting Ts")
    windows.counts <- cbind(windows.counts, "T"=overlapSums(nome.grange.ol, windows, values(nome.grange.ol)[[paste0(c(case, control)[i], ".cov")]])-windows.counts[,"C"])
    windows.counts <- cbind(windows.counts,  "Cov"= (windows.counts[,1]+windows.counts[,2]))
    colnames(windows.counts) = paste0(c(case, control)[i], ".", c("C", "T", "cov"))
    
    if(i == 1){
      mcols(windows)= as.data.frame(windows.counts)
    }else{
      mcols(windows)= cbind(mcols(windows), as.data.frame(windows.counts))
    }
  }
  
  ### filter windows for high(ish) counts
  
  #case = grep("low", samples$Sample, value = T)
  #control = grep("high", samples$Sample, value = T)
  
  cov_mat = as.matrix(mcols(windows)[,grep(".cov", colnames(mcols(windows)))])
  
  if(is.null(case) & is.null(control)){
    control_idx = match(paste0(control, ".cov"), colnames(cov_mat))
    case_idx = match(paste0(case, ".cov"), colnames(cov_mat))
    sample_idx = c(case_idx, control_idx)
  }else{
    sample_idx = grep(".cov", colnames(cov_mat))
  }
  
  filtered_rows = which(rowSums(cov_mat[,sample_idx] >= min_cov) >= min_sample)
  
  windows = windows[filtered_rows]
  return(windows)
}


