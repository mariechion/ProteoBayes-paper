simu_data = function(
    nb_peptide = 5,
    nb_sample = 5,
    list_mean_diff = c(0, 1, 5, 10),
    list_var = c(1, 1, 1, 1),
    range_peptide = c(0, 50), 
    multivariate = FALSE){
  
  nb_group = length(list_mean_diff)
  
  
  db = tibble(
    'Peptide' = rep(paste0('Peptide_', 1:nb_peptide), each= nb_group*nb_sample),
    'Group' = rep(rep(1:nb_group, each = nb_sample),nb_peptide),
    'Sample' = rep(rep( 1:nb_sample, nb_group*nb_peptide))
    ) 
  if(multivariate){
    db = db %>%
      dplyr::group_by(Peptide) %>%
      dplyr::mutate(
        'Mean' = runif(1, range_peptide[1], range_peptide[2]) +
          list_mean_diff[Group]) %>% 
      arrange(Group) %>% 
      group_by(Group) %>% 
      dplyr::mutate(
        'Output' = Mean + mvtnorm::rmvnorm(
          nb_sample,
          rep(0, nb_peptide), 
          diag(list_var[Group], nb_peptide) +
          rWishart(1, nb_peptide,diag(nb_peptide))[,, 1]
          ) %>% as.vector()
      ) %>%
      ungroup() 
  } else {
    db = db %>%
      dplyr::group_by(Peptide) %>%
      dplyr::mutate(
        'Mean' = runif(1, range_peptide[1], range_peptide[2]) +
          list_mean_diff[Group],
        'Output' = Mean + rnorm(n(), 0, list_var[Group])
      ) %>%
      ungroup() 
  }
  
  return(db)
}


eval <- function(
    nb_peptide = 5,
    nb_sample = 5,
    list_mean_diff = c(0, 1, 5, 10),
    list_var = c(1, 1, 1, 1),
    multivariate = FALSE,
    t_test = FALSE,
    limma = FALSE,
    msqrob = FALSE,
    proDA = FALSE,
    missing_ratio = 0,
    imputation = FALSE,
    mu_0 = NULL,
    lambda_0,
    beta_0,
    alpha_0){
  
  db = simu_data(
    nb_peptide = nb_peptide,
    nb_sample = nb_sample,
    list_mean_diff = list_mean_diff,
    list_var = list_var,
    multivariate = multivariate) 

    if(imputation){
      db = db %>%
        group_by(Peptide, Group) %>%
        mutate(Average = mean(Output)) %>%
        mutate(Missing = rbinom(nb_sample, 1, missing_ratio)) %>%
        mutate(Output = if_else(Missing == 1, Average, Output)) %>%
        dplyr::select(- c(Missing, Average)) %>%
        ungroup()
      
      ## MNAR case
      # db = db %>%
      #   group_by(Peptide, Group) %>% 
      #   mutate(Missing = (Output < quantile(Output, missing_ratio))) %>% 
      #   mutate(Average = mean(ifelse(Missing, NA, Output), na.rm=T)) %>% 
      #   mutate(Output = if_else(Missing == 1, Average, Output)) %>% 
      #   dplyr::select(- c(Missing, Average)) %>% 
      #   ungroup()
    } else
    {
      db = db %>%
        group_by(Peptide, Sample) %>%
        mutate(Missing = rbinom(1, 1, missing_ratio)) %>%
        dplyr::filter(Missing == 0) %>%
        dplyr::select(- c(Missing)) %>%
        ungroup()
      
      ## MNAR case
      # db = db %>% 
      #   group_by(Peptide, Group) %>% 
      #   mutate(Missing = (Output < quantile(Output, missing_ratio))) %>% 
      #   dplyr::filter(Missing == 0) %>% 
      #   dplyr::select(- c(Missing)) %>% 
      #   ungroup()
    }

  res = db %>%
    posterior_mean(mu_0 = mu_0, lambda_0 = lambda_0,
                   alpha_0 = alpha_0, beta_0 = beta_0) %>%
    identify_diff() %>%
    dplyr::filter(Group == 1) %>%
    left_join(db %>% 
                dplyr::select(c(Peptide, Group, Mean)) %>%
                dplyr::distinct() %>%
                dplyr::rename(Group2 = Group),
              by = c('Peptide', 'Group2')) %>%
    mutate('MSE' = (Mean - mean2)^2,
           'CIC' = ((Mean > mean2 + qt(0.025, df = df2)*sqrt(var2)) &
                    (Mean < mean2 + qt(0.975, df = df2)*sqrt(var2))) * 100,
           'CI_width' = 2 * 1.96*sqrt(var2),
           'Diff_mean' = mean2 - mean) %>% 
    mutate(across(c(Group, Group2), .fns = as.character)) %>% 
    mutate('Multivariate' = FALSE)
  
  if(multivariate){
    res = res %>%
      bind_rows(
        db %>% multi_posterior_mean(mu_0 = mu_0, lambda_0 = lambda_0,
                                    nu_0 = alpha_0, Sigma_0 = beta_0) %>% 
          identify_diff() %>%
          dplyr::filter(Group == 1) %>%
          left_join(db %>% 
                      dplyr::select(c(Peptide, Group, Mean)) %>%
                      dplyr::distinct() %>%
                      dplyr::rename(Group2 = Group),
                    by = c('Peptide', 'Group2')) %>%
          mutate('MSE' = (Mean - mean2)^2,
                 'CIC' = ((Mean > mean2 + qt(0.025, df = df2)*sqrt(var2)) &
                          (Mean < mean2 + qt(0.975, df = df2)*sqrt(var2))) * 100,
                 'CI_width' = 2 * 1.96*sqrt(var2),
                 'Diff_mean' = mean2 - mean,
                 'Group' = as.character(Group),
                 'Group2' = as.character(Group2)) 
      )
  }
    
  if(t_test){
    res <- res %>% 
      left_join(multi_t_test(db),
                by = c('Peptide', 'Group', 'Group2'))
  }
  
  if(limma){
    res <- res %>% 
      left_join(multi_limma(db),
                by = c('Peptide', 'Group', 'Group2'))
  }
  
  if(msqrob){
    res <- res %>% 
      left_join(multi_msqrob(db),
                by = c('Peptide', 'Group', 'Group2'))
  }
  
  if(proDA){
    res <- res %>% 
      left_join(multi_proDA(db),
                by = c('Peptide'))
  }

  return(res)
}

eval_multi <- function(
    mean_ref, 
    cov_ref,
    mean_compare,
    cov_compare,
    lambda_0,
    nu_0,
    sigma_0 = NULL,
    nb_rep = 100,
    nb_sample = 5){
  
  res = tibble(mean_diff = c(), 
             MSE = c(), 
             CI_coverage = c()
            )
  for(i in 1:nb_rep)
  {
  ## Generate reference dataset
  db_ref = mvtnorm::rmvnorm(nb_sample, mean = mean_ref, sigma = cov_ref) %>% 
    as_tibble() %>%
    rename('Peptide_1' = V1,
           'Peptide_2' = V2,
           'Peptide_3' = V3) %>%
    mutate(Sample = row_number()) %>% 
    pivot_longer(-Sample, names_to = 'Peptide', values_to = 'Output') %>%
    mutate('Group' = 0)

  db_compare = mvtnorm::rmvnorm(
      nb_sample,
      mean = mean_compare,
      sigma = cov_compare) %>% 
    as_tibble() %>%
    rename('Peptide_1' = V1,
           'Peptide_2' = V2,
           'Peptide_3' = V3) %>%
    mutate(Sample = row_number()) %>% 
    pivot_longer(-Sample, names_to = 'Peptide', values_to = 'Output')  %>%
    mutate('Group' = 1) %>%
    bind_rows(db_ref)

  # post = multi_posterior_mean(db_compare, lambda_0 = lambda_0, nu_0 = nu_0, Sigma_0 = sigma_0) 
    post = posterior_mean(db_compare, lambda_0 = 1e-10, alpha_0 = 0.01, beta_0 = 0.3) 

  # mean_diff = post %>%
  #   multi_identify_diff(plot = FALSE) %>% 
  #   pluck('Diff_mean') %>%
  #   dplyr::filter(Group == 0, Group2 != 0) %>%
  #   pull(Diff_mean) %>%
  #   mean()

  mean_diff = post %>% 
    identify_diff() %>% 
    filter(Group == 0) %>%
    pull(Diff_mean) %>%
    mean()

  # mse = post %>%
  #   multi_identify_diff(plot = FALSE) %>% 
  #   pluck('Diff_mean') %>%
  #   dplyr::filter(Group == 0, Group2 != 0) %>%
  #   pull(Mean) %>%
  #   `^`(2) %>%
  #   mean() 

  mse = post %>% 
    identify_diff() %>% 
    dplyr::filter(Group == 0) %>% 
    mutate('Squared_error' = (mean(mean_compare) + Diff_mean)^2) %>%
    pull(Squared_error) %>%
    mean()

  ## Compute the posterior covariance matrix for the reference group
  dim = n_distinct(post$Peptide)   
    
  post_mean_ref = post %>%
    filter(Group == 0) %>% 
    pull(mu) %>%
    unique()
    
  is_in_CI = post %>%
      identify_diff() %>%
      filter(Group == 0) %>%
      mutate(
        'CI_inf' = mean + sqrt(var) * qt(0.025, df),
        'CI_sup' = mean + sqrt(var) * qt(0.975, df)
      ) %>% 
      mutate(
        'is_in_CI' = (mean(mean_ref) >= CI_inf) & (mean(mean_ref) <= CI_sup)
      ) %>%
    pull(is_in_CI) %>%
    min()

  # post_cov_ref = post %>%
  #   filter(Group == 0) %>% 
  #   mutate(Cov = Sigma / (lambda * (nu - dim + 1 ) ) ) %>%
  #   pull(Cov) %>%
  #   matrix(nrow = dim, ncol = dim)
   
  # ## Compute the 95% credible interval coverage for the reference group
  # is_in_CI = multi_CI(
  #   data = mean_ref,
  #   mean = post_mean_ref,
  #   cov = post_cov_ref,
  #   df = (unique(post$nu) - dim + 1) ) %>%  
  #   as.vector()

  res = res %>%
    bind_rows(
      tibble('mean_diff' = mean_diff,
             'MSE' = mse,
             'CI_coverage' = is_in_CI * 100)
    )
  }

  return(res)
} 

summarise_eval <- function(eval){
  eval %>%
    group_by(Group, Group2, Multivariate) %>%
    summarise(across(where(is.double),
                     .fns = list('Mean' = mean, 'Sd' = sd)),
                     .groups = 'drop') %>%
    mutate('MSE_Mean' = sqrt(MSE_Mean), 'MSE_Sd' = sqrt(MSE_Sd)) %>%
    mutate(across(where(is.double), ~ round(.x, 2))) %>%
    return()
}

multi_t_test <- function(data){
  
  data %>%
    tidyr::expand_grid('Group2' = unique(Group)) %>%
    dplyr::left_join(
      data %>%
        dplyr::select(- Mean) %>%
        dplyr::rename('Group2' = Group,
                      'Output2' = Output),
      by = c("Peptide", 'Group2', 'Sample')) %>%
    dplyr::filter(Group == 1, Group2 != 1) %>%
    dplyr::group_by(Peptide, Group, Group2) %>%
    dplyr::reframe('p_value' = tryCatch(t.test(Output,Output2)$p.value,
                                      error = function(e){return(NA)}),
                   'Diff_t_test' = tryCatch(t.test(Output,Output2)$estimate[2] -
                                     t.test(Output,Output2)$estimate[1],
                                     error = function(e){return(NA)}),
                   'Signif' = (p_value < 0.05))  %>% 
    mutate(across(c(Group, Group2), .fns = as.character)) %>% 
    return()
}

multi_proDA <- function(data){
  
  db = data %>% 
    dplyr::select(-Mean) %>% 
    tidyr::pivot_wider(names_from = c(Group, Sample),
                       values_from = Output,
                       names_prefix = 'Condition_') %>% 
    column_to_rownames('Peptide') %>% 
    as.matrix()
  
    design = str_sub(colnames(db), end = -3) %>% factor()
  
    fit = proDA(db, design = design, 
                max_iter = 100)
    
    test_diff(fit, Condition_1 - Condition_2) %>% 
      dplyr::mutate(diff = abs(diff)) %>% 
      dplyr::rename(Peptide = name, 
             Diff_proDA = diff, 
             p_value_proDA = pval) %>%
      dplyr::select(Peptide,Diff_proDA, p_value_proDA) %>%
      return()
}

multi_msqrob <- function(data){

  db_msqrob <- data %>% 
    dplyr::select(-Mean) %>% 
    # Merge Group and Sample columns to denote biological samples analysed
    unite(col = "Cond_Rep", c(Group, Sample), sep = "_") %>% 
    # Reshape data in wide format
    spread(key = "Cond_Rep", value = "Output")
  
  ## Create quantitative data matrix to match DAPAR requirements
  qdata <- db_msqrob %>%  
    column_to_rownames("Peptide") %>% 
    as.matrix()
  
  ## Create design dataframe to match DAPAR requirements
  metadata <- data.frame(Sample.name = colnames(qdata)) %>% 
    separate(Sample.name, into = c("Condition", "Rep"), sep = "_") %>% 
    mutate(Condition = factor(Condition))
  
  rownames(metadata) <- colnames(qdata)
  
  se <- SummarizedExperiment(
    assays = list(intensity = qdata),
    colData = metadata,
    rowData = data.frame(Peptide = rownames(qdata))
  )
  
  fit <- msqrob(se, ~ Condition, modelColumnName = "msqrob_protein")
  
  extract_msqrob <- function(fit) {
    
    models <- rowData(fit)$msqrob_protein
    
    res_list <- lapply(names(models), function(p) {
      
      m <- models[[p]]
      
      coef <- m@params$coefficients
      vcov <- m@params$vcovUnscaled * m@varPosterior
      df   <- m@dfPosterior
      
      # enlever intercept
      idx <- grep("^Condition", names(coef))
      
      beta <- coef[idx]
      se <- sqrt(diag(vcov)[idx])
      
      tstat <- beta / se
      pval <- 2 * pt(-abs(tstat), df = df)
      
      data.frame(
        Peptide = p,
        Condition = names(beta),
        logFC = as.numeric(beta),
        pvalue = as.numeric(pval)
      )
    })
    
    do.call(rbind, res_list)
  }
  
  extract_msqrob(fit) %>% 
    rename(Diff_msqrob = logFC, p_value_msqrob = pvalue, Group2 = Condition) %>%
    mutate(Group2 = str_sub(Group2, start = -1)) %>%
    mutate(Group = '1') %>% 
    return()
}


multi_limma <- function(data){
  # Note, limma and DAPAR use data in wide format.
  db_limma <- data %>% 
    dplyr::select(-Mean) %>% 
    # Merge Group and Sample columns to denote biological samples analysed
    unite(col = "Cond_Rep", c(Group,Sample), sep = "_") %>% 
    # Reshape data in wide format
    spread(key = "Cond_Rep", value = "Output")
  
  ## Create quantitative data matrix to match DAPAR requirements
  qdata <- db_limma %>%  
    column_to_rownames(var = "Peptide") 
  
  ## Create design dataframe to match DAPAR requirements
  metadata <- data.frame(Sample.name = colnames(qdata)) %>% 
    separate(Sample.name, sep = "_", into = c("Condition",NA), remove = F) %>% 
    mutate(Bio.Rep = 1:ncol(qdata))
  
  qdata <- qdata %>% 
    select(metadata$Sample.name)
  
  ## Moderated t-test - OnevsOne setting
  res_limma <- DAPAR::limmaCompleteTest(qData = as.matrix(qdata),
                                        sTab = metadata,
                                        comp.type = "OnevsOne")
  
  P_Value <- res_limma$P_Value %>%
    rownames_to_column(var = "Peptide") %>%
    pivot_longer(-Peptide,
                 names_to = "Comparison",
                 values_to = "LM_p_value") %>%
    separate(col = Comparison,
             into = c("Group", NA, "Group2", NA),
             sep = "_")
  
  return(P_Value)
}


ci_coverage <- function(
    posterior,
    CI_level_seq = seq(0.05, 1, 0.05)
    ){
  
    posterior %>%
      dplyr::mutate('sigma' = sqrt(.data$beta / (.data$lambda * .data$alpha)),
                    'df' =  2 * .data$alpha) %>%
      uncount(length(CI_level_seq)) %>% 
      mutate(CI_level = rep(CI_level_seq, nrow(posterior))) %>% 
      dplyr::reframe(
        .data$Peptide, .data$Group, .data$mu, .data$CI_level,
        'CI_inf' = extraDistr::qlst(
          p = CI_level/2,
          df =.data$df,
          mu = .data$mu,
          sigma = .data$sigma,
          lower.tail = T),
        'CI_sup' =  extraDistr::qlst(
          p = CI_level/2,
          df =.data$df,
          mu = .data$mu,
          sigma = .data$sigma,
          lower.tail = F)
      ) %>% 
    return()
}

data_preprocessing <- function(data, output_str_id, prop_NA, nb_group, nb_rep,
                               maxquant, normalize){
  if(maxquant){
    db_temp <- data %>% 
      ## Select columns of interest
      dplyr::select(Sequence, Leading.razor.protein, starts_with(output_str_id)) %>% 
      ## Rename columns to match ProteoBayes requirements
      dplyr::rename(Peptide = Sequence, 
             Protein = Leading.razor.protein) %>%
      ## Remove reverse and contaminant proteins, remove iRT. 
      dplyr::filter(!str_detect(Protein, "CON") & !str_detect(Protein, "REV") &
               !str_detect(Protein, "iRT")) %>% 
      ## Replace 0 intensity values by NA
      dplyr::mutate(across(starts_with(output_str_id), ~ if_else(.x == 0, 
                                                         true = NA, 
                                                         false = log2(.x))))
    if(normalize){
      db <- db_temp %>%
        dplyr::select(-Protein) %>% 
        column_to_rownames(var = "Peptide") %>%
        as.matrix() %>% 
        preprocessCore::normalize.quantiles(copy = F) %>% 
        as_tibble(rownames = "Peptide") %>% 
        left_join(x = db_temp %>% select("Peptide", "Protein"),
                  by = "Peptide") %>% 
        pivot_longer(-c("Peptide","Protein"), 
                     names_to = c("Group", "Sample"), names_sep = "_",
                     values_to = "Output") %>% 
        mutate(Group = str_replace(Group, output_str_id, ""))
    }
    else{
      db <- db_temp %>%
        pivot_longer(-c("Peptide","Protein"), 
                     names_to = c("Group", "Sample"), names_sep = "_",
                     values_to = "Output") %>% 
        mutate(Group = str_replace(Group, output_str_id, ""))
    }
      
  }
  else{db <- data}
  
  # Get quantified peptides: 
  ## *max_NA* missing values in each of the *nb_group* groups.
  max_NA = round(nb_rep * prop_NA)
  qtfd_pept <- db %>% 
    group_by(Peptide, Group) %>% 
    summarise(Count_NA = sum(is.na(Output))) %>% 
    mutate(Is_Qtfd = Count_NA <= max_NA) %>%  
    group_by(Peptide) %>% 
    summarise(Qtfd = sum(Is_Qtfd)) %>% 
    filter(Qtfd == nb_group) %>% 
    pull(Peptide)
  
  db <- db %>% 
    ## Keep only quantified peptides
    ## and remove missing values
    filter(Peptide %in% qtfd_pept &
             !is.na(Output))
  
  return(db)
}

PB_DiffAna <- function(data, multi, mu_0, lambda_0, beta_0, alpha_0){
  if (!multi) {
    out <- ProteoBayes::posterior_mean(data = data,
                          mu_0 = mu_0,
                          lambda_0 = lambda_0, 
                          beta_0 = beta_0, 
                          alpha_0 = alpha_0)
    
    diff_out <- ProteoBayes::identify_diff(out)
  }
  
  return(diff_out)
}

LM_DiffAna <- function(data, group_labels, nb_rep, alpha, FDR){
  # Note, limma and DAPAR use data in wide format.
  db_limma <- data %>% 
    # Merge Group and Sample columns to denote biological samples analysed
    unite(col = "Cond_Rep", c(Group,Sample), sep = "_") %>% 
    # Reshape data in wide format
    spread(key = "Cond_Rep", value = "Output")
  
  ## Create quantitative data matrix to match DAPAR requirements
  qdata <- db_limma %>%  
    column_to_rownames(var = "Peptide") %>% 
    select(-Protein)
  
  ## Create design dataframe to match DAPAR requirements
  metadata <- data.frame(Sample.name = colnames(qdata)) %>% 
    separate(Sample.name, sep = "_", into = c("Condition",NA), remove = F) %>% 
    mutate(Bio.Rep = 1:ncol(qdata)) %>% 
    arrange(factor(Condition, levels = group_labels))
  
  qdata <- qdata %>% 
    select(metadata$Sample.name)
  
  ## Moderated t-test - OnevsOne setting
  res_limma <- DAPAR::limmaCompleteTest(qData = as.matrix(qdata),
                                        sTab = metadata,
                                        comp.type = "OnevsOne")
  
  P_Value <- res_limma$P_Value %>%
    rownames_to_column(var = "Peptide")
  
  if(!is.null(FDR)){
    P_Value <- P_Value %>%
      mutate(across(-Peptide,
                    ~ cp4p::adjust.p(p = ., alpha = FDR)$adjp %>%
                      select(adjusted.p) %>% pull))
  }
  
  P_Value <- P_Value %>%
    pivot_longer(-Peptide,
                 names_to = "Comparison",
                 values_to = "pval") %>%
    separate(col = Comparison,
             into = c("Group", NA, "Group2", NA),
             sep = "_") %>%
    mutate(Signif = pval < alpha)
  
  FC <- res_limma$logFC %>%
    rownames_to_column(var = "Peptide") %>%
    pivot_longer(-Peptide,
                 names_to = "Comparison",
                 values_to = "log2FC") %>%
    separate(col = Comparison,
             into = c("Group", NA, "Group2", NA),
             sep = "_")
  
  P_Value %>%
    left_join(y = FC, by = c("Peptide", "Group", "Group2")) %>%
    return()
}

CombineDA <- function(data, group_labels, fmol_labels, diff_str_id,
                      PB_res, LM_res, summary){
  ## Experimental condition/design to calculate true FC
  exp_cond <- tibble(Group = group_labels,
                     fmol = fmol_labels) %>%
    expand_grid(Group2 = unique(Group)) %>%
    filter(Group != Group2) %>%
    left_join(tibble(Group2 = group_labels,
                     fmol2 = fmol_labels),
              by = "Group2") %>%
    mutate(True_log2FC = log2(fmol/fmol2), .keep = "unused")
  
  # Create a "reference mean"
  ref_pept <- data %>%
    arrange(Peptide) %>%
    left_join(tibble(Group = group_labels,
                     fmol = fmol_labels),
              by = "Group") %>%
    mutate(log2FC = log2(max(fmol)/fmol),
           Mean = if_else(str_detect(Protein, diff_str_id),
                          true = Output + log2FC,
                          false = Output)) %>%
    group_by(Peptide, Protein) %>%
    mutate(Mean = if_else(str_detect(Protein, diff_str_id),
                          true = mean(Mean, na.rm = T) - log2FC,
                          false = mean(Mean, na.rm = T))) %>%
    group_by(Peptide, Protein, Group, Mean) %>%
    summarise(Avg = mean(Output)) %>%
    ungroup %>%
    mutate(Diff = Mean - Avg, Case = str_detect(Protein, diff_str_id))
  
  # Merge result tables
  db_results <- PB_res %>%
    filter(Group2 == tail(group_labels, n=1)) %>%
    # Add limma results
    left_join(y = LM_res %>%
                filter(Group2 == tail(group_labels, n=1)) %>%
                rename(LM_log2FC = log2FC),
              by = join_by(Peptide, Group, Group2)) %>% 
    # Add Protein Column
    left_join(y = data %>%
                select(Peptide,Protein) %>%
                distinct,
               by = "Peptide") %>%
    relocate(Protein, .after = "Peptide") %>%   
    # Add ground truth column
    mutate(Truth = if_else(str_detect(Protein, diff_str_id),
                           true = T, false = F)) %>%
    # Add fmol equivalent to groups
    left_join(y = exp_cond, by = join_by(Group, Group2)) %>%
    # Add FC from ProteoBayes and true FC
    mutate(PB_log2FC = mu - mean2,
           True_log2FC = if_else(Truth, true = True_log2FC, false = 0)) %>%
    # Add ref mean
    left_join(ref_pept %>% select(Peptide, Protein, Group, Mean),
              by = join_by(Peptide, Protein, Group == Group))
  
  db_eval <- db_results %>%
    mutate('MSE' = (Mean - mu)^2,
           'CIC' = ((Mean > CI_inf) & (Mean < CI_sup)) * 100,
           'Diff_mean' = PB_log2FC,
           'Diff_LM' = LM_log2FC,
           'CI_width' = CI_sup - CI_inf,
           'Distinct' = (Distinct == Truth)*100,
           'Signif' = (Signif == Truth)*100,
           "MSE_PBT" = (PB_log2FC - True_log2FC)^2,
           "MSE_LMT" = (LM_log2FC - True_log2FC)^2) %>%
    mutate(True_diff_mean = round(True_log2FC, digits = 2),
           .after = 'Truth', .keep = "unused")

  if (summary){
    db_eval <- db_eval %>%
      group_by(Group, Group2, Truth, True_diff_mean) %>%
      summarise(across(c(MSE, CIC, Diff_mean, Diff_LM,
                         pval, CI_width, Distinct, Signif,
                         MSE_PBT, MSE_LMT),
                       .fns = list('Mean' = ~mean(.x,na.rm = T),
                                   'Sd' = ~sd(.x,na.rm=T))),
                PB_na = sum(is.na(PB_log2FC)),
                LM_na = sum(is.na(pval)),
                nb_pept = n(),
                .groups = 'drop') %>%
      mutate('MSE_Mean' = sqrt(MSE_Mean), 'MSE_Sd' = sqrt(MSE_Sd),
             'MSE_PBT_Mean' = sqrt(MSE_PBT_Mean), 'MSE_PBT_Sd' = sqrt(MSE_PBT_Sd),
             'MSE_LMT_Mean' = sqrt(MSE_LMT_Mean), 'MSE_LMT_Sd' = sqrt(MSE_LMT_Sd)) %>%
      mutate(across(MSE_Mean:MSE_LMT_Sd, ~ round(.x, 2))) %>%
      reframe(Group, Group2, Truth, True_diff_mean, nb_pept, PB_na, LM_na,
              'PB_diff_mean' =  paste0(Diff_mean_Mean, ' (', Diff_mean_Sd, ')'),
              'CI_width' =  paste0(CI_width_Mean, ' (', CI_width_Sd, ')'),
              'Distinct' = paste0(Distinct_Mean, ' (', Distinct_Sd, ')'),
              'LM_diff_mean' =  paste0(Diff_LM_Mean, ' (', Diff_LM_Sd, ')'),
              'p_value' =  paste0(pval_Mean, ' (', pval_Sd, ')'),
              'Signif' = paste0(Signif_Mean, ' (', Signif_Sd, ')'),
              'RMSE' = paste0(MSE_Mean, ' (', MSE_Sd, ')'),
              'CIC' =  paste0(CIC_Mean, ' (', CIC_Sd, ')'),
              'RMSE_PBtoT' = paste0(MSE_PBT_Mean, ' (', MSE_PBT_Sd, ')'),
              'RMSE_LMtoT' = paste0(MSE_LMT_Mean, ' (', MSE_LMT_Sd, ')'))
  }
  
  return(db_eval)  
}

real_data_eval <- function(data, type, maxquant = T, normalize,
                           prop_NA = 0.2,
                           multi = F,
                           mu_0 = NULL, 
                           lambda_0 = 1e-10, beta_0 = 1, alpha_0 = 1,
                           alpha = 0.05, FDR = NULL,
                           summary = T){
  # Set parameters depending on the data used
  if (type == "ARATH"){
    nb_rep = 3
    output_str_id = "LFQ.intensity."
    group_labels = c("Point1", "Point2", "Point3",
                     "Point4", "Point5", "Point6",
                     "Point7")
    nb_group = length(group_labels)
    fmol_labels = c(0.05, 0.25, 0.5, 1.25, 2.5, 5, 10)
    diff_str_id = "UPS"
  }
  if (type == "YST"){
    nb_rep = 3
    output_str_id = "Intensity.TG."
    group_labels = c("0.5fmol", "1fmol", "2.5fmol",
                     "5fmol", "10fmol", "25fmol")
    nb_group = length(group_labels)
    fmol_labels = c(0.5, 1, 2.5, 5, 10, 25)
    diff_str_id = "ups"
  }
  if (type == "YST_B"){
    nb_rep = 4
    output_str_id = "Intensity."
    group_labels = c("10amol", "50amol", "100amol", "250amol", "500amol",
                      "1fmol", "5fmol", "10fmol", "25fmol","50fmol")
    nb_group = length(group_labels)
    fmol_labels = c(0.01,0.05,0.1,0.25,0.5, 1, 5, 10, 25, 50)
    diff_str_id = "ups"
  }
  if (type == "MOUSE"){
    nb_rep = 5
    output_str_id = NULL
    group_labels = c("S1", "S2", "S3", "S4", "S5")
    nb_group = length(group_labels)
    fmol_labels = c(0.75, 0.83, 1.07, 2.04, 7.54)
    diff_str_id = "_UPS"
    
    data_temp <- data %>% 
      select(R.Condition, R.Replicate, PG.ProteinAccessions, EG.StrippedSequence,
             FG.MS1PeakArea) %>% 
      rename(Peptide = EG.StrippedSequence,
             Protein = PG.ProteinAccessions,
             Group = R.Condition,
             Sample = R.Replicate,
             Intensity = FG.MS1PeakArea) %>% 
      group_by(Group,Sample, Protein, Peptide) %>% 
      summarise(Output = mean(Intensity)) %>% 
      mutate(Output = if_else(Output <= 1, NA, log2(Output))) %>% 
      ungroup()
    
    if(normalize){
      data <- data_temp %>% 
        pivot_wider(names_from = c("Group","Sample"), values_from = "Output") %>% 
        select(-Protein) %>% 
        column_to_rownames(var = "Peptide") %>%
        as.matrix() %>% 
        preprocessCore::normalize.quantiles(copy = F) %>% 
        as_tibble(rownames = "Peptide") %>% 
        left_join(x = data_temp %>% select("Peptide", "Protein") %>% distinct,
                  by = "Peptide") %>% 
        pivot_longer(-c("Peptide","Protein"), 
                     names_to = c("Group", "Sample"), names_sep = "_",
                     values_to = "Output")
    }
    else{
      data <- data_temp
    }
  }
  
  # Data preprocessing to ProteoBayes format
  db <- data_preprocessing(data = data, 
                           output_str_id = output_str_id,
                           prop_NA = prop_NA, 
                           nb_group = nb_group,
                           nb_rep = nb_rep,
                           maxquant = maxquant,
                           normalize = normalize)
  # ProteoBayes Differential Analysis
  PB_res <- PB_DiffAna(data = db, multi = multi,
                       mu_0 = mu_0, lambda_0 = lambda_0,
                       alpha_0 = alpha_0, beta_0 = beta_0)

  # Limma Differential Analysis
  LM_res <- LM_DiffAna(data = db, group_labels = group_labels, nb_rep = nb_rep,
                       alpha = alpha, FDR = FDR)

  # Combine results
  results <- CombineDA(data = db, group_labels = group_labels, 
                       fmol_labels = fmol_labels, 
                       diff_str_id = diff_str_id,
                       PB_res = PB_res, LM_res = LM_res,
                       summary = summary)
  
  
  
  return(list(
    data = db,
    ProteoBayes = PB_res,
    DAPAR = LM_res,
    results = results
  ))
  
}

sensitivity_uni <- function(
    nb_peptide = 1000, 
    nb_sample = 5, 
    param, 
    grid_param,
    alpha_0 = 0.01 , 
    beta_0 = 0.3,
    lambda_0 = 1e-10,
    list_mean_diff = c(0, 1, 10, 10),
    list_var = c(1, 10, 1, 100),
    CI_level = 0.05){

  nb_sample = 5
  db = simu_data(nb_peptide = nb_peptide,
                 nb_sample = nb_sample,
                 list_mean_diff = list_mean_diff,
                 list_var = list_var)
  
  
  list_post = tibble::tibble()
  
  if(param == 'lambda'){
    for(i in grid_param){
      
      list_post = list_post %>% 
        bind_rows(posterior_mean(db, 
                                 alpha_0 = alpha_0,
                                 beta_0 = beta_0,
                                 lambda_0 = i) %>% 
                    mutate(Param = i)
                  ) 
    }
  } else if(param == 'alpha'){
    for(i in grid_param){
      
      list_post = list_post %>% 
        bind_rows(posterior_mean(db, 
                                 alpha_0 = i,
                                 beta_0 = beta_0,
                                 lambda_0 = lambda_0) %>% 
                    mutate(Param = i)) 
    }
  } else if(param == 'beta'){
    for(i in grid_param){
      
      list_post = list_post %>% 
        bind_rows(posterior_mean(db, 
                                 alpha_0 = alpha_0,
                                 beta_0 = i,
                                 lambda_0 = lambda_0) %>% 
                    mutate(Param = i)) 
    }
  }
    
  list_post %>% 
    left_join(db %>% select(c(Peptide, Group, Mean)) %>% distinct(),
              by = c('Peptide', 'Group')) %>% 
      mutate(
        'NLL' = - extraDistr::dlst(
          x = Mean, 
          df = 2 * alpha,
          mu = mu,
          sigma = sqrt(beta / (lambda * alpha)), 
          log = TRUE),
        'CI_inf' = extraDistr::qlst(
          p = CI_level/2, 
          df = 2 * alpha,
          mu = mu,
          sigma = sqrt(beta / (lambda * alpha)),
          lower.tail = T),
        'CI_sup' = extraDistr::qlst(
          p = CI_level/2, 
          df = 2 * alpha,
          mu = mu,
          sigma = sqrt(beta / (lambda * alpha)),
          lower.tail = F),
        'CIC' = ((Mean > CI_inf) & (Mean <  CI_sup)) * 100) %>% 
      return()
}

sensitivity_multi <- function(
    nb_rep = 100, 
    nb_sample = 5, 
    param, 
    grid_param,
    lambda_0 = 1e-10, 
    cov_diag = FALSE){

  dim = 3
  true_mean = rep(0, dim)
  nu_0 = 2
  cov = matrix(c(1, 0.7, 0.2, 0.7, 1, 0.5, 0.2, 0.5, 1), nrow = dim, ncol = dim)
  
  list_res = tibble::tibble()
  
  for(i in 1:nb_rep)
  {
    db = mvtnorm::rmvnorm(nb_sample,
                          mean = true_mean,
                          sigma = cov) %>% 
      as_tibble() %>%
      rename('Peptide_1' = V1,
             'Peptide_2' = V2,
             'Peptide_3' = V3) %>%
      mutate(Sample = row_number()) %>% 
      pivot_longer(-Sample, names_to = 'Peptide', values_to = 'Output') %>%
      mutate('Group' = 1)
  
  if(param == 'Sigma'){
    
    list_post = tibble::tibble()
    
    for(j in grid_param){
      post = multi_posterior_mean(
        db, 
        lambda_0 = 1e-10,
        nu_0 = nu_0, 
        Sigma_0 = diag(cov) + j * (cov- diag(cov)) ) 
      
      post_cov = post %>%
        mutate(Cov = Sigma / (lambda * (nu - dim + 1) ) ) %>%
        pull(Cov) %>%
        matrix(nrow = dim, ncol = dim)
      
      post_mu = post %>% 
        dplyr::select(Peptide, mu) %>%
        dplyr::distinct() %>% 
        pull(mu)
      
      post_nu = unique(post$nu) - dim + 1
      post_lambda = post$lambda %>% unique()
      
      # ## Compute the 95% credible interval coverage for the reference group
      is_in_CI = multi_CI(
        data = true_mean,
        mean = post_mu,
        cov = post_cov,
        df = post_nu) %>%
        as.vector()
      
      NLL = - LaplacesDemon::dmvt(
        x = true_mean, 
        df = post_nu,
        mu = post_mu,
        S = post_cov / (post_lambda * post_nu),
        log = TRUE)
      
      list_post = list_post %>%
        bind_rows(tibble(
          'Param' = j,
          'NLL' = NLL,
          'CIC' = is_in_CI * 100))
      
    }
  } else if(param == 'nu'){
    
    list_post = tibble::tibble()
    for(j in grid_param){

      if(cov_diag){
        cov_0 = diag(dim)
      } else{
        cov_0 = cov
      }
      
      post = multi_posterior_mean(
        db, 
        lambda_0 = 1e-10,
        nu_0 = j, 
        Sigma_0 = (j - dim - 1) * cov_0)
      
      post_cov = post %>%
        mutate(Cov = Sigma / (lambda * (nu - dim + 1) ) ) %>%
        pull(Cov) %>%
        matrix(nrow = dim, ncol = dim)
      
      post_mu = post %>% 
        dplyr::select(Peptide, mu) %>%
        dplyr::distinct() %>% 
        pull(mu)
      
      post_nu = unique(post$nu) - dim + 1
      post_lambda = post$lambda %>% unique()
      
      # ## Compute the 95% credible interval coverage for the reference group
      is_in_CI = multi_CI(
        data = true_mean,
        mean = post_mu,
        cov = post_cov,
        df = post_nu) %>%
        as.vector()
      
      NLL = - LaplacesDemon::dmvt(
        x = true_mean, 
        df = post_nu,
        mu = post_mu,
        S = post_cov / (post_lambda * post_nu),
        log = TRUE)
      
      list_post = list_post %>%
        bind_rows(tibble(
          'Rep' = i,
          'Param' = j,
          'NLL' = NLL,
          'CIC' = is_in_CI * 100))
      }
    }
    list_res = list_res %>%
    bind_rows(list_post %>%
                mutate('Rep' = i))
  }
  
    return(list_res)
}


