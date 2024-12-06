
data_preprocessing <- function(data, output_str_id, prop_NA, nb_group, nb_rep,
                               maxquant, normalize){
  if(maxquant){
    db_temp <- data %>% 
      ## Select columns of interest
      select(Sequence, Leading.razor.protein, starts_with(output_str_id)) %>% 
      ## Rename columns to match ProteoBayes requirements
      rename(Peptide = Sequence, 
             Protein = Leading.razor.protein) %>%
      ## Remove reverse and contaminant proteins, remove iRT. 
      filter(!str_detect(Protein, "CON") & !str_detect(Protein, "REV") &
               !str_detect(Protein, "iRT")) %>% 
      ## Replace 0 intensity values by NA
      mutate(across(starts_with(output_str_id), ~ if_else(.x == 0, 
                                                         true = NA, 
                                                         false = log2(.x))))
    if(normalize){
      db <- db_temp %>%
        select(-Protein) %>% 
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
                      PB_res, LM_res){
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
    mutate(Mean = if_else(str_detect(Protein, diff_str_id),
                          true = Output + log2(max(fmol)/fmol),
                          false = Output)) %>%
    group_by(Peptide, Protein) %>%
    mutate(Mean = if_else(str_detect(Protein, diff_str_id),
                          true = mean(Mean, na.rm = T) - log2(max(fmol)/fmol),
                          false = mean(Mean, na.rm = T))) %>%
    group_by(Peptide, Protein, Group, Mean) %>%
    summarise(Avg = mean(Output)) %>%
    ungroup %>%
    mutate(Diff = Mean - Avg, Case = str_detect(Protein, diff_str_id))
  
  # Merge result tables
  db_results <- PB_res %>%
    filter(Group2 == tail(group_labels, n=1)) %>%
    # Add Protein Column
    left_join(y = .,
              x = data %>%
                select(Peptide,Protein) %>%
                distinct,
              by = "Peptide") %>%
    # Add limma results
    left_join(y = LM_res %>%
                filter(Group2 == tail(group_labels, n=1)) %>%
                rename(LM_log2FC = log2FC),
              by = join_by(Peptide, Group, Group2)) %>%
    # Add ground truth column
    mutate(Truth = if_else(str_detect(Protein, diff_str_id),
                           true = T, false = F)) %>%
    # Add fmol equivalent to groups
    left_join(y = exp_cond, by = join_by(Group, Group2)) %>%
    # Add FC from ProteoBayes and true FC
    mutate(PB_log2FC = mu - mu2,
           True_log2FC = if_else(Truth, true = True_log2FC, false = 0)) %>%
    # Add ref mean
    left_join(ref_pept %>% select(Peptide, Protein, Group, Mean),
              by = join_by(Peptide, Protein, Group))
  
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
    group_by(Group, Group2, Truth) %>%
    summarise(across(c(MSE, CIC, Diff_mean, Diff_LM,
                       pval, CI_width, Distinct, Signif,
                       MSE_PBT, MSE_LMT),
                     .fns = list('Mean' = ~mean(.x,na.rm = T),
                                 'Sd' = ~sd(.x,na.rm=T))),
              .groups = 'drop') %>%
    mutate('MSE_Mean' = sqrt(MSE_Mean), 'MSE_Sd' = sqrt(MSE_Sd),
           'MSE_PBT_Mean' = sqrt(MSE_PBT_Mean), 'MSE_PBT_Sd' = sqrt(MSE_PBT_Sd),
           'MSE_LMT_Mean' = sqrt(MSE_LMT_Mean), 'MSE_LMT_Sd' = sqrt(MSE_LMT_Sd)) %>%
    mutate(across(MSE_Mean:MSE_LMT_Sd, ~ round(.x, 2))) %>%
    reframe(Group, Group2, Truth,
            'PB_diff_mean' =  paste0(Diff_mean_Mean, ' (', Diff_mean_Sd, ')'),
            'CI_width' =  paste0(CI_width_Mean, ' (', CI_width_Sd, ')'),
            'Distinct' = paste0(Distinct_Mean, ' (', Distinct_Sd, ')'),
            'LM_diff_mean' =  paste0(Diff_LM_Mean, ' (', Diff_LM_Sd, ')'),
            'p_value' =  paste0(pval_Mean, ' (', pval_Sd, ')'),
            'Signif' = paste0(Signif_Mean, ' (', Signif_Sd, ')'),
            'RMSE' = paste0(MSE_Mean, ' (', MSE_Sd, ')'),
            'CIC' =  paste0(CIC_Mean, ' (', CIC_Sd, ')'),
            'RMSE_PBtoT' = paste0(MSE_PBT_Mean, ' (', MSE_PBT_Sd, ')'),
            'RMSE_LMtoT' = paste0(MSE_LMT_Mean, ' (', MSE_LMT_Sd, ')')) %>%
    left_join(db_results %>%
                select(Group, Group2, Truth, True_log2FC) %>%
                distinct,
              join_by(Group, Group2, Truth)) %>%
    mutate(True_diff_mean = round(True_log2FC, digits = 2),
           .after = 'Truth', .keep = "unused")
  
  return(list(DiffAna = db_eval %>%
                select(Group, Group2, Truth, Distinct, Signif),
              DiffMean = db_eval %>%
                select(Group, Group2, Truth,
                       True_diff_mean, PB_diff_mean, CI_width, LM_diff_mean,
                       RMSE_PBtoT, RMSE_LMtoT),
              EstimQual = db_eval %>%
                select(Group, Group2, Truth, True_diff_mean, RMSE, CIC)))
}

real_data_eval <- function(data, type, maxquant = T, normalize,
                           prop_NA = 0.2,
                           multi = F,
                           mu_0 = NULL, 
                           lambda_0 = 1e-10, beta_0 = 1, alpha_0 = 1,
                           alpha = 0.05, FDR = NULL){
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
                       PB_res = PB_res, LM_res = LM_res)
  
  return(list(
    data = db,
    ProteoBayes = PB_res,
    DAPAR = LM_res,
    results = results
  ))
  
}
