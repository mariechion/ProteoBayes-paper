# -- LIBRAIRIES -- # 
library(tidyverse)
library(ProteoBayes)
library(cp4p)
library(mi4p)

# -- FUNCTIONS -- #
PB_preprocess <- function(data, nb_group, max_NA, dapar = F){
  ## Data preprocessing
  db <- data %>% 
    ## Select columns of interest
    select(Sequence, Leading.razor.protein, starts_with("Intensity.")) %>% 
    ## Replace 0 intensity values by NA
    mutate(across(starts_with("Intensity."), ~ if_else(.x == 0, 
                                                       true = NA, 
                                                       false = .x))) %>% 
    ## Reshape data in long format
    gather(key = "Condition", value = "Intensity", 
           -c(Sequence, Leading.razor.protein)) %>% 
    ## Transform the key into Group and Sample
    separate(col = Condition, into = c("Group","Sample"), sep = "_") %>% 
    mutate(Group = str_replace(Group, "Intensity.TG.", "")) %>% 
    ## Rename columns to match ProteoBayes requirements
    rename(Peptide = Sequence, 
           Protein = Leading.razor.protein,
           Output = Intensity) %>%
    ## Remove reverse and contaminant proteins, remove iRT. 
    filter(!str_detect(Protein, "CON") & !str_detect(Protein, "REV") &
             !str_detect(Protein, "iRT"))
  
  ## Get quantified peptides: 
  ## *max_NA* missing values in each of the *nb_group* groups.
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
             !is.na(Output)) %>% 
    ## Log2-transform intensity values
    mutate(Output = log2(Output))
  
  return(db)
}

PB_DiffAna <- function(data, multi = F, 
                       mu_0 = NULL, lambda_0 = 1, beta_0 = 1, alpha_0 = 1){
  if (!multi) {
    out <- posterior_mean(data = db_YST,
                          mu_0 = mu_0,
                          lambda_0 = lambda_0, 
                          beta_0 = beta_0, 
                          alpha_0 = alpha_0)
    
    diff_out <- identify_diff(out)
  }
  
  return(diff_out)
}

LM_DiffAna <- function(data, alpha = 0.05, FDR = NULL){
  # Note, limma and DAPAR use data in wide format.
  db_YST_limma <- db_YST %>% 
    # Merge Group and Sample columns to denote biological samples analysed
    unite(col = "Cond_Rep", c(Group,Sample),sep = "_") %>% 
    # Reshape data in wide format
    spread(key = "Cond_Rep", value = "Output")
  
  # Create quantitative data matrix to match DAPAR requirements
  qYST <- db_YST_limma %>%
    # Tibble to dataframe to add row names to be propagated in results tables
    as.data.frame(row.names = pull(db_YST_limma %>% select(Peptide))) %>% 
    select(-c(Peptide,Protein)) 
  
  # Create design dataframe to match DAPAR requirements
  metadata <- data.frame(Sample.name = colnames(qYST),
                         Condition = rep(c("0.5fmol", "1fmol", "2.5fmol",
                                           "5fmol", "10fmol", "25fmol"), 
                                         rep(3,6)),
                         Bio.Rep = 1:ncol(qYST))
  
  # Moderated t-test - OnevsOne setting
  #res_DAPAR <- limmaCompleteTest(qData = as.matrix(qARATH),
  #                               sTab = metadata,
  #                               comp.type = "OnevsOne")
  # Issue with the formatting of the output of DAPAR function, 
  # Therefore, using the one in mi4p.
  res_limma <- mi4p::limmaCompleteTest.mod(qData = as.matrix(qYST),
                                           sTab = metadata,
                                           comp.type = "OnevsOne")$res.l
  
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

CombineDA <- function(data, PB_res, LM_res){
  
  # Experimental condition/design to calculate true FC
  exp_cond <- tibble(Group = c("0.5fmol", "1fmol", "2.5fmol",
                               "5fmol", "10fmol", "25fmol"),
                     fmol = c(0.5, 1, 2.5, 5, 10, 25)) %>% 
    expand_grid(Group2 = unique(Group)) %>% 
    filter(Group != Group2) %>% 
    left_join(tibble(Group2 = c("0.5fmol", "1fmol", "2.5fmol",
                                "5fmol", "10fmol", "25fmol"),
                     fmol2 = c(0.5, 1, 2.5, 5, 10, 25)), 
              by = "Group2") %>% 
    mutate(True_log2FC = log2(fmol/fmol2), .keep = "unused")
  
  # Create a "reference mean" 
  ref_pept <- data %>% 
    arrange(Peptide) %>% 
    mutate(fmol = as.numeric(str_replace(Group, "fmol","")),
           Mean = if_else(str_detect(Protein, "ups"), 
                          true = Output + log2(10/fmol), 
                          false = Output)) %>% 
    group_by(Peptide, Protein) %>% 
    mutate(Mean = if_else(str_detect(Protein, "ups"), 
                          true = mean(Mean, na.rm = T) - log2(10/fmol),
                          false = mean(Mean, na.rm = T))) %>% 
    group_by(Peptide, Protein, Group, Mean) %>% 
    summarise(Avg = mean(Output)) %>% 
    ungroup %>% 
    mutate(Diff = Mean - Avg, Case = str_detect(Protein, "UPS"))
  
  # Merge result tables
  db_results <- PB_res %>%
    filter(Group2 == "25fmol") %>% 
    # Add Protein Column
    left_join(y = .,
              x = data %>%
                select(Peptide,Protein) %>%
                distinct,
              by = "Peptide") %>% 
    # Add limma results
    left_join(y = LM_res %>% 
                 filter(Group2 == "25fmol") %>% 
                 rename(LM_log2FC = log2FC),
               by = join_by(Peptide, Group, Group2)) %>% 
    # Add ground truth column
    mutate(Truth = if_else(str_detect(Protein, "ups"), 
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
            'CI_width' =  paste0(CIC_width_Mean, ' (', CIC_width_Sd, ')'),
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
                select(Group, Group2, Truth, RMSE, CIC)))
  
} 


# -- DATA ANALYSIS -- #

# Set random generator
set.seed(17)

# Load peptide-level data
peptides <- read.delim("DATA/Yeast_UPS/peptides.txt")

# Preprocess data for ProteoBayes setting
db_YST <- PB_preprocess(peptides, max_NA = 2, nb_group = 6)

# ProteoBayes - Univariate setting
diff_PB <- PB_DiffAna(data = db_YST, 
                      mu_0 = db_YST %>% 
                        group_by(Group) %>% 
                        mutate(mu_0 = mean(Output)) %>% pull(mu_0))

diff_PB <- PB_DiffAna(db_YST)

diff_PB <- posterior_mean(db_YST, lambda_0 = 2)

db_YST %>%
  posterior_mean(lambda_0 = 2) %>%
  sample_distrib() %>% 
  plot_distrib(group1 = "1fmol", group2 = "25fmol", 
               #peptide = "AAFTECCQAADK") 
               peptide = "ADGLAVIGVLMK") 

diff_PB <- PB_DiffAna(db_YST, 
                      lambda_0 = 2)

diff_PB <- PB_DiffAna(db_YST, 
                      lambda_0 = 0.5)

diff_PB <- PB_DiffAna(db_YST, 
                      lambda_0 = 0.1)

diff_PB <- PB_DiffAna(db_YST, 
                      lambda_0 = 0.01)

diff_PB <- PB_DiffAna(db_YST, 
                      lambda_0 = 0.001)

diff_PB <- PB_DiffAna(db_YST, 
                      alpha_0 = 2)

diff_PB <- PB_DiffAna(db_YST, 
                      alpha_0 = 100)

diff_PB <- PB_DiffAna(db_YST,
                      lambda_0 = 0.01,
                      alpha_0 = 0.1)

diff_PB <- PB_DiffAna(db_YST,
                      lambda_0 = 0.01,
                      alpha_0 = 0.1,
                      beta_0 = 0.1)

diff_PB <- PB_DiffAna(db_YST, 
                      beta_0 = 100)


diff_PB <- PB_DiffAna(db_YST, 
                      alpha_0 = 5)

# DAPAR
diff_LM <- LM_DiffAna(db_YST, FDR = 0.01)

# Combine results
db_eval_YST <- CombineDA(db_YST, diff_PB, diff_LM)

db_eval_YST$DiffAna %>% view()
db_eval_YST$DiffMean %>% view()
db_eval_YST$PBPerf %>% view()



# --- BULK and DRAFT --- #

Start_time <- Sys.time()
res_uni_YST <- posterior_mean(data = db_YST,
                              mu_0 = db_YST %>% 
                                group_by(Group) %>% 
                                mutate(mu_0 = mean(Output)) %>% pull(mu_0),
                              lambda_0 = 2)
End_time <- Sys.time()
duration_uni_YST <- End_time - Start_time

diff_uni_YST <- identify_diff(res_uni_YST)


## Method of moments to estimate alpha and beta
# db_YST %>% 
#   group_by(Group, Peptide) %>% 
#   summarise(sigma2_inv = 1/(sd(Output)^2)) %>% 
#   group_by(Group) %>% 
#   summarise(alpha_0 = mean(sigma2_inv, 
#                            na.rm = T)^2/(sd(sigma2_inv, na.rm = T)^2),
#             beta_0 = mean(sigma2_inv, na.rm = T)/(sd(sigma2_inv, na.rm = T)^2))
  

# Moderated t-test - using DAPAR
# Note, limma and DAPAR use data in wide format.
db_YST_limma <- db_YST %>% 
  # Merge Group and Sample columns to denote biological samples analysed
  unite(col = "Cond_Rep", c(Group,Sample),sep = "_") %>% 
  # Reshape data in wide format
  spread(key = "Cond_Rep", value = "Output")

# Create quantitative data matrix to match DAPAR requirements
qYST <- db_YST_limma %>%
  # Tibble to dataframe to add row names to be propagated in results tables
  as.data.frame(row.names = pull(db_YST_limma %>% select(Peptide))) %>% 
  select(-c(Peptide,Protein)) 

# Create design dataframe to match DAPAR requirements
metadata <- data.frame(Sample.name = colnames(qYST),
                       Condition = rep(c("0.5fmol", "1fmol", "2.5fmol",
                                         "5fmol", "10fmol", "25fmol"), rep(3,6)),
                       Bio.Rep = 1:ncol(qYST))
                       
# Moderated t-test - OnevsOne setting
Start_time <- Sys.time()
#res_DAPAR <- limmaCompleteTest(qData = as.matrix(qARATH),
#                               sTab = metadata,
#                               comp.type = "OnevsOne")
# Issue with the formatting of the output of DAPAR function, 
# Therefore, using the one in mi4p.
res_limma <- mi4p::limmaCompleteTest.mod(qData = as.matrix(qYST),
                              sTab = metadata,
                              comp.type = "OnevsOne")$res.l
End_time <- Sys.time()
duration_limma <- End_time - Start_time

res_limma_adj <- res_limma$P_Value %>% 
  rownames_to_column() %>% 
  as_tibble %>% 
  mutate(across(starts_with("Point"), 
                ~ cp4p::adjust.p(p = ., alpha = 0.01)$adjp %>% 
                  select(adjusted.p) %>% pull))

res_limma_FC <- res_limma$logFC %>% 
  rownames_to_column(var = "Peptide") %>% 
  gather(key = "Comparison", 
         value = "log2FC", -Peptide) %>% 
  separate(col = Comparison, 
           into = c("Comparison", NA), 
           sep = "_logFC")

# -- Results -- #

exp_cond <- tibble(Group = c("0.5fmol", "1fmol", "2.5fmol",
                   "5fmol", "10fmol", "25fmol"),
                   fmol = c(0.5, 1, 2.5, 5, 10, 25)) %>% 
  expand_grid(Group2 = unique(Group)) %>% 
  filter(Group != Group2) %>% 
  left_join(tibble(Group2 = c("0.5fmol", "1fmol", "2.5fmol",
                              "5fmol", "10fmol", "25fmol"),
                   fmol2 = c(0.5, 1, 2.5, 5, 10, 25)), 
            by = "Group2")

ref_pept <- db_YST %>% 
  arrange(Peptide) %>% 
  mutate(fmol = as.numeric(str_replace(Group, "fmol","")),
         Mean = if_else(str_detect(Protein, "ups"), 
                          true = Output + log2(10/fmol), 
                          false = Output)) %>% 
  group_by(Peptide, Protein) %>% 
  mutate(Mean = if_else(str_detect(Protein, "ups"), 
                        true = mean(Mean, na.rm = T) - log2(10/fmol),
                        false = mean(Mean, na.rm = T))) %>% 
  group_by(Peptide, Protein, Group, Mean) %>% 
  summarise(Avg = mean(Output)) %>% 
  ungroup %>% 
  mutate(Diff = Mean - Avg, Case = str_detect(Protein, "UPS"))

ggplot(data = ref_pept, aes(x = Group, y = Diff, fill = Case)) +
  geom_boxplot()


db_results <- diff_uni_YST %>%
  rename(ProteoBayes = Distinct) %>% 
  filter(Group2 == "25fmol") %>% 
  # Add Protein Column
  left_join(y = .,
            x = db_YST %>%
              select(Peptide,Protein) %>%
              distinct,
            by = "Peptide") %>% 
  # Add pval from limma
  right_join(y = res_limma_adj %>% 
              pivot_longer(cols = - rowname, 
                           names_to = "Comparison", 
                           values_to = "LM_pval") %>% 
              separate(col = Comparison, 
                       into = c("Group", NA, "Group2", NA), 
                       sep = "_"),
            by = join_by(Peptide == rowname, Group, Group2)) %>% 
  # Add column to flag differential peptides as in limma
  # and ground truth
  mutate(Limma = if_else(LM_pval < 0.05, true = T, false = F, missing = NA),
         Truth = if_else(str_detect(Protein, "ups"), true = T, false = F)
  ) %>% 
  filter(!is.na(ProteoBayes) & !is.na(Limma)) %>% 
  # Add FC from limma
  left_join(y = res_limma$logFC %>% 
              rownames_to_column(var = "Peptide") %>% 
              pivot_longer(cols = - Peptide, 
                           names_to = "Comparison", 
                           values_to = "LM_log2FC") %>% 
              separate(col = Comparison, 
                       into = c("Group", NA, "Group2", NA), 
                       sep = "_"),
            by = join_by(Peptide, Group, Group2)) %>% 
  # Add fmol equivalent to groups
  left_join(y = exp_cond, by = join_by(Group, Group2)) %>% 
  # Add FC from ProteoBayes and true FC
  mutate(PB_log2FC = mu - mu2,
         True_log2FC = if_else(Truth, true = log2(fmol/fmol2), false = 0)) %>% 
  # Add ref mean
  left_join(ref_pept %>% select(Peptide, Protein, Group, Mean), 
            by = join_by(Peptide, Protein, Group))

db_eval <- db_results %>% 
  mutate('MSE' = (Mean - mu)^2,
         'CIC' = ((Mean > CI_inf) & (Mean < CI_sup)) * 100,
         'Diff_mean' = PB_log2FC,
         'Diff_LM' = LM_log2FC,
         'Diff_LM2' = LM_log2FC/log(2),
         'CIC_width' = CI_sup - CI_inf, 
         'Distinct' = (ProteoBayes == Truth)*100,
         'Signif' = (Limma == Truth)*100,
         "MSE_PBT" = (PB_log2FC - True_log2FC)^2,
         "MSE_LMT" = (LM_log2FC - True_log2FC)^2) %>% 
  group_by(Group, Group2, Truth) %>% 
  summarise(across(c(MSE, CIC, Diff_mean, Diff_LM, Diff_LM2,
                     LM_pval, CIC_width, Distinct, Signif,
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
          'CIC_width' =  paste0(CIC_width_Mean, ' (', CIC_width_Sd, ')'),
          'Distinct' = paste0(Distinct_Mean, ' (', Distinct_Sd, ')'),
          'LM_diff_mean' =  paste0(Diff_LM_Mean, ' (', Diff_LM_Sd, ')'),
          'LM_diff_mean_2' =  paste0(Diff_LM2_Mean, ' (', Diff_LM2_Sd, ')'),
          'p_value' =  paste0(LM_pval_Mean, ' (', LM_pval_Sd, ')'),
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


# Perfomance indicators for ProteoBayes
perf_PB <- db_results %>% 
  as.data.frame() %>% 
  mutate(Truth = as.factor(Truth),
         ProteoBayes = as.factor(ProteoBayes)) %>% 
  group_by(Group,Group2) %>% 
  #group_by(Comparison) %>%
  group_modify(.f = ~ conf_mat(data = .x,
                               truth = Truth,
                               estimate = ProteoBayes) %>%
                 tidy %>%
                 mutate(name = recode(name,
                                      "cell_1_1" = "TN",
                                      "cell_2_1" = "FP",
                                      "cell_1_2" = "FN",
                                      "cell_2_2" = "TP"))) %>%
  spread(key = "name", value = "value") %>%
  mutate(Sensitivity = round(TP/(TP + FN)*100, digits = 1),
         Specificity = round(TN/(TN + FP)*100, digits = 1),
         Accuracy = round((TP + TN)/(TP + TN + FP + FN)*100, digits = 1),
         F1Score = round(2*TP/(2*TP+FP+FN)*100, digits = 1))

# Perfomance indicators for Limma
perf_LM <- db_results %>% 
  as.data.frame() %>% 
  mutate(Truth = as.factor(Truth),
         Limma = as.factor(Limma)) %>% 
  group_by(Group, Group2) %>%
  group_modify(.f = ~ conf_mat(data = .x,
                               truth = Truth,
                               estimate = Limma) %>%
                 tidy %>%
                 mutate(name = recode(name,
                                      "cell_1_1" = "TN",
                                      "cell_2_1" = "FP",
                                      "cell_1_2" = "FN",
                                      "cell_2_2" = "TP"))) %>%
  spread(key = "name", value = "value") %>%
  mutate(Sensitivity = round(TP/(TP + FN)*100, digits = 1),
         Specificity = round(TN/(TN + FP)*100, digits = 1),
         Accuracy = round((TP + TN)/(TP + TN + FP + FN)*100, digits = 1),
         F1Score = round(2*TP/(2*TP+FP+FN)*100, digits = 1))

db_perf <- bind_rows(ProteoBayes = perf_PB, Limma = perf_LM, .id = "Method") %>% 
  arrange(Comparison) %>% relocate(Comparison)

db_ARATH %>% 
  mutate(Is_diff = if_else(str_detect(Protein, "UPS"), true = T, false = F)) %>% 
  group_by(Group, Is_diff) %>% 
  summarise(Count = n())

# Save results
# save(res_uni_YST,
#      file = "DATA/Yeast_UPS/Yeast_UPS_PB_res")
# save(diff_uni_YST,
#      file = "DATA/Yeast_UPS/Yeast_UPS_PB_diff")
# save(res_limma,
#      file = "DATA/Yeast_UPS/Yeast_UPS_LM_res")
# save(db_perf,
#      file = "DATA/Yeast_UPS/Yeast_UPS_performance")

# Load results
load(file = "DATA/Arabido_UPS/Arabido_UPS_PB_res")
load(file = "DATA/Arabido_UPS/Arabido_UPS_PB_diff")
load(file = "DATA/Arabido_UPS/Arabido_UPS_LM_res")
load(file = "DATA/Arabido_UPS/Arabido_UPS_performance")


sample_post <- sample_distrib(res_uni_ARATH)

test_post <- sample_post %>% slice(1:100) %>%  
  group_by(Peptide) %>% 
  reframe(pivot_wider(., names_from = Group, values_from = Sample)) %>% rbind

# ----------------------------------------------------------- #
load(file = "Arabido_UPS_1of3inall_ImpMLE_long")
load(file = "Arabido_UPS_1of3inall_ImpMLE_res_dapar")
head(data.lg)

db_ARATH <- data.lg %>% 
  rename(Draw = Imp.Draw, Peptide = ID, Output = Intensity)
db_ARATH %>% head

# - Univariate approach - #

start_uni <- Sys.time()
res_uni_ARATH <- posterior_mean(db_ARATH)
end_uni <- Sys.time()
duration_uni <- end_uni - start_uni #15 sec env

start_sample_uni <- Sys.time()
sample_uni_ARATH <- sample_distrib(res_uni_ARATH)
end_sample_uni <- Sys.time()
duration_sample_uni <- end_sample_uni - start_sample_uni #38 sec env

# Add Protein to results
res_uni_ARATH <- right_join(x = res_uni_ARATH,
                            y = db_ARATH %>% 
                              select(Peptide,Protein) %>% 
                              distinct,
                            by = "Peptide")

sample_uni_ARATH <- right_join(x = sample_uni_ARATH,
                         y = db_ARATH %>% 
                           select(Peptide,Protein) %>% 
                           distinct,
                         by = "Peptide")

# - Multivariate approach - #
start_multi <- Sys.time()
res_multi_ARATH <- db_ARATH %>% 
  filter(Protein %in% (db_ARATH %>% 
                          select(Protein) %>% 
                          distinct %>% 
                          pull(Protein))) %>% 
  group_by(Protein) %>% 
  group_modify(.f = ~ multi_posterior_mean(.x, nu_0 = 30,
                                           vectorised = T))
end_multi <- Sys.time()
duration_multi <- end_multi - start_multi #56 min -> 9min

start_sample_multi <- Sys.time()
sample_multi_ARATH <- sample_distrib(res_multi_ARATH)
end_sample_multi <- Sys.time()
duration_sample_multi <- end_sample_multi - start_sample_multi

# Error in sample_distrib(res_multi_ARATH) : 
# The 'nu' parameter is too small compared to the number of Peptides 
# (nu < P - 1). Consider changing prior value for 'nu' or decreasing the number 
# of Peptides.

# nu = nu + 50 (mutate)

# Apply multivariate framework on each protein to benefit from intra-protein 
# correlation.
# However, when there is only 1 peptide in the protein group, there is no 
# sense in using a multivariate framework. 
# BUT we have imputed data --> How to manage the univariate framework ? 
# --> Uses more
