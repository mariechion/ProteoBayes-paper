# -- LIBRAIRIES -- # 
library(tidyverse)
library(ProteoBayes)
library(cp4p)
library(DAPAR)

# -- FUNCTIONS -- # 
source("Experiments/Real_data/utils.R")

# -- REAL DATASETS -- #
db_ARATH <- read.delim("REAL_DATA_XP/Arabido_UPS/peptides.txt")
db_YST <- read.delim("REAL_DATA_XP/Yeast_UPS/peptides.txt")
db_MOUSE <- read_delim("REAL_DATA_XP/Mouse_UPS/Spike-in-biol-var-OT-SN-Report.txt", 
                       delim = "\t", escape_double = FALSE, trim_ws = TRUE)
db_YST_B <- read.delim("REAL_DATA_XP/Yeast_UPS_B/peptides.txt")

# -- EXPERIMENTS -- #
normalize = T
mu_0 = NULL
lambda_0 = 1e-10
alpha_0 = 0.01
beta_0 = 0.3
alpha = 0.05
FDR = NULL

# Experiment 1 : Comparison with limma
prop_NA = 0.2
multi = F

set.seed(17)

res_ARATH <- real_data_eval(data = db_ARATH, type = "ARATH", 
                            maxquant = T, 
                            normalize = normalize,
                            prop_NA = prop_NA,
                            multi = multi,
                            mu_0 = mu_0,
                            lambda_0 = lambda_0,
                            alpha_0 = alpha_0,
                            beta_0 = beta_0,
                            alpha = alpha,
                            FDR = FDR,
                            summary = F)

res_YST <- real_data_eval(data = db_YST, type = "YST", 
                          maxquant = T,
                          normalize = normalize,
                          prop_NA = prop_NA,
                          multi = multi,
                          mu_0 = mu_0,
                          lambda_0 = lambda_0,
                          alpha_0 = alpha_0,
                          beta_0 = beta_0,
                          alpha = alpha,
                          FDR = FDR,
                          summary = F)

res_MOUSE <- real_data_eval(data = db_MOUSE, type = "MOUSE", 
                            maxquant = F,
                            normalize = normalize,
                            prop_NA = prop_NA,
                            multi = multi,
                            mu_0 = mu_0,
                            lambda_0 = lambda_0,
                            alpha_0 = alpha_0,
                            beta_0 = beta_0,
                            alpha = alpha,
                            FDR = FDR,
                            summary = F)

res_YST_B <- real_data_eval(data = db_YST_B, type = "YST_B", 
                            maxquant = T, 
                            normalize = normalize,
                            prop_NA = prop_NA,
                            multi = multi,
                            mu_0 = mu_0,
                            lambda_0 = lambda_0,
                            alpha_0 = alpha_0,
                            beta_0 = beta_0,
                            alpha = alpha,
                            FDR = FDR,
                            summary = F)

write_csv(x = res_ARATH$results, file = "REAL_DATA_XP/Exp1_res_ARATH.csv")
write_csv(x = res_YST$results, file = "REAL_DATA_XP/Exp1_res_YST.csv")
write_csv(x = res_MOUSE$results, file = "REAL_DATA_XP/Exp1_res_MOUSE.csv")
write_csv(x = res_YST_B$results, file = "REAL_DATA_XP/Exp1_res_YST_B.csv")

## Output mean difference comparison table

bind_rows(
  res_YST$results %>% 
    group_by(True_diff_mean) %>% 
    summarise(Count = n(),
              PB_log2FC_mean = round(mean(PB_log2FC, na.rm = T), digits = 2),
              PB_log2FC_sd = round(sd(PB_log2FC, na.rm = T), digits = 2),
              LM_log2FC_mean = round(mean(LM_log2FC, na.rm = T), digits = 2),
              LM_log2FC_sd = round(sd(LM_log2FC, na.rm = T), digits = 2)) %>% 
    mutate(PB_Diff_Mean = paste0(PB_log2FC_mean," ", "(", PB_log2FC_sd, ")"),
           LM_Diff_Mean = paste0(LM_log2FC_mean," ", "(", LM_log2FC_sd, ")"),
           .keep = "unused") %>% 
    mutate(Experiment = "Muller2016", .before = True_diff_mean),
  res_MOUSE$results %>% 
    group_by(True_diff_mean) %>% 
    summarise(Count = n(),
              PB_log2FC_mean = round(mean(PB_log2FC, na.rm = T), digits = 2),
              PB_log2FC_sd = round(sd(PB_log2FC, na.rm = T), digits = 2),
              LM_log2FC_mean = round(mean(LM_log2FC, na.rm = T), digits = 2),
              LM_log2FC_sd = round(sd(LM_log2FC, na.rm = T), digits = 2)) %>% 
    mutate(PB_Diff_Mean = paste0(PB_log2FC_mean," ", "(", PB_log2FC_sd, ")"),
           LM_Diff_Mean = paste0(LM_log2FC_mean," ", "(", LM_log2FC_sd, ")"),
           .keep = "unused") %>% 
    mutate(Experiment = "Huang2020", .before = True_diff_mean),
  res_YST_B$results %>% 
    group_by(True_diff_mean) %>% 
    summarise(Count = n(),
              PB_log2FC_mean = round(mean(PB_log2FC, na.rm = T), digits = 2),
              PB_log2FC_sd = round(sd(PB_log2FC, na.rm = T), digits = 2),
              LM_log2FC_mean = round(mean(LM_log2FC, na.rm = T), digits = 2),
              LM_log2FC_sd = round(sd(LM_log2FC, na.rm = T), digits = 2)) %>% 
    mutate(PB_Diff_Mean = paste0(PB_log2FC_mean," ", "(", PB_log2FC_sd, ")"),
           LM_Diff_Mean = paste0(LM_log2FC_mean," ", "(", LM_log2FC_sd, ")"),
           .keep = "unused") %>% 
    mutate(Experiment = "Bouyssie2020", .before = True_diff_mean),
  res_ARATH$results %>% 
    group_by(True_diff_mean) %>% 
    summarise(Count = n(),
              PB_log2FC_mean = round(mean(PB_log2FC, na.rm = T), digits = 2),
              PB_log2FC_sd = round(sd(PB_log2FC, na.rm = T), digits = 2),
              LM_log2FC_mean = round(mean(LM_log2FC, na.rm = T), digits = 2),
              LM_log2FC_sd = round(sd(LM_log2FC, na.rm = T), digits = 2)) %>% 
    mutate(PB_Diff_Mean = paste0(PB_log2FC_mean," ", "(", PB_log2FC_sd, ")"),
           LM_Diff_Mean = paste0(LM_log2FC_mean," ", "(", LM_log2FC_sd, ")"),
           .keep = "unused") %>% 
    mutate(Experiment = "Chion2022", .before = True_diff_mean)
) %>% write_csv(file = "REAL_DATA_XP/Exp1_summary.csv")



# Experiment 2 : Evaluation of ProteoBayes performance
prop_NA = 1
multi = F

set.seed(17)

res_PB_ARATH <- real_data_eval(data = db_ARATH, type = "ARATH", 
                               maxquant = T, 
                               normalize = normalize,
                               prop_NA = prop_NA,
                               multi = multi,
                               mu_0 = mu_0,
                               lambda_0 = lambda_0,
                               alpha_0 = alpha_0,
                               beta_0 = beta_0, 
                               alpha = alpha, 
                               FDR = FDR,
                               summary = T)

res_PB_YST <- real_data_eval(data = db_YST, type = "YST", 
                             maxquant = T,
                             normalize = normalize,
                             prop_NA = prop_NA,
                             multi = multi,
                             mu_0 = mu_0,
                             lambda_0 = lambda_0,
                             alpha_0 = alpha_0,
                             beta_0 = beta_0,
                             alpha = alpha,
                             FDR = FDR,
                             summary = T)

res_PB_MOUSE <- real_data_eval(data = db_MOUSE, type = "MOUSE", 
                               maxquant = F,
                               normalize = normalize,
                               prop_NA = prop_NA,
                               multi = multi,
                               mu_0 = mu_0,
                               lambda_0 = lambda_0,
                               alpha_0 = alpha_0,
                               beta_0 = beta_0,
                               alpha = alpha,
                               FDR = FDR,
                               summary = T)


res_PB_YST_B <- real_data_eval(data = db_YST_B, type = "YST_B", 
                               maxquant = T, 
                               normalize = normalize,
                               prop_NA = prop_NA,
                               multi = multi,
                               mu_0 = mu_0,
                               lambda_0 = lambda_0,
                               alpha_0 = alpha_0,
                               beta_0 = beta_0,
                               alpha = alpha,
                               FDR = FDR, 
                               summary = T)

write_csv(x = res_PB_ARATH$results, file = "REAL_DATA_XP/Exp2_res_ARATH.csv")
write_csv(x = res_PB_YST$results, file = "REAL_DATA_XP/Exp2_res_YST.csv")
write_csv(x = res_PB_MOUSE$results, file = "REAL_DATA_XP/Exp2_res_MOUSE.csv")
write_csv(x = res_PB_YST_B$results, file = "REAL_DATA_XP/Exp2_res_YST_B.csv")

## Output ProteoBayes Evaluation table

res_Mean1 <- bind_rows(
  res_PB_YST$results %>% 
    filter(!is.na(Distinct)) %>% 
    group_by(True_diff_mean) %>% 
    summarise(MSE_mean = round(mean(sqrt(MSE), na.rm = T), digits = 2),
              MSE_sd = round(sd(sqrt(MSE), na.rm = T), digits = 2),
              CIC_mean = round(mean(CIC, na.rm = T), digits = 2),
              CIC_sd = round(sd(CIC, na.rm = T), digits = 2),
              Count = n()) %>% 
    mutate(RMSE = paste0(MSE_mean," ", "(", MSE_sd, ")"),
           CIC = paste0(CIC_mean," ", "(", CIC_sd, ")"),
           .keep = "unused") %>% 
    mutate(Experiment = "Muller2016", .before = True_diff_mean),
  res_PB_MOUSE$results %>% 
    filter(!is.na(Distinct)) %>% 
    group_by(True_diff_mean) %>% 
    summarise(MSE_mean = round(mean(sqrt(MSE), na.rm = T), digits = 2),
              MSE_sd = round(sd(sqrt(MSE), na.rm = T), digits = 2),
              CIC_mean = round(mean(CIC, na.rm = T), digits = 2),
              CIC_sd = round(sd(CIC, na.rm = T), digits = 2),
              Count = n()) %>% 
    mutate(RMSE = paste0(MSE_mean," ", "(", MSE_sd, ")"),
           CIC = paste0(CIC_mean," ", "(", CIC_sd, ")"),
           .keep = "unused") %>% 
    mutate(Experiment = "Huang2020", .before = True_diff_mean),
  res_PB_YST_B$results %>% 
    filter(!is.na(Distinct)) %>% 
    group_by(True_diff_mean) %>% 
    summarise(MSE_mean = round(mean(sqrt(MSE), na.rm = T), digits = 2),
              MSE_sd = round(sd(sqrt(MSE), na.rm = T), digits = 2),
              CIC_mean = round(mean(CIC, na.rm = T), digits = 2),
              CIC_sd = round(sd(CIC, na.rm = T), digits = 2),
              Count = n()) %>% 
    mutate(RMSE = paste0(MSE_mean," ", "(", MSE_sd, ")"),
           CIC = paste0(CIC_mean," ", "(", CIC_sd, ")"),
           .keep = "unused") %>% 
    mutate(Experiment = "Bouyssie2020", .before = True_diff_mean),
  res_PB_ARATH$results %>% 
    filter(!is.na(Distinct)) %>% 
    group_by(True_diff_mean) %>% 
    summarise(MSE_mean = round(mean(sqrt(MSE), na.rm = T), digits = 2),
              MSE_sd = round(sd(sqrt(MSE), na.rm = T), digits = 2),
              CIC_mean = round(mean(CIC, na.rm = T), digits = 2),
              CIC_sd = round(sd(CIC, na.rm = T), digits = 2),
              Count = n()) %>% 
    mutate(RMSE = paste0(MSE_mean," ", "(", MSE_sd, ")"),
           CIC = paste0(CIC_mean," ", "(", CIC_sd, ")"),
           .keep = "unused") %>% 
    mutate(Experiment = "Chion2022", .before = True_diff_mean)
) %>% filter(!is.na(True_diff_mean)) %>% 
  rename("True" = "True_diff_mean") %>%  write_csv(file = "REAL_DATA_XP/Exp2_summary.csv")



# Highlighted experiment: Muller 2016

res_YST$results %>% 
  group_by(True_diff_mean) %>% 
  summarise(Count = n(),
            LM_log2FC_mean = round(mean(LM_log2FC, na.rm = T), digits = 2),
            LM_log2FC_sd = round(sd(LM_log2FC, na.rm = T), digits = 2),
            PB_log2FC_mean = round(mean(PB_log2FC, na.rm = T), digits = 2),
            PB_log2FC_sd = round(sd(PB_log2FC, na.rm = T), digits = 2),
            CI_width_mean = round(mean(CI_width, na.rm = T), digits = 2),
            CI_width_sd = round(sd(CI_width, na.rm = T), digits = 2),
            MSE_mean = round(mean(sqrt(MSE), na.rm = T), digits = 2),
            MSE_sd = round(sd(sqrt(MSE), na.rm = T), digits = 2),
            CIC_mean = round(mean(CIC, na.rm = T), digits = 2),
            CIC_sd = round(sd(CIC, na.rm = T), digits = 2)) %>% 
  mutate(LM_Diff_Mean = paste0(LM_log2FC_mean," ", "(", LM_log2FC_sd, ")"),
         PB_Diff_Mean = paste0(PB_log2FC_mean," ", "(", PB_log2FC_sd, ")"),
         CI_width =  paste0(CI_width_mean, ' (', CI_width_sd, ')'),
         RMSE = paste0(MSE_mean," ", "(", MSE_sd, ")"),
         CIC = paste0(CIC_mean, " ", "(", CIC_sd, ")"),
         .keep = "unused") 

  