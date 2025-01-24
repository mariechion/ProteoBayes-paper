# -- LIBRAIRIES -- # 
library(tidyverse)
library(ProteoBayes)
library(cp4p)
library(DAPAR)

# -- FUNCTIONS -- # 
source("REAL_DATA_XP/Functions.R")

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

## Output mean difference comparison table

bind_rows(
  res_YST$results %>% 
    group_by(True_diff_mean) %>% 
    summarise(PB_log2FC_mean = round(mean(PB_log2FC, na.rm = T), digits = 2),
              PB_log2FC_sd = round(sd(PB_log2FC, na.rm = T), digits = 2),
              LM_log2FC_mean = round(mean(LM_log2FC, na.rm = T), digits = 2),
              LM_log2FC_sd = round(sd(LM_log2FC, na.rm = T), digits = 2)) %>% 
    mutate(PB_Diff_Mean = paste0(PB_log2FC_mean," ", "(", PB_log2FC_sd, ")"),
           LM_Diff_Mean = paste0(LM_log2FC_mean," ", "(", LM_log2FC_sd, ")"),
           .keep = "unused") %>% 
    mutate(Experiment = "Muller2016", .before = True_diff_mean),
  res_MOUSE$results %>% 
    group_by(True_diff_mean) %>% 
    summarise(PB_log2FC_mean = round(mean(PB_log2FC, na.rm = T), digits = 2),
              PB_log2FC_sd = round(sd(PB_log2FC, na.rm = T), digits = 2),
              LM_log2FC_mean = round(mean(LM_log2FC, na.rm = T), digits = 2),
              LM_log2FC_sd = round(sd(LM_log2FC, na.rm = T), digits = 2)) %>% 
    mutate(PB_Diff_Mean = paste0(PB_log2FC_mean," ", "(", PB_log2FC_sd, ")"),
           LM_Diff_Mean = paste0(LM_log2FC_mean," ", "(", LM_log2FC_sd, ")"),
           .keep = "unused") %>% 
    mutate(Experiment = "Huang2020", .before = True_diff_mean),
  res_YST_B$results %>% 
    group_by(True_diff_mean) %>% 
    summarise(PB_log2FC_mean = round(mean(PB_log2FC, na.rm = T), digits = 2),
              PB_log2FC_sd = round(sd(PB_log2FC, na.rm = T), digits = 2),
              LM_log2FC_mean = round(mean(LM_log2FC, na.rm = T), digits = 2),
              LM_log2FC_sd = round(sd(LM_log2FC, na.rm = T), digits = 2)) %>% 
    mutate(PB_Diff_Mean = paste0(PB_log2FC_mean," ", "(", PB_log2FC_sd, ")"),
           LM_Diff_Mean = paste0(LM_log2FC_mean," ", "(", LM_log2FC_sd, ")"),
           .keep = "unused") %>% 
    mutate(Experiment = "Bouyssie2020", .before = True_diff_mean),
  res_ARATH$results %>% 
    group_by(True_diff_mean) %>% 
    summarise(PB_log2FC_mean = round(mean(PB_log2FC, na.rm = T), digits = 2),
              PB_log2FC_sd = round(sd(PB_log2FC, na.rm = T), digits = 2),
              LM_log2FC_mean = round(mean(LM_log2FC, na.rm = T), digits = 2),
              LM_log2FC_sd = round(sd(LM_log2FC, na.rm = T), digits = 2)) %>% 
    mutate(PB_Diff_Mean = paste0(PB_log2FC_mean," ", "(", PB_log2FC_sd, ")"),
           LM_Diff_Mean = paste0(LM_log2FC_mean," ", "(", LM_log2FC_sd, ")"),
           .keep = "unused") %>% 
    mutate(Experiment = "Chion2022", .before = True_diff_mean)
) 



# Experiment 2 : Evaluation of ProteoBayes performance
prop_NA = 0.2
multi = F
source("REAL_DATA_XP/Functions.R")
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
                               summary = F)


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
                             summary = F)

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
                               summary = F)


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
                               summary = F)


## Output ProteoBayes Evaluation table

res_Mean1 <- bind_rows(
  res_PB_YST$results %>% 
    # filter(!is.na(Distinct)) %>% 
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
  rename("True" = "True_diff_mean") %>%  view()


# Check results for Bouyssie2020 using ProteoBayes native functions
# (not the ones from Functions.R)

data_B20 <- db_YST_B %>% 
  select(Sequence, Leading.razor.protein, starts_with("Intensity.")) %>% 
  mutate(across(starts_with("Intensity."), ~ if_else(.x == 0, 
                                                      true = NA, 
                                                      false = log2(.x))))

# Normalise dataset
data_B20[,-c(1,2)] <- preprocessCore::normalize.quantiles(as.matrix(data_B20[,-c(1,2)]),
                                                     copy = F)

# Transform into ProteoBayes form
data_B20_PB <- data_B20 %>% 
  pivot_longer(-c("Sequence", "Leading.razor.protein"), 
               names_to = c("Group", "Sample"), names_sep = "_",
               values_to = "Output") %>% 
  mutate(Group = str_replace(Group, "Intensity.", "")) %>% 
  rename("Peptide" = "Sequence", "Protein" = "Leading.razor.protein")

# Perform ProteoBayes
set.seed(17)
out_B20_PB <- posterior_mean(data = data_B20_PB,
                             lambda_0 = 1e-10,
                             alpha_0 = 0.01,
                             beta_0 = 0.3)

diff_B20_PB <- identify_diff(out_B20_PB) 

# Create results table
res_B20_PB <- diff_B20_PB %>% 
  # Add Protein Name
  left_join(y = data_B20_PB %>% 
              select(Peptide, Protein) %>% 
              distinct(),
            by = "Peptide") %>% 
  relocate(Protein, .after = Peptide) %>% 
  # Add flag to see which peptide is supposed to be differentially expressed
  mutate(Truth = if_else(str_detect(Protein, "ups"),
                         true = T, false = F))



# Experimental condition/design to calculate true FC
group_labels = c("10amol", "50amol", "100amol", "250amol", "500amol",
                 "1fmol", "5fmol", "10fmol", "25fmol","50fmol")
fmol_labels = c(0.01,0.05,0.1,0.25,0.5, 1, 5, 10, 25, 50)

exp_cond <- tibble(Group = group_labels,
                   fmol = fmol_labels) %>%
  expand_grid(Group2 = unique(Group)) %>%
  filter(Group != Group2) %>%
  left_join(tibble(Group2 = group_labels,
                   fmol2 = fmol_labels),
            by = "Group2") %>%
  mutate(True_log2FC = log2(fmol/fmol2), .keep = "unused")

# Add reference mean 

ref_pept <- data_B20_PB %>%
  arrange(Peptide) %>%
  left_join(tibble(Group = group_labels,
                   fmol = fmol_labels),
            by = "Group") %>%
  mutate(True_log2FC = log2(max(fmol)/fmol),
         Mean = if_else(str_detect(Protein, "ups"),
                        true = Output + True_log2FC,
                        false = Output)) %>% 
  group_by(Peptide, Protein) %>%
  mutate(Mean = if_else(str_detect(Protein, "ups"),
                        true = mean(Mean, na.rm = T) - True_log2FC,
                        false = mean(Mean, na.rm = T))) %>% 
  group_by(Peptide, Protein, Group, Mean) %>%
  summarise(Avg = mean(Output)) %>%
  ungroup %>%
  mutate(Diff = Mean - Avg, Case = str_detect(Protein, "ups"))

db_res_B20 <- res_B20_PB %>%
  filter(Group2 == tail(group_labels, n=1)) %>%
  # Add ref mean
  left_join(ref_pept %>% 
              select(Peptide, Protein, Group, Mean) %>% 
              distinct,
            by = join_by(Peptide, Protein, Group2 == Group)) %>% 
  # Add fmol equivalent to groups
  left_join(y = exp_cond, by = join_by(Group, Group2)) %>%
  # Add FC from ProteoBayes and true FC
  mutate(PB_log2FC = mu2 - mu,
         True_log2FC = if_else(Truth, true = -True_log2FC, false = 0))
  
db_res_B20 %>% 
  mutate('MSE' = (Mean - mu2)^2,
         'CIC' = ((Mean > CI_inf2) & (Mean < CI_sup2)) * 100) %>% 
  group_by(Group, Group2, Truth) %>%
  summarise(across(c(MSE, CIC),
                   .fns = list('Mean' = ~mean(.x,na.rm = T),
                               'Sd' = ~sd(.x,na.rm=T))),
            .groups = 'drop')
  
