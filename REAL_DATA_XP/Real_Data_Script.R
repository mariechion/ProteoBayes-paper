# -- LIBRAIRIES -- # 
library(tidyverse)
library(ProteoBayes)
library(cp4p)
library(DAPAR)

# -- FUNCTIONS -- # 
source("REAL_DATA_XP/Functions.R")


# -- EXPERIMENTS -- #
prop_NA = 0.2
normalize = T
multi = F
mu_0 = NULL
lambda_0 = 1e-10
alpha_0 = 0.01
beta_0 = 0.3
alpha = 0.05
FDR = NULL

set.seed(17)

db_ARATH <- read.delim("REAL_DATA_XP/Arabido_UPS/peptides.txt")
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
                            FDR = FDR)

ARATH_DiffMean <- res_ARATH$results$DiffMean
ARATH_EstimQual <- res_ARATH$results$EstimQual

db_YST <- read.delim("REAL_DATA_XP/Yeast_UPS/peptides.txt")
res_YST <- real_data_eval(data = db_YST, type = "YST", maxquant = T,
                          normalize = normalize,
                          prop_NA = prop_NA,
                          multi = multi,
                          mu_0 = mu_0,
                          lambda_0 = lambda_0,
                          alpha_0 = alpha_0,
                          beta_0 = beta_0,
                          alpha = alpha,
                          FDR = FDR)

YST_DiffMean <- res_YST$results$DiffMean
YST_EstimQual <- res_YST$results$EstimQual

db_MOUSE <- read_delim("REAL_DATA_XP/Mouse_UPS/Spike-in-biol-var-OT-SN-Report.txt", 
                       delim = "\t", escape_double = FALSE, trim_ws = TRUE)
res_MOUSE <- real_data_eval(data = db_MOUSE, type = "MOUSE", maxquant = F,
                            normalize = normalize,
                            prop_NA = prop_NA,
                            multi = multi,
                            mu_0 = mu_0,
                            lambda_0 = lambda_0,
                            alpha_0 = alpha_0,
                            beta_0 = beta_0,
                            alpha = alpha,
                            FDR = FDR)

MOUSE_DiffMean <- res_MOUSE$results$DiffMean
MOUSE_EstimQual <- res_MOUSE$results$EstimQual


db_YST_B <- read.delim("REAL_DATA_XP/Yeast_UPS_B/peptides.txt")
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
                            FDR = FDR)

YST_B_DiffMean <- res_YST_B$results$DiffMean
YST_B_EstimQual <- res_YST_B$results$EstimQual

# -- COMBINED TABLES -- #

diff_mean <- rbind(YST_DiffMean %>% 
                     mutate(Experiment = "Muller_2016"),
                   YST_B_DiffMean %>% 
                     mutate(Experiment = "Bouyssie_2020"), 
                   MOUSE_DiffMean %>% 
                     mutate(Experiment = "Huang_2020"),
                   ARATH_DiffMean %>% 
                     mutate(Experiment = "Chion_2022")) %>% 
  relocate(Experiment)

#save(diff_mean, file = "REAL_DATA_XP/Diff_Mean_Results")
#write_csv(diff_mean, file = "REAL_DATA_XP/Diff_Mean_Results.csv")

estim_qual <- rbind(YST_EstimQual %>% 
                      mutate(Experiment = "Muller_2016"),
                   YST_B_EstimQual %>% 
                     mutate(Experiment = "Bouyssie_2020"), 
                   MOUSE_EstimQual %>% 
                     mutate(Experiment = "Huang_2020"),
                   ARATH_EstimQual %>% 
                     mutate(Experiment = "Chion_2022")) %>% 
  relocate(Experiment)

estim_qual %>% 
  separate(col = CIC, into = c("CIC", NA), sep = " ") %>% 
  mutate(CIC = as.numeric(CIC)) %>% 
  group_by(Experiment, Group) %>% 
  summarise(CIC_mean = mean(CIC)) %>% arrange(CIC_mean) %>% print(n = Inf)

#save(estim_qual, file = "REAL_DATA_XP/Estim_Qual_Results")
#write_csv(estim_qual, file = "REAL_DATA_XP/Estim_Qual_Results.csv")

combined_results <- diff_mean %>% 
  left_join(y = estim_qual, 
            by = c("Experiment", "Group", "Group2", "Truth", "True_diff_mean")) %>% 
  select(-RMSE_PBtoT, -RMSE_LMtoT) %>% 
  relocate(LM_diff_mean, .after = "True_diff_mean") %>% 
  mutate(Group = case_match(Group,
                            "Point1" ~ "0.05 fmol",
                            "Point2" ~ "0.25 fmol",
                            "Point3" ~ "0.5 fmol",
                            "Point4" ~ "1.25 fmol",
                            "Point5" ~ "2.5 fmol",
                            "Point6" ~ "5 fmol",
                            "S1" ~ "0.75 amol",
                            "S2" ~ "0.83 amol",
                            "S3" ~ "1.07 amol",
                            "S4" ~ "2.04 amol",
                            .default = Group),
         Group2 = case_match(Group2,
                             "Point7" ~ "10 fmol",
                             "S5" ~ "7.5 fmol",
                             .default = Group2),
         Type = case_when(
           Truth ~ "UPS",
           !Truth & Experiment %in% c("Muller_2016", "Bouyssie_2020") ~ "YEAST",
           !Truth & Experiment == "Huang_2020" ~ "MOUSE",
           !Truth & Experiment == "Chion_2022" ~ "ARATH",
         ), .before = "Truth") %>% 
  select(-Truth)


# save(combined_results, file = "REAL_DATA_XP/Combined_Results")
# write_csv(combined_results, file = "REAL_DATA_XP/Combined_Results.csv")
