# -- LIBRAIRIES -- # 
library(tidyverse)
library(ProteoBayes)
library(cp4p)
library(DAPAR)

# -- FUNCTIONS -- # 
source("REAL_DATA_XP/Functions.R")


# -- EXPERIMENTS -- #
set.seed(17)

prop_NA = 0.2
multi = F
mu_0 = NULL
lambda_0 = 1e-10
beta_0 = 0.01
alpha_0 = 0.01
alpha = 0.05
FDR = NULL

db_ARATH <- read.delim("REAL_DATA_XP/Arabido_UPS/peptides.txt")
res_ARATH <- real_data_eval(data = db_ARATH, type = "ARATH", maxquant = T)

ARATH_DiffMean <- res_ARATH$results$DiffMean
ARATH_EstimQual <- res_ARATH$results$EstimQual

db_YST <- read.delim("REAL_DATA_XP/Yeast_UPS/peptides.txt")
res_YST <- real_data_eval(data = db_YST, type = "YST", maxquant = T)

YST_DiffMean <- res_YST$results$DiffMean
YST_EstimQual <- res_YST$results$EstimQual

db_MOUSE <- read_delim("REAL_DATA_XP/Mouse_UPS/Spike-in-biol-var-OT-SN-Report.txt", 
                       delim = "\t", escape_double = FALSE, trim_ws = TRUE)
res_MOUSE <- real_data_eval(data = db_MOUSE, type = "MOUSE", maxquant = F)

MOUSE_DiffMean <- res_MOUSE$results$DiffMean
MOUSE_EstimQual <- res_MOUSE$results$EstimQual


db_YST_B <- read.delim("~/Documents/GitHub/ProteoBayes_draft/REAL_DATA_XP/Yeast_UPS_B/peptides.txt")
res_YST_B <- real_data_eval(data = db_YST_B, type = "YST_B", maxquant = T)

YST_B_DiffMean <- res_YST_B$results$DiffMean
YST_B_EstimQual <- res_YST_B$results$EstimQual

# -- COMBINED TABLES -- #

diff_mean <- rbind(ARATH_DiffMean %>% 
                     mutate(Experiment = "Chion_2022"), 
                   YST_DiffMean %>% 
                     mutate(Experiment = "Muller_2016"), 
                   MOUSE_DiffMean %>% 
                     mutate(Experiment = "Huang_2020"),
                   YST_B_DiffMean %>% 
                     mutate(Experiment = "Bouyssie_2020")) %>% 
  relocate(Experiment)

save(diff_mean, file = "REAL_DATA_XP/Diff_Mean_Results")
write_csv(diff_mean, file = "REAL_DATA_XP/Diff_Mean_Results.csv")

estim_qual <- rbind(ARATH_EstimQual %>% 
                     mutate(Experiment = "Chion_2022"), 
                   YST_EstimQual %>% 
                     mutate(Experiment = "Muller_2016"), 
                   MOUSE_EstimQual %>% 
                     mutate(Experiment = "Huang_2020"),
                   YST_B_EstimQual %>% 
                     mutate(Experiment = "Bouyssie_2020")) %>% 
  relocate(Experiment)

save(estim_qual, file = "REAL_DATA_XP/Estim_Qual_Results")
write_csv(estim_qual, file = "REAL_DATA_XP/Estim_Qual_Results.csv")

