# -- LIBRAIRIES -- # 
library(tidyverse)
library(ProteoBayes)
library(cp4p)
library(DAPAR)

# -- FUNCTIONS -- # 
source("REAL_DATA_XP/Functions.R")


# -- EXPERIMENTS -- #
set.seed(17)

db_ARATH <- read.delim("REAL_DATA_XP/Arabido_UPS/peptides.txt")
res_ARATH <- real_data_eval(data = db_ARATH, type = "ARATH",
                            prop_NA = 0.2,
                            multi = F,
                            mu_0 = NULL, 
                            lambda_0 = 1e-10, beta_0 = 0.1, alpha_0 = 0.1,
                            alpha = 0.05, FDR = NULL)

ARATH_DiffMean <- res_ARATH$results$DiffMean
ARATH_EstimQual <- res_ARATH$results$EstimQual

db_YST <- read.delim("REAL_DATA_XP/Yeast_UPS/peptides.txt")
res_YST <- real_data_eval(data = db_YST, type = "YST",
                          prop_NA = 0.2,
                          multi = F,
                          mu_0 = NULL, 
                          lambda_0 = 1e-10, beta_0 = 0.1, alpha_0 = 0.1,
                          alpha = 0.05, FDR = NULL)

YST_DiffMean <- res_YST$results$DiffMean
YST_EstimQual <- res_YST$results$EstimQual

db_MOUSE <- read_delim("REAL_DATA_XP/DIAHuang2020/Spike-in-biol-var-OT-SN-Report.txt", 
                       delim = "\t", escape_double = FALSE, trim_ws = TRUE)

res_MOUSE <- real_data_eval(data = db_MOUSE, type = "MOUSE",
                            maxquant = F,
                            prop_NA = 0.2,
                            multi = F,
                            mu_0 = NULL, 
                            lambda_0 = 1e-10, beta_0 = 0.1, alpha_0 = 0.1,
                            alpha = 0.05, FDR = NULL)


MOUSE_DiffMean <- res_MOUSE$results$DiffMean
MOUSE_EstimQual <- res_MOUSE$results$EstimQual