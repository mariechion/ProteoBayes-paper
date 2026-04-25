library(tidyverse)
library(ProteoBayes)
library(mvtnorm)
# library(proDA)
# library(msqrob2)
# library(SummarizedExperiment)
# library(DAPAR)

source("R/utils.R")


#### Simulation study  ####

#### Experiment 1: Evaluation of posteriors for different effect sizes ####
set.seed(42)

res1 = eval(
  nb_peptide = 200,
  nb_sample = 1000,
  list_mean_diff = c(0, 1, 5, 10, 1, 1, 1),
  list_var = c(1, 1, 1, 1, 5, 10, 20),
  lambda_0 = 1e-10,
  t_test = FALSE,
  limma = FALSE,
  msqrob = TRUE,
  proDA = FALSE,
  alpha_0 = 0.01, 
  beta_0 = 0.3
)

sum_res1 = summarise_eval(res1)

#### Experiment 2: Evaluation of posteriors for different variances ####
set.seed(1)

res2 = eval(
  nb_peptide = 1000,
  nb_sample = 5,
  list_mean_diff = c(0, 1, 1, 1, 1),
  list_var = c(1, 1, 5, 10, 20),
  lambda_0 = 1e-10,
  alpha_0 = 0.01,
  beta_0 = 0.3
)

summarise_eval(res2)

#### Experiment 3: Differences between univariate and multivariate versions ####

set.seed(42)

res3_loop = c()
for(i in 1:100){
  res3_10 = eval(
    nb_peptide = 10,
    nb_sample = 100,
    list_mean_diff = c(0, 1, 1),
    list_var = c(1, 1, 10),
    lambda_0 = 1e-10,
    alpha_0 = 10,
    beta_0 = 10,
    multivariate = TRUE
  ) %>% 
    mutate("Nb_peptide" = 10)
  
  res3_100 = eval(
    nb_peptide = 100,
    nb_sample = 5,
    list_mean_diff = c(0, 1, 1, 1, 1),
    list_var = c(1, 1, 1, 10, 10),
    lambda_0 = 1e-10,
    alpha_0 = 100,
    beta_0 = 100,
    multivariate = TRUE
  ) %>% 
    mutate("Nb_peptide" = 100)
  
  res3_loop = res3_loop %>%
    bind_rows(res3_10) %>%
    bind_rows(res3_100)
}

sum_res3 = res3_loop %>% 
  group_by(Group, Group2, Multivariate, Nb_peptide) %>% 
  summarise(across(c(Diff_mean, MSE, CIC, CI_width, Diff_mean), 
                 .fns = mean),
            .groups= 'keep')

#write_csv(sum_res3, 'Results_simu/comparison_uni_multi.csv')
sum_res3 = read_csv('Results_simu/comparison_uni_multi.csv')

#### Experiment 4: Evaluation of running times #####
set.seed(1)

floop = function(n){
  tib = c()
  for(i in c(10, 100, 1000, 10000)){
    
    data = simu_data(nb_peptide = i)
    
    t1 = Sys.time()
    dummy = posterior_mean(data)
    dummy2 = identify_diff(dummy)
    t2 = Sys.time()
    dummy = multi_t_test(data)
    t3 = Sys.time()
    dummy = multi_posterior_mean(data)
    dummy2 = multi_identify_diff(dummy)
    t4 = Sys.time()
    dummy = multi_limma(data)
    t5 = Sys.time()
    dummy = multi_proDA(data)
    t6 = Sys.time()
    dummy = multi_msqrob(data)
    t7 = Sys.time()
    
    tib = bind_rows(
      tib, 
      tibble('Nb_peptide' = i, 
             'Time' = (t2 - t1),
             'Time_t_test' = (t3 - t2),
             'Time_multi' = (t4 - t3),
             'Time_lima' = (t5 - t4),
             'Time_proDA' = (t6 - t5),
             'Time_msqrob' = (t7 - t6))
    )
  }
  return(tib)
}
res4 = lapply(1:10, floop) %>%
  bind_rows() %>%
  group_by(Nb_peptide) %>% 
  summarise(across(everything(), list(mean = mean, sd = sd)))

#write_csv(res4, 'running_time.csv')
  
#### Experiment 5: Evaluation of the uncertainty bias coming from imputation ####

floop2 = function(i)
{
  no_imput = eval(
    nb_peptide = 100,
    nb_sample = 10,
    list_mean_diff = c(0, 1),
    list_var = c(1, 1),
    multivariate = FALSE,
    imputation = F,
    missing_ratio = i, 
    t_test = F,
    limma = T,
    proDA = T,
    lambda_0 = 1e-10,
    alpha_0 = 0.01, 
    beta_0 = 0.3
  ) %>% 
    mutate(Missing_ratio = i) %>% 
    mutate(Imputation = F)
  
  imput = eval(
    nb_peptide = 100,
    nb_sample = 10,
    list_mean_diff = c(0, 1),
    list_var = c(1, 1),
    multivariate = FALSE,
    imputation = T,
    missing_ratio = i, 
    t_test = F,
    limma = T,
    proDA = T,
    lambda_0 = 1e-10,
    alpha_0 = 0.01, 
    beta_0 = 0.3
  ) %>% 
    mutate(Missing_ratio = i) %>% 
    mutate(Imputation = T)

  no_imput %>% 
    bind_rows(imput) %>% 
    return()
}
res5 = c(0, 0.2, 0.5, 0.8) %>% 
  lapply(floop2) %>% 
  bind_rows()

sum_res5 = res5 %>%
  drop_na() %>%  
  mutate(Proba_differential = 1- Proba_differential) %>% 
  group_by(Missing_ratio, Imputation) %>%
  summarise(across(where(is.double), 
                   .fns = list('Mean' = mean, 'Sd' = sd)),
            .groups = 'drop') %>%
  mutate('MSE_Mean' = sqrt(MSE_Mean), 'MSE_Sd' = sqrt(MSE_Sd)) %>%
  mutate(across(where(is.double), ~ round(.x, 2))) %>% 
  mutate(across(ends_with("_Mean"), 
                ~ sprintf("%.2f (%.2f)", ., get(sub("_Mean$", "_Sd", cur_column()))),
                .names = "{sub('_Mean$', '', .col)}")) %>%
  select(-ends_with("_Mean"), -ends_with("_Sd"))

#### Experiment 6: Evaluation of the effect size and uncertainty quantification in the multivariate case ####

set.seed(42)

mean_0 = c(0, 0, 0)
mean_1 = c(1, 1, 1)
mean_10 = c(10, 10, 10)

cov_1 = diag(c(1,1,1))
cov = matrix(c(1, 0.7, 0.2, 0.7, 1, 0.5, 0.2, 0.5, 1) , nrow = 3, ncol = 3)
cov_high = matrix(c(10, 0.7, 0.2, 0.7, 10, 0.5, 0.2, 0.5, 10) , nrow = 3, ncol = 3)

list_compare = list(
  'Mean_1_Cov_diag' = list(mean_1,  cov_1),
  'Mean_1_Cov_correlated' = list(mean_1, cov),
  'Mean_10_Cov_correlated' = list(mean_10, cov),
  'Mean_1_Cov_high_correlated' = list(mean_1, cov_high)
)

res = tibble(mean_diff = c(), MSE = c(), CI_coverage = c())

for(i in names(list_compare)){
  res = eval_multi(
    mean_ref = mean_0,
    cov_ref = cov,
    mean_compare = list_compare[[i]][[1]],
    cov_compare = list_compare[[i]][[2]],
    lambda_0 = 1e-10,
    nu_0 = 3,
    nb_rep = 100,
    nb_sample = 100
  ) %>%
  mutate(
    'Distrib_compare' = i) %>% 
  bind_rows(res)
}

sum_res = res %>% group_by(Distrib_compare) %>%
  summarise(across(where(is.double), 
                   .fns = list('Mean' = mean, 'Sd' = sd)),
            .groups = 'drop') %>%
  mutate('MSE_Mean' = sqrt(MSE_Mean), 'MSE_Sd' = sqrt(MSE_Sd)) %>%
  mutate(across(where(is.double), ~ round(.x, 4)))

# write_csv(sum_res, 'Results_simu/summary_uni_on_multi_simu_evaluation_100samples.csv')

####

#### Experiment 7: Sensitivity analysis HPs univariate ####

lambda = sensitivity_uni(
    nb_peptide = 1000, 
    nb_sample = 5, 
    param = "lambda", 
    grid_param = c(1e-3, 1e-2, 1e-1, 1, 10, 100, 1000),
    alpha_0 = 0.01, 
    beta_0 = 0.3,
    lambda_0 = 1e-10, 
    list_mean_diff = c(0, 1, 10, 10),
    list_var = c(1, 10, 1, 100),
    CI_level = 0.05)

sum_lambda = lambda %>% 
  dplyr::select(NLL, CIC, Param) %>% 
  group_by(Param) %>% 
  summarise(across(where(is.double), 
                   .fns = list('Mean' = mean, 'Sd' = sd)),
            .groups = 'drop') %>%
  mutate(across(where(is.double), ~ round(.x, 2)))

#write_csv(sum_lambda, 'Results/Simulations/sensitivity_lambda.csv')

alpha = sensitivity_uni(
  nb_peptide = 1000, 
  nb_sample = 5, 
  param = "alpha", 
  grid_param = c(1e-3, 1e-2, 1e-1, 1, 10, 100, 1000),
  alpha_0 = 0.01, 
  beta_0 = 0.3,
  lambda_0 = 1e-10, 
  list_mean_diff = c(0, 1, 10, 10),
  list_var = c(1, 10, 1, 100),
  CI_level = 0.05)

sum_alpha = alpha %>% 
  dplyr::select(NLL, CIC, Param) %>% 
  group_by(Param) %>% 
  summarise(across(where(is.double), 
                   .fns = list('Mean' = mean, 'Sd' = sd)),
            .groups = 'drop') %>%
  mutate(across(where(is.double), ~ round(.x, 2)))

#write_csv(sum_alpha, 'Results/Simulations/sensitivity_alpha.csv')

beta = sensitivity_uni(
  nb_peptide = 1000, 
  nb_sample = 5, 
  param = "beta", 
  grid_param = c(1e-3, 1e-2, 1e-1, 1, 10, 100, 1000),
  alpha_0 = 0.01, 
  beta_0 = 0.3,
  lambda_0 = 1e-10, 
  list_mean_diff = c(0, 1, 10, 10),
  list_var = c(1, 10, 1, 100),
  CI_level = 0.05)

sum_beta = beta %>% 
  dplyr::select(NLL, CIC, Param) %>% 
  group_by(Param) %>% 
  summarise(across(where(is.double),  
                   .fns = list('Mean' = mean, 'Sd' = sd)),
            .groups = 'drop') %>%
  mutate(across(where(is.double), ~ round(.x, 2)))

#write_csv(sum_beta, 'Results/Simulations/sensitivity_beta.csv')

  

#### Experiment 7: Sensitivity analysis HPs multivariate ####

Sigma_5 = sensitivity_multi(
    nb_rep = 300, 
    nb_sample = 5, 
    param = "Sigma", 
    grid_param = c(1e-3, 5e-2, 1e-2, 5e-1, 1e-1, 0.5, 1),
    lambda_0 = 1e-10) %>% 
  mutate(Nb_sample = 5)

Sigma_10 = sensitivity_multi(
  nb_rep = 300, 
  nb_sample = 10, 
  param = "Sigma", 
  grid_param = c(1e-3, 5e-2, 1e-2, 5e-1, 1e-1, 0.5, 1),
  lambda_0 = 1e-10) %>% 
  mutate(Nb_sample = 10)

Sigma_50 = sensitivity_multi(
  nb_rep = 300, 
  nb_sample = 50, 
  param = "Sigma", 
  grid_param = c(1e-3, 5e-2, 1e-2, 5e-1, 1e-1, 0.5, 1),
  lambda_0 = 1e-10) %>% 
  mutate(Nb_sample = 50)

Sigma_100 = sensitivity_multi(
  nb_rep = 300, 
  nb_sample = 100, 
  param = "Sigma", 
  grid_param = c(1e-3, 5e-2, 1e-2, 5e-1, 1e-1, 0.5, 1),
  lambda_0 = 1e-10) %>% 
  mutate(Nb_sample = 100)

Sigma = Sigma_5 %>% 
  bind_rows(Sigma_10) %>% 
  bind_rows(Sigma_50) %>% 
  bind_rows(Sigma_100)

sum_Sigma = Sigma %>% 
  mutate(NLL = NLL / Nb_sample) %>%
  dplyr::select(NLL, CIC, Param, Nb_sample) %>% 
  group_by(Param, Nb_sample) %>% 
  summarise(across(where(is.numeric), 
                   .fns = list('Mean' = mean, 'Sd' = sd)),
            .groups = 'drop') %>%
  mutate(across(where(is.double), ~ round(.x, 2)))

#write_csv(sum_Sigma, 'Results/Simulations/sensitivity_Sigma.csv')

nu_diag = sensitivity_multi(
  nb_rep = 500, 
  nb_sample = 5, 
  param = "nu", 
  grid_param = c(5, 10, 20, 50, 100),
  lambda_0 = 1e-10, 
  cov_diag = TRUE) %>% 
  mutate(Cov_diag = TRUE)

nu_cov = sensitivity_multi(
  nb_rep = 500, 
  nb_sample = 5, 
  param = "nu", 
  grid_param = c(5, 10, 20, 50, 100),
  lambda_0 = 1e-10, 
  cov_diag = FALSE) %>% 
  mutate(Cov_diag = FALSE)


nu = nu_diag %>% 
  bind_rows(nu_cov)

sum_nu = nu %>% 
  mutate(NLL = NLL) %>%
  dplyr::select(NLL, CIC, Param, Cov_diag) %>% 
  group_by(Param, Cov_diag) %>% 
  summarise(across(where(is.numeric), 
                   .fns = list('Mean' = mean, 'Sd' = sd)),
            .groups = 'drop') %>%
  mutate(across(where(is.double), ~ round(.x, 2)))

#write_csv(sum_nu, 'Results/Simulations/sensitivity_nu.csv')

