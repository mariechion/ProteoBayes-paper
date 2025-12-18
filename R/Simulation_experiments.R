library(tidyverse)
library(ProteoBayes)
library(mvtnorm)
# library(DAPAR)

#### Simulation study  ####

#### Experiment 1: Evaluation of posteriors for different effect sizes ####
set.seed(42)

res1 = eval(
  nb_peptide = 1000,
  nb_sample = 5,
  list_mean_diff = c(0, 0, 1, 5, 10),
  list_var = c(1, 1, 1, 1, 1),
  lambda_0 = 1e-10,
  alpha_0 = 0.01, 
  beta_0 = 0.3
  )

summarise_eval(res1)

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
    t2 = Sys.time()
    dummy = multi_t_test(data) 
    t3 = Sys.time()
    dummy = multi_posterior_mean(data)
    t4 = Sys.time()
    dummy = multi_limma(data) 
    t5 = Sys.time()
    
    tib = bind_rows(
      tib, 
      tibble('Nb_peptide' = i, 
             'Time' = (t2 - t1),
             'Time_t_test' = (t3 - t2),
             'Time_multi' = (t4 - t3), 
             'Time_lima' = (t5 - t4))
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
    nb_peptide = 10000,
    nb_sample = 10,
    list_mean_diff = c(0, 1),
    list_var = c(1, 1),
    multivariate = FALSE,
    imputation = F,
    missing_ratio = i, 
    t_test = T,
    limma = T,
    lambda_0 = 1e-10,
    alpha_0 = 0.01, 
    beta_0 = 0.3
  ) %>% 
    mutate(Missing_ratio = i) %>% 
    mutate(Imputation = F)
  
  imput = eval(
    nb_peptide = 10000,
    nb_sample = 10,
    list_mean_diff = c(0, 1),
    list_var = c(1, 1),
    multivariate = FALSE,
    imputation = T,
    missing_ratio = i, 
    t_test = T,
    limma = T,
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
  dplyr::select(MSE:Imputation) %>% 
  group_by(Missing_ratio, Imputation) %>%
  summarise(across(where(is.double), 
                   .fns = list('Mean' = mean, 'Sd' = sd)),
            .groups = 'drop') %>%
  mutate('MSE_Mean' = sqrt(MSE_Mean), 'MSE_Sd' = sqrt(MSE_Sd)) %>%
  mutate(across(where(is.double), ~ round(.x, 2)))

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

  
