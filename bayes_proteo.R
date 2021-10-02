# Load required packages
library(tidyverse)
library(mvtnorm)
library(ggplot2)
library(Matrix)

post_mean_diff = function(
  data,
  mu_0, 
  lambda_0,
  Sigma_0,
  nu_0
){
  ## data: A tibble or data frame containing imputed data setsfor all groups.
  ##  Required columns: ID, Output, Group, Draw
  ## mu_0: A vector, corresponding to the prior mean 
  ## lamba_0: A number, corresponding to the prior covariance scaling parameter
  ## Sigma_0: A matrix, corresponding to the prior covariance parameter
  ## nu_0: A number, corresponding to the prior degrees of freedom
  ## return: A vector providing the empirical posterior distribution of the 
  ##   mean's difference between the desired groups.
  t1 = Sys.time()
  ## Loop over the groups
  floop_k = function(k){
    
    ## Extract the adequate group
    data_k = data %>% filter(Group == k)
    ## Collect all the different groups
    list_draw = data_k$Draw %>% unique()
    n_draw = data_k$Draw %>% n_distinct()
    
    list_mat = lapply(list_draw, floop_d, k = k) 

    ((1/n_draw) * Reduce('+', list_mat)) %>% 
      as_tibble() %>% 
      mutate(Group = k, .before = 1) %>% 
    return()
  }
  
  ## Loop over the draws
  floop_d = function(d, k){
    t2 = Sys.time()
    paste0('Group n°',k , ' - Draw n° ', d, ' - ', t2 - t1) %>% print()
    ## Extract the adequate draws
    data_k_d = data %>%
      filter(Group == k,  Draw == d)
    
    N_k = data_k_d$Sample %>% n_distinct()
    list_sample = data_k_d$Sample %>% unique()
    
    ## Compute the mean 1/N sum_1^N{y_n}
    mean_yn_k = data_k_d %>% 
      group_by(ID) %>% 
      summarise(Output = mean(Output)) %>% 
      pull(Output)
    
    ##Compute the mean 1/N sum_1^N{(y_n - \bar{y_n}^)2}
    cov_yn = 0
    for(n in list_sample)
    {
      yn_k = data_k_d %>%
        filter(Sample == n) %>%
        pull(Output)
      
      centred_y = (yn_k - mean_yn_k) 
      
      cov_yn = cov_yn + tcrossprod(centred_y)
    }
    
    centred_mean = mean_yn_k - mu_0
    ## Compute the updated posterior hyper-parameters
    mu_N = (lambda_0 * mu_0 + N_k * mean_yn_k) / (lambda_0 + N_k) 
    lambda_N = lambda_0 + N_k
    Sigma_N = Sigma_0 + cov_yn + 
      (lambda_0 * N_k) / lambda_N * tcrossprod(centred_mean)
    nu_N = nu_0 + N_k
    ## Draw from the adequate T-distribution
    rmvt(n = 10, sigma = Sigma_N / (nu_N * lambda_N),
                    df = nu_N, delta = mu_N) %>% 
      return()
  }
  
  ## Collect all the different groups
  list_group = data$Group %>% unique()
  
  lapply(list_group, floop_k) %>% 
    bind_rows() %>% 
    return()
}

post_mean_diff_uni = function(
  data,
  mu_0, 
  lambda_0,
  beta_0,
  alpha_0
){
  ## data: A tibble or data frame containing imputed data setsfor all groups.
  ##  Required columns: ID, Output, Group, Draw
  ## mu_0: A vector, corresponding to the prior mean 
  ## lamba_0: A number, corresponding to the prior covariance scaling parameter
  ## beta_0: A matrix, corresponding to the prior covariance parameter
  ## alpha_0: A number, corresponding to the prior degrees of freedom
  ## return: A vector providing the empirical posterior distribution of the 
  ##   mean's difference between the desired groups.
  t1 = Sys.time()
  ## Loop over the groups
  floop_k = function(k){
    t2 = Sys.time()
    paste0('Group n°',k , ' Time : ', t2 - t1) %>% print()
    
    ## Extract the adequate group
    data_k = data %>% filter(Group == k)
    ## Extract the number of samples
    N_k = data_k$Sample %>% n_distinct()
    ## Extract the list of samples
    list_sample = data_k$Sample %>% unique()
    
    ## Compute the mean 1/N sum_1^N{y_n}
    mean_yn_k = data_k %>% 
      group_by(ID) %>% 
      summarise(Output = mean(Output)) %>% 
      pull(Output)
    
    ##Compute the mean 1/N sum_1^N{(y_n - \bar{y_n}^)2}
    cov_yn = 0
    for(n in list_sample)
    {
      yn_k = data_k %>%
        filter(Sample == n) %>%
        pull(Output)
      
      cov_yn = cov_yn + (yn_k - mean_yn_k)^2
    }
    
    centred_mean = 
    ## Compute the updated posterior hyper-parameters
    mu_N = (lambda_0 * mu_0 + N_k * mean_yn_k) / (lambda_0 + N_k) 
    lambda_N = lambda_0 + N_k
    beta_N = beta_0 + cov_yn + 
      (lambda_0 * N_k) / lambda_N * (mean_yn_k - mu_0)^2
    alpha_N = alpha_0 + N_k / 2
    ## Draw from the adequate T-distribution
    tibble(
      'Group' = k,
      'Mean' = mu_N + sqrt(beta_N) * rt(n = 10000, df = alpha_N)
      )%>%
      return()
  }
  
  ## Collect all the different groups
  list_group = data$Group %>% unique()
  
  lapply(list_group, floop_k) %>% 
    bind_rows() %>% 
    return()
}


plot_dif = function(emp_dist, groups, peptide){
  if(length(groups) == 2){
    db1 = emp_dist %>% filter(Group == groups[[1]]) %>% select(- Group)
    db2 = emp_dist %>% filter(Group == groups[[2]]) %>% select(- Group)
    db = tibble(Dif = (db1 - db2)[, peptide])
    bar = 0
  } else if(length(groups) == 1){
    db1 = emp_dist %>% filter(Group == groups[[1]]) %>% select(- Group)
    db = tibble(Dif = db1 %>% pull(peptide))
    bar = mean(db$Dif)
  }
  ggplot(db) +
    geom_density(aes(x = Dif), fill = "#00B2EE") %>% return() +
    geom_vline(xintercept = bar, color = 'red') +
    theme_classic()
}
