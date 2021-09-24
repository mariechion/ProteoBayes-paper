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
  
  ## Loop over the groups 
  floop_k = function(k){
    ## Extract the adequate group
    data_k = data %>% filter(Group == k)
    ## Collect all the different groups
    list_Draw = data_k$Draw %>% unique()
    
    sapply(list_Draw, floop_d, k = k) %>% 
      group_by(Group) %>% 
      summarise_all(Mean_emp_dist = mean())
    return()
  }
  
  ## Loop over the draws
  floop_d = function(d, k){
    ## Extract the adequate draws
    data_k_d = data %>%
      filter(Group == k,  Draw = d) %>% 
      Select(- Group)
    
    ## Compute the sum_1^N{y_n}
    sum_yn = data_k_d %>% 
      group_by(Draw) %>% 
      summarise(Output = sum(Output)) 
    
    ## Compute the updated posterior hyper-parameters
    mu_N = NA
    lambda_N = NA
    Sigma_N = NA
    nu_N = NA
    ## Draw from the adequate T-distribution
    emp_dist = rt(..., size = 10^6)
    
    emp_dist %>% 
      as_tibble() %>% 
      mutate(Draw = k)%>% 
      return()
  }
  
  ## Collect all the different groups
  list_group = data$Group %>% unique()
  
  sapply(list_group, floop_k) %>% 
    return()
}