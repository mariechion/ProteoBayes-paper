library(tidyverse)
library(ProteoBayes)
library(gganimate)

#### Useful functions#####

simu_data = function(
    nb_peptide = 5,
    nb_sample = 5,
    list_mean_diff = c(0, 1, 5, 10),
    list_var = c(1, 1, 1, 1),
    list_cov = c(0.1, 0.1, 0.1, 0.1),
    range_peptide = c(0, 50), 
    multivariate = FALSE){
  
  nb_group = length(list_mean_diff)
  
  
  db = tibble(
    'Peptide' = rep(paste0('Peptide_', 1:nb_peptide), each= nb_group*nb_sample),
    'Group' = rep(rep(1:nb_group, each = nb_sample),nb_peptide),
    'Sample' = rep(rep( 1:nb_sample, nb_group*nb_peptide))
    ) 
  if(multivariate){
    db = db %>%
      dplyr::group_by(Peptide) %>%
      dplyr::mutate(
        'Mean' = runif(1, range_peptide[1], range_peptide[2]) +
          list_mean_diff[Group]) %>% 
      arrange(Group) %>% 
      group_by(Group) %>% 
      dplyr::mutate(
        'Output' = Mean + mvtnorm::rmvnorm(
          nb_sample,
          rep(0, nb_peptide), 
          diag(list_var[Group], nb_peptide) +
            matrix(list_cov[Group], nrow = nb_peptide, ncol = nb_peptide)
          ) %>% as.vector()
      ) %>%
      ungroup() 
  } else {
    db = db %>%
      dplyr::group_by(Peptide) %>%
      dplyr::mutate(
        'Mean' = runif(1, range_peptide[1], range_peptide[2]) +
          list_mean_diff[Group],
        'Output' = Mean + rnorm(n(), 0, list_var[Group])
      ) %>%
      ungroup() 
  }
  
  return(db)
}


eval <- function(
    nb_peptide = 5,
    nb_sample = 5,
    list_mean_diff = c(0, 1, 5, 10),
    list_var = c(1, 1, 1, 1),
    list_cov = c(0.1, 0.1, 0.1, 0.1),
    multivariate = FALSE,
    missing_ratio = 0
    ){
  
  db = simu_data(
    nb_peptide = nb_peptide,
    nb_sample = nb_sample,
    list_mean_diff = list_mean_diff,
    list_var = list_var,
    list_cov = list_cov,
    multivariate = multivariate) %>% 
    group_by(Peptide, Group) %>% 
    mutate(Average = mean(Output)) %>% 
    mutate(Missing = rbinom(1, 1, missing_ratio)) %>% 
    mutate(Output = if_else(Missing == 1, Average, Output)) %>% 
    dplyr::select(- c(Missing, Average))
  
  # Note, limma and DAPAR use data in wide format.
  db_limma <- db %>% 
    select(-Mean) %>% 
    # Merge Group and Sample columns to denote biological samples analysed
    unite(col = "Cond_Rep", c(Group,Sample), sep = "_") %>% 
    # Reshape data in wide format
    spread(key = "Cond_Rep", value = "Output")
  
  ## Create quantitative data matrix to match DAPAR requirements
  qdata <- db_limma %>%  
    column_to_rownames(var = "Peptide") 
  
  ## Create design dataframe to match DAPAR requirements
  metadata <- data.frame(Sample.name = colnames(qdata)) %>% 
    separate(Sample.name, sep = "_", into = c("Condition",NA), remove = F) %>% 
    mutate(Bio.Rep = 1:ncol(qdata))
  
  qdata <- qdata %>% 
    select(metadata$Sample.name)
  
  ## Moderated t-test - OnevsOne setting
  res_limma <- DAPAR::limmaCompleteTest(qData = as.matrix(qdata),
                                        sTab = metadata,
                                        comp.type = "OnevsOne")
  
  P_Value <- res_limma$P_Value %>%
    rownames_to_column(var = "Peptide") %>%
    pivot_longer(-Peptide,
                 names_to = "Comparison",
                 values_to = "p_value") %>%
    separate(col = Comparison,
             into = c("Group", NA, "Group2", NA),
             sep = "_") %>%
    mutate(Signif = p_value < alpha)
  
  res = db %>%
    posterior_mean() %>%
    identify_diff() %>%
    dplyr::filter(Group == 1) %>%
    left_join(db %>% select(c('Peptide', 'Group', 'Mean')) %>%
                              distinct() %>%
                              rename('Group2' = Group),
                            by = c('Peptide', 'Group2')) %>%
    mutate('MSE' = (Mean - mu2)^2,
           'CIC' = ((Mean > CI_inf2) & (Mean < CI_sup2)) * 100,
           'CI_width' = CI_sup2 - CI_inf2,
           'Diff_mean' = abs(mu2 - mu)) %>% 
    mutate(across(c(Group, Group2), .fns = as.character)) %>% 
    left_join(y = P_Value, by = c("Peptide", "Group", "Group2")) %>% 
    mutate('Multivariate' = FALSE)
  
  if(multivariate){
    res = res %>%
      bind_rows(
        db %>% multi_posterior_mean() %>% 
          identify_diff() %>%
          dplyr::filter(Group == 1) %>%
          left_join(db %>% select(c('Peptide', 'Group', 'Mean')) %>%
                      distinct() %>%
                      rename('Group2' = Group),
                    by = c('Peptide', 'Group2')) %>%
          mutate('MSE' = (Mean - mu2)^2,
                 'CIC' = ((Mean > CI_inf2) & (Mean < CI_sup2)) * 100,
                 'Diff_mean' = abs(mu2 - mu),
                 'CI_width' = CI_sup2 - CI_inf2,
                 'Multivariate' = TRUE) 
      )
  }
    
    # res <- res %>% mutate('p_value' = 0, Signif = 0) %>% 
    #   # left_join(multi_t_test(db), by = c('Peptide', 'Group', 'Group2')) %>% 
    #   return()
  
  return(res)
}

summarise_eval <- function(eval){
  eval %>%
    group_by(Group, Group2, Multivariate) %>%
    summarise(across(c(MSE, CIC, Diff_mean, CI_width, Distinct,
                       p_value, Signif),
                     .fns = list('Mean' = mean, 'Sd' = sd)),
                     .groups = 'drop') %>%
    mutate('MSE_Mean' = sqrt(MSE_Mean), 'MSE_Sd' = sqrt(MSE_Sd)) %>%
    mutate(across(MSE_Mean:Signif_Sd, ~ round(.x, 2))) %>%
    reframe(Group, Group2, Multivariate,
            'RMSE' = paste0(MSE_Mean, ' (', MSE_Sd, ')'),
            'CIC' =  paste0(CIC_Mean, ' (', CIC_Sd, ')'),
            'Diff_mean' =  paste0(Diff_mean_Mean, ' (', Diff_mean_Sd, ')'),
            'CI_width' =  paste0(CI_width_Mean, ' (', CI_width_Sd, ')'),
            'Distinct' =  paste0(Distinct_Mean*100, ' (', Distinct_Sd*100, ')'),
            'p_value' =  paste0(p_value_Mean, ' (', p_value_Sd, ')')) %>%
    return()
}

multi_t_test <- function(data){
  data %>%
    tidyr::expand_grid('Group2' = unique(Group)) %>%
    dplyr::left_join(
      data %>%
        dplyr::select(- Mean) %>%
        dplyr::rename('Group2' = Group,
                      'Output2' = Output),
      by = c("Peptide", 'Group2', 'Sample')) %>%
    dplyr::filter(Group == 1, Group2 != 1) %>%
    dplyr::group_by(Peptide, Group, Group2) %>%
    dplyr::reframe('p_value' = t.test(Output, Output2)$p.value,
                   'Signif' = (p_value < 0.05)) %>% 
    return()
}

#### Simulation study  ####

## Experiment 1: Evaluation of posteriors for different effect sizes
set.seed(42)

res1 = eval(
  nb_peptide = 1000,
  nb_sample = 5,
  list_mean_diff = c(0, 0, 1, 5, 10),
  list_var = c(1, 1, 1, 1, 1)
  )

co = summarise_eval(res1)

## Experiment 2: Evaluation of posteriors for different variances 
set.seed(1)

res2 = eval(
  nb_peptide = 1000,
  nb_sample = 5,
  list_mean_diff = c(0, 1, 1, 1, 1),
  list_var = c(1, 1, 5, 10, 20)
)

summarise_eval(res2)

## Experiment 3: Differences between univariate and multivariate versions

set.seed(42)

res3_loop = c()
for(i in 1:100)
{
  res3_10 = eval(
    nb_peptide = 10,
    nb_sample = 5,
    list_mean_diff = c(0, 1, 1, 1, 1),
    list_var = c(1, 1, 1, 10, 10),
    list_cov = c(0, 0.1, 1 ,1,  10),
    multivariate = TRUE
  ) %>% 
    mutate("Nb_peptide" = 10)
  
  res3_100 = eval(
    nb_peptide = 100,
    nb_sample = 5,
    list_mean_diff = c(0, 1, 1, 1, 1),
    list_var = c(1, 1, 1, 10, 10),
    list_cov = c(0, 0.1, 1 ,1,  10),
    multivariate = TRUE
  ) %>% 
    mutate("Nb_peptide" = 100)
  
  res3_loop = res3_loop %>%
    bind_rows(res3_10) %>%
    bind_rows(res3_100)
}

co = res3_loop %>% 
  group_by(Peptide, Group, Group2, Multivariate, Nb_peptide) %>% 
  summarise(across(c(MSE, CIC, CIC_width, Diff_mean, Distinct, p_value, Signif), 
                 .fns = mean),
            .groups= 'keep') %>% 
  summarise_eval()

## Experiment 4: Evaluation of running times
set.seed(1)
res4 = c()
for(i in c(10, 100, 1000)){
  
  data = simu_data(nb_peptide = i)
  
  t1 = Sys.time()
  dummy = posterior_mean(data)
  t2 = Sys.time()
  dummy = multi_t_test(data) 
  t3 = Sys.time()
  #dummy = multi_posterior_mean(data)
  t4 = Sys.time()
  
   res4 = bind_rows(
     res4, 
     tibble('Nb_peptide' = i, 
            'Time_t_test' = (t2 - t1) %>% round(2),
            'Time' = (t3 - t2) %>% round(2),
            'Time_multi' = (t4 - t3) %>% round(2))
     ) 
}

## Experiment 5: Evaluation of the uncertainty bias coming from imputation 

res5 = eval(
  nb_peptide = 1000,
  nb_sample = 100,
  list_mean_diff = c(0, 1),
  list_var = c(1, 1),
  multivariate = FALSE,
  missing_ratio = 0.8
)

sum_res5 = summarise_eval(res5)


#### GIF ProteoBayes visualisation ####
set.seed(1)
nb_sample = 5
data = simu_data(nb_peptide = 1, 
                 nb_sample = nb_sample,
                 list_mean_diff = c(-2, 5),
                 list_var = c(1, 3),
                 range_peptide = c(10,10))

## Prepare the data and posterior distributions for animation
data_anim <- tibble::tibble()
post_anim <- tibble::tibble()
for (j in 1:n_distinct(data$Sample)) {
  
  data_anim = data %>%
    dplyr::filter(Sample %in% 1:j) %>% 
    dplyr::mutate('Index' = j) %>% 
    mutate('Group' = as.factor(Group)) %>%
    dplyr::bind_rows(data_anim)
  
  ## Extract the sample of the 'j' first data points
  post_anim <- data %>%
    dplyr::filter(Sample %in% 1:j) %>% 
    posterior_mean(mu_0 = 10, lambda_0 = 1, alpha_0 = 1, beta_0 = 1) %>% 
    dplyr::mutate('Index' = j) %>%
    dplyr::bind_rows(post_anim)
}

## Define the prior distribution
prior = tibble('Peptide' = 'Peptide_1', 'Group' = 1:2 , 'mu' = 10, 
               'lambda' = 1, 'alpha' = 1, 'beta' = 1, 'Index' = 0)

## Sample empirical distributions from the posteriors
samples = post_anim %>%
  bind_rows(prior) %>% 
  group_by(Group, Index) %>% 
  sample_distrib(nb_sample = 100000) %>%
  mutate('Group' = as.factor(Group))

## Sample empirical distributions from the priors
samples_prior = prior %>%
  uncount(6) %>%
  mutate(Index = rep(0:5, 2)) %>% 
  group_by(Group, Index) %>% 
  sample_distrib(nb_sample = 100000)

## Create the animation
gg = ggplot(samples) +
  geom_density(data=samples_prior, aes(x = Sample), fill = 'black', alpha=0.3) +
  geom_density(aes(x = Sample, fill = Group), alpha = 0.6)  + 
  geom_point(data = data_anim, aes(x = Output, y = 0, col = Group), size = 3) + 
  theme_classic() + xlim(2, 21) + xlab('Output') + ylab('Density') + 
  gganimate::transition_states(Index)


gg_anim = animate(gg, height = 1800, width = 3200, res = 300)


anim_save('anim_ProteoBayes_wide2.gif', gg_anim) 
