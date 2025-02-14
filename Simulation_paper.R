library(tidyverse)
library(ProteoBayes)
library(gganimate)
library(DAPAR)

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
    t_test = FALSE,
    limma = FALSE,
    missing_ratio = 0,
    imputation = F,
    mu_0 = NULL,
    lambda_0,
    beta_0,
    alpha_0){
  
  db = simu_data(
    nb_peptide = nb_peptide,
    nb_sample = nb_sample,
    list_mean_diff = list_mean_diff,
    list_var = list_var,
    list_cov = list_cov,
    multivariate = multivariate) 
  
    if(imputation){
      db = db %>% 
        group_by(Peptide, Group) %>% 
        mutate(Average = mean(Output)) %>% 
        mutate(Missing = rbinom(nb_sample, 1, missing_ratio)) %>% 
        mutate(Output = if_else(Missing == 1, Average, Output)) %>% 
        dplyr::select(- c(Missing, Average)) %>% 
        ungroup()
    } else
    {
      db = db %>% 
        group_by(Peptide, Sample) %>% 
        mutate(Missing = rbinom(1, 1, missing_ratio)) %>% 
        dplyr::filter(Missing == 0) %>% 
        dplyr::select(- c(Missing)) %>% 
        ungroup()
    }
  
  res = db %>%
    posterior_mean(mu_0 = mu_0, lambda_0 = lambda_0,
                   alpha_0 = alpha_0, beta_0 = beta_0) %>%
    identify_diff() %>%
    dplyr::filter(Group == 1) %>%
    left_join(db %>% select(c('Peptide', 'Group', 'Mean')) %>%
                              distinct() %>%
                              rename('Group2' = Group),
                            by = c('Peptide', 'Group2')) %>%
    mutate('MSE' = (Mean - mu2)^2,
           'CIC' = ((Mean > CI_inf2) & (Mean < CI_sup2)) * 100,
           'CI_width' = CI_sup2 - CI_inf2,
           'Diff_mean' = mu2 - mu) %>% 
    mutate(across(c(Group, Group2), .fns = as.character)) %>% 
    mutate('Multivariate' = FALSE)
  
  if(multivariate){
    res = res %>%
      bind_rows(
        db %>% multi_posterior_mean(mu_0 = mu_0, lambda_0 = lambda_0,
                                    nu_0 = alpha_0, Sigma_0 = beta_0) %>% 
          identify_diff() %>%
          dplyr::filter(Group == 1) %>%
          left_join(db %>% select(c('Peptide', 'Group', 'Mean')) %>%
                      distinct() %>%
                      rename('Group2' = Group),
                    by = c('Peptide', 'Group2')) %>%
          mutate('MSE' = (Mean - mu2)^2,
                 'CIC' = ((Mean > CI_inf2) & (Mean < CI_sup2)) * 100,
                 'CI_width' = CI_sup2 - CI_inf2,
                 'Diff_mean' = mu2 - mu,
                 'Multivariate' = TRUE, 
                 'Group' = as.character(Group),
                 'Group2' = as.character(Group2)) 
      )
  }
    
  if(t_test){
    res <- res %>% 
      left_join(multi_t_test(db),
                by = c('Peptide', 'Group', 'Group2'))
  }
  
  if(limma){
    res <- res %>% 
      left_join(multi_limma(db),
                by = c('Peptide', 'Group', 'Group2'))
  }

  return(res)
}

summarise_eval <- function(eval){
  eval %>%
    group_by(Group, Group2, Multivariate) %>%
    summarise(across(where(is.double),
                     .fns = list('Mean' = mean, 'Sd' = sd)),
                     .groups = 'drop') %>%
    mutate('MSE_Mean' = sqrt(MSE_Mean), 'MSE_Sd' = sqrt(MSE_Sd)) %>%
    mutate(across(where(is.double), ~ round(.x, 2))) %>%
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
    dplyr::reframe('p_value' = tryCatch(t.test(Output,Output2)$p.value,
                                      error = function(e){return(NA)}),
                   'Signif' = (p_value < 0.05))  %>% 
    mutate(across(c(Group, Group2), .fns = as.character)) %>% 
    return()
}

multi_limma <- function(data){
  # Note, limma and DAPAR use data in wide format.
  db_limma <- data %>% 
    dplyr::select(-Mean) %>% 
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
                 values_to = "LM_p_value") %>%
    separate(col = Comparison,
             into = c("Group", NA, "Group2", NA),
             sep = "_")
  
  return(P_Value)
}


ci_coverage <- function(
    posterior,
    CI_level_seq = seq(0.05, 1, 0.05)
    ){
  
    posterior %>%
      dplyr::mutate('sigma' = sqrt(.data$beta / (.data$lambda * .data$alpha)),
                    'df' =  2 * .data$alpha) %>%
      uncount(length(CI_level_seq)) %>% 
      mutate(CI_level = rep(CI_level_seq, nrow(posterior))) %>% 
      dplyr::reframe(
        .data$Peptide, .data$Group, .data$mu, .data$CI_level,
        'CI_inf' = extraDistr::qlst(
          p = CI_level/2,
          df =.data$df,
          mu = .data$mu,
          sigma = .data$sigma,
          lower.tail = T),
        'CI_sup' =  extraDistr::qlst(
          p = CI_level/2,
          df =.data$df,
          mu = .data$mu,
          sigma = .data$sigma,
          lower.tail = F)
      ) %>% 
    return()
}

#### Simulation study  ####
## Experiment 1: Evaluation of posteriors for different effect sizes
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

## Experiment 2: Evaluation of posteriors for different variances 
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

## Experiment 3: Differences between univariate and multivariate versions

set.seed(42)

res3_loop = c()
for(i in 1:1000){
  res3_10 = eval(
    nb_peptide = 10,
    nb_sample = 5,
    list_mean_diff = c(0, 1, 1, 1, 1),
    list_var = c(1, 1, 1, 10, 10),
    list_cov = c(0, 0.1, 1 , 1, 10),
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
    list_cov = c(0, 0.1, 1 ,1,  10),
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

## Experiment 4: Evaluation of running times
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
  
## Experiment 5: Evaluation of the uncertainty bias coming from imputation 


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

####
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

#### Additional visualisations ####
## Pre-experiment: Role of alpha_0 and beta_0 on CI_Width and CIC
## and impact of the number of samples.

db_sim <- map(.x = c(3,5,10,100),
              .f = ~ simu_data(nb_peptide = 1000,
                               nb_sample = .x,
                               list_mean_diff = c(0,1,5,10,1,1),
                               list_var = c(1,1,1,1,5,10)) %>% 
                mutate(nb_samples = .x)) %>% 
  bind_rows()

exp_sim <- tibble(alpha0 = rep(rep(c(0.01, 0.1, 1), rep(3,3)),4),
                  beta0 = rep(c(0.01,0.1,1),12)) %>% 
  rowwise() %>% 
  reframe(db_sim %>%
            group_by(nb_samples) %>%
            ProteoBayes::posterior_mean(lambda_0 = 1e-10, 
                                        alpha_0 = alpha0,
                                        beta_0 = beta0) %>% 
            ProteoBayes::identify_diff() %>% 
            filter(Group == 1) %>% 
            left_join(db_sim %>% 
                        filter(nb_samples == nb_samples) %>% 
                        select(c('nb_samples','Peptide', 'Group', 'Mean')) %>%
                        distinct() %>%
                        rename('Group2' = Group),
                      by = c('Peptide', 'Group2')) %>%
            mutate('MSE' = (Mean - mu2)^2,
                   'CIC' = ((Mean > CI_inf2) & (Mean < CI_sup2)) * 100,
                   'CI_width' = CI_sup2 - CI_inf2,
                   'Diff_mean' = mu2 - mu) %>% 
            mutate(alpha_0 = alpha0, beta_0 = beta0))

### Instead of wanting a straightforward function, make it annoyingly manually.

db_sim <- map(.x = c(3,5,10,100),
              .f = ~ simu_data(nb_peptide = 1000,
                               nb_sample = .x,
                               list_mean_diff = c(0,1,5,10,1,1),
                               list_var = c(1,1,1,1,5,10)) %>% 
                mutate(nb_samples = .x))

exp_sim <- lapply(db_sim, function(db){
  tibble(alpha0 = rep(c(0.01, 0.1, 1), rep(3,3)),
         beta0 = rep(c(0.01,0.1,1),3)) %>% 
    rowwise() %>% 
    reframe(db %>%
              ProteoBayes::posterior_mean(lambda_0 = 1e-10, 
                                          alpha_0 = alpha0,
                                          beta_0 = beta0) %>% 
              ProteoBayes::identify_diff() %>% 
              filter(Group == 1) %>% 
              mutate(alpha_0 = alpha0, beta_0 = beta0)) %>% 
    left_join(db %>% 
                select(c('Peptide', 'Group', 'Mean')) %>%
                distinct() %>%
                rename('Group2' = Group),
              by = c('Peptide', 'Group2')) %>%
    mutate('MSE' = (Mean - mu2)^2,
           'CIC' = ((Mean > CI_inf2) & (Mean < CI_sup2)) * 100,
           'CI_width' = CI_sup2 - CI_inf2,
           'Diff_mean' = mu2 - mu)
  
}) 

exp_sim %>% 
  bind_rows(.id = "nb_samples") %>% 
  mutate(nb_samples = case_match(nb_samples, 
                                 '1' ~ '3', '2' ~ '5', '3' ~ '10', '4' ~ '100'),
         Group2 = case_match(Group2,
                             2 ~ "N(1,1)",
                             3 ~ "N(5,1)",
                             4 ~ "N(10,1)",
                             5 ~ "N(1,5)",
                             6 ~ "N(1,10)")) %>% 
  select(nb_samples, alpha_0, beta_0, Group2, 
         Peptide, MSE, CIC, CI_width, Diff_mean) %>% 
  mutate(nb_samples = factor(nb_samples, levels = c('3','5','10','100'))) %>% 
  group_by(nb_samples, alpha_0, beta_0, Group2) %>% 
  summarise(CIC_mean = mean(CIC)) %>% 
  filter(Group2 == "N(5,1)") %>% 
  ggplot() +
  geom_raster(aes(x = alpha_0, y = beta_0, fill = CIC_mean),
              interpolate = T) +
  scale_x_continuous(transform = "log10") +
  scale_y_continuous(transform = "log10") +
  # geom_boxplot(aes(x = nb_samples, 
  #                  y = CI_width, 
  #                  fill = Group2)) +
  # facet_grid(beta_0 ~ alpha_0, labeller = "label_both") +
  theme_minimal()


## Heatmaps

db_heatmap <- simu_data(nb_peptide = 1000,
                        nb_sample = 3,
                        list_mean_diff = c(1,5),
                        list_var = c(1,1))

grid_hp = seq(-2,2, 0.25)

exp_heatmap <- expand_grid(alpha_0 = grid_hp,
                           beta_0 = grid_hp) %>% 
  mutate(across(alpha_0:beta_0, ~10^(.))) %>% 
  rowwise() %>% 
  reframe(db_heatmap %>%
            ProteoBayes::posterior_mean(lambda_0 = 1e-10, 
                                        alpha_0 = alpha_0,
                                        beta_0 = beta_0) %>%
            ci_coverage() %>%
            mutate(alpha_0 = alpha_0, beta_0 = beta_0)) %>% 
  left_join(db_heatmap %>% 
              select(c('Peptide', 'Group', 'Mean')) %>%
              distinct(),
            by = c('Peptide', 'Group')) %>%
  mutate('CIC' = ((Mean > CI_inf) & (Mean < CI_sup)) * 100) %>% 
  group_by(alpha_0, beta_0, CI_level) %>% 
  summarise(Coverage = mean(CIC)) %>% 
  mutate(CI_level = 100*(1-CI_level)) %>% 
  mutate(Error = abs(Coverage - CI_level)) %>% 
  group_by(alpha_0, beta_0) %>% 
  summarise(CIC_error = mean(Error))

ggplot(exp_heatmap) +
  geom_raster(aes(x = alpha_0, y = beta_0, 
                  fill = CIC_error),
              interpolate = T) + 
  scale_x_continuous(transform = "log10") +
  scale_y_continuous(transform = "log10") +
  scale_fill_gradientn(colours = c(
    "white",
    "#FDE0DD",
    "#FCC5C0",
    "#FA9FB5",
    "#F768A1",
    "#DD3497",
    "#AE017E",
    "#7A0177"), trans = 'reverse') +
  theme_classic()
         
## Calibration of uncertainty graphical check

ci_valid = db_heatmap %>%
  ProteoBayes::posterior_mean(lambda_0 = 1e-10, 
                              alpha_0 = 0.01,
                              beta_0 = 0.3) %>%
  ci_coverage(seq(0.01, 1, 0.01)) %>% 
  left_join(db_heatmap %>% 
              select(c('Peptide', 'Group', 'Mean')) %>%
              distinct(),
            by = c('Peptide', 'Group')) %>%
  mutate('CIC' = ((Mean > CI_inf) & (Mean < CI_sup))) %>% 
  group_by(CI_level) %>% 
  summarise(Coverage = mean(CIC)) %>% 
  mutate(CI_level = (1-CI_level))

ggplot(ci_valid) + 
  geom_line(aes(x = CI_level, y = Coverage), col = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  theme_classic()


# exp_heatmap %>% 
#   mutate(Group2 = case_match(Group2,
#                              2 ~ "N(5,1)")) %>% 
#   group_by(alpha_0, beta_0, Group2) %>% 
#   summarise(CIC_mean = mean(CIC)) %>%
#   mutate(Confidence = case_when(CIC_mean < 94 ~ "< 94",
#                                 CIC_mean > 96 ~ "> 96", 
#                                 .default = "94-96")) %>% 
#   filter(Group2 == "N(5,1)") %>% 
#   ggplot() +
#   geom_raster(aes(x = alpha_0, y = beta_0, 
#                   fill = CIC_mean),
#               interpolate = T) +
#   scale_x_continuous(transform = "log10") +
#   scale_y_continuous(transform = "log10") +
#   scale_fill_gradientn(colors = c(rep("red",70),"green","blue"), 
#                        limits = c(20,100)) +
#   theme_minimal()
#### Graph of results ####

full_res = read_csv("REAL_DATA_XP/Exp2_summary.csv") %>%
  separate(RMSE, c('RMSE', NA), " ") %>%
  separate(CIC, c('CIC', NA), " ") %>% 
  mutate(True = - True, RMSE = as.numeric(RMSE), CIC = as.numeric(CIC))

gg1 = ggplot(full_res, aes(x = True, y = CIC, col = Experiment, shape = Experiment)) +  
  geom_point(position=position_dodge(0.05)) + 
  geom_hline(yintercept = 95, linetype = 'dashed') +
  ylim(0, 100) + xlab('Mean difference') + ylab('95% Credible Interval Coverage') +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 95, 100), limits = c(0,100)) +
  theme_classic()

gg2 = ggplot(full_res, aes(x = True, y = RMSE, col = Experiment, shape = Experiment))+  
  geom_point(position=position_dodge(0.05)) + 
  xlab('Mean difference') + ylab('RMSE') +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5)) +
  theme_classic()


ggsave('FIGURES/CIC_all_exp.png', gg1, width = 4, height = 3, dpi = 600)

