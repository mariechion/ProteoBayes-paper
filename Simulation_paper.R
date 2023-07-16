library(tidyverse)
library(ProteoBayes)

#### Useful functions#####

simu_data = function(
    nb_peptide = 5,
    nb_sample = 5,
    list_mean_diff = c(0, 1, 5, 10),
    list_var = c(1, 1, 1, 1),
    range_peptide = c(0, 50)){
  nb_group = length(list_mean_diff)
  tibble(
    'Peptide' = rep(paste0('Peptide_', 1:nb_peptide), each= nb_group*nb_sample),
    'Group' = rep(rep(1:nb_group, each = nb_sample),nb_peptide),
    'Sample' = rep(rep( 1:nb_sample, nb_group*nb_peptide))
    ) %>%
    dplyr::group_by(.data$Peptide) %>%
    dplyr::mutate(
      'Mean' = runif(1, range_peptide[1], range_peptide[2]) +
        list_mean_diff[.data$Group],
      'Output' = .data$Mean + rnorm(n(), 0, list_var[.data$Group])
      ) %>%
    ungroup() %>%
  return()
}


eval <- function(
    nb_peptide = 5,
    nb_sample = 5,
    list_mean_diff = c(0, 1, 5, 10),
    list_var = c(1, 1, 1, 1)){
  
  db = simu_data(
    nb_peptide = nb_peptide,
    nb_sample = nb_peptide,
    list_mean_diff = list_mean_diff,
    list_var = list_var)
  
  t_test =  db %>%
    tidyr::expand_grid('Group2' = unique(.data$Group)) %>%
    dplyr::left_join(
      db %>%
        dplyr::select(- Mean) %>%
        dplyr::rename('Group2' = .data$Group,
                      'Output2' = .data$Output),
        by = c("Peptide", 'Group2', 'Sample')) %>%
    dplyr::filter(Group == 1, Group2 != 1) %>%
    dplyr::group_by(Peptide, Group, Group2) %>%
    dplyr::reframe('p_value' = t.test(Output, Output2)$p.value,
                    'Signif' = (p_value < 0.05))
  db %>%
    posterior_mean() %>%
    identify_diff() %>%
    dplyr::filter(Group == 1) %>%
    left_join(db %>% select(c('Peptide', 'Group', 'Mean')) %>%
                              distinct() %>%
                              rename('Group2' = Group),
                            by = c('Peptide', 'Group2')) %>%
    mutate('MSE' = (Mean - mu2)^2,
                         'CIC' = ((Mean > CI_inf2) & (Mean < CI_sup2)) * 100,
                         'Diff_mean' = abs(mu2 - mu)) %>%
    left_join(t_test, by = c('Peptide', 'Group', 'Group2')) %>%
    return()
}

summarise_eval <- function(eval){
  eval %>%
    group_by(Group, Group2) %>%
    summarise(across(c(MSE, CIC, Diff_mean, Distinct, p_value, Signif),
                     .fns = list('Mean' = mean, 'Sd' = sd))) %>%
    mutate('MSE_Mean' = sqrt(MSE_Mean), 'MSE_Sd' = sqrt(MSE_Sd)) %>%
    mutate(across(MSE_Mean:Signif_Sd, ~ round(.x, 2))) %>%
    reframe(Group2,
            'RMSE' = paste0(MSE_Mean, ' (', MSE_Sd, ')'),
            'CIC' =  paste0(CIC_Mean, ' (', CIC_Sd, ')'),
            'Diff_mean' =  paste0(Diff_mean_Mean, ' (', Diff_mean_Sd, ')'),
            'Distinct' =  paste0(Distinct_Mean*100, ' (', Distinct_Sd*100, ')'),
            'p_value' =  paste0(p_value_Mean, ' (', p_value_Sd, ')'),
            'Signif' =  paste0(Signif_Mean*100, ' (', Signif_Sd*100, ')')) %>%
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

summarise_eval(res1)

## Experiment 2: Evaluation of posteriors for different variances 

set.seed(42)

res2 = eval(
  nb_peptide = 1000,
  nb_sample = 100,
  list_mean_diff = c(0, 1, 1, 1, 1),
  list_var = c(1, 1, 5, 10, 20)
  )

summarise_eval(res2)