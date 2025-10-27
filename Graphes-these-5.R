 library(tidyverse)
library(ProteoBayes)
 
load("Arabido_raw_data_lg")
load("Arabido_imp_data_lg")

# source("bayes_proteo.R")

pept_ups = "AALEELVK"
pept_arath = "EVQELAQEAAER"
# db_pept_ups = db.raw %>% filter(ID == pept_ups) %>% drop_na()
# db_pept_arath = db.raw %>% filter(ID == pept_arath) %>% drop_na()

db_prot_ups = data.lg %>% 
  dplyr::rename(Draw = Imp.Draw, Output = Intensity) %>% 
  filter(Protein == "P12081ups|SYHC_HUMAN_UPS") %>% 
  drop_na()

db_prot_arath = data.lg  %>% 
  dplyr::rename(Draw = Imp.Draw, Output = Intensity) %>%
  filter(Protein == "sp|F4I893|ILA_ARATH") %>%
   drop_na()

## Graph comparison with and without correlations
dim_prot_ups = db_prot_ups$ID %>% n_distinct()
mean_ups = db_prot_ups$Output %>% mean()
mean_arath = db_prot_arath$Output %>% mean()

#### New version of multivariate ProteoBayes ####

db = db_prot_arath %>%
  select(-Protein) %>%
  rename(Peptide = ID) %>% 
  mutate(Peptide = substr(Peptide, 1, 4)) %>% 
  mutate(Group = substr(Group, 6, 7)) %>% 
  filter(Group %in% c("1", "7"))

post_arath = multi_posterior_mean(
  data = db, 
  mu_0 = NULL,
  lambda_0 = 1e-10,
  nu_0 = 10)

diff = multi_identify_diff(post_arath)

## Graph of multivariate differences
gg = plot_multi_diff(diff)

# ggsave("FIGURES/multi_diff_diff.png", gg,
#        dpi = 600, width = 6, height = 6)

## Graph of marginal differences for each peptide
sample = sample_distrib(post_arath, n = 10000)
list_pept = db$Peptide %>% unique()
p1 = plot_distrib(sample, group1 = "1", group2 = "7", peptide = list_pept[1]) +
  xlim(c(-2,2)) + ylim(c(0,3.5)) +
  xlab(list_pept[1])
p2 = plot_distrib(sample, group1 = "1", group2 = "7", peptide = list_pept[2]) +
  xlim(c(-2,2)) + ylim(c(0,3.5)) + ylab("") +
  xlab(list_pept[2])
p3 = plot_distrib(sample, group1 = "1", group2 = "7", peptide = list_pept[3]) +
  xlim(c(-2,2)) + ylim(c(0,3.5)) + ylab("") +
  xlab(list_pept[3])
p4 = plot_distrib(sample, group1 = "1", group2 = "7", peptide = list_pept[4]) +
  xlim(c(-2,2)) + ylim(c(0,3.5)) +
  xlab(list_pept[4])
p5 = plot_distrib(sample, group1 = "1", group2 = "7", peptide = list_pept[5]) +
  xlim(c(-2,2)) + ylim(c(0,3.5)) + ylab("") +
  xlab(list_pept[5])
p6 = plot_distrib(sample, group1 = "1", group2 = "7", peptide = list_pept[6]) +
  xlim(c(-2,2)) + ylim(c(0,3.5)) + ylab("") +
  xlab(list_pept[6])
p7 = plot_distrib(sample, group1 = "1", group2 = "7", peptide = list_pept[7]) +
  xlim(c(-2,2)) + ylim(c(0,3.5)) +
  xlab(list_pept[7])
p8 = plot_distrib(sample, group1 = "1", group2 = "7", peptide = list_pept[8]) +
  xlim(c(-2,2)) + ylim(c(0,3.5)) + ylab("") +
  xlab(list_pept[8])
p9 = plot_distrib(sample, group1 = "1", group2 = "7", peptide = list_pept[9]) +
  xlim(c(-2,2)) + ylim(c(0,3.5)) + ylab("") +
  xlab(list_pept[9])

gg = cowplot::plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, ncol = 3, label_x = 'Bla') 

# ggsave('FIGURES/multivariate_marginals.png', gg, dpi = 600, width = 6400, height = 3600, units = "px")
