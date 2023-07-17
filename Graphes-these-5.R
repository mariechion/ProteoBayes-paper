load("Arabido_raw_data_lg")
load("Arabido_imp_data_lg")

source("bayes_proteo.R")

pept_ups = "AALEELVK"
pept_arath = "EVQELAQEAAER"
# db_pept_ups = db.raw %>% filter(ID == pept_ups) %>% drop_na()
# db_pept_arath = db.raw %>% filter(ID == pept_arath) %>% drop_na()

db_prot_ups = data.lg %>% 
  dplyr::rename(Draw = Imp.Draw, Output = Intensity) %>% 
  filter(Protein == "P12081ups|SYHC_HUMAN_UPS") %>% drop_na()
db_prot_arath = data.lg  %>% 
  dplyr::rename(Draw = Imp.Draw, Output = Intensity) %>%
  filter(Protein == "sp|F4I893|ILA_ARATH") %>% drop_na()

## Graph comparison with and without correlations
dim_prot_ups = db_prot_ups$ID %>% n_distinct()
mean_ups = db_prot_ups$Output %>% mean()
mean_arath = db_prot_arath$Output %>% mean()

res_graph5_ups_multi = post_mean_diff(
  data = db_prot_ups,
  mu_0 = rep(mean_ups, dim_prot_ups), 
  lambda_0 = 1,
  Sigma_0 = diag(1, ncol = dim_prot_ups, nrow = dim_prot_ups),
  nu_0 = 10
)

dim_prot_arath = db_prot_arath$ID %>% n_distinct()

res_graph5_arath_multi = post_mean_diff(
  data = db_prot_arath,
  mu_0 = rep(mean_arath, dim_prot_arath), 
  lambda_0 = 1,
  Sigma_0 = diag(1, ncol = dim_prot_arath, nrow = dim_prot_arath),
  nu_0 = 10
)

list_pept = db_prot_ups$ID %>% unique()
gg5_1 = plot_dif(res_graph5_ups_multi, c('Point1', 'Point7'), peptide = list_pept[1]) +
  xlim(c(-7,6)) + ylim(c(0,1.8)) +
  xlab(list_pept[1])
gg5_2 = plot_dif(res_graph5_ups_multi, c('Point1', 'Point7'), peptide = list_pept[2]) +
  xlim(c(-7,6)) + ylim(c(0,1.8)) +
  xlab(list_pept[2])
gg5_3 = plot_dif(res_graph5_ups_multi, c('Point1', 'Point7'), peptide = list_pept[3]) +
  xlim(c(-7,6)) + ylim(c(0,1.8)) +
  xlab(list_pept[3])
gg5_4 = plot_dif(res_graph5_ups_multi, c('Point1', 'Point7'), peptide = list_pept[4]) +
  xlim(c(-7,6)) + ylim(c(0,1.8)) +
  xlab(list_pept[4])
gg5_5 = plot_dif(res_graph5_ups_multi, c('Point1', 'Point7'), peptide = list_pept[5]) +
  xlim(c(-7,6)) + ylim(c(0,1.8)) +
  xlab(list_pept[5])
gg5_6 = plot_dif(res_graph5_ups_multi, c('Point1', 'Point7'), peptide = list_pept[6]) +
  xlim(c(-7,6)) + ylim(c(0,1.8)) +
  xlab(list_pept[6])
gg5_7 = plot_dif(res_graph5_ups_multi, c('Point1', 'Point7'), peptide = list_pept[7]) +
  xlim(c(-7,6)) + ylim(c(0,1.8)) +
  xlab(list_pept[7])
gg5_8 = plot_dif(res_graph5_ups_multi, c('Point1', 'Point7'), peptide = list_pept[8]) +
  xlim(c(-7,6)) + ylim(c(0,1.8)) +
  xlab(list_pept[8])
gg5_9 = plot_dif(res_graph5_ups_multi, c('Point1', 'Point7'), peptide = list_pept[9]) +
  xlim(c(-7,6)) + ylim(c(0,2)) +
  xlab(list_pept[9])

png('FIGURES/ch5_graph5.png', res = 600, width = 4800, height = 2700, units = "px")
cowplot::plot_grid(gg5_1, gg5_3, gg5_5, gg5_7, gg5_9,
                   gg5_2, gg5_4, gg5_6, gg5_8, nrow = 3, ncol = 3)
dev.off()
