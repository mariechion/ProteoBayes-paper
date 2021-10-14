load("Arabido_raw_data_lg")
load("Arabido_imp_data_lg")

source("bayes_proteo.R")

pept_ups = "AALEELVK" 
pept_arath = "EVQELAQEAAER" 
db_pept_ups = db.raw %>% filter(ID == pept_ups) %>% drop_na()
db_pept_arath = db.raw %>% filter(ID == pept_arath) %>% drop_na()

db_pept_ups_imp = data.lg %>% 
  dplyr::rename(Draw = Imp.Draw, Output = Intensity) %>% 
  filter(ID == pept_ups) %>% 
  group_by(Sample) %>%
  summarise(ID = max(ID), Group = max(Group), Output = mean(Output))

db_pept_arath_imp = data.lg %>% 
  dplyr::rename(Draw = Imp.Draw, Output = Intensity) %>% 
  filter(ID == pept_arath) %>% 
  group_by(Sample) %>%
  summarise(ID = max(ID), Group = max(Group), Output = mean(Output))

## Graph comparison with and without correlations
res_graph3_ups_uni = post_mean_diff_uni(
  data = db_pept_ups,
  mu_0 = 20, 
  lambda_0 = 1,
  beta_0 = 1,
  alpha_0 = 2
)
res_graph3_ups_uni_imp = post_mean_diff_uni(
  data = db_pept_ups_imp,
  mu_0 = 20, 
  lambda_0 = 1,
  beta_0 = 1,
  alpha_0 = 1
)

res_graph3_arath_uni = post_mean_diff_uni(
  data = db_pept_arath,
  mu_0 = 20, 
  lambda_0 = 1,
  beta_0 = 1,
  alpha_0 = 2
)
res_graph3_arath_uni_imp = post_mean_diff_uni(
  data = db_pept_arath_imp,
  mu_0 = 20, 
  lambda_0 = 1,
  beta_0 = 1,
  alpha_0 = 1
)

gg3_1 = plot_dif(res_graph3_ups_uni, c('Point1', 'Point7'), peptide = pept_ups) +
  xlim(c(-1,10)) + ylim(c(0,0.4))
  xlab(TeX('$\\mu_1 - \\mu_7$ - ups'))
gg3_2 = plot_dif(res_graph3_ups_uni_imp, c('Point1', 'Point7'), peptide = pept_ups) +
  xlim(c(-1,10)) + ylim(c(0,0.4)) 
  xlab(TeX('$\\mu_1 - \\mu_7$- ups - imp'))
gg3_3 = plot_dif(res_graph3_arath_uni, c('Point1', 'Point7'), peptide = pept_arath) +
  xlim(c(-2,3)) + ylim(c(0,0.9))
  xlab(TeX('$\\mu_1 - \\mu_7$- arath'))
gg3_4 = plot_dif(res_graph3_arath_uni_imp, c('Point1', 'Point7'), peptide = pept_arath) +
  xlim(c(-2,3)) + ylim(c(0,0.9)) 
  xlab(TeX('$\\mu_1 - \\mu_7$'))

#png('FIGURES/ch5_graph3.png', res = 600, width = 4800, height = 3600, units = "px")
cowplot::plot_grid(gg3_1, gg3_3, gg3_2, gg3_4, nrow = 2, ncol = 2)
#dev.off()
