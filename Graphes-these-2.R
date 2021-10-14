load("Arabido_raw_data_lg")
load("Arabido_imp_data_lg")

source("bayes_proteo.R")

pept_ups = "AALEELVK" 
pept_arath = "ALADPNTDVR" 
db_pept_ups = db.raw %>% filter(ID == pept_ups) %>% drop_na()
db_pept_arath = db.raw %>% filter(ID == pept_arath) %>% drop_na()

db_prot_ups = data.lg %>% 
  dplyr::rename(Draw = Imp.Draw, Output = Intensity) %>% 
  filter(Protein == "P12081ups|SYHC_HUMAN_UPS") %>% drop_na()
db_prot_arath = data.lg  %>% 
  dplyr::rename(Draw = Imp.Draw, Output = Intensity) %>%
  filter(Protein == "sp|F4I893|ILA_ARATH") %>% drop_na()

## Graph comparison with and without correlations
dim_prot_ups = db_prot_ups$ID %>% n_distinct()
mean_ups = db_pept_ups$Output %>% mean()
mean_arath = db_pept_arath$Output %>% mean()

res_graph2_ups_uni = post_mean_diff_uni(
  data = db_pept_ups,
  mu_0 = mean_ups, 
  lambda_0 = 1,
  beta_0 = 0.5,
  alpha_0 = 5
)

res_graph2_ups_multi = post_mean_diff(
  data = db_prot_ups,
  mu_0 = rep(mean_ups, dim_prot_ups), 
  lambda_0 = 1,
  Sigma_0 = diag(1, ncol = dim_prot_ups, nrow = dim_prot_ups),
  nu_0 = 10
) %>% filter(ID == pept_ups)

res_graph2_arath_uni = post_mean_diff_uni(
  data = db_pept_arath,
  mu_0 = mean_arath, 
  lambda_0 = 1,
  beta_0 = 0.5,
  alpha_0 = 5
)

dim_prot_arath = db_prot_arath$ID %>% n_distinct()

res_graph2_arath_multi = post_mean_diff(
  data = db_prot_arath,
  mu_0 = rep(mean_arath, dim_prot_arath), 
  lambda_0 = 1,
  Sigma_0 = diag(1, ncol = dim_prot_arath, nrow = dim_prot_arath),
  nu_0 = 10
) %>% filter(ID == pept_arath)

gg2_1 = plot_dif(res_graph2_ups_uni, c('Point1', 'Point7'), peptide = pept_ups) +
  xlim(c(-1,7))
gg2_2 = plot_dif(res_graph2_ups_multi, c('Point1', 'Point7'), peptide = pept_ups) +
  xlim(c(-1,7))
# gg2_3 = plot_dif(res_graph2_arath_uni, c('Point4', 'Point7'), peptide = pept_arath) +
#   xlim(c(-5,5)) 
# gg2_4 = plot_dif(res_graph2_arath_multi, c('Point4', 'Point7'), peptide = pept_arath) +
#    xlim(c(-5,5))

png('FIGURES/ch5_graph2.png', res = 600, width = 4800, height = 3600, units = "px")
cowplot::plot_grid(gg2_1, gg2_2, nrow = 2, ncol = 1)
dev.off()

