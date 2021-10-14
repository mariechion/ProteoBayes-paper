load("Arabido_raw_data_lg")
load("Arabido_imp_data_lg")

source("bayes_proteo.R")

pept_ups = "AALEELVK"
db_pept_ups = db.raw %>% filter(ID == pept_ups) %>% drop_na()
mean = db_pept_ups$Output %>% mean()

## Graph comparison with and without correlations
res_ups_uni = post_mean_diff_uni(
  data = db_pept_ups,
  mu_0 = mean,
  lambda_0 = 1,
  beta_0 = 1,
  alpha_0 = 2
)

gg4_1 = plot_dif(res_ups_uni, c('Point1', 'Point2'), peptide = pept_ups) +
  xlim(c(-7.5,2.5)) + ylim(c(0,0.7))
  gg4_2 = plot_dif(res_ups_uni, c('Point1', 'Point4'), peptide = pept_ups) +
  xlim(c(-7.5,2.5)) + ylim(c(0,0.7))
gg4_3 = plot_dif(res_ups_uni, c('Point1', 'Point5'), peptide = pept_ups) +
  xlim(c(-7.5,2.5)) + ylim(c(0,0.7)) 
gg4_4 = plot_dif(res_ups_uni, c('Point1', 'Point7'), peptide = pept_ups) +
  xlim(c(-7.5,2.5)) + ylim(c(0,0.7))

png('FIGURES/ch5_graph4.png', res = 600, width = 4800, height = 3600, units = "px")
cowplot::plot_grid(gg4_1, gg4_2, gg4_4, nrow = 3, ncol = 1)
dev.off()
