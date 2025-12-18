library(tidyverse)
library(ProteoBayes)
library(gganimate)
library(mvtnorm)
library(patchwork)
library(cowplot)

source("R/utils.R")
 
load("Data/Arabido_UPS/Arabido_raw_data_lg")
load("Data/Arabido_UPS/Arabido_imp_data_lg")

pept_ups = "AALEELVK"
pept_arath = "EVQELAQEAAER"

db_pept_ups = db.raw %>% 
  filter(ID == pept_ups) %>% 
  drop_na() %>% 
  select(-Protein) %>%
  rename(Peptide = ID)

db_pept_arath = db.raw %>%
  filter(ID == pept_arath) %>% 
  drop_na() %>% 
  select(-Protein) %>%
  rename(Peptide = ID)

db_prot_ups = data.lg %>% 
  dplyr::rename(Draw = Imp.Draw, Output = Intensity) %>% 
  filter(Protein == "P12081ups|SYHC_HUMAN_UPS") %>% 
  drop_na() %>% 
  select(-Protein) %>%
  rename(Peptide = ID)

db_prot_arath = data.lg  %>% 
  dplyr::rename(Draw = Imp.Draw, Output = Intensity) %>%
  filter(Protein == "sp|F4I893|ILA_ARATH") %>%
   drop_na() %>% 
  select(-Protein) %>%
  rename(Peptide = ID)

#### Illustration of Univariate Inference ####
db_ups = db_prot_ups %>%
  select(-Protein) %>%
  rename(Peptide = ID)

res_ups = posterior_mean(
  data = db_ups,
  mu_0 = 25, 
  lambda_0 = 1,
  beta_0 = 1,
  alpha_0 = 2
)

samp_ups = sample_distrib(res_ups, n = 10000)

gg_ups <- plot_distrib(samp_ups, group1 = 'Point1', group2 = 'Point7') 

db_arath = db_prot_arath %>%
  select(-Protein) %>%  
  rename(Peptide = ID)

res_arath = posterior_mean(
  data = db_arath,
  mu_0 = 25,
  lambda_0 = 1,
  beta_0 = 1,
  alpha_0 = 2
)

samp_arath = sample_distrib(res_arath, n = 10000)

gg_arath <- plot_distrib(samp_arath, group1 = 'Point1', group2 = 'Point7')

#########################################################################

#### The mirage of imputation ####

## Impute missing data with the mean 
db_pept_arath_imp = data.lg %>% 
  dplyr::rename(Draw = Imp.Draw, Output = Intensity) %>% 
  filter(ID == pept_arath) %>% 
  group_by(Sample) %>%
  summarise(ID = max(ID), Group = max(Group), Output = mean(Output)) %>% 
  rename(Peptide = ID)

res_raw = posterior_mean(
  data = db_pept_arath,
  mu_0 = 20, 
  lambda_0 = 1,
  beta_0 = 1,
  alpha_0 = 2
)

res_imp = posterior_mean(
  data = db_pept_arath_imp,
  mu_0 = 20, 
  lambda_0 = 1,
  beta_0 = 1,
  alpha_0 = 2
)

samp_raw = sample_distrib(res_raw, n = 10000)
samp_imp = sample_distrib(res_imp, n = 10000)

gg_raw = plot_distrib(samp_raw, 'Point1', group2 = 'Point4') +
  xlim(c(-5,5)) + ylim(c(0,0.9))
gg_imp = plot_distrib(samp_imp, 'Point1', group2 = 'Point4')  +
  xlim(c(-5,5)) + ylim(c(0,0.9))

#########################################################################

#### Illustration of the effect size ####

res_effect = posterior_mean(
  data = db_pept_ups,
  mu_0 = db_pept_ups$Output %>% mean(), 
  lambda_0 = 1,
  beta_0 = 1,
  alpha_0 = 2
)

samp_effect = sample_distrib(res_effect, n = 10000)

gg_effect_1_2 = plot_distrib(samp_effect, group1 = 'Point1', group2 = 'Point2') +
  xlim(c(-7.5,2.5)) + ylim(c(0,0.7))

gg_effect_1_4 = plot_distrib(samp_effect, group1 = 'Point1', group2 = 'Point4') +
  xlim(c(-7.5,2.5)) + ylim(c(0,0.7))

gg_effect_1_7 = plot_distrib(samp_effect, group1 = 'Point1', group2 = 'Point7') +
  xlim(c(-7.5,2.5)) + ylim(c(0,0.7))

#########################################################################



#########################################################################

#### Difficulty of Multivariate Inference ####

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

## Graph of marginal differences for each peptide ##
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

gg = cowplot::plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, ncol = 3, label_x = '') 

# ggsave('Figures/multivariate_marginals.png', gg, dpi = 600, width = 6400, height = 3600, units = "px")

## Multivariate differential inference ##

diff = multi_identify_diff(post_arath)

gg = plot_multi_diff(diff)

# ggsave("Figures/multi_diff_diff.png", gg, dpi = 600, width = 6, height = 6)

#########################################################################

#### Summary of real data experiments ####

full_res = read_csv("Results/Real_Data/Bouyssie2020.csv") %>%
  mutate(Experiment = 'Bouyssie2020') %>%
  bind_rows(read_csv("Results/Real_Data/Chion2022.csv") %>% 
    mutate(Experiment = 'Chion2022')) %>%
  bind_rows(read_csv("Results/Real_Data/Huang2020.csv") %>% 
    mutate(Experiment = 'Huang2020')) %>%
  bind_rows(read_csv("Results/Real_Data/Muller2016.csv") %>% 
    mutate(Experiment = 'Muller2016')) %>%
  separate(RMSE, c('RMSE', NA), " ") %>%
  separate(CIC, c('CIC', NA), " ") %>% 
  mutate(Mean_diff = - True_diff_mean, RMSE = as.numeric(RMSE), CIC = as.numeric(CIC))

gg_CIC = ggplot(full_res,
             aes(x = Mean_diff, y = CIC, col = Experiment, shape = Experiment)) +  
  geom_point(position=position_dodge(0.05)) + 
  geom_hline(yintercept = 95, linetype = 'dashed') +
  xlab('Mean difference') + ylab('95% Credible Interval Coverage') +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 95, 100), limits = c(0,100)) +
  theme_classic()

gg_RMSE = ggplot(full_res, 
             aes(x = Mean_diff, y = RMSE, col = Experiment, shape = Experiment))+  
  geom_point(position=position_dodge(0.05)) + 
  xlab('Mean difference') + ylab('RMSE') +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5)) +
  theme_classic()

#########################################################################


#### Role of alpha_0 and beta_0 on CI_Width and CIC ####

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

gg_heatmap <- ggplot(exp_heatmap) +
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

#########################################################################

##### Calibration of uncertainty graphical check ####

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

#########################################################################

#### Draw 2D overlapping Gaussians ####

mu1    <- c(0, 0)
Sigma1 <- matrix(c(1.0, 0.6, 0.6, 1.2), 2)
mu2    <- c(1.8, 1.0)
Sigma2 <- matrix(c(1.0,-0.4,-0.4, 0.8), 2)
probs  <- seq(0.25, 0.95, length.out = 5)     # 5 niveaux d'ellipses
levs   <- qchisq(probs, df = 2)

# --- 2D grid + densities ---
rngx <- range(mu1[1] + c(-4,4)*sqrt(Sigma1[1,1]), mu2[1] + c(-4,4)*sqrt(Sigma2[1,1]))
rngy <- range(mu1[2] + c(-4,4)*sqrt(Sigma1[2,2]), mu2[2] + c(-4,4)*sqrt(Sigma2[2,2]))
x <- seq(rngx[1], rngx[2], length.out = 220)
y <- seq(rngy[1], rngy[2], length.out = 220)
G <- expand.grid(x = x, y = y)
X <- cbind(as.numeric(G$x), as.numeric(G$y))
G$d21 <- mahalanobis(X, mu1, Sigma1)
G$d22 <- mahalanobis(X, mu2, Sigma2)
G$f1  <- dmvnorm(X, mean = mu1, sigma = Sigma1)
G$f2  <- dmvnorm(X, mean = mu2, sigma = Sigma2)
G$fo  <- pmin(G$f1, G$f2)

# --- Central panel ---
p_main <- ggplot(G, aes(x, y)) +
  geom_contour(aes(z = d21), breaks = levs, linewidth = 0.4, col = '#366bd5ff') +
  geom_contour(aes(z = d22), breaks = levs, linetype = "dashed", col = 'black') +
  geom_point(data = data.frame(x = mu1[1], y = mu1[2]), size = 2, col = '#366bd5ff') +
  geom_point(data = data.frame(x = mu2[1], y = mu2[2]), size = 2) +
  coord_equal(xlim = rngx, ylim = rngy) +
  labs(x = "x", y = "y") +
  theme_classic()

# --- Marginal X (top) ---
dx <- data.frame(x = x,
                 g1 = dnorm(x, mu1[1], sqrt(Sigma1[1,1])),
                 g2 = dnorm(x, mu2[1], sqrt(Sigma2[1,1])))
p_top <- ggplot(dx, aes(x)) +
  geom_ribbon(aes(ymin = 0, ymax = pmin(g1, g2)), alpha = 0.7, fill = "#F8B9C5") +
  geom_line(aes(y = g1),  col = '#366bd5ff') +
  geom_line(aes(y = g2), linetype = "dashed") +
  coord_cartesian(xlim = rngx) +
  theme_void() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = margin(0,0,2,0))

# --- Marginal Y (left) ---
dy <- data.frame(y = y,
                 g1 = dnorm(y, mu1[2], sqrt(Sigma1[2,2])),
                 g2 = dnorm(y, mu2[2], sqrt(Sigma2[2,2])))
p_left <- ggplot(dy, aes(y = y)) +
  geom_ribbon(aes(xmin = 0, xmax = pmin(g1, g2), y = y), alpha = 0.7,
fill = "#F8B9C5", orientation = "y") +
  geom_path(aes(x = g1, y = y),  col = '#366bd5ff') +
  geom_path(aes(x = g2, y = y), linetype = "dashed") +
  scale_x_reverse() +
  coord_cartesian(ylim = rngy, xlim = c(0, 0.8)) +
  theme_void() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = margin(0,0,0,2))

#  Precise placement
gg_2D_gaussians = ggdraw() +
  draw_plot(p_left, 0.00, 0.10, 0.28, 0.75) +   # left
  draw_plot(p_main, 0.15, 0.00, 0.85, 0.85) +   # center
  draw_plot(p_top,  0.35, 0.85, 0.50, 0.15)     # top

#################################################################################

#### GIF ProteoBayes illustration ####
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


# anim_save('anim_ProteoBayes_wide2.gif', gg_anim) 
#########################################################################

