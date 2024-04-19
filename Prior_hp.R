## Illustrate the impact of prior hyperparameters values
library(tidyverse)
library(ProteoBayes)
library(extraDistr)

# Sampling from a Gaussian-inverse-Gamma distribution
rnorm_inv_gamma <- function(n, mu_0, lambda_0, alpha_0, beta_0){
  if (length(n) > 1) 
    n <- length(n)
  sigma2 <- rinvgamma(n = n, alpha = alpha_0, beta = beta_0)
  mu <- rnorm(n, mean = mu_0, sd = sqrt(sigma2/lambda_0))
  data.frame(mu = mu, sigma2 = sigma2)
}

# Compute the PDF -> useless?
pdf_norminvgamma <- function(mu, sigma2, mu_0, lambda_0, alpha_0, beta_0){
  return(sqrt(lambda_0/(2*pi))*(beta_0^alpha_0/gamma(alpha_0))*
           (1/sigma2)^(alpha_0+3/2)*
           exp(-(2*beta_0 + lambda_0*(mu-mu_0)^2)/(2*sigma2)))
}


# Create tibble of hyperparameters

db_hp <- tibble(experiment = as.character(1:4), mu_0 = 0, 
                lambda_0 = c(1,2,5,10),
                beta_0 = 1,
                alpha_0 = 1)

obs_hp <- db_hp %>% 
  group_by(experiment) %>% 
  reframe(rnorm_inv_gamma(n = 10000, mu_0 = mu_0, lambda_0 = lambda_0,
                          beta_0 = beta_0, alpha_0 = alpha_0)) %>% 
  left_join(y = db_hp, by = "experiment") 

ggplot(obs_hp, aes(y = mu, x = lambda_0, col = experiment)) +
  geom_violin()

