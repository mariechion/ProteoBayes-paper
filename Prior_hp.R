## Illustrate the impact of prior hyperparameters values
library(tidyverse)
library(ProteoBayes)
library(extraDistr)
library(gganimate)

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

#### Arthur

nb_sample = 1000
lambda = 1
alpha = 1
beta = 1

sample_prior = function(nb_sample, alpha, beta, lambda, param, seq_param)
{
  floop = function(i){
    alpha = if_else(param == 'alpha', i, alpha)
    beta = if_else(param == 'beta', i, beta)
    lambda = if_else(param == 'lambda', i, lambda)

    sigma = rinvgamma(n = nb_sample, alpha = alpha, beta = beta)
    mu = rnorm(nb_sample, mean = rep(0, nb_sample), sd = sqrt((1/lambda) * sigma))
    
    param_value = case_when(
      param == 'alpha' ~ alpha, 
      param == 'beta' ~ beta,
      param == 'lambda' ~ lambda
    )
    
    tibble(Mean = mu, Param = i) %>% 
      return()
  }
  lapply(seq_param, floop) %>% 
    bind_rows() %>%
    return()
}


text_param = 'lambda'
db = sample_prior(nb_sample, alpha, beta, lambda, text_param, seq(0.1, 10, 1))


gif = ggplot(db) + 
  geom_density(aes(x = Mean), col = "#DB15C1", fill = "#FA9FB5") + 
  geom_text(aes(x= -10, y= 0.3,
                label=paste("Value of",text_param , ": ", Param)), 
            hjust=0)+
  geom_vline(aes(xintercept = 0), col = "black") +
  theme_classic() + xlim(c(-10, 10)) + 
  transition_states(Param,
                    transition_length = 2,
                    state_length = 1)

anim_save('Priors_study/gif_lambda.gif', gif, height = 5, width = 8, 
          units = 'in', res = 300)


 

