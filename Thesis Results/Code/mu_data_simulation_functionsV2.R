########################
# Joshua D. Ingram
# 5/13/2021
# Mu Data Simulation Functions V2
########################

library(poweRlaw)
library(tidyverse)
library(readr)
library(pracma)
source("D:/main/Projects/Solar-Flare-Project/Thesis Results/Code/mu_likelihood_functionsV2.R")

pi_bs <- function(H, theta, H0){
  
  # parameters for bounded sigmoid model
  alpha <- theta[1]
  kappa <- theta[2]
  beta <- theta[3]
  gamma <- theta[4]
  xi <- theta[5]
  
  # bounded sigmoid pi function
  f <- (H - H0)/(xi - H)
  pi <- (gamma * f^beta)/(1 + gamma * f^beta)
  pi[H < H0] <- 0
  pi[H > xi] <- 1
  
  return(pi)
  
}

# simulates data from bounded sigmoid model for sample size n, with total=1mil being original amount drawn (n < total)
rmu_bs <- function(n, theta, H0, total = 10^4.6){
  
  # parameters for bounded sigmoid model
  alpha <- theta[1]
  kappa <- theta[2]
  beta <- theta[3]
  gamma <- theta[4]
  xi <- theta[5]
  
  # randomly drawing n values from powerlaw and putting in log_10 scale
  pow_sim <- log(rplcon(n = total, xmin = 10^H0, alpha = alpha), 10)
  
  # probability of detecting data from power-law
  prob_detect <- pi_bs(pow_sim, theta, H0)
  
  # detection (0 if not detected, 1 if detected) based on probabilities
  detection <- rbinom(n = total, size = 1, prob = prob_detect)
  
  # dataframe of simulated values
  sim_data <- data.frame(pow_sim, prob_detect, detection)
  
  # simulated that will be fit to (what would be observed)
  rmu_sim <- sim_data %>% filter(detection == 1) %>% select(pow_sim)
  rmu_sim <- sample(rmu_sim[,1], n)
  
  # creating binned data
  # rmui_sim <- ggplot(data = rmu_sim, aes(x= pow_law)) + geom_histogram(bins = bins)
  # rmui_sim <- ggplot_build(rmui_sim)$data[[1]]
  # rmui_sim <- data.frame(count = rmui_sim$count, lb = rmui_sim$xmin, ub = rmui_sim$xmax, mid = rmui_sim$x)
  
  return(rmu_sim)
  
}

pi_4pl <- function(H, theta, H0){
  
  # parameters for 4-parameter logistic model
  alpha <- theta[1]
  kappa <- theta[2]
  b <- theta[3]
  g <- theta[4]
  c <- 1
  d <- 0
  
  # 4-parameter logistic pi function
  pi <- c + (d - c)/(1 + exp(b*(H-g)))
  
  return(pi)
  
}

# simulates data from bounded sigmoid model for sample size n
rmu_4pl <- function(n, theta, H0, total = 10^3.7){
  
  # parameters for 4-parameter logistic model
  alpha <- theta[1]
  kappa <- theta[2]
  b <- theta[3]
  g <- theta[4]
  c <- 1
  d <- 0
  
  # randomly drawing n values from powerlaw and putting in log_10 scale
  pow_sim <- log(rplcon(n = total, xmin = 10^H0, alpha = alpha), 10)
  
  # probability of detecting data from power-law
  prob_detect <- pi_4pl(pow_sim, theta, H0)
  
  # detection (0 if not detected, 1 if detected) based on probabilities
  detection <- rbinom(n = total, size = 1, prob = prob_detect)
  
  # dataframe of simulated values
  sim_data <- data.frame(pow_sim, prob_detect, detection)
  
  # simulated that will be fit to (what would be observed)
  rmu_sim <- as.vector(sim_data %>% filter(detection == 1) %>% select(pow_sim))
  rmu_sim <- sample(rmu_sim[,1], n)
  
  # creating binned data
  # rmui_sim <- ggplot(data = rmu_sim, aes(x= pow_law)) + geom_histogram(bins = bins)
  # rmui_sim <- ggplot_build(rmui_sim)$data[[1]]
  # rmui_sim <- data.frame(count = rmui_sim$count, lb = rmui_sim$xmin, ub = rmui_sim$xmax, mid = rmui_sim$x)
  
  return(rmu_sim)
  
}