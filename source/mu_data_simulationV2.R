########################
# Joshua D. Ingram
# 5/13/2021
# Mu Data Simulations Creation
########################

library(poweRlaw)
library(tidyverse)
library(readr)
library(pracma)
library(progress)
source("D:/main/Projects/Solar-Flare-Project/Thesis Results/Code/mu_likelihood_functionsV2.R")
source("D:/main/Projects/Solar-Flare-Project/Thesis Results/Code/mu_data_simulation_functionsV2.R")

# sample sizes and bin counts
sample_sizes <- c(250, 500, 900)
bin_counts <- c(30, 50, 75)

# number of simulations per dataset to conduct
nsims <- 500

# set H0 value
H0 <- 24

# parameter sets
# for bs model - theta(alpha, kappa, beta, gamma, xi)
theta_bs1 <- c(1.8, 2, 2.5, 0.25, 27)
theta_bs2 <- c(2.2, 2, 2.5, 0.25, 26)
theta_bs3 <- c(1.4, 2, 2.5, 0.25, 27)
theta_bs4 <- c(1.8, 2, 2.0, 0.15, 26)
theta_bs5 <- c(2.2, 2, 3.4, 0.50, 27)

# for 4pl model - theta(alpha, kappa, b, g)
theta_4pl1 <- c(1.8, 23, 08, 24.5)
theta_4pl2 <- c(2.2, 23, 08, 24.5)
theta_4pl3 <- c(1.4, 23, 08, 24.5)
theta_4pl4 <- c(1.8, 23, 08, 24.4)
theta_4pl5 <- c(2.2, 23, 09, 24.7)

thetas_bs <- data.frame(alpha = numeric(0), kappa = numeric(0), beta = numeric(0), gamma = numeric(0), xi = numeric(0))
thetas_4pl <- data.frame(alpha = numeric(0), kappa = numeric(0), b = numeric(0), g = numeric(0))

thetas_bs <- rbind(theta_bs1, theta_bs2, theta_bs3, theta_bs4, theta_bs5)
colnames(thetas_bs) <- c("alpha", "kappa", "beta", "gamma", "xi")
thetas_4pl <- rbind(theta_4pl1, theta_4pl2, theta_4pl3, theta_4pl4, theta_4pl5)
colnames(thetas_4pl) <- c("alpha", "kappa", "b", "g")

# # sample sizes and bin counts
# sample_sizes <- c(500, 600)
# bin_counts <- c(40, 45, 50)
# 
# # number of simulations per dataset to conduct
# nsims <- 100
# 
# # set H0 value
# H0 <- 24
# 
# thetas_bs <- rbind(theta_bs1, theta_bs2)
# colnames(thetas_bs) <- c("alpha", "kappa", "beta", "gamma", "xi")
# thetas_4pl <- rbind(theta_4pl1, theta_4pl2)
# colnames(thetas_4pl) <- c("alpha", "kappa", "b", "g")

total <- nrow(thetas_bs) * length(sample_sizes) * length(bin_counts) * nsims
# starting up progress bar
progress <- progress_bar$new(total = total, 
                             format = " Creating Simulated Data [:bar] :percent Elapsed Time: :elapsed")

# for each theta, sample size, and bin count create a simulated data set
for (parms in 1:nrow(thetas_bs)){
  
  for (n in 1:length(sample_sizes)){
    
    for (bin in 1:length(bin_counts)){
      
      for (i in 1:nsims){
        
        sim_data_bs <- rmu_bs(sample_sizes[n], thetas_bs[parms,], H0 = H0)
        sim_data_bs <- data.frame(values = sim_data_bs)
        # creating binned data
        sim_bin_bs <- ggplot(data = sim_data_bs, aes(x=values)) + geom_histogram(bins = bin_counts[bin])
        sim_bin_bs <- ggplot_build(sim_bin_bs)$data[[1]]
        sim_bin_bs <- data.frame(count = sim_bin_bs$count, lb = sim_bin_bs$xmin, ub = sim_bin_bs$xmax, mid = sim_bin_bs$x)
        
        file_name_bs <- paste0("sim_bs_", "s", sample_sizes[n], "bin", bin_counts[bin], 
                               "H", H0, "a", as.numeric(thetas_bs[parms,])[1], "k", as.numeric(thetas_bs[parms,])[2], "b", as.numeric(thetas_bs[parms,])[3],
                               "g", as.numeric(thetas_bs[parms,])[4], "x", as.numeric(thetas_bs[parms,])[5], "i", i)
        
        write.csv(sim_bin_bs, 
                  file = paste0("D:/main/Datasets/Solar_REU/thesis_simulations/simulations_bs/", file_name_bs, ".csv"), 
                  row.names = FALSE)
        
        # # simulated data for 4pl model
        # sim_data_4pl <- rmu_bs(sample_sizes[n], thetas_bs[parms,], H0 = H0)
        # sim_data_4pl <- data.frame(values = sim_data_4pl)
        # # creating binned data
        # sim_bin_4pl <- ggplot(data = sim_data_4pl, aes(x=values)) + geom_histogram(bins = bin_counts[bin])
        # sim_bin_4pl <- ggplot_build(sim_bin_4pl)$data[[1]]
        # sim_bin_4pl <- data.frame(count = sim_bin_4pl$count, lb = sim_bin_4pl$xmin, ub = sim_bin_4pl$xmax, mid = sim_bin_4pl$x)
        # 
        # file_name_4pl <- paste0("sim_4pl_", "s", sample_sizes[n], "bin", bin_counts[bin], 
        #                         "H", H0, "a", as.numeric(thetas_4pl[parms,])[1], "k", as.numeric(thetas_4pl[parms,])[2], "b", as.numeric(thetas_4pl[parms,])[3],
        #                         "g", as.numeric(thetas_4pl[parms,])[4], "i", i)
        # 
        # write.csv(sim_bin_4pl, 
        #           file = paste0("D:/main/Datasets/Solar_REU/thesis_simulations/simulations_4pl/", file_name_4pl, ".csv"), 
        #           row.names = FALSE)
        
        # progress bar tick
        progress$tick()
        
      }
        
    }
    
  }
  
}
