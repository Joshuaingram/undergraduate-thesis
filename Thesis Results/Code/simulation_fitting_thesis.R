########################
# Joshua D. Ingram
# 5/21/2021
# Mu Data Simulations Results
########################

library(pracma)
library(tidyverse)
library(plotly)
library(progress)
source("D:/main/Projects/Solar-Flare-Project/Thesis Results/Code/mu_likelihood_functionsV2.R")
source("D:/main/Projects/Solar-Flare-Project/Thesis Results/Code/mu_data_simulation_functionsV2.R")

# directory for simulation data and all file names
directory <- "D:/main/Datasets/Solar_REU/thesis_simulations/simulations_bs/"
file_names <- dir(directory, pattern = ".csv")

# starting up progress bar
progress <- progress_bar$new(total = length(file_names), format = " Processing Simulations [:bar] :percent Elapsed Time: :elapsed")

sim_df <- data.frame(matrix(ncol = 24, nrow = 0))
colnames(sim_df) <- c("id", "n", "bins", "H0", "alpha" , "kappa", "beta", "gamma", "xi", "index", 
                      "convergence", "message", "runtime", "alpha_hat", "kappa_hat", "beta_hat", 
                      "gamma_hat", "xi_hat", "alpha_se", "kappa_se", "beta_se", "gamma_se", "xi_se",
                      "error")

write.csv(sim_df, "D:/main/Datasets/Solar_REU/thesis_simulations/results/bs_simulations_results.csv",
          row.names = FALSE,
          quote = FALSE)

for (i in 1:length(file_names)){
  
  skip <- FALSE
  
  # reading in data
  tryCatch({
    
    data <- read.csv(paste(directory, file_names[i], sep = ""), sep = ",")
    
  },
  error = function(e) { skip <<- TRUE}
  )
  
  if (skip == TRUE) {
    next
  }
  
  # getting all parameter values for simulations from file name 
  # sim id
  nchr <- nchar(file_names[i])
  id <- substr(file_names[i], 1, nchr-4)
  nchr <- nchar(id)
  # sample size
  n <- as.numeric(substr(file_names[i], 9, 11))
  # number of bins
  bins <- as.numeric(substr(file_names[i], 15, 16))
  # H min
  H0 <- as.numeric(substr(file_names[i], 18, 19))
  # alpha
  alpha <- as.numeric(substr(file_names[i], 21, 23))
  # kappa
  kappa <- as.numeric(substr(file_names[i], 25, 25))
  # beta
  beta <- as.numeric(substr(file_names[i], 27, 29))
  # gamma
  gamma <- as.numeric(substr(file_names[i], 31, 34))
  # xi
  xi <- as.numeric(substr(file_names[i], 36, 37))
  # index
  index <- as.numeric(substr(file_names[i], 39, nchr))
  
  error <- FALSE
  
  info <- c(id, n, bins, H0, alpha, kappa, beta, gamma, xi, index)
  
  # tries running fit... if error then moves on
  tryCatch(
    {
      start_time <- Sys.time()
      
      # initial value for optim
      theta_0 <- c()
      
      theta0 = c(alpha, kappa, beta, gamma, max(data$lb))
      
      ML_results <- optim(par = theta0, fn = ll_bs, df = data, H0 = 24,
                           lower = c(1.001, 0, 1, 0.0001, min(bins_data$ub)), upper = c(3, 12, 10, 10, 32),
                           hessian = TRUE, control = list(fnscale = -1), method = "L-BFGS-B")
      
      se_fit <- sqrt(diag(solve(-ML_results$hessian)))
      
      end_time <- Sys.time()
      
      # elapsed time ML fit
      ML_time <- difftime(end_time, start_time, units = "secs")[[1]]
    },
    
    error = function(e) {error <<- TRUE}
    
  )
    
    # progress bar tick
    progress$tick()
    
    if (error == FALSE){
      
      ML_info <- c(ML_results$convergence, ML_results$message, ML_time,  ML_results$par[1],  ML_results$par[2], ML_results$par[3], ML_results$par[4], ML_results$par[5], se_fit, 0)
      
    }else {
      
      ML_info <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 1)
      
    }
    
    sim_info <- c(info, ML_info)
    
    sim_df <- t(data.frame(sim_info))
    colnames(sim_df) <- c("id", "n", "bins", "H0", "alpha" , "kappa", "beta", "gamma", "xi", "index", 
                          "convergence", "message", "runtime", "alpha_hat", "kappa_hat", "beta_hat", 
                          "gamma_hat", "xi_hat", "alpha_se", "kappa_se", "beta_se", "gamma_se", "xi_se",
                          "error")
    
    
    write.table(sim_df, "D:/main/Datasets/Solar_REU/thesis_simulations/results/bs_simulations_results.csv",
                append = TRUE,
                sep = ",",
                col.names = FALSE,
                row.names = FALSE,
                quote = FALSE)
  
}

