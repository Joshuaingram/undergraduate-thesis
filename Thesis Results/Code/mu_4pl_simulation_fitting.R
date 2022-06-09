########################
# Joshua D. Ingram
# 5/13/2021
# Mu Data Simulation Fitting 4pl
########################

library(poweRlaw)
library(tidyverse)
library(readr)
library(pracma)
library(progress)
source("D:/main/Projects/Solar-Flare-Project/Thesis Results/Code/mu_likelihood_functionsV2.R")

# directory for simulation data and all file names
directory <- "D:/main/Datasets/Solar_REU/thesis_simulations/simulations_4pl/"
file_names <- dir(directory, pattern = ".csv")

fit_results <- data.frame(matrix(ncol = 19, nrow = 0))
colnames(fit_results) <- c("id", "n", "bin_count", "ealpha", "ekappa", "eb", "eg", "alpha", "kappa", "b", "g", "convergence", "error", "runtime", "H0", "alphase", "kappase", "bse", "gse")

# saving dataframe as a .csv file
write.csv(fit_results, "D:/main/Datasets/Solar_REU/thesis_simulations/results/mu_4pl_sim_results.csv", row.names = FALSE)

# starting up progress bar
progress <- progress_bar$new(total = length(file_names), format = " Processing Simulations [:bar] :percent Elapsed Time: :elapsed")

# for loop that loads in each dataset, finds ML fit, and adds results to a dataframe
for (i in 1:length(file_names)){
  
  skip <- FALSE
  
  # reading in data
  tryCatch({
    
    data <- read.csv(paste(directory, file_names[i], sep = ""), sep = ",")
    
  },
  error = function(e) { skip <<- TRUE})
  
  if (skip == TRUE) {
    next
  }
  
  # simulation ID
  id <- gsub('.{4}$', '', file_names[i])
  # sample size
  n <- as.numeric(substr(id, 10,12))
  # bin count
  bin_count <- as.numeric(substr(id, 16, 17))
  # H0
  H0 <- as.numeric(substr(id, 19, 20))
  # alpha
  alpha <- as.numeric(substr(id, 22, 24))
  # kappa
  kappa <- as.numeric(substr(id, 26, 27))
  # b
  b <- as.numeric(substr(id, 29, 29))
  # g
  g <- as.numeric(substr(id, 31, 34))
  
  error <- FALSE
  
  # tries running fit... if error then moves on
  tryCatch(
    {
      start_time <- Sys.time()
      
      ll_4pl <- function(theta, H0, df)
      
      # initial value for optim
      theta0 <- c(alpha, kappa, b, g)
      
      results <- optim(par = theta0, fn = ll_4pl, H0 = H0, df = data,
                       method = "L-BFGS-B",
                       lower = c(1.4, 15, 4, 20), 
                       hessian = TRUE)
      
      end_time <- Sys.time()
      
      # elapsed time ML fit
      runtime <- difftime(end_time, start_time, units = "secs")[[1]]
    },
    
    error = function(e) {error <<- TRUE}
    
  )
  
  # progress bar tick
  progress$tick()
  
  if (error == FALSE){
    
    se <- sqrt(diag(solve(results$hessian)))
    
    c("id", "n", "bin_count", "ealpha", "ekappa", "eb", "eg", "alpha", "kappa", "b", "g", "convergence", "error", "runtime", "H0", "alphase", "kappase", "bse", "gse")
    
    sim_4pl_results <- c(id, n, bin_count, results$par[1],  results$par[2], results$par[3], results$par[4], 
                         alpha, kappa,  b, g, results$convergence, 0, runtime, H0, se[1], se[2], se[3], se[4])
    
  }else {
    
    sim_4pl_results <- c(id, n, bin_count, NA, NA, NA, NA, alpha, kappa,  b, g, NA, 1, NA, H0, NA, NA, NA, NA)
    
  }
  
  # write.table(sim_4pl_results, 
  #             file = "D:/main/Datasets/Solar_REU/thesis_simulations/results/mu_4pl_sim_results.csv",
  #             append = TRUE,
  #             sep = ",",
  #             row.names = FALSE,
  #             col.names = FALSE,
  #             quote = FALSE)
  
  print(sim_4pl_results)
  
}

