########################
# Joshua D. Ingram
# 5/13/2021
# Mu and Likelihood Functions V2
########################

library(pracma)

# mu(H) function for bounded sigmoid probability model. Used to take numeric integral for mu_i
mu_bs <- function(H, theta, H0){
  
  # parameters for bounded sigmoid model
  alpha <- theta[1]
  kappa <- theta[2]
  beta <- theta[3]
  gamma <- theta[4]
  xi <- theta[5]
  
  # lambda function
  # C <- H0^(alpha - 1) * (alpha - 1) * log(10)
  # lambda <- C * 10^(-H*(alpha - 1) + kappa)
  lambda <- (10^H0)^(alpha - 1)*(alpha - 1) * (10^H)^(-alpha) * log(10) * 10^H * 10^kappa
  
  # bounded sigmoid pi function
  f <- (H - H0)/(xi - H)
  pi <- (gamma * f^beta)/(1 + gamma * f^beta)
  pi[H <= H0] <- 0
  pi[H >= xi] <- 1
  pi[is.infinite(pi)] <- 1
  
  # mu function
  mu <- lambda * pi
  
  return(mu)
  
}

# mu(H) function for 4-parameter logistic probability model. Used to take numeric integral for mu_i
mu_4pl <- function(H, theta, H0){
  
  # parameters for 4-parameter logistic model
  alpha <- theta[1]
  kappa <- theta[2]
  b <- theta[3]
  g <- theta[4]
  c <- 1
  d <- 0
  
  # lambda function
  
  #C <- H0^(alpha - 1) * (alpha - 1) * log(10)
  #lambda <- C * 10^(-H*(alpha - 1) + kappa)
  lambda <- (alpha - 1)*(10^H0)^(alpha - 1)*log(10)*10^(-H*(alpha - 1) + kappa)
  
  # 4-parameter logistic pi function
  pi <- c + (d - c)/(1 + exp(b*(H-g)))
  
  # mu function
  mu <- lambda * pi
  
  return(mu)
  
}

# log likelihood function for bounded sigmoid model that we minimize (because we take negative of the sum) using optim function
ll_bs <- function(theta, H0, df, print = FALSE){
  
  if (print == TRUE){
    
    print(theta)
    
  }
  
  # bin counts
  y <- df[,1]
  # bounds of bins
  Hmin <- df[,2]
  Hmax <- df[,3]
  bounds <- data.frame(min = Hmin, max = Hmax)
  
  # mui function after taking numeric integral
  mui <- apply(bounds, 1, function(x) integral(fun = mu_bs, method = "Kron", xmin = x[1], xmax = x[2], theta = theta, H0 = H0)/ diff(x))
  mui <- mui[mui > 0]
  ind <- which((mui > 0))
  
  # negative of function proportional to log-likelihood
  y <- y[ind]
  ll <- sum(y*log(mui) - mui)
  
  return(ll)
  
}

# log likelihood function for 4-parameter logistic model that we minimize (because we take negative of the sum) using optim function
ll_4pl <- function(theta, H0, df, print = FALSE){
  
  if (print == TRUE){
    
    print(theta)
    
  }
  
  # bin counts
  y <- df[,1]
  # bounds of bins
  Hmin <- df[,2]
  Hmax <- df[,3]
  bounds <- data.frame(min = Hmin, max = Hmax)
  
  # mui function after taking numeric integral
  mui <- apply(bounds, 1, function(x) integral(fun = mu_4pl, xmin = x[1], xmax = x[2], theta = theta, H0 = H0)/ diff(x))
  
  # negative of function proportional to log-likelihood
  ll <- sum(y*log(mui) - mui)
  
  return(ll)
  
}

# mu(H) function for 4-parameter logistic probability model. Used to take numeric integral for mu_i
mu_test <- function(H, theta, H0){
  
  # parameters for 4-parameter logistic model
  alpha <- theta[1]
  kappa <- theta[2]
  b <- theta[3]
  g <- theta[4]
  c <- 1
  d <- 0
  
  # lambda function
  # C <- *H0^(alpha - 1) * (alpha - 1) * log(10)
  lambda <- kappa^2 * 10^(-H*(alpha - 1))
  
  # 4-parameter logistic pi function
  pi <- c + (d - c)/(1 + exp(b*(H-g)))
  pi[is.infinite(pi)] <- 1
  
  # mu function
  mu <- lambda * pi
  
  return(mu)
  
}
# log likelihood function for 4-parameter logistic model that we minimize (because we take negative of the sum) using optim function
ll_test <- function(theta, H0, df){
  
  # bin counts
  y <- df[,1]
  # bounds of bins
  Hmin <- df[,2]
  Hmax <- df[,3]
  bounds <- data.frame(min = Hmin, max = Hmax)
  
  print(theta)
  
  # mui function after taking numeric integral
  mui <- apply(bounds, 1, function(x) integral(fun = mu_test, xmin = x[1], xmax = x[2], theta = theta, H0 = H0)/ diff(x))
  mui <- mui[mui > 0]
  ind <- which((mui > 0))
  
  # negative of function proportional to log-likelihood
  y <- y[ind]
  ll <- -sum(y*log(mui) - mui)
  
  return(ll)
  
}

# mu(H) function for bounded sigmoid probability model. Used to take numeric integral for mu_i
mu_test2 <- function(H, theta, H0){
  
  # parameters for bounded sigmoid model
  alpha <- theta[1]
  kappa <- theta[2]
  beta <- theta[3]
  gamma <- theta[4]
  xi <- theta[5]
  
  # lambda function
  C <- H0^(alpha - 1) * (alpha - 1) * log(10)
  lambda <- 10^(-H*(alpha - 1) + kappa)
  
  # bounded sigmoid pi function
  f <- (H - H0)/(xi - H)
  pi <- (gamma * f^beta)/(1 + gamma * f^beta)
  pi[H <= H0] <- 0
  pi[H >= xi] <- 1
  pi[is.infinite(pi)] <- 1
  
  # mu function
  mu <- lambda * pi
  
  return(mu)
  
}

# log likelihood function for bounded sigmoid model that we minimize (because we take negative of the sum) using optim function
ll_test2 <- function(theta, H0, df){
  
  # bin counts
  y <- df[,1]
  # bounds of bins
  Hmin <- df[,2]
  Hmax <- df[,3]
  bounds <- data.frame(min = Hmin, max = Hmax)
  print(theta)
  
  # mui function after taking numeric integral
  mui <- apply(bounds, 1, function(x) integral(fun = mu_test2, xmin = x[1], xmax = x[2], theta = theta, H0 = H0)/ diff(x))
  mui <- mui[mui > 0]
  ind <- which((mui > 0))
  
  # negative of function proportional to log-likelihood
  y <- y[ind]
  ll <- -sum(y*log(mui) - mui)
  
  return(ll)
  
}