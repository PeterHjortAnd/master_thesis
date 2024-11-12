library(MASS)
library(tidyverse)
library(complex)
library(knitr)

## Functions used throughout the document
Wstat <- function(L, K, p, t = 1, X) {
  # Initialize W to 0
  W_value <- 0
  
  # Precompute constants
  Kseq <- seq(0, t, length.out = K+1)
  Lseq <- seq(0, t, length.out = L+1)
  L_K <- L / K
  
  # Precompute differences in X for all L partitions
  deltaX_L_array <- abs(diff(X))^p
  
  # Iterate over the K indices
  for (i in 1:K) {
    # Calculate the difference in X for the K partition
    deltaX_K <- abs(X[1 + i * L_K] - X[1 + (i - 1) * L_K])^p
    
    # Precompute the denominator for the L partition over the i-th K partition
    deltaX_L <- sum(deltaX_L_array[(1 + (i - 1) * L_K):(i * L_K)])
    
    # Accumulate the sum
    W_value <- W_value + deltaX_K / deltaX_L * (Kseq[i + 1] - Kseq[i])
  }
  
  return(W_value)
} #Improved Wstat fct
mqdelta <- function(cap_delta, q, X) {
  # Pre-allocate the result vector
  m_qdelta <- numeric(cap_delta)
  
  # Loop over cap_delta values
  for (i in 1:cap_delta) {
    # Extract elements with the specified step
    delta_elements <- X[seq(i, length(X), by = cap_delta)]
    
    # Calculate log differences and take absolute value raised to the power q
    vol_term <- abs(diff(log(delta_elements)))^q
    
    # Calculate the mean of vol_term and store in m_qdelta
    m_qdelta[i] <- mean(vol_term)
  }
  
  # Calculate the final mean and return its logarithm
  final_m <- mean(m_qdelta)
  return(log(final_m))
} # Optimized version
# For fBm spectral sim
a_plus <- function(j, lambda) 2 * pi * j + lambda
a_minus <- function(j, lambda) 2 * pi * j - lambda
B3_lambda_H <- function(lambda, H) {
  j_values <- 1:3
  term1 <- (a_plus(j_values, lambda))^(-2 * H - 1)
  term2 <- (a_minus(j_values, lambda))^(-2 * H - 1)
  sum_result <- sum(term1 + term2)
  
  # Additional terms for j=3 and j=4
  additional_terms <- sum((a_plus(c(3, 4), lambda))^(-2 * H) + (a_minus(c(3, 4), lambda))^(-2 * H)) / (8 * H * pi)
  
  return(sum_result + additional_terms)
}
f_lambda <- function(lambda, H, gamma_term = gamma(2 * H + 1), sine_term = sin(pi * H)) {
  if (lambda == 0) return(Inf)  # Handle singularity
  B3_term <- B3_lambda_H(lambda, H)
  2 * sine_term * gamma_term * (1 - cos(lambda)) * (abs(lambda)^(-2 * H - 1) + B3_term)
}
compute_a_k <- function(l, H) {
  f_tk <- sapply(t_k[1:(l + 1)], f_lambda, H = H)
  
  # Initialize a_k as a vector of complex numbers
  a_k <- complex(real = rep(0, 2 * l), imaginary = rep(0, 2 * l))
  
  # Vectorized calculation for a_k values
  indices <- 1:(l - 1)
  a_k[indices + 1] <- 0.5 * (U_0[indices] + 1i * U_1[indices]) * sqrt(f_tk[indices + 1] / l)
  a_k[l + 1] <- U_0[l] * sqrt(f_tk[l + 1] / l)
  indices_2 <- (l+1):(2*l-1)
  a_k[indices_2+1] <- 0.5 * (U_0[2 * l - indices_2] - 1i * U_1[2 * l - indices_2]) * sqrt(f_tk[2*l-indices_2+1] / l)
  
  return(a_k)
}


# Simulation from RSFV
## Parameters
T <- 4900 #2000 in volis rep
delta <- 1/5000 #1/20000 in volis rep
L <- T
K <- sqrt(L)
n <- round(T/delta, digits = 0)
H <- 0.5 # 0.14
nu <- 0.2 # 0.3 in vol is rough replication
m <- -7 #-5 in vol is rough replication
alpha <- 5*10^(-4)
l <- ifelse(n*3>=30000, n*3, 30000)
dt <- T / n 
time <- seq(0, T, by = dt)

set.seed(21)
U_0 <- rnorm(l)  # U^(0) for k = 0, ..., l-1
U_1 <- rnorm(l)  # U^(1) for k = 0, ..., l-1
t_k <- 0:l
for (k in 0:l){
  t_k[k+1] <- pi * k/l
}
U <- rnorm(n)

n_plus_1 <- n + 1
S_0 <- 1
roughness_fct <- function(X){
  objective_function <- function(p) {
    (Wstat(L = L, K = K, p, t=1 , X = X) - 1)^2  # Squared difference to minimize
  }
  return(1/optimize(objective_function, interval = c(1,30))$minimum)
}
Delta_values <- 1:50
log_Delta <- log(Delta_values)
q_values <- c(0.5, 1, 1.5, 2, 3)
smoothness_fct <- function(X){
  log_mq_list <- lapply(q_values, function(q) {
    sapply(Delta_values, function(delta) mqdelta(delta, q, X = X))  # Calculate log_mq for each delta
  })
  log_delta <- c(log_Delta, log_Delta, log_Delta, log_Delta, log_Delta)  # Vector of log(Δ) values
  log_m <- c(log_mq_list[[1]], log_mq_list[[2]], log_mq_list[[3]], log_mq_list[[4]], log_mq_list[[5]])      # Vector of log(m(q, Δ)) values
  group <- c(rep(1,length(Delta_values)), rep(2,length(Delta_values)), rep(3,length(Delta_values)), rep(4,length(Delta_values)), rep(5,length(Delta_values)))
  data <- data.frame(log_delta, log_m, group)
  data_split <- split(data, data$group)
  models <- lapply(data_split, function(group_data) {
    lm(log_m ~ log_delta, data = group_data)
  })
  epsilon_q <- numeric()
  for (i in 1:length(q_values)){
    epsilon_q[i] <- models[[i]]$coefficients[[2]]
  }
  fit_epsilon <- lm(epsilon_q ~ q_values)
  return(fit_epsilon$coefficients[[2]])
}

output_given_H <- function(H){
  
  d_fbm <- Re(fft(compute_a_k(l, H))[1:n])*(T/n)^H
  
  X <- numeric(n_plus_1)
  X[1] <- m
  for (i in 2:n){
    X[i] = nu*(W[i]-W[i-1])+alpha*delta*(m-X[i-1])+X[i-1]
  } #RFSV simulation
  vol <- exp(X) #Simulated volatility
  
  S <- numeric(n_plus_1)
  S[1] <- S_0
  for (i in 2:n) {
    S[i] <- S[i-1] + S[i-1] * vol[i-1] * sqrt(delta) * U[i-1]
  } #Efficient price simulation

  # Realized volatility
  log_returns <- diff(log(abs(S)))
  squared_log_returns <- log_returns^2
  # Realized volatility calculation (vectorized)
  real_vol <- sqrt(rowSums(matrix(squared_log_returns, nrow = T, ncol = 1/delta, byrow = TRUE)))
  RV <- real_vol
  IV <- sigma[seq(1, n, by = 1/delta)]
  
  plot(S[1:10000]-S_0, type = "l", ylab = "Price S_t", xlab = "", xaxt = "n")
  plot(RV[1:10000], type = "l", ylab = "Realized volatility", xlab = "", xaxt = "n")
  plot(IV[1:10000], type = "l", ylab = "Instantaneous volatility", xlab = "", xaxt = "n")  
  
  H_IV <- roughness_fct(IV)
  H_RV <- roughness_fct(RV)
  H_IV_lres <- smoothness_fct(IV)
  H_RV_lres <- smoothness_fct(RV)
  
  return(list(H_IV = H_IV, H_RV = H_RV, H_IV_lres = H_IV_lres, H_RV_lres = H_RV_lres))
}

# Display the final table using knitr::kable
kable(results, col.names = c("H", "Instantaneous volatility", "Realized volatility", "IV (log-res)", "RV (log-res)"), 
      align = "c", caption = "Volatility Table")

plot(results$H, results$Instantaneous_volatility, type = "l", col = "blue", xlab = "", 
     ylab = "", xlim = c(0,0.8), ylim = c(0,0.8))
abline(v = seq(0, 0.8, by = 0.1), col = "gray", lty = "dotted", lwd = 0.5) # Vertical lines
abline(h = seq(0, 0.8, by = 0.1), col = "gray", lty = "dotted", lwd = 0.5) # Horizontal lines
lines(results$H, results$Realized_volatility, col = "red", lwd = 2)
lines(results$H, results$IV_logres, col = "blue", lwd = 2, lty = 2)
lines(results$H, results$RV_logres, col = "red", lwd = 2, lty = 2)

legend("topleft", legend = c('Realized Vol', 'Instantaneous Vol'), col = c('red', 'blue'), lty = 1, cex = 0.7)


