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
f_lambda <- function(lambda, H, start_term = 2 * gamma(2 * H + 1) * sin(pi * H)) {
  if (lambda == 0) return(Inf)  # Handle singularity
  B3_term <- B3_lambda_H(lambda, H)
  start_term * (1 - cos(lambda)) * (abs(lambda)^(-2 * H - 1) + B3_term)
}
compute_a_k <- function(l, H) {
  #f_tk <- sapply(t_k[1:(l_plus_one)], f_lambda, H = H)
  
  # Find the row index corresponding to the given H in H_values
  i <- which(H_values == H)
  # Check if H is in H_values
  if (length(i) == 0) {
    stop("The provided H value is not in the precomputed H_values.")
  }
  # Extract the precomputed f_tk for this H value
  f_tk <- f_tk_mat[i, ]
  
  # Initialize a_k as a vector of complex numbers
  a_k <- complex(real = rep(0, two_l), imaginary = rep(0, two_l))
  
  # Vectorized calculation for a_k values
  a_k[indi_plus_1] <- ak_constant * sqrt(f_tk[indi_plus_1] / l)
  a_k[l_plus_one] <- U_0[l] * sqrt(f_tk[l_plus_one] / l)
  a_k[indi2_plus_1] <- ak_constant_2 * sqrt(f_tk[indi_2l] / l)
  
  return(a_k)
}
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
rough_expo <- function(k, y_t, n_max) {
  # k: Current level
  # y_t: Observations at the finest resolution (n_max)
  # n_max: Maximum level of discretization (finest resolution)
  
  var_theta <- numeric(2^k)  # Initialize the ϑ_k,j array
  
  # Number of points in y_t corresponds to the finest resolution
  num_points <- 2^(n_max + 2) + 1
  
  # Time step at the finest resolution
  dt_fine <- 1 / (2^(n_max + 2))
  
  # Scaling factor for time points to match the current k
  dt_current <- 1 / (2^(k + 2))
  scale_factor <- dt_current / dt_fine  # Should equal 2^(n_max - k)
  
  for (j in 0:(2^k - 1)) {
    # Compute indices for current resolution
    idx1 <- round(4 * j * scale_factor) + 1
    idx2 <- round((4 * j + 1) * scale_factor) + 1
    idx3 <- round((4 * j + 3) * scale_factor) + 1
    idx4 <- round((4 * j + 4) * scale_factor) + 1
    
    # Ensure indices are within bounds
    if (idx4 > num_points) {
      stop("Index out of bounds: Ensure y_t has enough points for the maximum level n_max.")
    }
    
    # Compute var_theta for the current j
    var_theta[j + 1] <- 2^(3 * k / 2 + 3) * (
      y_t[idx1] - 2 * y_t[idx2] + 2 * y_t[idx3] - y_t[idx4]
    )
  }
  
  # Compute the sum of squared θ terms
  sum_theta_terms <- sum(var_theta^2)
  
  # Compute R_k(y_t)
  r_hat <- 1 - 1/k * log2(sqrt(sum_theta_terms))
  return(r_hat)
} # For sequential scale estimator
roughness_fct <- function(X, threshold = 0.1) {
  objective_function <- function(H_hat) {
    (Wstat(L = L, K = K, 1 / H_hat, t = 1, X = X) - 1)^2
  }
  
  # Perform optimization
  result <- optimize(objective_function, interval = c(0.0001, 0.9999))
  
  # Return roughness estimate if objective is within threshold, else 0
  return(ifelse(result$objective <= threshold, result$minimum, 0))
} # For p-th variation estimation
m <- 3
alpha_all <- c(1,1,1,1)
scale_est_fct <- function(X){
  objective_function <- function(lambda) {
    sum <- 0
    for (k in (small_n - m):small_n) {
      r_k <- rough_expo(k, lambda * X, n_max = small_n)
      r_k_minus_1 <- rough_expo(k - 1, lambda * X, n_max = small_n)
      sum <- sum + alpha_all[small_n - k + 1] * (r_k - r_k_minus_1)^2
    }
    return(sum)
  }
  lambda_opt <- optimize(objective_function, interval = c(0.001,300000000))$minimum
  est <- rough_expo(k = small_n, y_t = lambda_opt*X, n_max = small_n)
  return(est)
} # Total sequential scale estimator
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
} # Log regression estimation


# Simulation from RSFV
## Parameters
L <- 150 * 150
K <- 150
window <- 300
n <- window * (L+1)
T <- 5 
dt <- T/n

#H <- 0.5 # 0.14
nu <- 0.3 # 0.3 in vol is rough replication
x_0 <- -5 #-5 in vol is rough replication
alpha <- 5*10^-4 # 5*10^(-4) in rep


## Used in both settings
l <- ifelse(n*3>=30000, n*3, 30000)
dt <- T / n 
#time <- seq(0, T, by = dt)

set.seed(21)
U_0 <- rnorm(l)  # U^(0) for k = 0, ..., l-1
U_1 <- rnorm(l)  # U^(1) for k = 0, ..., l-1
t_k <- 0:l
for (k in 0:l){
  t_k[k+1] <- pi * k/l
}
U <- rnorm(n)

n_plus_1 <- n + 1
window_length <- floor(n_plus_1 / window)
sqrt_term <- sqrt(n_plus_1 / (T * window))
S_0 <- 1


## Constants for fBm sim
indices_1 <- 1:(l - 1)
ak_constant <- 0.5 * (U_0[indices_1] + 1i * U_1[indices_1])
indices_2 <- (l+1):(2*l-1)
ak_constant_2 <- 0.5 * (U_0[2 * l - indices_2] - 1i * U_1[2 * l - indices_2])
indi_plus_1 <- indices_1 + 1
indi2_plus_1 <- indices_2+1
indi_2l <- 2*l-indices_2+1
l_plus_one <- l +1 
two_l <- 2 * l


# Sequential scale estimator
small_n <- 9
capital_N <- small_n + 13
number_of_points <- 2^capital_N # Excluding T=0
T_new <- 1
dt_new <- T_new/number_of_points
k <- 1:(2^(small_n+2))  # Sequence for k
indices <- 2^(capital_N - small_n - 2) * k  # Compute the desired indices


output_given_H <- function(H){
  
  d_fbm <- Re(fft(compute_a_k(l, H))[1:n])*(dt)^H
  
  X <- numeric(n_plus_1)
  X[1] <- x_0
  for (i in 1:n){
    X[i+1] = nu*(d_fbm[i])+alpha*dt*(x_0-X[i])+X[i]
  } #RFSV simulation
  vol <- exp(X) #Simulated volatility
  
  S <- numeric(n_plus_1)
  S[1] <- S_0
  for (i in 1:n) {
    S[i+1] <- S[i] + S[i] * vol[i] * sqrt(dt) * U[i]
  } #Efficient price simulation
  
  # Realized volatility
  log_returns <- diff(log(S))
  squared_log_returns <- log_returns^2
  
  # Calculate Realized Volatility (RV)
  real_vol <- sqrt(rowSums(matrix(squared_log_returns, nrow = window_length, ncol = window, byrow = TRUE)))
  RV <- real_vol*sqrt_term
  IV <- vol[seq(1, n, by = window)]
  
  #plot(S[1:10000]-S_0, type = "l", ylab = "Price S_t", xlab = "", xaxt = "n")
  #plot(RV[1:10000], type = "l", ylab = "Realized volatility", xlab = "", xaxt = "n")
  #plot(IV[1:10000], type = "l", ylab = "Instantaneous volatility", xlab = "", xaxt = "n")  
  
  H_IV <- roughness_fct(IV)
  H_RV <- roughness_fct(RV)
  H_IV_lres <- smoothness_fct(IV)
  H_RV_lres <- smoothness_fct(RV)
  
  dfbm_new <- (dt_new/dt)^H*d_fbm[1:number_of_points]
  
  X <- numeric(number_of_points + 1)
  X[1] <- x_0
  for (i in 1:number_of_points){
    X[i+1] = nu*(dfbm_new[i])+alpha*dt_new*(x_0-X[i])+X[i]
  } #RFSV simulation
  vol <- exp(X) #Simulated volatility
  
  S <- numeric(number_of_points + 1)
  S[1] <- S_0
  for (i in 1:number_of_points) {
    S[i+1] <- S[i] + S[i] * vol[i] * sqrt(dt_new) * U[i]
  }
  
  log_returns <- diff(log(S))
  squared_log_returns <- log_returns^2
  cumulative_sum <- cumsum(squared_log_returns)
  
  # Extract the elements
  result <- cumulative_sum[indices]
  RV_scale <- c(0, result)
  
  H_RV_scale <- scale_est_fct(RV_scale)
  
  sum_terms <- cumsum(exp(2*X[2:(number_of_points+1)]))
  result <- sum_terms[indices]
  IV_scale <- c(0, 2^-capital_N * result)
  
  H_IV_scale <- scale_est_fct(IV_scale)
  
  return(list(H_IV = H_IV, H_RV = H_RV, H_IV_lres = H_IV_lres, H_RV_lres = H_RV_lres, H_IV_scale = H_IV_scale, H_RV_scale = H_RV_scale))
}

H_values <- seq(0.10, 0.80, by = 0.10)
f_tk_mat <- matrix(0, nrow = length(H_values), ncol = l + 1)
# Loop over H values with indexing
for (i in seq_along(H_values)) {
  H <- H_values[i]
  f_tk_mat[i, ] <- sapply(t_k[1:(l + 1)], f_lambda, H = H)
}

store_results <- list()
for (H in H_values) {
  store_results[[paste0("H_", H)]] <- output_given_H(H)
}
results <- data.frame(H = numeric(), Instantaneous_volatility = numeric(), Realized_volatility = numeric(), IV_logres = numeric(), RV_logres = numeric(),
                      IV_scale_vol = numeric(), RV_scale_vol = numeric())
for (H in H_values) {
  IV_H_est <- store_results[[paste0("H_", H)]]$H_IV
  RV_H_est <- store_results[[paste0("H_", H)]]$H_RV
  IV_H_est_logres <- store_results[[paste0("H_", H)]]$H_IV_lres
  RV_H_est_logres <- store_results[[paste0("H_", H)]]$H_RV_lres
  IV_H_est_scale <- store_results[[paste0("H_", H)]]$H_IV_scale
  RV_H_est_scale <- store_results[[paste0("H_", H)]]$H_RV_scale
  # Append the results to the data frame
  results <- rbind(results, data.frame(H = H, 
                                       Instantaneous_volatility = IV_H_est, 
                                       Realized_volatility = RV_H_est, 
                                       IV_logres = IV_H_est_logres, 
                                       RV_logres = RV_H_est_logres,
                                       IV_scale_vol = IV_H_est_scale,
                                       RV_scale_vol = RV_H_est_scale))
}
# Display the final table using knitr::kable
kable(results, col.names = c("H", "Instantaneous volatility", "Realized volatility", "IV (log-res)", "RV (log-res)", "IV (Scale)", "RV (Scale)"), 
      align = "c", caption = "Volatility Table")

plot(results$H, results$Instantaneous_volatility, type = "l", col = "blue", xlab = "", 
     ylab = "", xlim = c(0,0.8), ylim = c(0,0.8), lwd = 1)
abline(v = seq(0, 0.8, by = 0.1), col = "gray", lty = "dotted", lwd = 0.5) # Vertical lines
abline(h = seq(0, 0.8, by = 0.1), col = "gray", lty = "dotted", lwd = 0.5) # Horizontal lines
lines(results$H, results$Realized_volatility, col = "red", lwd = 1)
lines(results$H, results$IV_logres, col = "blue", lwd = 1, lty = 2)
lines(results$H, results$RV_logres, col = "red", lwd = 1, lty = 2)
lines(results$H, results$IV_scale_vol, col = "blue", lwd = 1, lty = 3)
lines(results$H, results$RV_scale_vol, col = "red", lwd = 1, lty = 3)

legend("topleft", legend = c('RV (p-th)', 'IV (p-th)', 'RV (log reg)', 
                             'IV (log reg)', 'RV (scale)', 'IV (scale)'), col = c('red', 'blue', 'red', 'blue','red','blue'), 
       lty = c(1,1,2,2,3,3) , cex = 0.7)


## Simulate 100 times

simoutput_given_H <- function(H){
  
  d_fbm <- Re(fft(compute_a_k(l, H))[1:n])*(dt)^H
  
  X <- numeric(n_plus_1)
  X[1] <- x_0
  for (i in 1:n){
    X[i+1] = nu*(d_fbm[i])+alpha*dt*(x_0-X[i])+X[i]
  } #RFSV simulation
  vol <- exp(X) #Simulated volatility
  
  S <- numeric(n_plus_1)
  S[1] <- S_0
  for (i in 1:n) {
    S[i+1] <- S[i] + S[i] * vol[i] * sqrt(dt) * U[i]
  } #Efficient price simulation
  
  # Realized volatility
  log_returns <- diff(log(abs(S)))
  squared_log_returns <- log_returns^2
  
  # Calculate Realized Volatility (RV)
  real_vol <- sqrt(rowSums(matrix(squared_log_returns, nrow = window_length, ncol = window, byrow = TRUE)))
  RV <- real_vol*sqrt_term
  IV <- vol[seq(1, n, by = window)]
  
  IV_H <- smoothness_fct(IV)
  RV_H <- smoothness_fct(RV)
  
  return(list(IV_H = IV_H, RV_H = RV_H))
}
sim_vol_table <- function(H_values){
  volatility_table <- data.frame(
    'Instantaneous volatility' = numeric(0),
    'Realized volatility' = numeric(0)
  )
  store_results <- list()
  for (H in H_values) {
    store_results[[paste0("H_", H)]] <- simoutput_given_H(H)
    new_row <- data.frame(
      'Instantaneous volatility' = store_results[[paste0("H_", H)]]$IV_H,
      'Realized volatility' = store_results[[paste0("H_", H)]]$RV_H
    )
    volatility_table <- rbind(volatility_table, new_row)
  }
  return(volatility_table)
}

simtim <- 10
simulation_results <- vector("list", simtim)
for (i in 1:simtim){
  U_0 <- rnorm(l)  # U^(0) for k = 0, ..., l-1
  U_1 <- rnorm(l)  # U^(1) for k = 0, ..., l-1
  ak_constant <- 0.5 * (U_0[indices_1] + 1i * U_1[indices_1])
  ak_constant_2 <- 0.5 * (U_0[2 * l - indices_2] - 1i * U_1[2 * l - indices_2])
  simulation_results[[i]] <- sim_vol_table(H_values)
}
print(simulation_results[[3]])

plot(H_values, simulation_results[[1]]$Instantaneous.volatility, type = "l", col = "blue", xlab = "", 
     ylab = "", xlim = c(0,0.8), ylim = c(0,0.8))
abline(v = seq(0, 0.8, by = 0.1), col = "gray", lty = "dotted", lwd = 0.5) # Vertical lines
abline(h = seq(0, 0.8, by = 0.1), col = "gray", lty = "dotted", lwd = 0.5) # Horizontal lines
legend("topleft", legend = c('Realized Vol', 'Instantaneous Vol'), col = c('red', 'blue'), lty = 1, cex = 0.7)

for (i in 1:simtim){
  lines(H_values, simulation_results[[i]]$Realized.volatility, col = "red", lwd = 1)
}
for (i in 2:simtim){
  lines(H_values, simulation_results[[i]]$Instantaneous.volatility, col = "blue", lwd = 1)
}




