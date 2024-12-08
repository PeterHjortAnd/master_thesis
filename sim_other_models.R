library(MASS)
library(tidyverse)
library(complex)
library(knitr)
library(data.table)
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


### Example 5
# Parameters
L = 500*500
K = 500
window <- 300
n <- window*(L+1)
#n <- window*10000
T <- 1
dt <- T / n 
time <- seq(0, T, by = dt)
S0 <- 1


# Simulate two independent Brownian motions B_t and W_t
set.seed(1)
dB <- rnorm(n, mean = 0, sd = sqrt(dt))
dW <- rnorm(n, mean = 0, sd = sqrt(dt))
W <- c(0, cumsum(dW))
sigma_vol <- abs(W)
#exp_W <- exp(sigma_vol)

# Simulate the process S_t using Euler-Maruyama
S <- numeric(n+1)
S[1] <- S0

for (i in 1:n) {
  S[i + 1] <- S[i] + sigma_vol[i] * S[i] * dB[i]
}

# Plot the results
plot(time[1:100000], S[1:100000], type = "l", col = "blue", lwd = 2, 
     xlab = "Time", ylab = "S_t", main = "Simulated Process S_t")

# Realized volatility
log_returns <- diff(log(S))
squared_log_returns <- log_returns^2
real_vol <- numeric()

for (i in 1:(floor((n+1)/window))){
  real_vol[i] <- sqrt(sum(squared_log_returns[(1+(i-1)*window):(i*window)]))
}
daily_vol <- real_vol*sqrt((n+1)/(window*T)) #sqrt(((n+1)/T)/window)
sigma_vol_win <- sigma_vol[seq(1, n, by = window)]#*sqrt((n+1)/T)#*sqrt(1/T)


# Plot
plot(daily_vol[1:10000], type = "l", col = "black", lwd = 2, 
     xlab = "Time",xaxt = "n", ylab = "volatility",  main = "")
lines(sigma_vol_win[1:10000], col = "red", lwd = 2)

IV_minus_RV <- sigma_vol_win-daily_vol
plot(IV_minus_RV[1:10000], type = "l", col = "black", lwd = 2, 
     xlab = "", xaxt = "n", ylab = "IV-RV", main = "")
log_IV_minus_RV <- log(sigma_vol_win)-log(daily_vol)
plot(log_IV_minus_RV[1:10000], type = "l", col = "black", lwd = 2, 
     xlab = "", xaxt = "n", ylab = "log(IV)-log(RV)", main = "")

X <- sigma_vol_win
X <- daily_vol

objective_function <- function(p) {
  (Wstat(L = L, K = K, p, t=1 , X = X) - 1)^2  # Squared difference to minimize
}
p_inv_sols <- 1/optimize(objective_function, interval = c(1,30))$minimum
p_inv_sols

inv_p_values <- seq(p_inv_sols-0.09, p_inv_sols+0.09, length.out = 1000)  # Avoid p=0 to prevent division by zero
p_values <- 1 / inv_p_values
W_values <- sapply(p_values, function(p) Wstat(L = L, K = K, p, t = 1, X = X))
 
plot(inv_p_values, W_values, log= "y", type = "l", col = "black", lwd = 1,
     main = paste("Estimated H=",round(p_inv_sols, digits = 4)),
     xlab = "H=1/p",
     ylab = expression(W(L, K, pi, p, t == 1, X)), yaxt = "n")  #

y_ticks <- c(10^-0.4, 10^-0.2, 10^0, 10^0.2, 10^0.4, 10^0.6)
y_ticks <- c(10^-0.5, 10^0, 10^0.5, 10^1, 10^1.5)
#y_ticks <- c(10^0, 10^5, 10^10, 10^15)
axis(2, at = y_ticks, labels = expression(10^-0.4, 10^-0.2, 10^0, 10^0.2, 10^0.4, 10^0.6))
axis(2, at = y_ticks, labels = expression(10^-0.5, 10^0, 10^0.5, 10^1, 10^1.5))
#axis(2, at = y_ticks, labels = expression(10^0, 10^5, 10^10, 10^15))
abline(h = 1, col = "blue", lwd = 1, lty = 1)  #Estimating \hat{p}
abline(v = p_inv_sols, col = "blue", lwd = 1, lty = 1)

## W against different values of K
roughness_fct <- function(K){
  objective_function <- function(p) {
    (Wstat(L = L, K = K, p, t=1 , X = X) - 1)^2  # Squared difference to minimize
  }
  return(1/optimize(objective_function, interval = c(1,30))$minimum)
}
K_values <- 10:990 
roughness_index <- sapply(K_values, function(K) roughness_fct(K))  # apply your function

plot(K_values, roughness_index, type = "l",  # 'l' for lines
     xlab = "K", ylab = "Estimated Roughness Index", 
     main = "")
abline(h = roughness_fct(sqrt(L)), col = "blue", lwd = 1, lty = 1) 
abline(v = sqrt(L), col = "blue", lwd = 1, lty = 1)

# Sequential scale estimator
small_n <- 17
capital_N <- small_n + 7
number_of_points <- 2^capital_N # Excluding T=0
T_new <- 1
dt_new <- T_new/number_of_points
time <- seq(0, T_new, by = dt_new)
#Rescale Brownian motions
dB_new <- sqrt(dt_new)*1/sqrt(dt)*dB[1:number_of_points]
dW_new <- sqrt(dt_new)*1/sqrt(dt)*dW[1:number_of_points]

W <- c(0, cumsum(dW_new))
sigma_vol <- abs(W)

# Simulate the process S_t using Euler-Maruyama
S <- numeric(number_of_points+1)
S[1] <- S0

for (i in 1:number_of_points) {
  S[i + 1] <- S[i] + sigma_vol[i] * S[i] * dB_new[i]
}

plot(time[1:100000], S[1:100000], type = "l", col = "blue", lwd = 2, 
     xlab = "Time", ylab = "S_t", main = "Rescaled Process")

log_returns <- diff(log(S))
squared_log_returns <- log_returns^2
cumulative_sum <- cumsum(squared_log_returns)

# Extract every 16th point, including the starting point 0
#y_t <- c(0, cumulative_sum[seq(2^4, length(cumulative_sum), by = 2^(capital_N-small_n-2))]) # Realized variance

k <- 1:(2^(small_n+2))  # Sequence for k
indices <- 2^(capital_N - small_n - 2) * k  # Compute the desired indices
# Extract the elements
result <- cumulative_sum[indices]
RV_scale <- c(0, result)
y_t <- RV_scale

sum_terms <- cumsum((sigma_vol[2:(n+1)])^2)
k <- 1:(2^(small_n+2))  # Sequence for k
indices <- 2^(capital_N - small_n - 2) * k  # Compute the desired indices
# Extract the elements
result <- sum_terms[indices]
IV_scale <- c(0, 2^-capital_N * result)
y_t <- IV_scale

plot(RV_scale[1:10000], type = "l", col = "black", lwd = 2, 
     xlab = "Time",xaxt = "n", ylab = "Integrated variance",  main = "")
lines(IV_scale[1:10000], col = "red", lwd = 2)
IV_minus_RV <- IV_scale-RV_scale
plot(IV_minus_RV[1:10000], type = "l", col = "black", lwd = 2, 
     xlab = "", xaxt = "n", ylab = "Difference", main = "")
log_IV_minus_RV <- log(IV_scale)-log(RV_scale)
plot(log_IV_minus_RV[1:10000], type = "l", col = "black", lwd = 2, 
     xlab = "", xaxt = "n", ylab = "Log difference", main = "")

plot(RV_scale[seq(1, length(RV_scale), by = 1000)], type = "l", col = "black", lwd = 2, 
     xlab = "Time",xaxt = "n", ylab = "Integrated variance",  main = "")
lines(IV_scale[seq(1, length(IV_scale), by = 1000)], col = "red", lwd = 2)
IV_minus_RV <- IV_scale-RV_scale
plot(IV_minus_RV[1:10000], type = "l", col = "black", lwd = 2, 
     xlab = "", xaxt = "n", ylab = "Estimation error", main = "")
log_IV_minus_RV <- log(IV_scale)-log(RV_scale)
plot(log_IV_minus_RV[1:10000], type = "l", col = "black", lwd = 2, 
     xlab = "", xaxt = "n", ylab = "Log estimation error", main = "")

# Estimate roughness from realized variance
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
}

m <- 3
alpha_all <- c(1,1,1,1)
obj_funct <- function(lambda) {
  sum <- 0
  for (k in (small_n - m):small_n) {
    r_k <- rough_expo(k, lambda * y_t, n_max = small_n)
    r_k_minus_1 <- rough_expo(k - 1, lambda * y_t, n_max = small_n)
    sum <- sum + alpha_all[small_n - k + 1] * (r_k - r_k_minus_1)^2
  }
  return(sum)
}

lambda_opt <- optimize(obj_funct, interval = c(0.001,30))$minimum
scale_est <- rough_expo(k = small_n, y_t = lambda_opt*y_t, n_max = small_n)
scale_est
rough_expo(k = small_n, y_t = y_t, n_max = small_n)

roughness_fct <- function(X, threshold = 0.1) {
  objective_function <- function(H_hat) {
    (Wstat(L = L, K = K, 1 / H_hat, t = 1, X = X) - 1)^2
  }
  
  # Perform optimization
  result <- optimize(objective_function, interval = c(0.0001, 0.9999))
  
  # Return roughness estimate if objective is within threshold, else 0
  return(ifelse(result$objective <= threshold, result$minimum, 0))
}
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
  lambda_opt <- optimize(objective_function, interval = c(0.001,30))$minimum
  est <- rough_expo(k = small_n, y_t = lambda_opt*X, n_max = small_n)
  return(est)
}

simtim <- 10
H_IV <- numeric(simtim)
H_RV <- numeric(simtim)
H_IV_scale <- numeric(simtim)
H_RV_scale <- numeric(simtim)
H_IV_expo <- numeric(simtim)
H_RV_expo <- numeric(simtim)
# Precompute constants
sqrt_dt <- sqrt(dt)
sqrt_dt_new <- sqrt(dt_new)
n_plus_1 <- n + 1
window_length <- floor(n_plus_1 / window)
sqrt_term <- sqrt(n_plus_1 / (T * window))
for (j in 1:simtim){
  dB <- rnorm(n, mean = 0, sd = sqrt_dt)
  dW <- rnorm(n, mean = 0, sd = sqrt_dt)
  W <- numeric(n_plus_1)
  W[2:(n+1)] <- cumsum(dW)
  sigma_vol <- abs(W)
  
  # Simulate the process S_t using Euler-Maruyama
  S <- numeric(n_plus_1)
  S[1] <- S0
  
  for (i in 1:n) {
    S[i + 1] <- S[i] + sigma_vol[i] * S[i] * dB[i]
  }
  
  log_returns <- diff(log(abs(S)))
  squared_log_returns <- log_returns^2
  # Realized volatility calculation (vectorized)
  real_vol <- sqrt(rowSums(matrix(squared_log_returns, nrow = window_length, ncol = window, byrow = TRUE)))
  RV <- real_vol * sqrt_term
  IV <- sigma_vol[seq(1, n, by = window)]
  
  H_IV[j] <- roughness_fct(IV)
  H_RV[j] <- roughness_fct(RV)
  
  #Rescale Brownian motions
  dB_new <- sqrt_dt_new*1/sqrt_dt*dB[1:number_of_points]
  dW_new <- sqrt_dt_new*1/sqrt_dt*dW[1:number_of_points]
  
  W <- c(0, cumsum(dW_new))
  sigma_vol <- abs(W)
  S <- numeric(number_of_points+1)
  S[1] <- S0
  
  for (i in 1:number_of_points) {
    S[i + 1] <- S[i] + sigma_vol[i] * S[i] * dB_new[i]
  }
  
  log_returns <- diff(log(S))
  squared_log_returns <- log_returns^2
  cumulative_sum <- cumsum(squared_log_returns)
  
  k <- 1:(2^(small_n+2))  # Sequence for k
  indices <- 2^(capital_N - small_n - 2) * k  # Compute the desired indices
  # Extract the elements
  result <- cumulative_sum[indices]
  RV_scale <- c(0, result)

  H_RV_scale[j] <- scale_est_fct(RV_scale)
  H_RV_expo[j] <- rough_expo(k = small_n, y_t = RV_scale, n_max = small_n)  
  
  sum_terms <- cumsum((sigma_vol[2:(n+1)])^2)
  result <- sum_terms[indices]
  IV_scale <- c(0, 2^-capital_N * result)
  
  H_IV_scale[j] <- scale_est_fct(IV_scale)
  H_IV_expo[j] <- rough_expo(k = small_n, y_t = IV_scale, n_max = small_n)   

  
}

plot(density(H_IV), 
     main = "",
     xlab = "Estimated roughness",
     ylab = "Density",
     xlim = c(min(c(H_IV,H_IV_scale, H_IV_expo)),max(c(H_IV,H_IV_scale, H_IV_expo))),
     col = "red", 
     lwd = 2)  # Density plot of roughness estimate
lines(density(H_IV_scale), col = 'blue', lwd =2)
lines(density(H_IV_expo), col = 'green', lwd =2)

plot(density(H_RV), 
     main = "",
     xlab = "Estimated roughness",
     ylab = "Density",
     xlim = c(min(c(H_RV,H_RV_scale, H_RV_expo)),max(c(H_RV,H_RV_scale, H_RV_expo))),
     col = "red", 
     lwd = 2)  # Density plot of roughness estimate
lines(density(H_RV_scale), col = 'blue', lwd =2)
lines(density(H_RV_expo), col = 'green', lwd =2)


### Example 6 OU-SV model
# Parameters
L = 300*300
K = 300
#n <- L+1
window <- 300
n <- window*(L+1)
#n <- 10000
T <- 5 # Use T=5 for this example
dt <- T / n 
time <- seq(0, T, by = dt)

sigma_0 <- 1
Y_0 <- 0
gamma <- 1
theta <- 1
S0 <- 1

set.seed(12)
# Simulate two independent Brownian motions
dB_1 <- rnorm(n, mean = 0, sd = sqrt(dt))
dB_2 <- rnorm(n, mean = 0, sd = sqrt(dt))

Y <- numeric()
Y[1] <- Y_0
for (i in 1:n){
  Y[i+1] <- Y[i] -gamma*Y[i]*dt + theta*dB_2[i]
}

sigma <- numeric()
sigma[1] <- sigma_0
for (i in 2:(n+1)){
  sigma[i] <- sigma_0*exp(Y[i])
}

S_or <- numeric()
S_or[1] <- S0
for (i in 1:n){
  S_or[i+1] <- S_or[i]+S_or[i]*sigma[i]*dB_1[i]
}
min(S_or)

plot(time[1:10000], S_or[1:10000], type = "l", col = "blue", lwd = 2, 
     xlab = "Time", ylab = "S_t", main = "Simulated Process S_t")

# Realized volatility
log_returns <- diff(log(S_or))
squared_log_returns <- log_returns^2
real_vol <- numeric()

for (i in 1:(floor((n+1)/window))){
  real_vol[i] <- sqrt(sum(squared_log_returns[(1+(i-1)*window):(i*window)]))
}
daily_vol <- real_vol*sqrt(((n+1)/T)/window)
sigma_vol_win <- sigma[seq(1, n, by = window)]

plot(daily_vol[1:10000], type = "l", col = "black", lwd = 2, 
     xlab = "Time",xaxt = "n", ylab = "volatility",  main = "")
lines(sigma_vol_win[1:10000], col = "red", lwd = 2)

IV_minus_RV <- sigma_vol_win[1:10000]-daily_vol[1:10000]
plot(IV_minus_RV[1:10000], type = "l", col = "black", lwd = 2, 
     xlab = "", xaxt = "n", ylab = "IV-RV", main = "")
log_IV_minus_RV <- log(sigma_vol_win[1:10000])-log(daily_vol[1:10000])
plot(log_IV_minus_RV[1:10000], type = "l", col = "black", lwd = 2, 
     xlab = "", xaxt = "n", ylab = "log(IV)-log(RV)", main = "")

X <- sigma_vol_win
X <- daily_vol
objective_function <- function(p) {
  (Wstat(L = L, K = K, p, t=T , X = X) - T)^2  # Squared difference to minimize
}
p_inv_sols <- 1/optimize(objective_function, interval = c(1,30))$minimum
p_inv_sols

Delta_values <- 1:50
log_Delta <- log(Delta_values)
q_values <- c(0.5, 1, 1.5, 2, 3)  # Different q values
# Calculate log(m_q(Δ)) for each q and Δ
log_mq_list <- lapply(q_values, function(q) {
  sapply(Delta_values, function(delta) mqdelta(delta, q, X))  # Calculate log_mq for each delta
})
colors <- c("black", "green", "red", "blue", "purple")
# Plot the first q value
plot(log_Delta, log_mq_list[[1]], type = "p", col = colors[1], pch = 21,
     xlab = expression(log(Delta)), 
     ylab = expression(log(m(q,Delta))),
     ylim = c(min(sapply(log_mq_list, min)), max(sapply(log_mq_list, max))),
     xlim = c(min(log_Delta), max(log_Delta))
)
# Add lines for the remaining q values
for (i in 2:length(q_values)) {
  points(log_Delta, log_mq_list[[i]], col = colors[i], pch = 21)  # Add points
}
# Add legend
legend("bottomright", legend = paste("q =", q_values), col = colors, pch = 21, cex = 0.7)
# Data
log_delta <- c(log_Delta, log_Delta, log_Delta, log_Delta, log_Delta)  # Vector of log(Δ) values
log_m <- c(log_mq_list[[1]], log_mq_list[[2]], log_mq_list[[3]], log_mq_list[[4]], log_mq_list[[5]])      # Vector of log(m(q, Δ)) values
group <- c(rep(1,length(Delta_values)), rep(2,length(Delta_values)), rep(3,length(Delta_values)), rep(4,length(Delta_values)), rep(5,length(Delta_values)))      # Vector indicating the group (1 to 5 for each set)
# Combine into a data frame
data <- data.frame(log_delta, log_m, group)
# Split data by group
data_split <- split(data, data$group)
# Perform linear regression for each group
models <- lapply(data_split, function(group_data) {
  lm(log_m ~ log_delta, data = group_data)
})
# Add regression lines to plot
for (i in 1:length(q_values)){
  abline(models[[i]], col = colors[i])
}
epsilon_q <- numeric()
for (i in 1:length(q_values)){
  epsilon_q[i] <- models[[i]]$coefficients[[2]]
}
fit_epsilon <- lm(epsilon_q ~ q_values)
est_smooth_H <- fit_epsilon$coefficients[[2]]
plot(q_values, epsilon_q, type = "p", pch = 21, main = paste("Estimated H=",round(est_smooth_H, digits = 4)), 
     xlab = 'q', ylab = expression(zeta[q]), xlim = c(0,max(q_values)), ylim = c(0,max(epsilon_q)))
abline(fit_epsilon)

# Sequential scale estimator
small_n <- 16
capital_N <- small_n + 6
number_of_points <- 2^capital_N # Excluding T=0
T <- 1
dt_new <- T/number_of_points
#Rescale Brownian motions
dB1_new <- sqrt(dt_new)*1/sqrt(dt)*dB_1[1:number_of_points]
dB2_new <- sqrt(dt_new)*1/sqrt(dt)*dB_2[1:number_of_points]

Y <- numeric()
Y[1] <- Y_0
for (i in 1:number_of_points){
  Y[i+1] <- Y[i] -gamma*Y[i]*dt_new + theta*dB2_new[i]
}

sigma <- numeric()
sigma[1] <- sigma_0
for (i in 2:(number_of_points+1)){
  sigma[i] <- sigma_0*exp(Y[i])
}

S_or <- numeric()
S_or[1] <- S0
for (i in 1:number_of_points){
  S_or[i+1] <- S_or[i]+S_or[i]*sigma[i]*dB1_new[i]
}
min(S_or)

plot(time[1:10000], S_or[1:10000], type = "l", col = "blue", lwd = 2, 
     xlab = "Time", ylab = "S_t", main = "Simulated Process S_t")

# Realized volatility
log_returns <- diff(log(S_or))
squared_log_returns <- log_returns^2
cumulative_sum <- cumsum(squared_log_returns)
# Extract every 16th point, including the starting point 0
y_t <- c(0, cumulative_sum[seq(2^4, length(cumulative_sum), by = 2^(capital_N-small_n-2))]) # Realized variance

k <- 1:(2^(small_n+2))  # Sequence for k
indices <- 2^4 * k  # Compute the desired indices
# Extract the elements
result <- cumulative_sum[indices]
RV_scale <- c(0, result)
y_t <- RV_scale

sum_terms <- cumsum(exp(2*Y[2:(n+1)]))
k <- 1:(2^(small_n+2))  # Sequence for k
indices <- 2^4 * k  # Compute the desired indices
# Extract the elements
result <- sum_terms[indices]
IV_scale <- c(0, 2^-capital_N * result)
y_t <- IV_scale

plot(RV_scale[1:10000], type = "l", col = "black", lwd = 2, 
     xlab = "Time",xaxt = "n", ylab = "volatility",  main = "")
lines(IV_scale[1:10000], col = "red", lwd = 2)

plot(RV_scale[seq(1, length(RV_scale), by = 100)], type = "l", col = "black", lwd = 2, 
     xlab = "Time",xaxt = "n", ylab = "volatility",  main = "")
lines(IV_scale[seq(1, length(IV_scale), by = 100)], col = "red", lwd = 2)


lambda_opt <- optimize(obj_funct, interval = c(0.001,30))$minimum
scale_est <- rough_expo(k = small_n, y_t = lambda_opt*y_t, n_max = small_n)
scale_est
rough_expo(k = small_n, y_t = y_t, n_max = small_n)

## Density plot
simtim <- 2500
RV_check <- matrix(0,simtim, L+1)
#IV_check <- matrix(0,simtim, L+1)
IV_new <- matrix(0,simtim, L+1)
# Precompute constants
sqrt_dt <- sqrt(dt)
n_plus_1 <- n + 1
window_length <- floor(n_plus_1 / window)
sqrt_term <- sqrt(n_plus_1 / (T * window))
for (k in 1:simtim){
  dB_1 <- rnorm(n, mean = 0, sd = sqrt_dt)
  dB_2 <- rnorm(n, mean = 0, sd = sqrt_dt)
  
  Y <- numeric(n_plus_1)
  Y[1] <- Y_0
  for (i in 1:n){
    Y[i+1] <- Y[i] -gamma*Y[i]*dt + theta*dB_2[i]
  }
  
  sigma <- sigma_0 * exp(Y)
  
  S <- numeric(n_plus_1)
  S[1] <- 100
  for (i in 1:n){
    S[i+1] <- S[i]+S[i]*sigma[i]*dB_1[i]
  }
  # Realized volatility
  log_returns <- diff(log(abs(S)))
  squared_log_returns <- log_returns^2
  
  # Realized volatility calculation (vectorized)
  real_vol <- sqrt(rowSums(matrix(squared_log_returns, nrow = window_length, ncol = window, byrow = TRUE)))
  daily_vol <- real_vol * sqrt_term
  
  # Average sigma in each window (vectorized)
  #sigma_vol_win <- rowMeans(matrix(sigma[2:n_plus_1], nrow = window_length, ncol = window, byrow = TRUE))
  sigma_vol_new <- sigma[seq(1, n, by = window)]
  # Store results in RV_check and IV_check
  RV_check[k, ] <- daily_vol
  #IV_check[k, ] <- sigma_vol_win
  IV_new[k, ] <- sigma_vol_new
}

X = RV_check
#X = IV_check
X = IV_new
p_inv_sols <- numeric(simtim)
for (i in 1:simtim){
  objective_function <- function(p) {
    (Wstat(L = L, K = K, p, t=T , X = X[i,]) - T)^2  # Squared difference to minimize
  }
  p_inv_sols[i] <- 1/optimize(objective_function, interval = c(1,30))$minimum
}
plot(density(p_inv_sols), 
     main = "",
     xlab = "H=1/p",
     ylab = "Density",
     col = "red", 
     lwd = 2)

quartiles <- as.numeric(quantile(p_inv_sols, probs = c(0.25, 0.75)))
#Only run one of the lines below
#table_dens <- matrix(0,2,6)
table_dens[1,] <- c(min(p_inv_sols),quartiles[1], median(p_inv_sols), mean(p_inv_sols), quartiles[2], max(p_inv_sols))
table_dens[2,] <- c(min(p_inv_sols),quartiles[1], median(p_inv_sols), mean(p_inv_sols), quartiles[2], max(p_inv_sols))
table_dens


### A fractional Ornstein-Uhlenbeck model (Example 7)
L = 30*30
K = 30
window <- 300
n <- window*(L+1)
l <- ifelse(n*3>=30000, n*3, 30000)
T <- 5
dt <- T / n 
#time <- seq(0, T, by = dt)

set.seed(1)
U_0 <- rnorm(l)  # U^(0) for k = 0, ..., l-1
U_1 <- rnorm(l)  # U^(1) for k = 0, ..., l-1

t_k <- 0:l
for (k in 0:l){
  t_k[k+1] <- pi * k/l
}

indices <- 1:(l - 1)
ak_constant <- 0.5 * (U_0[indices] + 1i * U_1[indices])
indices_2 <- (l+1):(2*l-1)
ak_constant_2 <- 0.5 * (U_0[2 * l - indices_2] - 1i * U_1[2 * l - indices_2])
indi_plus_1 <- indices + 1
indi2_plus_1 <- indices_2+1
indi_2l <- 2*l-indices_2+1
l_plus_one <- l +1 
two_l <- 2 * l

sigma_0 <- 1
Y_0 <- 0
gamma <- 1
theta <- 1
S_0 <- 1

sqrt_dt <- sqrt(dt)
n_plus_1 <- n + 1
window_length <- floor(n_plus_1 / window)
sqrt_term <- sqrt(n_plus_1 / (T * window))
dB <- rnorm(n, mean = 0, sd = sqrt_dt)
Y_decay_factor <- 1 - gamma * dt

roughness_fct <- function(X, threshold = 0.1) {
  objective_function <- function(H_hat) {
    (Wstat(L = L, K = K, 1 / H_hat, t = 1, X = X) - 1)^2
  }
  
  # Perform optimization
  result <- optimize(objective_function, interval = c(0.0001, 0.9999))
  
  # Return roughness estimate if objective is within threshold, else 0
  return(ifelse(result$objective <= threshold, result$minimum, 0))
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

H <- H_values[7]
H
d_fbm <- Re(fft(compute_a_k(l, H))[1:n])*(dt)^H

Y <- numeric(n_plus_1)
Y[1] <- Y_0
for (i in seq_len(n)) {
  Y[i + 1] <- Y[i] * Y_decay_factor + theta * d_fbm[i]
}
sigma <- sigma_0 * exp(Y)
S <- numeric(n_plus_1)
S[1] <- S_0
for (i in seq_len(n)){
  S[i+1] <- S[i]*(1+sigma[i]*dB[i])
}
# Realized volatility
log_returns <- diff(log(abs(S)))
squared_log_returns <- log_returns^2
# Realized volatility calculation (vectorized)
real_vol <- sqrt(rowSums(matrix(squared_log_returns, nrow = window_length, ncol = window, byrow = TRUE)))
RV <- real_vol * sqrt_term
IV <- sigma[seq(1, n, by = window)]

plot(S[seq(1, length(S), by = 100)], type = "l", ylab = "Price S_t", xlab = "", xaxt = "n")
plot(S[1:10000]-S_0, type = "l", ylab = "Price S_t", xlab = "", xaxt = "n")
plot(RV[1:10000], type = "l", ylab = "Realized volatility", xlab = "", xaxt = "n")
plot(IV[1:10000], type = "l", ylab = "Instantaneous volatility", xlab = "", xaxt = "n")  

H_IV <- roughness_fct(IV)
H_RV <- roughness_fct(RV)
H_IV_lres <- smoothness_fct(IV)
H_RV_lres <- smoothness_fct(RV)
print(c(H_IV, H_IV_lres, H_RV, H_RV_lres))

# Sequential scale estimator
small_n <- 12
capital_N <- small_n + 6
number_of_points <- 2^capital_N # Excluding T=0
T <- 1
dt_new <- T/number_of_points
Y_decay_factor_new <- 1 - gamma * dt_new

dfbm_new <- (dt_new/dt)^H*d_fbm[1:number_of_points]
dB_new <- sqrt(dt_new)*1/sqrt(dt)*dB[1:number_of_points]

Y <- numeric(number_of_points+1)
Y[1] <- Y_0
for (i in seq_len(number_of_points)) {
  Y[i + 1] <- Y[i] * Y_decay_factor_new + theta * dfbm_new[i]
}
sigma <- sigma_0 * exp(Y)
S <- numeric(number_of_points+1)
S[1] <- S_0
for (i in seq_len(number_of_points)){
  S[i+1] <- S[i]*(1+sigma[i]*dB_new[i])
}
# Realized volatility
log_returns <- diff(log(abs(S)))
squared_log_returns <- log_returns^2
cumulative_sum <- cumsum(squared_log_returns)
#y_t <- c(0, cumulative_sum[seq(2^4, length(cumulative_sum), by = 2^(capital_N-small_n-2))]) # Realized variance

k <- 1:(2^(small_n+2))  # Sequence for k
indices <- 2^4 * k  # Compute the desired indices
# Extract the elements
result <- cumulative_sum[indices]
RV_scale <- c(0, result)
y_t <- RV_scale

sum_terms <- cumsum(exp(2*Y[2:(n+1)]))
k <- 1:(2^(small_n+2))  # Sequence for k
indices <- 2^4 * k  # Compute the desired indices
# Extract the elements
result <- sum_terms[indices]
IV_scale <- c(0, 2^-capital_N * result)
y_t <- IV_scale

plot(RV_scale[1:10000], type = "l", col = "black", lwd = 2, 
     xlab = "Time",xaxt = "n", ylab = "volatility",  main = "")
lines(IV_scale[1:10000], col = "red", lwd = 2)

plot(RV_scale[seq(1, length(RV_scale), by = 100)], type = "l", col = "black", lwd = 2, 
     xlab = "Time",xaxt = "n", ylab = "volatility",  main = "")
lines(IV_scale[seq(1, length(IV_scale), by = 100)], col = "red", lwd = 2)

lambda_opt <- optimize(obj_funct, interval = c(0.001,30))$minimum
scale_est <- rough_expo(k = small_n, y_t = lambda_opt*y_t, n_max = small_n)
scale_est
rough_expo(k = small_n, y_t = y_t, n_max = small_n)

output_given_H <- function(H){
  
  d_fbm <- Re(fft(compute_a_k(l, H))[1:n])*(dt)^H

  Y <- numeric(n_plus_1)
  Y[1] <- Y_0
  for (i in seq_len(n)) {
    Y[i + 1] <- Y[i] * Y_decay_factor + theta * d_fbm[i]
  }
  sigma <- sigma_0 * exp(Y)
  S <- numeric(n_plus_1)
  S[1] <- S_0
  for (i in seq_len(n)){
    S[i+1] <- S[i]*(1+sigma[i]*dB[i])
  }
  # Realized volatility
  log_returns <- diff(log(abs(S)))
  squared_log_returns <- log_returns^2
  # Realized volatility calculation (vectorized)
  real_vol <- sqrt(rowSums(matrix(squared_log_returns, nrow = window_length, ncol = window, byrow = TRUE)))
  RV <- real_vol * sqrt_term
  IV <- sigma[seq(1, n, by = window)]
  
  plot(S[1:10000]-S_0, type = "l", ylab = "Price S_t", xlab = "", xaxt = "n")
  plot(RV[1:10000], type = "l", ylab = "Realized volatility", xlab = "", xaxt = "n")
  plot(IV[1:10000], type = "l", ylab = "Instantaneous volatility", xlab = "", xaxt = "n")  
  
  H_IV <- roughness_fct(IV)
  H_RV <- roughness_fct(RV)
  H_IV_lres <- smoothness_fct(IV)
  H_RV_lres <- smoothness_fct(RV)

  
  return(list(H_IV = H_IV, H_RV = H_RV, H_IV_lres = H_IV_lres, H_RV_lres = H_RV_lres))
}

H_values <- seq(0.10, 0.80, by = 0.10)
f_tk_mat <- matrix(0, nrow = length(H_values), ncol = l_plus_one)
# Loop over H values with indexing
for (i in seq_along(H_values)) {
  H <- H_values[i]
  f_tk_mat[i, ] <- sapply(t_k[1:l_plus_one], f_lambda, H = H)
}

store_results <- list()
for (H in H_values) {
  store_results[[paste0("H_", H)]] <- output_given_H(H)
}
results <- data.frame(H = numeric(), Instantaneous_volatility = numeric(), Realized_volatility = numeric(), IV_logres = numeric(), RV_logres = numeric(),
                      scale_vol = numeric())
for (H in H_values) {
  IV_H_est <- store_results[[paste0("H_", H)]]$H_IV
  RV_H_est <- store_results[[paste0("H_", H)]]$H_RV
  IV_H_est_logres <- store_results[[paste0("H_", H)]]$H_IV_lres
  RV_H_est_logres <- store_results[[paste0("H_", H)]]$H_RV_lres
  RV_H_est_scale <- store_results[[paste0("H_", H)]]$scale_est
  # Append the results to the data frame
  results <- rbind(results, data.frame(H = H, 
                                       Instantaneous_volatility = IV_H_est, 
                                       Realized_volatility = RV_H_est, 
                                       IV_logres = IV_H_est_logres, 
                                       RV_logres = RV_H_est_logres,
                                       Scale_est = RV_H_est_scale))
}
# Display the final table using knitr::kable
kable(results, col.names = c("H", "Instantaneous volatility", "Realized volatility", "IV (log-res)", "RV (log-res)", "RV (Scale)"), 
      align = "c", caption = "Volatility Table")

plot(results$H, results$Instantaneous_volatility, type = "l", col = "blue", xlab = "", 
     ylab = "", xlim = c(0,0.8), ylim = c(0,0.8))
abline(v = seq(0, 0.8, by = 0.1), col = "gray", lty = "dotted", lwd = 0.5) # Vertical lines
abline(h = seq(0, 0.8, by = 0.1), col = "gray", lty = "dotted", lwd = 0.5) # Horizontal lines
lines(results$H, results$Realized_volatility, col = "red", lwd = 2)
lines(results$H, results$IV_logres, col = "blue", lwd = 2, lty = 2)
lines(results$H, results$RV_logres, col = "red", lwd = 2, lty = 2)

legend("topleft", legend = c('Realized Vol', 'Instantaneous Vol'), col = c('red', 'blue'), lty = 1, cex = 0.7)



## Simulate 100 times

simoutput_given_H <- function(H){
  
  d_fbm <- Re(fft(compute_a_k(l, H))[1:n])*(dt)^H
  
  Y <- numeric(n_plus_1)
  Y[1] <- Y_0
  for (i in 1:n){
    Y[i+1] <- Y[i] -gamma*Y[i]*dt + theta*d_fbm[i]
  }
  sigma <- sigma_0 * exp(Y)
  S <- numeric(n_plus_1)
  S[1] <- S_0
  for (i in 1:n){
    S[i+1] <- S[i]+S[i]*sigma[i]*dB[i]
  }
  # Realized volatility
  log_returns <- diff(log(abs(S)))
  squared_log_returns <- log_returns^2
  # Realized volatility calculation (vectorized)
  real_vol <- sqrt(rowSums(matrix(squared_log_returns, nrow = window_length, ncol = window, byrow = TRUE)))
  RV <- real_vol * sqrt_term
  IV <- sigma[seq(1, n, by = window)]
  
  IV_H <- roughness_fct(IV)
  RV_H <- roughness_fct(RV)
  
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

simtim <- 100
simulation_results <- vector("list", simtim)
for (i in 1:simtim){
  U_0 <- rnorm(l)  # U^(0) for k = 0, ..., l-1
  U_1 <- rnorm(l)  # U^(1) for k = 0, ..., l-1
  ak_constant <- 0.5 * (U_0[indices] + 1i * U_1[indices])
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


# Initialize matrices to store volatilities across runs
realized_vol_matrix <- matrix(NA, nrow = length(H_values), ncol = simtim)
instantaneous_vol_matrix <- matrix(NA, nrow = length(H_values), ncol = simtim)
# Fill the matrices with realized and instantaneous volatilities for each run
for (i in seq_along(simulation_results)) {
  realized_vol_matrix[, i] <- simulation_results[[i]]$Realized.volatility
  instantaneous_vol_matrix[, i] <- simulation_results[[i]]$Instantaneous.volatility
}

# Calculate the mean, lower bound, and upper bound for realized and instantaneous volatilities
mean_realized_vol <- rowMeans(realized_vol_matrix)
lower_bound_realized_vol <- apply(realized_vol_matrix, 1, quantile, probs = 0.125)
upper_bound_realized_vol <- apply(realized_vol_matrix, 1, quantile, probs = 0.875)

mean_instantaneous_vol <- rowMeans(instantaneous_vol_matrix)
lower_bound_instantaneous_vol <- apply(instantaneous_vol_matrix, 1, quantile, probs = 0.125)
upper_bound_instantaneous_vol <- apply(instantaneous_vol_matrix, 1, quantile, probs = 0.875)

# Plotting
lines(H_values, mean_realized_vol, type = "l", col = "black", lwd = 2)
lines(H_values, lower_bound_realized_vol, col = "black", lty = 2)
lines(H_values, upper_bound_realized_vol, col = "black", lty = 2)

lines(H_values, mean_instantaneous_vol, col = "black", lwd = 2)
lines(H_values, lower_bound_instantaneous_vol, col = "black", lty = 2)
lines(H_values, upper_bound_instantaneous_vol, col = "black", lty = 2)



