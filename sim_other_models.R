library(MASS)
library(tidyverse)
library(complex)
## Functions used throughout the document
set.seed(1)
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

# Parameters
L = 500*500
K = 500
n <- L+1
window <- 300
n <- window*(L+1)
T <- 10
dt <- T / n 
time <- seq(0, T, by = dt)

# Simulate two independent Brownian motions B_t and W_t
dB <- rnorm(n, mean = 0, sd = sqrt(dt))
B <- numeric(n+1)
B[2:(n+1)] <- cumsum(dB)
dW <- rnorm(n, mean = 0, sd = sqrt(dt))
W <- numeric(n+1)
W[2:(n+1)] <- cumsum(dW)


# Simulate the process S_t using Euler-Maruyama
S <- numeric(n+1)
S[1] <- 100

sigma_vol <- abs(W)
for (i in 1:n) {
  S[i + 1] <- S[i] + sigma_vol[i] * S[i] * dB[i]
}

# Plot the results
plot(time[1:(20000*window)], S[1:(20000*window)], type = "l", col = "blue", lwd = 2, 
     xlab = "Time", ylab = "S_t", main = "Simulated Process S_t")

# Realized volatility
log_returns <- diff(log(S))
squared_log_returns <- log_returns^2
real_vol <- numeric(1)
window <- 300

#for (i in 1:window){
#  real_vol[i] <- sqrt(sum(squared_log_returns[1:i]))
#}
#for (i in window:(n)){
#  real_vol[i] <- sqrt(sum(squared_log_returns[(i-window+1):(i)]))
#}
#daily_vol <- real_vol*sqrt(n/window)
#sigma_vol <- sigma_vol[1:n]

for (i in 1:(floor((n+1)/window))){
  real_vol[i] <- sqrt(sum(squared_log_returns[(1+(i-1)*window):(i*window)]))
}
daily_vol <- real_vol*sqrt(((n+1)/T)/window)
sigma_vol_win <- numeric(1)
for (i in 1:(floor((n+1)/window))){
  sigma_vol_win[i] <- mean(sigma_vol[(1+(i-1)*window):(i*window)])
}
max(sigma_vol_win)
max(daily_vol)
min(S)

plot(daily_vol[10000:20000], type = "l", col = "black", lwd = 2, 
     xlab = "Time",xaxt = "n", ylab = "volatility",  main = "")
lines(sigma_vol_win[10000:20000], col = "red", lwd = 2)

IV_minus_RV <- sigma_vol_win-daily_vol
plot(IV_minus_RV, type = "l", col = "black", lwd = 2, 
     xlab = "", xaxt = "n", ylab = "IV-RV", main = "")
log_IV_minus_RV <- log(sigma_vol_win)-log(daily_vol)
plot(log_IV_minus_RV, type = "l", col = "black", lwd = 2, 
     xlab = "", xaxt = "n", ylab = "log(IV)-log(RV)", main = "")

#X <- sigma_vol_win
X <- daily_vol
objective_function <- function(p) {
  (Wstat(L = L, K = K, p, t=1 , X = X) - T)^2  # Squared difference to minimize
}
p_inv_sols <- 1/optimize(objective_function, interval = c(1,20))$minimum
p_inv_sols

inv_p_values <- seq(p_inv_sols-0.09, p_inv_sols+0.09, length.out = 1000)  # Avoid p=0 to prevent division by zero
p_values <- 1 / inv_p_values

W_values <- sapply(p_values, function(p) Wstat(L = L, K = K, p, t = 1, X = X))
 
plot(inv_p_values, W_values, log= "y", type = "l", col = "black", lwd = 1,
     main = "",
     xlab = "1/p",
     ylab = "W", yaxt = "n")  #

#y_ticks <- c(10^-0.4, 10^-0.2, 10^0, 10^0.2, 10^0.4, 10^0.6)
y_ticks <- c(10^-0.5, 10^0, 10^0.5, 10^1, 10^1.5)
#axis(2, at = y_ticks, labels = expression(10^-0.4, 10^-0.2, 10^0, 10^0.2, 10^0.4, 10^0.6))
axis(2, at = y_ticks, labels = expression(10^-0.5, 10^0, 10^0.5, 10^1, 10^1.5))
abline(h = 1, col = "blue", lwd = 1, lty = 1)  #Estimating \hat{p}
abline(v = p_inv_sols, col = "blue", lwd = 1, lty = 1)

## W against different values of K
roughness_fct <- function(K){
  objective_function <- function(p) {
    (Wstat(L = L, K = K, p, t=1 , X = X, tol = 10^-3) - T)^2  # Squared difference to minimize
  }
  return(1/optimize(objective_function, interval = c(1,15), tol = 10^-3)$minimum)
}
K_values <- 2:500 
roughness_index <- sapply(K_values, function(K) roughness_fct(K))  # apply your function

plot(K_values, roughness_index, type = "l",  # 'l' for lines
     xlab = "K", ylab = "Estimated Roughness Index", 
     main = "")
abline(h = roughness_fct(sqrt(L)), col = "blue", lwd = 1, lty = 1) 
abline(v = sqrt(L), col = "blue", lwd = 1, lty = 1)

#Scaling analysis


### Example 6 OU-SV model
# Parameters
L = 300*300
K = 300
n <- L+1
window <- 300
n <- window*(L+1)
T <- 10
dt <- T / n 
time <- seq(0, T, by = dt)

sigma_0 <- 1
Y_0 <- 0
gamma <- 1
theta <- 1

# Simulate two independent Brownian motions B_t and W_t
dB_1 <- rnorm(n, mean = 0, sd = sqrt(dt))
B_1 <- numeric(n+1)
B_1[2:(n+1)] <- cumsum(dB_1)
dB_2 <- rnorm(n, mean = 0, sd = sqrt(dt))
B_2 <- numeric(n+1)
B_2[2:(n+1)] <- cumsum(dB_2)

Y <- numeric(1)
Y[1] <- Y_0
for (i in 1:n){
  Y[i+1] <- Y[i] -gamma*Y[i]*dt + theta*dB_2[i]
}

sigma <- numeric(1)
sigma[1] <- sigma_0
for (i in 2:(n+1)){
  sigma[i] <- sigma_0*exp(Y[i])
}

S <- numeric(1)
S[1] <- 100
for (i in 1:n){
  S[i+1] <- S[i]+S[i]*sigma[i]*dB_1[i]
}

# Realized volatility
log_returns <- diff(log(S))
squared_log_returns <- log_returns^2
real_vol <- numeric(1)

for (i in 1:(floor((n+1)/window))){
  real_vol[i] <- sqrt(sum(squared_log_returns[(1+(i-1)*window):(i*window)]))
}
daily_vol <- real_vol*sqrt(((n+1)/T)/window)
sigma_vol_win <- numeric(1)
for (i in 1:(floor((n+1)/window))){
  sigma_vol_win[i] <- mean(sigma[(1+(i-1)*window):(i*window)])
}
max(sigma_vol_win)
max(daily_vol)
max(S)
plot(daily_vol, type = "l", col = "black", lwd = 2, 
     xlab = "Time",xaxt = "n", ylab = "volatility",  main = "")
lines(sigma_vol_win, col = "red", lwd = 2)

IV_minus_RV <- sigma_vol_win-daily_vol
plot(IV_minus_RV, type = "l", col = "black", lwd = 2, 
     xlab = "", xaxt = "n", ylab = "IV-RV", main = "")
log_IV_minus_RV <- log(sigma_vol_win)-log(daily_vol)
plot(log_IV_minus_RV, type = "l", col = "black", lwd = 2, 
     xlab = "", xaxt = "n", ylab = "log(IV)-log(RV)", main = "")

plot(time[1:10000], S[1:10000], type = "l", col = "blue", lwd = 2, 
    xlab = "Time", ylab = "S_t", main = "Simulated Process S_t")

X <- sigma_vol_win
X <- daily_vol
objective_function <- function(p) {
  (Wstat(L = L, K = K, p, t=T , X = X) - T)^2  # Squared difference to minimize
}
p_inv_sols <- 1/optimize(objective_function, interval = c(1,20))$minimum
p_inv_sols

simtim <- 100
RV_check <- matrix(0,simtim, L+1)
IV_check <- matrix(0,simtim, L+1)
for (k in 1:simtim){
  dB_1 <- rnorm(n, mean = 0, sd = sqrt(dt))
  B_1 <- numeric(n+1)
  B_1[2:(n+1)] <- cumsum(dB_1)
  dB_2 <- rnorm(n, mean = 0, sd = sqrt(dt))
  B_2 <- numeric(n+1)
  B_2[2:(n+1)] <- cumsum(dB_2)
  
  Y <- numeric(1)
  Y[1] <- Y_0
  for (i in 1:n){
    Y[i+1] <- Y[i] -gamma*Y[i]*dt + theta*dB_2[i]
  }
  
  sigma <- numeric(1)
  sigma[1] <- sigma_0
  for (i in 2:(n+1)){
    sigma[i] <- sigma_0*exp(Y[i])
  }
  
  S <- numeric(1)
  S[1] <- 100
  for (i in 1:n){
    S[i+1] <- S[i]+S[i]*sigma[i]*dB_1[i]
  }
  # Realized volatility
  log_returns <- diff(log(S))
  squared_log_returns <- log_returns^2
  real_vol <- numeric(1)
  for (i in 1:(floor((n+1)/window))){
    real_vol[i] <- sqrt(sum(squared_log_returns[(1+(i-1)*window):(i*window)]))
  }
  daily_vol <- real_vol*sqrt(((n+1)/T)/window)
  sigma_vol_win <- numeric(1)
  for (i in 1:(floor((n+1)/window))){
    sigma_vol_win[i] <- mean(sigma[(1+(i-1)*window):(i*window)])
  }
  RV_check[k,] <- daily_vol
  IV_check[k,] <- sigma_vol_win
}

simtim <- 10
RV_check <- matrix(0, simtim, L+1)
IV_check <- matrix(0, simtim, L+1)
# Precompute constants
sqrt_dt <- sqrt(dt)
n_plus_1 <- n + 1
window_length <- floor(n_plus_1 / window)
sqrt_term <- sqrt(n_plus_1 / (T * window))

for (k in 1:simtim) {
  # Generate Brownian motions
  dB_1 <- rnorm(n, mean = 0, sd = sqrt_dt)
  dB_2 <- rnorm(n, mean = 0, sd = sqrt_dt)
  
  # Cumulative sum of dB_1 and dB_2 for B_1 and B_2
  B_1 <- c(0, cumsum(dB_1))
  B_2 <- c(0, cumsum(dB_2))
  
  # Simulate Y using vectorized calculation
  Y <- numeric(n_plus_1)
  Y[1] <- Y_0
  Y[2:n_plus_1] <- Y_0 - gamma * cumsum(Y[1:n] * dt) + theta * B_2[2:n_plus_1]
  
  # Compute sigma using Y (vectorized)
  sigma <- sigma_0 * exp(Y)
  
  # Simulate S using vectorized calculation
  S <- numeric(n_plus_1)
  S[1] <- 100
  S[2:n_plus_1] <- S[1] * cumprod(1 + sigma[2:n_plus_1] * dB_1)
  
  # Calculate log returns and squared log returns (vectorized)
  log_returns <- diff(log(S))
  squared_log_returns <- log_returns^2
  
  # Realized volatility calculation (vectorized)
  real_vol <- sqrt(rowSums(matrix(squared_log_returns, nrow = window_length, ncol = window, byrow = TRUE)))
  daily_vol <- real_vol * sqrt_term
  
  # Average sigma in each window (vectorized)
  sigma_vol_win <- rowMeans(matrix(sigma[2:n_plus_1], nrow = window_length, ncol = window, byrow = TRUE))
  
  # Store results in RV_check and IV_check
  RV_check[k, ] <- daily_vol
  IV_check[k, ] <- sigma_vol_win
}

X = RV_check
X = IV_check

p_inv_sols <- numeric(simtim)
for (i in 1:simtim){
  objective_function <- function(p) {
    (Wstat(L = L, K = K, p, t=T , X = X[i,]) - T)^2  # Squared difference to minimize
  }
  p_inv_sols[i] <- 1/optimize(objective_function, interval = c(1,20))$minimum
}


plot(density(p_inv_sols), 
     main = "Density Plot of Data",
     xlab = "Value",
     ylab = "Density",
     col = "blue", 
     lwd = 2)
