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
a_plus <- function(j, lambda) {
  return(2 * pi * j + lambda)
}
a_minus <- function(j, lambda) {
  return(2 * pi * j - lambda)
}
B3_lambda_H <- function(lambda, H) {
  sum_result <- 0
  # First sum from j=1 to 3
  for (j in 1:3) {
    term1 <- (a_plus(j, lambda))^(-2 * H - 1)
    term2 <- (a_minus(j, lambda))^(-2 * H - 1)
    sum_result <- sum_result + term1 + term2
  }
  
  # Additional terms for j=3 and j=4
  term3_1 <- (a_plus(3, lambda))^(-2 * H)
  term3_2 <- (a_minus(3, lambda))^(-2 * H)
  term4_1 <- (a_plus(4, lambda))^(-2 * H)
  term4_2 <- (a_minus(4, lambda))^(-2 * H)
  
  additional_terms <- (term3_1 + term3_2 + term4_1 + term4_2) / (8 * H * pi)
  
  return(sum_result + additional_terms)
}
f_lambda <- function(lambda, H) {
  # Handle the singularity at lambda = 0
  if (lambda == 0) {
    return(Inf)  # In practice, you may want to set a very large number instead of Inf
  }
  
  gamma_term <- gamma(2 * H + 1)  # Gamma(2H + 1)
  sine_term <- sin(pi * H)  # sin(pi * H)
  cos_term <- 1 - cos(lambda)  # (1 - cos(lambda))
  
  # Compute the B3(lambda, H) term
  B3_term <- B3_lambda_H(lambda, H)
  
  # Compute the full f(lambda)
  f_result <- 2 * sine_term * gamma_term * cos_term * (abs(lambda)^(-2 * H - 1) + B3_term)
  
  return(f_result)
}
#U_0 <- rnorm(l)  # U^(0) for k = 0, ..., l-1
#U_1 <- rnorm(l)  # U^(1) for k = 0, ..., l-1
#t_k <- 0:l
#for (k in 0:l){
  #t_k[k+1] <- pi * k/l
#}
compute_a_k <- function(l, H) {
  # Generate i.i.d. standard normal random variables U^(0) and U^(1)

  
  f_tk <- 0:l
  for (k in 0:l){
    f_tk[k+1] <- f_lambda(t_k[k+1],H)
  }
  
  # Initialize a_k as a vector of complex numbers
  a_k <- complex(real = rep(0, 2 * l), imaginary = rep(0, 2 * l))
  
  # Loop to compute a_k for each k = 0, ..., 2l - 1
  for (k in 0:(2 * l - 1)) {
    if (k == 0) {
      # a_0 = 0
      a_k[k + 1] <- 0
    } else if (k >= 1 && k <= (l - 1)) {
      # a_k = 1/2 * (U^(0)_(k-1) + i * U^(1)_(k-1)) * sqrt(f(t_k) / l)
      a_k[k + 1] <- 0.5 * (U_0[k] + 1i * U_1[k]) * sqrt(f_tk[k+1] / l)
    } else if (k == l) {
      # a_l = U^(0)_(l-1) * sqrt(f(t_k) / l)
      a_k[k + 1] <- U_0[l] * sqrt(f_tk[l+1] / l)
    } else {
      # a_k = 1/2 * (U^(0)_(2l-k-1) - i * U^(1)_(2l-k-1)) * sqrt(f(t_k) / l)
      a_k[k + 1] <- 0.5 * (U_0[2 * l - k] - 1i * U_1[2 * l - k]) * sqrt(f_tk[2*l-k+1] / l)
    }
  }
  
  return(a_k)
}
# fbm_sim <- function(n, l, H, T){
  test <- Re(fft(compute_a_k(l, H))[1:n])
  X <- cumsum(test)
  X <- X*(T/n)^H
  return(X)
} # Final sim fBm function
mqdelta <- function(cap_delta,q, X){
  m_qdelta <- numeric()
  for (i in 1:cap_delta){
    delta_elements <- X[seq(i, length(X), by = cap_delta)]
    vol_term <- numeric()
    for (j in 1:(length(delta_elements)-1)){
      vol_term[j] <- abs(log(delta_elements[j+1])-log(delta_elements[j]))^q 
    }
    m_qdelta[i] <- 1/(length(delta_elements)-1)*sum(vol_term)
  }
  final_m <- mean(m_qdelta)
  return(log(final_m))
} # For log-regression smoothness

### Example 5
## First attempt to replicate (just 300 point each day)
# Parameters
L = 300*300
K = 300
#n <- L+1
window <- 300
n <- window*(L+1)
#n <- 10000
T <- 100
dt <- T / n 
time <- seq(0, T, by = dt)
S0 <- 1


# Second attempt Parameters
L = 300*300
K = 300
trade_hours <- 8
window <- 300
n <- (L+1)*window*trade_hours
T <- 1
#delta <- 1/300
#n <- round(T/delta, digits = 0)
dt <- T / n 
time <- seq(0, T, by = dt)#[1:20000]
S0 <- 1

# Simulate two independent Brownian motions B_t and W_t
set.seed(1)
dB <- rnorm(n, mean = 0, sd = sqrt(dt))
dW <- rnorm(n, mean = 0, sd = sqrt(dt))
W <- c(0, cumsum(dW))
sigma_vol <- abs(W)
exp_W <- exp(sigma_vol)


# Simulate the process S_t using Euler-Maruyama
S <- numeric(n+1)
S[1] <- S0

for (i in 1:n) {
  S[i + 1] <- S[i] + sigma_vol[i] * S[i] * dB[i]
}

# Plot the results
plot(time[1:20000], S[1:20000], type = "l", col = "blue", lwd = 2, 
     xlab = "Time", ylab = "S_t", main = "Simulated Process S_t")
min(S)

# First attempt Realized volatility
log_returns <- diff(log(S))
squared_log_returns <- log_returns^2
real_vol <- numeric()

for (i in 1:(floor((n+1)/window))){
  real_vol[i] <- sqrt(sum(squared_log_returns[(1+(i-1)*window):(i*window)]))
}
daily_vol <- real_vol*sqrt(((n+1)/T)/window)
#sigma_vol_win <- numeric()
#for (i in 1:(floor((n+1)/window))){
  #sigma_vol_win[i] <- mean(sigma_vol[(1+(i-1)*window):(i*window)])
#}
sigma_vol_win <- sigma_vol[seq(1, n, by = window)]

## For log-res
log_returns <- diff(log(exp_W))
squared_log_returns <- log_returns^2
real_vol <- numeric()

for (i in 1:(floor((n+1)/window))){
  real_vol[i] <- sqrt(sum(squared_log_returns[(1+(i-1)*window):(i*window)]))
}
daily_vol <- real_vol*sqrt(((n+1)/T)/window)

sigma_vol_win <- numeric()
for (i in 1:(floor((n+1)/window))){
  sigma_vol_win[i] <- exp_W[(1+(i-1)*window)]
}

# Second attempt Realized volatility
log_returns <- diff(log(S))
squared_log_returns <- log_returns^2
real_vol <- numeric()

for (i in 1:(floor((n+1)/(window)))){
  real_vol[i] <- sqrt(sum(squared_log_returns[(1+(i-1)*window):(i*window)]))
}
daily_vol <- numeric()
for (i in 1:(n/window/trade_hours)){
  daily_vol[i] <- mean(real_vol[(1+(i-1)*trade_hours):(i*trade_hours)])*sqrt(((n+1)/T)/window)
}

mean_vol <- numeric()
for (i in 1:(floor(n/window))){
  mean_vol[i] <- mean(sigma_vol[((1+(i-1)*window)):(i*window)])
}
sigma_vol_win <- numeric()
for (i in 1:(n/window/trade_hours)){
  sigma_vol_win[i] <- mean(mean_vol[(1+(i-1)*trade_hours):(i*trade_hours)])* sqrt(T^-1)
}

# First attempt plot
plot(daily_vol[1:10000], type = "l", col = "black", lwd = 2, 
     xlab = "Time",xaxt = "n", ylab = "volatility",  main = "")
lines(sigma_vol_win[1:10000], col = "red", lwd = 2)
min(sigma_vol_win)

IV_minus_RV <- sigma_vol_win-daily_vol
plot(IV_minus_RV[1:10000], type = "l", col = "black", lwd = 2, 
     xlab = "", xaxt = "n", ylab = "IV-RV", main = "")
log_IV_minus_RV <- log(sigma_vol_win)-log(daily_vol)
plot(log_IV_minus_RV[1:10000], type = "l", col = "black", lwd = 2, 
     xlab = "", xaxt = "n", ylab = "log(IV)-log(RV)", main = "")

# Second attempt plot
plot(daily_vol, type = "l", col = "black", lwd = 2, 
     xlab = "Time",xaxt = "n", ylab = "volatility",  main = "")
lines(sigma_vol_win, col = "red", lwd = 2)

IV_minus_RV <- sigma_vol_win-daily_vol
plot(IV_minus_RV, type = "l", col = "black", lwd = 2, 
     xlab = "", xaxt = "n", ylab = "IV-RV", main = "")
log_IV_minus_RV <- log(sigma_vol_win)-log(daily_vol)
plot(log_IV_minus_RV, type = "l", col = "black", lwd = 2, 
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
     xlab = "1/p",
     ylab = "W", yaxt = "n")  #

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



# Define a range of cap_delta (Δ) values and q values
Delta_values <- 1:50
log_Delta <- log(Delta_values)
q_values <- c(0.5, 1, 1.5, 2, 3)  # Different q values

#test_size <- L
#window <- 300
#dW <- rnorm(test_size, mean = 0, sd = 1)
#W <- c(0, cumsum(dW))
#abs_W <- abs(W)
#test <- numeric()
#for (i in 1:(test_size-window)){
  #test[i] <- mean(abs_W[i:(i-1+window)])
#}

#mean_test <- numeric()
#for (k in 1:(test_size/window)){
  #mean_test[k] <- mean(abs_W[((1+(k-1)*window)):(k*window)])
#}
#test2 <- numeric()
#for (i in 1:(test_size/window-24)){
  #test2[i] <- mean(mean_test[i:(i-1+24)])
#}
#test3 <- numeric()
#for (i in 1:(test_size/window/8)){
  #test3[i] <- mean(mean_test[(1+(i-1)*8):(i*8)])
#}

# Calculate log(m_q(Δ)) for each q and Δ
log_mq_list <- lapply(q_values, function(q) {
  sapply(Delta_values, function(delta) mqdelta(delta, q, X = X))  # Calculate log_mq for each delta
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

legend("bottomright", legend = paste("q =", q_values), col = colors, pch = 21, cex = 0.7)

# Data
log_delta <- c(log_Delta, log_Delta, log_Delta, log_Delta, log_Delta)  # Vector of log(Δ) values
log_m <- c(log_mq_list[[1]], log_mq_list[[2]], log_mq_list[[3]], log_mq_list[[4]], log_mq_list[[5]])      # Vector of log(m(q, Δ)) values
group <- c(rep(1,length(Delta_values)), rep(2,length(Delta_values)), rep(3,length(Delta_values)), rep(4,length(Delta_values)), rep(5,length(Delta_values)))      # Vector indicating the group (1 to 5 for each set)
# Combine into a data frame
data <- data.frame(log_delta, log_m, group)
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
est_smooth_H

plot(q_values, epsilon_q, type = "p", pch = 21, main = paste("Estimated H=",round(est_smooth_H, digits = 4)), 
     xlab = 'q', ylab = expression(zeta[q]), xlim = c(0,max(q_values)), ylim = c(0,max(epsilon_q)))
abline(fit_epsilon)



## Different approach to simulating the model
# Initialize variables
S_t <- numeric(n + 1)
S_t[1] <- S0
# Simulate independent Brownian motion W_t
W_t <- cumsum(c(0, sqrt(dt) * rnorm(n)))
# Time-dependent volatility sigma_t = |W_t|
sigma_t <- abs(W_t)
# Simulate S_t using Euler-Maruyama method
for (i in 1:n) {
  dB_t <- sqrt(dt) * rnorm(1)  # Brownian increment for B_t
  S_t[i + 1] <- S_t[i] + sigma_t[i] * S_t[i] * dB_t
}
# Plot the volatility process (sigma_t = |W_t|)
plot(time[1:20000], sigma_t[1:20000], type = "l", col = "red", xlab = "Time", ylab = expression(sigma[t]), main = expression(sigma[t] == abs(W[t])))
# Plot the S_t process
plot(time[1:20000], S_t[1:20000], type = "l", col = "blue", xlab = "Time", ylab = expression(S[t]), main = expression(S[t] ~ "with stochastic volatility"))

log_returns <- diff(log(S_t))
squared_log_returns <- log_returns^2
real_vol <- numeric()

for (i in 1:(floor((n+1)/window))){
  real_vol[i] <- sqrt(sum(squared_log_returns[(1+(i-1)*window):(i*window)]))
}
daily_vol <- real_vol*sqrt(((n+1)/T)/window)
sigma_vol_win <- numeric()
for (i in 1:(floor((n+1)/window))){
  sigma_vol_win[i] <- mean(sigma_t[(1+(i-1)*window):(i*window)])
}

simtim <- 100
RV_check <- matrix(0,simtim, L+1)
IV_check <- matrix(0,simtim, L+1)
# Precompute constants
sqrt_dt <- sqrt(dt)
n_plus_1 <- n + 1
window_length <- floor(n_plus_1 / window)
sqrt_term <- sqrt(n_plus_1 / (T * window))
for (k in 1:simtim){
  dB <- rnorm(n, mean = 0, sd = sqrt(dt))
  dW <- rnorm(n, mean = 0, sd = sqrt(dt))
  W <- numeric(n_plus_1)
  W[2:(n+1)] <- cumsum(dW)
  sigma_vol <- abs(W)
  
  # Simulate the process S_t using Euler-Maruyama
  S <- numeric(n_plus_1)
  S[1] <- S0
  
  for (i in 1:n) {
    S[i + 1] <- S[i] + sigma_vol[i] * S[i] * dB[i]
  }
  # Realized volatility
  log_returns <- diff(log(abs(S)))
  squared_log_returns <- log_returns^2
  
  # Realized volatility calculation (vectorized)
  real_vol <- sqrt(rowSums(matrix(squared_log_returns, nrow = window_length, ncol = window, byrow = TRUE)))
  daily_vol <- real_vol * sqrt_term
  
  # Average sigma in each window (vectorized)
  sigma_vol_win <- rowMeans(matrix(sigma_vol[2:n_plus_1], nrow = window_length, ncol = window, byrow = TRUE))
  
  # Store results in RV_check and IV_check
  RV_check[k, ] <- daily_vol
  IV_check[k, ] <- sigma_vol_win
}

X = RV_check
X = IV_check
p_inv_sols <- numeric(simtim)
for (i in 1:simtim){
  objective_function <- function(p) {
    (Wstat(L = L, K = K, p, t=1 , X = X[i,]) - 1)^2  # Squared difference to minimize
  }
  p_inv_sols[i] <- 1/optimize(objective_function, interval = c(1,30))$minimum
}
plot(density(p_inv_sols), 
     main = "Density Plot of Data",
     xlab = "H=1/p",
     ylab = "Density",
     col = "red", 
     lwd = 2)  # Density plot of roughness estimate
hist(p_inv_sols)
p_inv_sols

H_smooth_sols <- numeric(simtim)
for (i in 1:simtim){
  H_smooth_sols[i] <- smoothness_fct(X = X[i,])
}
hist(H_smooth_sols)

### Example 6 OU-SV model
# Parameters
L = 300*300
K = 300
#n <- L+1
window <- 300
n <- window*(L+1)
#n <- 10000
T <- 1
dt <- T / n 
time <- seq(0, T, by = dt)

sigma_0 <- 1
Y_0 <- 0
gamma <- 1
theta <- 1
S0 <- 1

set.seed(12)
# Simulate two independent Brownian motions B_t and W_t
dB_1 <- rnorm(n, mean = 0, sd = sqrt(dt))
#B_1 <- c(0, cumsum(dB_1))
dB_2 <- rnorm(n, mean = 0, sd = sqrt(dt))
#B_2 <- c(0, cumsum(dB_2))

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
#sigma_vol_win <- numeric()
#for (i in 1:(floor((n+1)/window))){
#  sigma_vol_win[i] <- mean(sigma[(1+(i-1)*window):(i*window)])
#}
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
p_inv_sols <- 1/optimize(objective_function, interval = c(1,20))$minimum
p_inv_sols

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
     main = "Density Plot of Data",
     xlab = "H=1/p",
     ylab = "Density",
     col = "red", 
     lwd = 2)




### A fractional Ornstein-Uhlenbeck model (Example 7)
set.seed(1)
L = 100*100
K = 100
window <- 300
n <- window*(L+1)
l <- ifelse(n*3>=30000, n*3, 30000)
T <- 1
dt <- T / n 
time <- seq(0, T, by = dt)

U_0 <- rnorm(l)  # U^(0) for k = 0, ..., l-1
U_1 <- rnorm(l)  # U^(1) for k = 0, ..., l-1
t_k <- 0:l
for (k in 0:l){
  t_k[k+1] <- pi * k/l
}

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

output_given_H <- function(H){
  
  d_fbm <- Re(fft(compute_a_k(l, H))[1:n])*(T/n)^H

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
  daily_vol <- real_vol * sqrt_term
  # Average sigma in each window (vectorized)
  #sigma_vol_win <- rowMeans(matrix(sigma[1:n], nrow = window_length, ncol = window, byrow = TRUE))
  sigma_vol_win <- sigma[seq(1, n, by = window)]
  # Average S in each window
  S_win <- rowMeans(matrix(S[2:n_plus_1], nrow = window_length, ncol = window, byrow = TRUE))

  IV <- sigma_vol_win
  RV <- daily_vol
  return(list(S_win = S_win, IV = IV, RV = RV))
}

H_values <- seq(0.10, 0.80, by = 0.10)
store_results <- list()
for (H in H_values) {
  store_results[[paste0("H_", H)]] <- output_given_H(H)
}

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

results <- data.frame(H = numeric(), Instantaneous_volatility = numeric(), Realized_volatility = numeric(), IV_logres = numeric(), RV_logres = numeric())
for (H in H_values) {
  IV_H_est <- roughness_fct(store_results[[paste0("H_", H)]]$IV)
  RV_H_est <- roughness_fct(store_results[[paste0("H_", H)]]$RV)
  IV_H_est_logres <- smoothness_fct(store_results[[paste0("H_", H)]]$IV)
  RV_H_est_logres <- smoothness_fct(store_results[[paste0("H_", H)]]$RV)
  # Append the results to the data frame
  results <- rbind(results, data.frame(H = H, 
                                       Instantaneous_volatility = IV_H_est, 
                                       Realized_volatility = RV_H_est, 
                                       IV_logres = IV_H_est_logres, 
                                       RV_logres = RV_H_est_logres))
}
# Display the final table using knitr::kable
kable(results, col.names = c("H", "Instantaneous volatility", "Realized volatility", "IV (log-res)", "RV (log-res)"), 
      align = "c", caption = "Volatility Table")

par(mfrow = c(2, 3))
# Plot the data for each H 
for (H in H_values) {
  # Simulate three vectors for H (use stored vectors from your list)
  plot(store_results[[paste0("H_", H)]]$S_win-S_0, type = "l", ylab = "Price S_t", xlab = "", xaxt = "n")
  plot(store_results[[paste0("H_", H)]]$RV, type = "l", ylab = "Realized volatility", xlab = "", xaxt = "n")
  plot(store_results[[paste0("H_", H)]]$IV, type = "l", ylab = "Instantaneous volatility", xlab = "", xaxt = "n")
}

par(mfrow = c(1, 1))
plot(results$H, results$Instantaneous_volatility, type = "l", col = "blue", xlab = "", 
     ylab = "", xlim = c(0,0.8), ylim = c(0,0.8))
abline(v = seq(0, 0.8, by = 0.1), col = "gray", lty = "dotted", lwd = 0.5) # Vertical lines
abline(h = seq(0, 0.8, by = 0.1), col = "gray", lty = "dotted", lwd = 0.5) # Horizontal lines
lines(results$H, results$Realized_volatility, col = "red", lwd = 2)

legend("topleft", legend = c('Realized Vol', 'Instantaneous Vol'), col = c('red', 'blue'), lty = 1, cex = 0.7)

results_logres <- data.frame(H = numeric(), Instantaneous_volatility = numeric(), Realized_volatility = numeric())
for (H in H_values) {
  IV_H_est <- smoothness_fct(store_results[[paste0("H_", H)]]$IV)
  RV_H_est <- smoothness_fct(store_results[[paste0("H_", H)]]$RV)
  # Append the results to the data frame
  results_logres <- rbind(results_logres, data.frame(H = H, 
                                       Instantaneous_volatility = IV_H_est, 
                                       Realized_volatility = RV_H_est))
}
# Display the final table using knitr::kable
kable(results_logres, col.names = c("H", "IV (log-res)", "RV (log-res)"), 
      align = "c", caption = "Volatility Table") # Only log-res smoothness estimates

results <- data.frame(H = numeric(), Instantaneous_volatility = numeric(), Realized_volatility = numeric())
for (H in H_values) {
  IV_H_est <- roughness_fct(store_results[[paste0("H_", H)]]$IV)
  RV_H_est <- roughness_fct(store_results[[paste0("H_", H)]]$RV)
  # Append the results to the data frame
  results <- rbind(results, data.frame(H = H, 
                                       Instantaneous_volatility = IV_H_est, 
                                       Realized_volatility = RV_H_est))
}
# Display the final table using knitr::kable
kable(results, col.names = c("H", "Instantaneous volatility", "Realized volatility"), 
      align = "c", caption = "Volatility Table") # Only p-th variation roughness estimates
