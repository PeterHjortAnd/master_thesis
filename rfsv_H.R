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

# Functions for scale estimator
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


# Simulation from RSFV
## Parameters
T <- 3601 #2000 in volis rep
delta <- 1/5000 #1/20000 in volis rep
per_day <- 1/delta
L <- T - 1
K <- sqrt(L)
n <- round(T*per_day, digits = 0)
H <- 0.5 # 0.14
nu <- 0.2 # 0.3 in vol is rough replication
m <- -7 #-5 in vol is rough replication
alpha <-1*10^-4 # 5*10^(-4) in rep

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
S_0 <- 1
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

## Constants for fBm sim
indices <- 1:(l - 1)
ak_constant <- 0.5 * (U_0[indices] + 1i * U_1[indices])
indices_2 <- (l+1):(2*l-1)
ak_constant_2 <- 0.5 * (U_0[2 * l - indices_2] - 1i * U_1[2 * l - indices_2])
indi_plus_1 <- indices + 1
indi2_plus_1 <- indices_2+1
indi_2l <- 2*l-indices_2+1
l_plus_one <- l +1 
two_l <- 2 * l


H <- H_values[1]
H
d_fbm <- Re(fft(compute_a_k(l, H))[1:n])*(dt)^H

X <- numeric(n_plus_1)
X[1] <- m
for (i in 1:n){
  X[i+1] = nu*(d_fbm[i])+alpha*delta*(m-X[i])+X[i]
} #RFSV simulation
vol <- exp(X) #Simulated volatility

S <- numeric(n_plus_1)
S[1] <- S_0
for (i in 1:n) {
  S[i+1] <- S[i] + S[i] * vol[i] * sqrt(delta) * U[i]
} #Efficient price simulation

# Realized volatility
log_returns <- diff(log(abs(S)))
squared_log_returns <- log_returns^2
# Realized volatility calculation (vectorized)
real_vol <- sqrt(rowSums(matrix(squared_log_returns, nrow = T, ncol = per_day, byrow = TRUE)))
RV <- real_vol
#IV <- vol[seq(1, n, by = 1/delta)]
IV <- rowMeans(matrix(vol[1:n], nrow = T, ncol = per_day, byrow = TRUE))

# Subset only the first 1000 points per day
first_1000_indices <- unlist(lapply(1:T, function(i) ((i - 1) * per_day + 1):((i - 1) * per_day + 1000)))
# Update log returns and squared log returns to include only the first 1000 points per day
squared_log_returns_subset <- squared_log_returns[first_1000_indices]
# Calculate Realized Volatility (RV)
real_vol <- sqrt(rowSums(matrix(squared_log_returns_subset, nrow = T, ncol = 1000, byrow = TRUE)))
RV <- real_vol*sqrt((per_day)/(1000))
# Update vol to include only the first 1000 points per day
vol_subset <- vol[first_1000_indices]
# Calculate Instantaneous Volatility (IV)
IV <- rowMeans(matrix(vol_subset, nrow = T, ncol = 1000, byrow = TRUE))


#plot(d_fbm[seq(1, length(d_fbm), by = 100)], type = "l", ylab = "d_fbm", xlab = "", xaxt = "n")
plot(vol[seq(1, length(vol), by = 100)], type = "l", ylab = "Vol", xlab = "", xaxt = "n")
plot(S[seq(1, length(S), by = 100)], type = "l", ylab = "Price S_t", xlab = "", xaxt = "n")
plot(S[1:10000]-S_0, type = "l", ylab = "Price S_t", xlab = "", xaxt = "n")
plot(RV, type = "l", ylab = "Realized volatility", xlab = "", xaxt = "n")
lines(IV, col = "red", lwd = 2)
#plot(IV, type = "l", ylab = "Instantaneous volatility", xlab = "", xaxt = "n")  

S[(1100000-5000):3490000]
S[seq(1, length(S), by = 10000)]


H_IV <- roughness_fct(IV)
H_RV <- roughness_fct(RV)
H_IV_lres <- smoothness_fct(IV)
H_RV_lres <- smoothness_fct(RV)

print(c(H_IV, H_IV_lres, H_RV, H_RV_lres))

log_mq_list <- lapply(q_values, function(q) {
  sapply(Delta_values, function(delta) mqdelta(delta, q, X = RV))  # Calculate log_mq for each delta
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

print(c(H_IV, H_IV_lres, H_RV, H_RV_lres))


# Sequential scale estimator
small_n <- 15
capital_N <- small_n + 6
number_of_points <- 2^capital_N # Excluding T=0
T <- 1
dt_new <- T/number_of_points

dfbm_new <- (dt_new/dt)^H*d_fbm[1:number_of_points]

X <- numeric(n_plus_1)
X[1] <- m
for (i in 1:n){
  X[i+1] = nu*(dfbm_new[i])+alpha*dt_new*(m-X[i])+X[i]
} #RFSV simulation
vol <- exp(X) #Simulated volatility

S <- numeric(n_plus_1)
S[1] <- S_0
for (i in 1:n) {
  S[i+1] <- S[i] + S[i] * vol[i] * sqrt(dt_new) * U[i]
} #Efficient price simulation
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

sum_terms <- cumsum(exp(2*X[2:(n+1)]))
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

lambda_opt <- optimize(obj_funct, interval = c(0.001,50))$minimum
scale_est <- rough_expo(k = small_n, y_t = lambda_opt*y_t, n_max = small_n)
scale_est
rough_expo(k = small_n, y_t = y_t, n_max = small_n)

# Subset only the first 1000 points per day
first_1000_indices <- unlist(lapply(1:T, function(i) ((i - 1) * per_day + 1):((i - 1) * per_day + 1000)))

output_given_H <- function(H){
  
  d_fbm <- Re(fft(compute_a_k(l, H))[1:n])*(dt)^H
  
  X <- numeric(n_plus_1)
  X[1] <- m
  for (i in 1:n){
    X[i+1] = nu*(d_fbm[i])+alpha*delta*(m-X[i])+X[i]
  } #RFSV simulation
  vol <- exp(X) #Simulated volatility
  
  S <- numeric(n_plus_1)
  S[1] <- S_0
  for (i in 1:n) {
    S[i+1] <- S[i] + S[i] * vol[i] * sqrt(delta) * U[i]
  } #Efficient price simulation

  # Realized volatility
  log_returns <- diff(log(abs(S)))
  squared_log_returns <- log_returns^2

  # Update log returns and squared log returns to include only the first 1000 points per day
  squared_log_returns_subset <- squared_log_returns[first_1000_indices]
  # Calculate Realized Volatility (RV)
  real_vol <- sqrt(rowSums(matrix(squared_log_returns_subset, nrow = T, ncol = 1000, byrow = TRUE)))
  RV <- real_vol*sqrt((per_day)/(1000))
  # Update vol to include only the first 1000 points per day
  vol_subset <- vol[first_1000_indices]
  # Calculate Instantaneous Volatility (IV)
  IV <- rowMeans(matrix(vol_subset, nrow = T, ncol = 1000, byrow = TRUE))
  
  plot(S[1:10000]-S_0, type = "l", ylab = "Price S_t", xlab = "", xaxt = "n")
  plot(RV, type = "l", ylab = "Realized volatility", xlab = "", xaxt = "n")
  plot(IV, type = "l", ylab = "Instantaneous volatility", xlab = "", xaxt = "n")  
  
  H_IV <- roughness_fct(IV)
  H_RV <- roughness_fct(RV)
  H_IV_lres <- smoothness_fct(IV)
  H_RV_lres <- smoothness_fct(RV)
  
  return(list(H_IV = H_IV, H_RV = H_RV, H_IV_lres = H_IV_lres, H_RV_lres = H_RV_lres))
}

H_values <- seq(0.10, 0.50, by = 0.10)
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
                      scale_vol = numeric())
for (H in H_values) {
  IV_H_est <- store_results[[paste0("H_", H)]]$H_IV
  RV_H_est <- store_results[[paste0("H_", H)]]$H_RV
  IV_H_est_logres <- store_results[[paste0("H_", H)]]$H_IV_lres
  RV_H_est_logres <- store_results[[paste0("H_", H)]]$H_RV_lres
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

plot(results$H, results$Instantaneous_volatility, type = "l", col = "blue", xlab = "", 
     ylab = "", xlim = c(0,0.6), ylim = c(0,0.6))
abline(v = seq(0, 0.6, by = 0.1), col = "gray", lty = "dotted", lwd = 0.5) # Vertical lines
abline(h = seq(0, 0.6, by = 0.1), col = "gray", lty = "dotted", lwd = 0.5) # Horizontal lines
lines(results$H, results$Realized_volatility, col = "red", lwd = 2)
lines(results$H, results$IV_logres, col = "blue", lwd = 2, lty = 2)
lines(results$H, results$RV_logres, col = "red", lwd = 2, lty = 2)

legend("topleft", legend = c('Realized Vol', 'Instantaneous Vol', 'RV (log reg)', 'IV (log reg)'), col = c('red', 'blue', 'red', 'blue'), lty = c(1,1,2,2) , cex = 0.7)



## Simulate 100 times

simoutput_given_H <- function(H){
  
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



