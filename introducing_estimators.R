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
## Simulate fBm by spectral simulation
# Define a+ and a- terms
a_plus <- function(j, lambda) 2 * pi * j + lambda
a_minus <- function(j, lambda) 2 * pi * j - lambda
# Define the B3(lambda, H) function
B3_lambda_H <- function(lambda, H) {
  j_values <- 1:3
  term1 <- (a_plus(j_values, lambda))^(-2 * H - 1)
  term2 <- (a_minus(j_values, lambda))^(-2 * H - 1)
  sum_result <- sum(term1 + term2)
  
  # Additional terms for j=3 and j=4
  additional_terms <- sum((a_plus(c(3, 4), lambda))^(-2 * H) + (a_minus(c(3, 4), lambda))^(-2 * H)) / (8 * H * pi)
  
  return(sum_result + additional_terms)
}
# Define the f(lambda) function using the new B3_lambda_H
f_lambda <- function(lambda, H, start_term = 2 * gamma(2 * H + 1) * sin(pi * H)) {
  if (lambda == 0) return(Inf)  # Handle singularity
  B3_term <- B3_lambda_H(lambda, H)
  start_term * (1 - cos(lambda)) * (abs(lambda)^(-2 * H - 1) + B3_term)
}
# Function to compute a_k sequence

compute_a_k <- function(l, H) {
  # Generate i.i.d. standard normal random variables U^(0) and U^(1)
  U_0 <- rnorm(l)  # U^(0) for k = 0, ..., l-1
  U_1 <- rnorm(l)  # U^(1) for k = 0, ..., l-1
  #t_k <- 0:l
  #f_tk <- 0:l
  #for (k in 0:l){
  #  t_k[k+1] <- pi * k/l
  #  f_tk[k+1] <- f_lambda(t_k[k+1],H)
  #}
  
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
# Final sim fBm function
fbm_sim <- function(n, l, H, T){
  d_fbm <- Re(fft(compute_a_k(l, H))[1:n])
  X <- c(0,cumsum(d_fbm))
  X <- X*(T/n)^H
  return(X)
}

T <- 10
n <- 1000
H <- 0.2
simtim <- 5
l <- ifelse(n*3>=30000, n*3, 30000)
dt <- T / n 
time <- seq(0, T, by = dt)

X_check <- matrix(0,simtim, n)

for (i in 1:simtim){
  X_check[i,] <- Re(fft(compute_a_k(l, H))[1:n])
} #FFT method

Y <- matrix(0,simtim,n+1)
for (i in 2:(n+1)){
  Y[,i] <- Y[,i-1]+X_check[,i-1]
}

Y <- Y*(T/n)^H

plot(time, Y[1,], type = "l", col = "blue", xlab = "time (s)", ylab = "y", ylim = c(min(Y),max(Y)))
for (i in 2:nrow(Y)) {
  lines(time, Y[i, ], col = rainbow(nrow(Y))[i], lwd = 2)
}

#Roughness estimator using simulations of fBm
set.seed(1)
L = 300*300
K = 300
n <- L
l <- ifelse(n*3>=30000, n*3, 30000)

H = 0.8 #Start from H=0,1 after set.seed to recreate the saved plots 
T= 1
X <- fbm_sim(n+1,l,H,T)

# Define a sequence of p values
inv_p_values <- seq(H-0.09, H+0.09, length.out = 1000)  # Avoid p=0 to prevent division by zero
p_values <- 1 / inv_p_values

# Compute W for each p
W_values <- sapply(p_values, function(p) Wstat(L = L, K = K, p, t = 1, X = X))

# Calculate log(W) to avoid log(0) issues
#log_W_values <- log(W_values)

# Log scale plot of W against 1/p 
plot(inv_p_values, W_values, log= "y", type = "l", col = "black", lwd = 1,
     main = "",
     xlab = "H=1/p",
     ylab = expression(W(L, K, pi, p, t == 1, X)), yaxt = "n")
# Define where to place the ticks (powers of 10)
y_ticks <- c(10^-0.3, 10^-0.15, 10^0, 10^0.15, 10^0.3) #H=0,8
y_ticks <- c(10^-0.4, 10^-0.2, 10^0, 10^0.2, 10^0.4, 10^0.6) #H=0,5
y_ticks <- c(10^-0.5, 10^0, 10^0.5, 10^1) #H=0,3
y_ticks <- c(10^0, 10^10, 10^20, 10^30) #H=0,1

# Add y-axis with labels as powers of 10
axis(2, at = y_ticks, labels = expression(10^-0.3, 10^-0.15, 10^0, 10^0.15, 10^0.3)) #H=0,8
axis(2, at = y_ticks, labels = expression(10^-0.4, 10^-0.2, 10^0, 10^0.2, 10^0.4, 10^0.6)) #H=0,5
axis(2, at = y_ticks, labels = expression(10^-0.5, 10^0, 10^0.5, 10^1)) #H=0,3
axis(2, at = y_ticks, labels = expression(10^0, 10^10, 10^20, 10^30)) #H=0,1

abline(h = 1, col = "blue", lwd = 1, lty = 1)  #Estimating \hat{p}
abline(v = H, col = "black", lwd = 1, lty = 2)  # True H

# Find the index where log(W) crosses y = 0
crossing_index <- which.min(abs(W_values-1))
x_at_y_0 <- inv_p_values[crossing_index]
abline(v = x_at_y_0, col = "blue", lwd = 1, lty = 1)


## Histogram of estimated roughness index
set.seed(1)
L <- 2000*2000
K <- 2000
n <- L
l <- ifelse(n*3>=30000, n*3, 30000)

simtim <- 150 #Number of times to simulate the fBm (>1)
H = 0.1
T= 1


p_inv_sols <- numeric(simtim)
for (i in 1:simtim){
  X <- fbm_sim(n,l,H,T)
  objective_function <- function(p) {
    (Wstat(L = L, K = K, p, t=1 , X = X) - T)^2  # Squared difference to minimize
  }
  p_inv_sols[i] <- 1/optimize(objective_function, interval = c(1,30))$minimum
} 

#dens <- density(p_inv_sols) # For K=2000 density
hist(p_inv_sols, main = paste("H=",H), xlab = "H=1/p", ylab = "Density") #,  ylim=range(dens$y) #add ylim for K=2000
#lines(dens, col = "blue", lwd = 2) # For K=2000 density plot

quartiles <- as.numeric(quantile(p_inv_sols, probs = c(0.25, 0.75)))

#Only run one of the lines below
#table_hist <- matrix(0,4,7)
table_hist[1,] <- c(H, min(p_inv_sols),quartiles[1], median(p_inv_sols), mean(p_inv_sols), quartiles[2], max(p_inv_sols)) #H=0,1
table_hist[2,] <- c(H, min(p_inv_sols),quartiles[1], median(p_inv_sols), mean(p_inv_sols), quartiles[2], max(p_inv_sols)) #H=0,3
table_hist[3,] <- c(H, min(p_inv_sols),quartiles[1], median(p_inv_sols), mean(p_inv_sols), quartiles[2], max(p_inv_sols)) #H=0,5
table_hist[4,] <- c(H, min(p_inv_sols),quartiles[1], median(p_inv_sols), mean(p_inv_sols), quartiles[2], max(p_inv_sols)) #H=0,8
table_hist
c(H, min(p_inv_sols),quartiles[1], median(p_inv_sols), mean(p_inv_sols), quartiles[2], max(p_inv_sols))


## Estimated H plotted against different values of K
set.seed(1)
L <- 300*300
n <- L
l <- ifelse(n*3>=30000, n*3, 30000)
H = 0.1
T= 1

X <- fbm_sim(n,l,H,T)

roughness_fct <- function(K){
  objective_function <- function(p) {
    (Wstat(L = L, K = K, p, t=1 , X = X) - T)^2  # Squared difference to minimize
  }
  return(1/optimize(objective_function, interval = c(1,30))$minimum)
}

K_values <- 2:1000 
roughness_index <- sapply(K_values, function(K) roughness_fct(K))  # apply your function

# Create the plot
plot(K_values, roughness_index, type = "l", xlab = "K", ylab = "Estimated roughness index", 
     main = "")
abline(h = H, col = "blue", lwd = 1, lty = 1) 
abline(v = sqrt(L), col = "blue", lwd = 1, lty = 1)


#### Roughness estimator via log regression
selected_data <- subset(oxfordmanrealizedvolatilityindices, Symbol == ".SPX")
oxford_real_vol <- selected_data$bv #Using oxford data
daily_vol_est <- sqrt(oxford_real_vol[1:3500])
data_vol <- sqrt(oxford_real_vol[1:3500])
plot(data_vol, type = "l", col = "blue", lwd = 2, 
     xlab = "", ylab = "", ylim = c(0, 0.06),  main = "")

mqdelta <- function(cap_delta,q){
  m_qdelta <- numeric()
  for (i in 1:cap_delta){
    delta_elements <- daily_vol_est[seq(i, length(daily_vol_est), by = cap_delta)]
    vol_term <- numeric()
    for (j in 1:(length(delta_elements)-1)){
      vol_term[j] <- (abs(log(delta_elements[j+1])-log(delta_elements[j])))^q 
    }
    m_qdelta[i] <- 1/(length(delta_elements)-1)*sum(vol_term)
  }
  final_m <- mean(m_qdelta)
  return(log(final_m))
}
# Define a range of cap_delta (Δ) values and q values
Delta_values <- 1:50
log_Delta <- log(Delta_values)
q_values <- c(0.5, 1, 1.5, 2, 3)  # Different q values

# Calculate log(m_q(Δ)) for each q and Δ
log_mq_list <- lapply(q_values, function(q) {
  sapply(Delta_values, function(delta) mqdelta(delta, q))  # Calculate log_mq for each delta
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

L <- 3500
K <- sqrt(L)

objective_function <- function(p) {
  (Wstat(L = L, K = K, p, t=1 , X = daily_vol_est) - 1)^2  # Squared difference to minimize
}
p_inv_sols <- 1/optimize(objective_function, interval = c(1,30))$minimum
p_inv_sols

# Define a sequence of p values
inv_p_values <- seq(p_inv_sols-0.09, p_inv_sols+0.09, length.out = 1000)  # Avoid p=0 to prevent division by zero
p_values <- 1 / inv_p_values

W_values <- sapply(p_values, function(p) Wstat(L = L, K = K, p, t = 1, X = daily_vol_est))
plot(inv_p_values, W_values, log= "y", type = "l", col = "black", lwd = 1,
     main = paste("Estimated H=",round(p_inv_sols, digits = 4)),
     xlab = "H=1/p",
     ylab = expression(W(L, K, pi, p, t == 1, X)), yaxt = "n")
# Define where to place the ticks (powers of 10)
y_ticks <- c(10^0, 10^1, 10^2, 10^3) 
axis(2, at = y_ticks, labels = expression(10^0, 10^1, 10^2, 10^3))
y_ticks <- c(10^-0.5, 10^0, 10^0.5, 10^1) #H=0,3
axis(2, at = y_ticks, labels = expression(10^-0.5, 10^0, 10^0.5, 10^1)) #H=0,3

abline(h = 1, col = "blue", lwd = 1, lty = 1)  #Estimating \hat{p}
crossing_index <- which.min(abs(W_values-1))
x_at_y_0 <- inv_p_values[crossing_index]
abline(v = x_at_y_0, col = "blue", lwd = 1, lty = 1)




#### Sequential scale estimator
small_n <- 13
capital_N <- small_n + 6
n <- 2^capital_N # Excluding T=0
T <- 1
H <- 0.7
l <- ifelse(n*3>=30000, n*3, 30000)
dt <- T/n

t_k <- 0:l
f_tk <- 0:l
for (k in 0:l){
  t_k[k+1] <- pi * k/l
  f_tk[k+1] <- f_lambda(t_k[k+1],H)
}

W <- fbm_sim(n,l,H,T)
sum_terms <- cumsum(W[2:(n+1)])
# Compute indices for 2^4 * k where k = 1, 2, ..., 2^(n+2)
k <- 1:(2^(small_n+2))  # Sequence for k
indices <- 2^4 * k  # Compute the desired indices
# Extract the elements
result <- sum_terms[indices]
y_t <- c(0, 2^-capital_N * result)

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

set.seed(123)
simtim <- 1000
expo_est_sols <- numeric(simtim)
scale_est_sols <- numeric(simtim)
for (i in 1:simtim){
  W <- fbm_sim(n,l,H,T)
  sum_terms <- cumsum(W[2:(n+1)])
  result <- sum_terms[indices]
  
  y_t <- c(0, 2^-capital_N * result)
  expo_est_sols[i] <- rough_expo(k = small_n, y_t = y_t, n_max = small_n)
  lambda_opt <- optimize(obj_funct, interval = c(0.001,30))$minimum
  scale_est_sols[i] <- rough_expo(k = small_n, y_t = lambda_opt*y_t, n_max = small_n)
} 


## Fractional OU
small_n <- 13
capital_N <- small_n + 6
n <- 2^capital_N # Excluding T=0
T <- 1
H <- 0.3
l <- ifelse(n*3>=30000, n*3, 30000)
dt <- T/n

t_k <- 0:l
f_tk <- 0:l
for (k in 0:l){
  t_k[k+1] <- pi * k/l
  f_tk[k+1] <- f_lambda(t_k[k+1],H)
}

d_fbm <- Re(fft(compute_a_k(l, H))[1:n])*(dt)^H
rho <- 0.2
mu <- 2
Y_0 <- 0

Y <- numeric(n + 1)
Y[1] <- Y_0
for (i in seq_len(n)) {
  Y[i + 1] <- Y[i] + rho * (mu - Y[i]) * dt + d_fbm[i]
}

sum_terms <- cumsum(exp(2*Y[2:(n+1)]))
# Compute indices for 2^4 * k where k = 1, 2, ..., 2^(n+2)
k <- 1:(2^(small_n+2))  # Sequence for k
indices <- 2^4 * k  # Compute the desired indices
# Extract the elements
result <- sum_terms[indices]
y_t <- c(0, 2^-capital_N * result)

lambda_opt <- optimize(obj_funct, interval = c(0.001,30))$minimum
scale_est <- rough_expo(k = small_n, y_t = lambda_opt*y_t, n_max = small_n)
scale_est
rough_expo(k = small_n, y_t = y_t, n_max = small_n)

set.seed(123)
simtim <- 100
expo_est_sols <- numeric(simtim)
scale_est_sols <- numeric(simtim)
for (i in 1:simtim){
  d_fbm <- Re(fft(compute_a_k(l, H))[1:n])*(dt)^H
  Y <- numeric(n + 1)
  Y[1] <- Y_0
  for (j in seq_len(n)) {
    Y[j + 1] <- Y[j] + rho * (mu - Y[j]) * dt + d_fbm[j]
  }
  
  sum_terms <- cumsum(exp(2*Y[2:(n+1)]))
  result <- sum_terms[indices]
  
  y_t <- c(0, 2^-capital_N * result)
  expo_est_sols[i] <- rough_expo(k = small_n, y_t = y_t, n_max = small_n)
  lambda_opt <- optimize(obj_funct, interval = c(0.001,30))$minimum
  scale_est_sols[i] <- rough_expo(k = small_n, y_t = lambda_opt*y_t, n_max = small_n)
} 

expo12 <- expo_est_sols
scale12 <- scale_est_sols
#expo13 <- expo_est_sols
#scale13 <- scale_est_sols
#expo14 <- expo_est_sols
#scale14 <- scale_est_sols
#expo15 <- expo_est_sols
#scale15 <- scale_est_sols
#expo16 <- expo_est_sols
#scale16 <- scale_est_sols

expo_data <- data.frame(
  value = c(expo12, expo13, expo14, expo15, expo16),
  group = rep(12:16, each = simtim)
)
scale_data <- data.frame(
  value = c(scale12, scale13, scale14, scale15, scale16),
  group = rep(12:16, each = simtim)
)
boxplot(expo_est_sols)

# Plot boxplot
boxplot(value ~ group, data = expo_data, 
        range = 0,
        col = "orange", # Fill color
        #border = "black", # Border color
        boxwex = 0.6, # Width of boxes
        #ylim = c(0.255, 0.275), # Y-axis limits
        #outline = FALSE,
        xlab = "", ylab = "", # Axis labels
        main = "", # Title
        medlwd = 1.5, # Adjust thickness of median line
        las = 1) # Rotate Y-axis labels for better readability

