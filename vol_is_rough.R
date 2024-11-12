library(MASS)
library(tidyverse)
library(complex)

## Functions used throughout the document
#simfbmonce <- function(n, H, T) {
  ph<- function(H, k) 1/2*(abs(k+1)^(2*H)+abs(k-1)^(2*H)-2*abs(k)^(2*H))  #autocovariance of fGn
  covmat <- matrix(data= NA, nrow =n, ncol = n)
  for (i in 1:n){
    for (j in 1:n){
      covmat[i,j]<- ph(H, abs(i-j))
    }
  }
  X <- mvrnorm(n = 1,mu = matrix(0, n, 1), Sigma = covmat)
  Y <- matrix(0,1,n)
  for (i in 2:n){
    Y[i] <- Y[i-1]+X[i]
  }
  Ychang <- Y*(T/n)^H
  return(Ychang)
} #Simulate fBm once (simtim=1)
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
compute_a_k <- function(l, H) {
  # Generate i.i.d. standard normal random variables U^(0) and U^(1)
  U_0 <- rnorm(l)  # U^(0) for k = 0, ..., l-1
  U_1 <- rnorm(l)  # U^(1) for k = 0, ..., l-1
  t_k <- 0:l
  f_tk <- 0:l
  for (k in 0:l){
    t_k[k+1] <- pi * k/l
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
fbm_sim <- function(n, l, H, T){
  d_fbm <- Re(fft(compute_a_k(l, H))[1:n])
  X <- c(0,cumsum(d_fbm))
  X <- X*(T/n)^H
  return(X)
} # Simulate fBm using a spectral method
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


set.seed(1)
T <- 4900 #2000 in volis rep
delta <- 1/5000 #1/20000 in volis rep
n <- round(T/delta, digits = 0)
H <- 0.5 # 0.14
nu <- 0.2 # 0.3 in vol is rough replication
m <- -7 #-5 in vol is rough replication
alpha <- 5*10^(-4)
l <- ifelse(n*3>=30000, n*3, 30000)
dt <- T / n 
time <- seq(0, T, by = dt)

W <- fbm_sim(n, l, H, T)
#save_points <- W

X <- numeric()
X[1] <- m
for (i in 2:n){
  X[i] = nu*(W[i]-W[i-1])+alpha*delta*(m-X[i-1])+X[i-1]
} #RFSV simulation
vol <- exp(X) #Simulated volatility

#plot(time, vol, type = "l", col = "blue", xlab = "", ylab = "")


U <- rnorm(n)
P <- numeric()
P[1] <- 1 # 2 seems to be working in vol is rough replication
for (i in 2:n) {
  P[i] <- P[i-1] + P[i-1] * vol[i-1] * sqrt(delta) * U[i-1]
} #Efficient price simulation

eta <- 0.25
tick_size <- 5*10^-6 # 5*10^-4

X_price <- numeric()
X_price[1] <- round(P[1], digits = 5) # 3 digits
for (i in 2:n){
  if(P[i]>(X_price[i-1]+(1/2+eta)*tick_size)){
    X_price[i] <- X_price[i-1]+tick_size
    while (P[i]>(X_price[i]+(1/2+eta)*tick_size)){
      X_price[i] <- X_price[i]+tick_size
    }
  } else if(P[i]<(X_price[i-1]-(1/2+eta)*tick_size)){
    X_price[i] <- X_price[i-1]-tick_size
    while (P[i]<(X_price[i]-(1/2+eta)*tick_size)){
      X_price[i] <- X_price[i]-tick_size
    }
  } else {
    X_price[i] <- X_price[i-1]
  }
} #Observed price simulation

plot(time[1:1000], P[1:1000], type = "l", col = "blue", xlab = "", ylab = "")
lines(time[1:1000], X_price[1:1000], col = "red", lwd = 2)

plot_P <- P[seq(1, length(P), by = 100)]
plot_time <- time[seq(1, length(P), by = 100)]
plot(plot_time, plot_P, type = "l", col = "blue", xlab = "", ylab = "")
min(P)

price_changes <- function(X){
  price <- numeric(2)
  for (i in 2:length(X)){
    price[i] <- X[i]-X[i-1]
  }
  return(price)
}
price_changes(X_price[1:1000])
# Function to estimate eta
estimate_eta <- function(price) {
  # Initialize counts
  N_c <- 0  # Number of continuations
  N_a <- -1  # Number of alternations
  last_change <- 0
  price_chan <- price_changes(price)
  # Loop through the price changes to count alternations and continuations
  for (i in 2:length(price_chan)) {
    if (price_chan[i]!=0){
      if (sign(price_chan[i]) == sign(last_change)) {
        # If the price change direction is the same as the previous one, it's a continuation
        N_c <- N_c + 1
      } else {
        # If the price change direction is different, it's an alternation
        N_a <- N_a + 1
      }
      last_change <- price_chan[i]
    }
  }
  # Compute eta
  if (N_a == 0) {
    warning("No alternations found. Eta cannot be computed.")
    return(NA)
  } else {
    eta_hat <- N_c / (2 * N_a)
    return(eta_hat)
  }
}
estimate_eta(X_price[1:20000])

hour_total_points <- 9167-8334
int_var_est <- numeric()
for (i in 1:T){
  X_eff_hat <- numeric(1)
  ten_to_11 <- X_price[(8334+(i-1)*1/delta):(9167+(i-1)*1/delta)]
  eta_hat <- estimate_eta(ten_to_11)
  X_eff_hat[1] <- ten_to_11[1]
  sum_term <- 0
  for (k in 2:833){
    X_eff_hat[k] <- ten_to_11[k] - sign(ten_to_11[k]-ten_to_11[k-1])*(1/2-eta_hat)*tick_size
    sum_term <- sum_term+((X_eff_hat[k]-X_eff_hat[k-1])/X_eff_hat[k-1])^2
  }
  int_var_est[i] <- sum_term
} #Make vol estimate 10 am to 11 am daily
daily_vol_est <- sqrt(int_var_est) * sqrt(n/hour_total_points) #Daily volatility estimate (normalized values)
daily_vol_est

int_var_est <- numeric()
for (i in 1:T){
  X_eff_hat <- numeric(1)
  ten_to_11 <- X_price[(8334+(i-1)*1/delta):(9167+(i-1)*1/delta)]
  change_points <- ten_to_11[c(TRUE, diff(ten_to_11) != 0)]  # Include the first element, then select changes
  eta_hat <- estimate_eta(ten_to_11)
  X_eff_hat[1] <- ten_to_11[1]
  sum_term <- 0
  for (k in 2:(length(change_points))){
    X_eff_hat[k] <- change_points[k] - sign(change_points[k]-change_points[k-1])*(1/2-eta_hat)*tick_size
    sum_term <- sum_term+((X_eff_hat[k]-X_eff_hat[k-1])/X_eff_hat[k-1])^2
  }
  int_var_est[i] <- sum_term
} #Make vol estimate 10 am to 11 am daily
daily_vol_est <- sqrt(int_var_est) * sqrt(n/hour_total_points) #Daily volatility estimate (only tau terms)
daily_vol_est

int_var_est <- numeric()
for (i in 1:T){
  ten_to_11 <- X_price[(8334+(i-1)*1/delta):(9167+(i-1)*1/delta)]
  log_returns <- diff(log(ten_to_11))
  squared_log_returns <- log_returns^2
  int_var_est[i] <- sqrt(sum(squared_log_returns))
} #Make vol estimate 10 am to 11 am daily
daily_vol_est <- int_var_est * sqrt(n/hour_total_points) #Daily volatility estimate using RV def from roughvol
daily_vol_est

int_var_est <- numeric(T)
for (i in 1:T){
  ten_to_11 <- P[(8334+(i-1)*1/delta):(9167+(i-1)*1/delta)]
  log_returns <- diff(log(ten_to_11))
  squared_log_returns <- log_returns^2
  int_var_est[i] <- sqrt(sum(squared_log_returns))
} #Make vol estimate 10 am to 11 am daily
daily_vol_est <- int_var_est*sqrt(n/hour_total_points) #Using efficient price directly and RV def from roughvol
daily_vol_est

mean_vol <- numeric()
for (i in 1:T){
  mean_vol[i] <- mean(vol[(8334+(i-1)*1/delta):(9167+(i-1)*1/delta)])
}
daily_vol_est <- mean_vol #Test using directly simulated volatility instead of estimating

mean_vol <- numeric()
for (i in 1:T){
  mean_vol[i] <- mean(vol[(1+(i-1)*1/delta):(i*1/delta)])
}
daily_vol_est <- mean_vol #Test using directly simulated volatility as mean of the whole day

int_var_est <- numeric()
for (i in 1:T){
  X_eff_hat <- numeric(1)
  ten_to_11 <- X_price[(1+(i-1)*1/delta):(i*1/delta)]
  change_points <- ten_to_11[c(TRUE, diff(ten_to_11) != 0)]  # Include the first element, then select changes
  eta_hat <- estimate_eta(ten_to_11)
  X_eff_hat[1] <- ten_to_11[1]
  sum_term <- 0
  for (k in 2:(length(change_points))){
    X_eff_hat[k] <- change_points[k] - sign(change_points[k]-change_points[k-1])*(1/2-eta_hat)*tick_size
    sum_term <- sum_term+((X_eff_hat[k]-X_eff_hat[k-1])/X_eff_hat[k-1])^2
  }
  int_var_est[i] <- sum_term
} 
daily_vol_est <- sqrt(int_var_est) * sqrt(n*delta) #Daily volatility estimate (only tau terms) using the whole day

selected_data <- subset(oxfordmanrealizedvolatilityindices, Symbol == ".SPX")
oxford_real_vol <- selected_data$bv #Using oxford data
daily_vol_est <- sqrt(oxford_real_vol[1:3500])
plot(oxford_real_vol, type = "l", col = "black", lwd = 2, 
     xlab = "Days",xaxt = "n", ylab = "oxford realized volatility",  main = "") ## RV5 seems to be the one plotted in roughvol
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

#0.1175238 # RV mixed with UZ estimate of H
#0.1175231 # UZ estimate of H

L <- 2000
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

