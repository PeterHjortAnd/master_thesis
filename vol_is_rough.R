library(MASS)
library(tidyverse)
library(complex)

## Functions used throughout the document
set.seed(1)
simfbmonce <- function(n, H, T) {
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

T <- 2000
delta <- 1/20000
n <- round(T/delta, digits = 0)
H <- 0.14
nu <- 0.3
m <- -5
alpha <- 5*10^(-4)
lchoice <- 90000000
l <- ifelse(lchoice>=n/2, lchoice, n/2)


#W <- simfbmonce(n,H,T)

#all_points <- fft(compute_a_k(l, H))
test <- Re(all_points[0:n])
X <- numeric(n)
X[1] <- test[1]
for (i in 2:n){
  X[i] <- X[i-1]+test[i]
}
W <- X*(T/n)^H #Simulated fBm

X <- numeric(n)
X[1] <- m
for (i in 2:n){
  X[i] = nu*(W[i]-W[i-1])+alpha*delta*(m-X[i-1])+X[i-1]
} #RSFV simulation
vol <- exp(X) #Simulated volatility

#time <- seq(0,T, length.out = n)
#plot(time, vol, type = "l", col = "blue", xlab = "", ylab = "")


U <- rnorm(n)
P <- numeric(2)
P[1] <- 2
for (i in 2:n) {
  P[i] <- P[i-1] + P[i-1] * vol[i-1] * sqrt(delta) * U[i-1]
} #Efficient price simulation

eta <- 0.25
tick_size <- 5*10^-4

X_price <- numeric(n)
X_price[1] <- round(P[1], digits = 3)
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

time <- seq(0,T, length.out = n)
plot(time[1:100000], P[1:100000], type = "l", col = "blue", xlab = "", ylab = "")
lines(time[1:100000], X_price[1:100000], col = "red", lwd = 2)

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
eta_estimates <- numeric(1)
for (i in 1:T){
  eta_estimates[i] <- estimate_eta(X_price[(8334+(i-1)*20000):(9167+(i-1)*20000)])
}

hour_total_points <- 9167-8334
int_var_est <- numeric(T)
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
daily_vol_est <- int_var_est * sqrt(n/hour_total_points) #Daily volatility estimate (normalized values)

int_var_est <- numeric(T)
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
daily_vol_est <- int_var_est * sqrt(n/hour_total_points) #Daily volatility estimate (only tau terms)
daily_vol_est

int_var_est <- numeric(T)
for (i in 1:T){
  ten_to_11 <- X_price[(8334+(i-1)*1/delta):(9167+(i-1)*1/delta)]
  change_points <- ten_to_11[c(TRUE, diff(ten_to_11) != 0)]  # Include the first element, then select changes
  log_returns <- diff(log(ten_to_11))
  squared_log_returns <- log_returns^2
  int_var_est[i] <- sqrt(sum(squared_log_returns))
} #Make vol estimate 10 am to 11 am daily
daily_vol_est <- int_var_est * sqrt(n/hour_total_points) #Daily volatility estimate using RV def from roughvol
daily_vol_est

int_var_est <- numeric(T)
for (i in 1:T){
  X_eff_hat <- numeric(1)
  ten_to_11 <- X_price[(8334+(i-1)*1/delta):(9167+(i-1)*1/delta)]
  change_points <- ten_to_11[c(TRUE, diff(ten_to_11) != 0)]  # Include the first element, then select changes
  eta_hat <- estimate_eta(ten_to_11)
  X_eff_hat[1] <- ten_to_11[1]
  sum_term <- 0
  for (k in 2:(length(change_points))){
    X_eff_hat[k] <- change_points[k] - sign(change_points[k]-change_points[k-1])*(1/2-eta_hat)*tick_size
  }
  log_returns <- diff(log(X_eff_hat))
  squared_log_returns <- log_returns^2
  int_var_est[i] <- sqrt(sum(squared_log_returns))
} #Make vol estimate 10 am to 11 am daily
daily_vol_est <- int_var_est * sqrt(n/hour_total_points) #Daily volatility estimate (only tau terms) with RV def from roughvol
daily_vol_est

int_var_est <- numeric(T)
for (i in 1:T){
  ten_to_11 <- P[(8334+(i-1)*1/delta):(9167+(i-1)*1/delta)]
  sum_term <- 0
  for (k in 2:length(ten_to_11)){
    sum_term <- sum_term+((ten_to_11[k]-ten_to_11[k-1])/ten_to_11[k-1])^2
  }
  int_var_est[i] <- sum_term
} #Make vol estimate 10 am to 11 am daily
daily_vol_est <- int_var_est*sqrt(n/hour_total_points) #Daily volatility using efficient price directly instead of estimating
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

mean_vol <- numeric(1)
for (i in 1:T){
  mean_vol[i] <- mean(vol[(8334+(i-1)*1/delta):(9167+(i-1)*1/delta)])
}
daily_vol_est <- mean_vol #Test using directly simulated volatility instead of estimating

mqdelta <- function(cap_delta,q){
  m_qdelta <- numeric(1)
  for (i in 1:cap_delta){
    delta_elements <- daily_vol_est[seq(i, length(daily_vol_est), by = cap_delta)]
    vol_term <- numeric(2)
    for (j in 2:(T/cap_delta)){
      vol_term[j] <- abs(log(delta_elements[j])-log(delta_elements[j-1]))^q 
    }
    m_qdelta[i] <- 1/(length(delta_elements-1))*sum(vol_term)
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

# Plot the first q value
plot(log_Delta, log_mq_list[[1]], type = "p", col = "blue", pch = 19,
     xlab = expression(log(Delta)), 
     ylab = expression(log(m(q,Delta))),
     ylim = c(min(sapply(log_mq_list, min)), max(sapply(log_mq_list, max))),
     xlim = c(min(log_Delta), max(log_Delta))
)

colors <- c("blue", "green", "red", "lightblue", "purple")
# Add lines for the remaining q values
for (i in 2:length(q_values)) {
  points(log_Delta, log_mq_list[[i]], col = colors[i], pch = i + 18)  # Add points
  #lines(log_Delta, log_mq_list[[i]], col = colors[i], lty = i)  # Add lines
}

# Add legend
legend("topright", legend = paste("q =", q_values), col = colors, pch = 19:23, lty = 1:length(q_values))



