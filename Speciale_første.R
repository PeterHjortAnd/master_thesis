library(MASS)
library(tidyverse)
library(complex)
## Functions used throughout the document
set.seed(1)
Wstat <- function(L, K, p, t=1, X) {
  # Initialize W to 0
  W_value <- 0
  Kseq <- seq(0,t, length.out = K+1)
  Lseq <- seq(0,t, length.out = L+1)
  
  # Iterate over the K indices
  for (i in 1:K) {
    # Calculate the difference in X for the K partition
    deltaX_K <- abs(X[1+i*L/K] - X[1+(i-1)*L/K])^p
    
    # Calculate the denominator
    deltaX_L <- 0
    for (j in (1+(i-1)*L/K):(i*L/K)) {
      # Calculate the difference in X for the L partition
      deltaX_L <- deltaX_L + abs(X[j+1] - X[j])^p
    }
    # Accumulate the sum
    W_value <- W_value + deltaX_K/deltaX_L * (Kseq[i+1]-Kseq[i])
  }
  
  return(W_value)
} #For estimation of roughness based on data
simfbm <- function(n, simtim, H, T) {
  ph<- function(H, k) 1/2*(abs(k+1)^(2*H)+abs(k-1)^(2*H)-2*abs(k)^(2*H))  #autocovariance of fGn
  covmat <- matrix(data= NA, nrow =n, ncol = n)
  for (i in 1:n){
    for (j in 1:n){
      covmat[i,j]<- ph(H, abs(i-j))
    }
  }
  X <- mvrnorm(n = simtim,mu = matrix(0, n, 1), Sigma = covmat)
  Y <- matrix(0,simtim,n)
  for (i in 2:n){
    Y[,i] <- Y[,i-1]+X[,i]
  }
  Ychang <- Y*(T/n)^H
  return(Ychang)
} #Naive way to simulate fBm
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

#Simulate fBm (naive)
n<- 1000
simtim <- 50 #Number of times to simulate the fBm (simtim>1)
H <- 0.5
T <- 10

X <- simfbm(n, simtim, H, T)
time <- seq(0,T, length.out = n)
time

plot(time, X[1,], type = "l", col = "blue", xlab = "time (s)", ylab = "y", ylim = c(min(X),max(X)))
for (i in 2:nrow(X)) {
  lines(time, X[i, ], col = rainbow(nrow(X))[i], lwd = 2)
}

## Simulate fBm by spectral simulation

# Define a+ and a- terms
a_plus <- function(j, lambda) {
  return(2 * pi * j + lambda)
}
a_minus <- function(j, lambda) {
  return(2 * pi * j - lambda)
}
# Define the B3(lambda, H) function
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
# Define the f(lambda) function using the new B3_lambda_H
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
# Function to compute a_k sequence

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
# Final sim fBm function
fbm_sim <- function(n, l, H, T){
  d_fbm <- Re(fft(compute_a_k(l, H))[1:n])
  X <- cumsum(d_fbm)
  X <- X*(T/n)^H
  return(X)
}

U_0 <- rnorm(l)  # U^(0) for k = 0, ..., l-1
U_1 <- rnorm(l)  # U^(1) for k = 0, ..., l-1
# Function to compute X^{(ell)}_n
compute_X_n <- function(n, l, H) {
  # Generate i.i.d. standard normal random variables U^(0) and U^(1)

  # Loop over k from 0 to l-1
  com_X_n <- 1:n
    for (i in 0:(n-1)){
      X_n <- 0
      for (k in 0:(l-1)) {
        t_k1 <- pi * (k + 1) / l  # t_k+1 = pi * (k+1) / l
        f_tk1 <- f_lambda(t_k1, H)  # Compute f(t_k+1) using B3(lambda, H)
    
        # Update the value of X_n based on the sum
        X_n <- X_n + sqrt(f_tk1 / l) * (cos(i * t_k1) * U_0[k+1] - sin(n * t_k1) * U_1[k+1])
      }
    com_X_n[i+1] <- X_n
  }
  
  return(com_X_n)
} #Spectral wo/ fft (very slow!)


test <- compute_X_n(n, l, H)
testny <- numeric(n)
testny[1] <- test[1]
for (i in 2:n){
  testny[i] <- testny[i-1]+test[i]
}
testny
testny <- testny*(T/n)^H

T <- 10
n <- 1000
H <- 0.2
simtim <- 5
l <- ifelse(n*3>=30000, n*3, 30000)
dt <- T / n 
time <- seq(0, T, by = dt)

X_fft <- Re(fft(compute_a_k(l, H))[1:n])
X_pure <- compute_X_n(n, l, H)
plot(X_fft, type = "l", col = "blue", xlab = "time (s)", ylab = "y")
lines(X_pure, col = "red", lwd = 2)

Y_fft <- c(0,cumsum(X_fft))*dt^H
Y_pure <- c(0,cumsum(X_pure))*dt^H
plot(time, Y_fft, type = "l", col = "blue", xlab = "time (s)", ylab = "y")
lines(time, Y_pure, col = "red", lwd = 2)

sum(X_pure[1:501])
Y_pure[990:1010]

X_check <- matrix(0,simtim, n)

for (i in 1:simtim){
  X_check[i,] <- Re(fft(compute_a_k(l, H))[1:n]) #(2*l/pi)^H*
} #FFT method
for (i in 1:simtim){
  X_check[i,] <- compute_X_n(n, l, H)
} #Pure X_n estimate
#time <- seq(0,T, length.out = n)

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
L = 300*300
K = 300
n <- L+1
lchoice <- 500000
l <- ifelse(lchoice>=n/2, lchoice, n/2)

H = 0.1
T= 1
#X <- simfbmonce(n,H,T)

X <- fbm_sim(n,l,H,T)


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
     xlab = "1/p",
     ylab = "W", yaxt = "n")
# Define where to place the ticks (powers of 10)
y_ticks <- c(10^-0.2, 10^-0.1, 10^0, 10^0.1, 10^0.2) #H=0,8
y_ticks <- c(10^-0.4, 10^-0.2, 10^0, 10^0.2, 10^0.4, 10^0.6) #H=0,5
y_ticks <- c(10^-0.5, 10^0, 10^0.5, 10^1) #H=0,3
y_ticks <- c(10^0, 10^5, 10^10, 10^15, 10^20, 10^25) #H=0,1

# Add y-axis with labels as powers of 10
axis(2, at = y_ticks, labels = expression(10^-0.2, 10^-0.1, 10^0, 10^0.1, 10^0.2)) #H=0,8
axis(2, at = y_ticks, labels = expression(10^-0.4, 10^-0.2, 10^0, 10^0.2, 10^0.4, 10^0.6)) #H=0,5
axis(2, at = y_ticks, labels = expression(10^-0.5, 10^0, 10^0.5, 10^1)) #H=0,3
axis(2, at = y_ticks, labels = expression(10^0, 10^5, 10^10, 10^15, 10^20, 10^25)) #H=0,1

abline(h = 1, col = "blue", lwd = 1, lty = 1)  #Estimating \hat{p}
abline(v = H, col = "black", lwd = 1, lty = 2)  # True H

# Find the index where log(W) crosses y = 0
crossing_index <- which.min(abs(W_values-1))
x_at_y_0 <- inv_p_values[crossing_index]
abline(v = x_at_y_0, col = "blue", lwd = 1, lty = 1)


## Histogram of estimated roughness index
L <- 2000*2000
K <- 2000
n <- L+1
lchoice <- 5000000
l <- ifelse(lchoice>=n/2, lchoice, n/2)

simtim <- 150 #Number of times to simulate the fBm (>1)
H = 0.1
T= 1

#X = simfbm(n, simtim, H, T) #Original simulation

X_check <- matrix(0,simtim, n)
for (i in 1:simtim){
  X_check[i,] <- Re(fft(compute_a_k(l, H))[1:n])
} #FFT method
Y <- matrix(0,simtim,n)
for (i in 2:n){
  Y[,i] <- Y[,i-1]+X_check[,i]
}
X <- Y*(T/n)^H

p_inv_sols <- numeric(simtim)
for (i in 1:simtim){
  objective_function <- function(p) {
    (Wstat(L = L, K = K, p, t=1 , X = X[i,]) - T)^2  # Squared difference to minimize
  }
  p_inv_sols[i] <- 1/optimize(objective_function, interval = c(1,15))$minimum
}
hist(p_inv_sols, main = paste("H=",H), xlab = "H=1/p", ylab = "Density")



##Estimated H plotted against different values of K
L <- 300*300

n <- L+1
lchoice <- 200000
l <- ifelse(lchoice>=n/2, lchoice, n/2)
H = 0.1
T= 1


test <- Re(fft(compute_a_k(l, H))[1:n])
X <- numeric(n)
X[1] <- test[1]
for (i in 2:n){
  X[i] <- X[i-1]+test[i]
}
X <- X*(T/n)^H #Simulated fBm (spectral method)


roughness_fct <- function(K){
  objective_function <- function(p) {
    (Wstat(L = L, K = K, p, t=1 , X = X) - T)^2  # Squared difference to minimize
  }
  return(1/optimize(objective_function, interval = c(1,20))$minimum)
}


K_values <- 2:500 
roughness_index <- sapply(K_values, function(K) roughness_fct(K))  # apply your function

# Create the plot
plot(K_values, roughness_index, type = "l",  # 'l' for lines
     xlab = "K", ylab = "Estimated Roughness Index", 
     main = "")

abline(h = roughness_fct(sqrt(L)), col = "blue", lwd = 1, lty = 1) 
abline(v = sqrt(L), col = "blue", lwd = 1, lty = 1)


