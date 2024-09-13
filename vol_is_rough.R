library(MASS)
library(tidyverse)

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

T <- 3500
delta <- 1/20
n <- T/delta
H <- 0.14
nu <- 0.3
m <- -5
alpha <- 5*10^(-4)
l <- n/2


#W <- simfbmonce(n,H,T)

test <- Re(fft(compute_a_k(l, H)))
X <- numeric(n)
X[1] <- test[1]
for (i in 2:n){
  X[i] <- X[i-1]+test[i]
}
W <- X*(T/n)^H

X <- numeric(n)
X[1] <- m
for (i in 2:n){
  X[i] = nu*(W[i]-W[i-1])+alpha*delta*(m-X[i-1])+X[i-1]
} #RSFV simulation
vol <- exp(X)

time <- seq(0,T, length.out = n)
plot(time, vol, type = "l", col = "blue", xlab = "", ylab = "")

U <- rnorm(T)
P <- numeric(T*1/delta)
P[1] <- 5
for (i in 1:T) {
  for (j in 1:(1/delta)) {
    idx <- (i - 1) * 1/delta + j  # Convert (n, delta) indices into a single index for P
    
    # Update P for the substep using the same U_n for all substeps in the partition
    P[idx + 1] <- P[idx] + P[idx] * X[idx] * sqrt(delta) * U[i]
  }
} #Price simulation

XDaily <- X[seq(1, length(X), by = 1/delta)]
#Estimating volatility
m <- function(q,Delta){
  volDelta <- exp(XDaily[seq(1, length(XDaily), by = Delta)])
  1/(T/Delta) *sum(abs(volDelta)^q)
}

q = 1
Delta_values <- 1:100
m_values <- sapply(Delta_values, function(Delta) m(q = q, Delta))
log_Delta <- log(Delta_values)
log_m <- log(m_values)

plot(log_Delta, log_m, xlab = "Log(Delta)", ylab = "Log(m(q,Delta))")



