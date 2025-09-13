 # Required libraries
  library(circular)  # For von Mises functions and datasets
library(openair)   # For mydata dataset
library(ggplot2)   # For advanced plots

# Section 1: Define functions for SSvM (proposed distribution)

# Density function for SSvM
dssvm <- function(theta, mu, kappa, eta) {
  if (abs(eta) >= 1) return(rep(0, length(theta)))
  (1 / (2 * pi * besselI(kappa, 0))) * exp(kappa * cos(theta - mu)) * (1 + eta * sin(theta - mu))
}

# Simulation from SSvM using rejection sampling
rssvm <- function(n, mu, kappa, eta) {
  if (abs(eta) >= 1) stop("|eta| must be < 1")
  M <- 1 + abs(eta)
  res <- numeric(n)
  i <- 1
  while (i <= n) {
    theta_prop <- rvonmises(1, mu = mu, kappa = kappa)
    u <- runif(1)
    if (u < (1 + eta * sin(theta_prop - mu)) / M) {
      res[i] <- theta_prop
      i <- i + 1
    }
  }
  res
}

# Log-likelihood for SSvM
loglik_ssvm <- function(par, data) {
  mu <- par[1]
  kappa <- par[2]
  eta <- par[3]
  ll <- sum(log(dssvm(data, mu, kappa, eta)))
  if (is.na(ll) || is.infinite(ll)) return(-1e10)
  ll
}

# MLE for SSvM
mle_ssvm <- function(data, start = c(mean.circular(data), 1, 0)) {
  optim(start, loglik_ssvm, data = data, method = "L-BFGS-B",
        lower = c(0, 0.001, -0.999), upper = c(2 * pi, 50, 0.999),
        control = list(fnscale = -1, maxit = 1000))$par
}

# Section 2: Define functions for AGvM (Kim et al.)

# Normalizing constant G0 for AGvM
g0_agvm <- function(k1, k2, maxj = 200) {
  k1 <- abs(k1)
  k2 <- abs(k2)  # Fix for negative k2
  i0k1 <- besselI(k1, 0)
  i0k2 <- besselI(k2, 0)
  sum_term <- 0
  for (j in 1:maxj) {
    i2j_k1 <- besselI(k1, 2 * j)
    ij_k2 <- besselI(k2, j)
    cos_term <- cos(j * pi / 2)
    sum_term <- sum_term + i2j_k1 * ij_k2 * cos_term
  }
  g0 <- i0k1 * i0k2 + 2 * sum_term
  if (is.na(g0) || g0 <= 0) g0 <- 1e-10  # Prevent NA or zero
  g0
}

# Density function for AGvM
dagvm <- function(theta, mu, k1, k2) {
  g0 <- g0_agvm(k1, k2)
  1 / (2 * pi * g0) * exp(k1 * cos(theta - mu) + k2 * sin(2 * (theta - mu)))
}

# Simulation from AGvM using rejection sampling
ragvm <- function(n, mu, k1, k2) {
  g0 <- g0_agvm(k1, k2)
  max_exp <- exp(abs(k1) + abs(k2))
  res <- numeric(n)
  i <- 1
  while (i <= n) {
    theta_prop <- runif(1, 0, 2 * pi)
    target_exp <- exp(k1 * cos(theta_prop - mu) + k2 * sin(2 * (theta_prop - mu)))
    u <- runif(1)
    if (u < target_exp / max_exp) {
      res[i] <- theta_prop
      i <- i + 1
    }
  }
  res
}

# Log-likelihood for AGvM
loglik_agvm <- function(par, data) {
  mu <- par[1]
  k1 <- par[2]
  k2 <- par[3]
  ll <- sum(log(dagvm(data, mu, k1, k2)))
  if (is.na(ll) || is.infinite(ll)) return(-1e10)
  ll
}

# MLE for AGvM
mle_agvm <- function(data, start = c(mean.circular(data), 1, 0)) {
  optim(start, loglik_agvm, data = data, method = "L-BFGS-B",
        lower = c(0, 0.001, -0.999), upper = c(2 * pi, 50, 0.999),
        control = list(fnscale = -1, maxit = 1000))$par
}

# Function to compute AIC and BIC
compute_criteria <- function(loglik, p, n) {
  aic <- -2 * loglik + 2 * p
  bic <- -2 * loglik + log(n) * p
  list(AIC = aic, BIC = bic)
}

# Section 3: Simulation Study

# 3.1 Simulate from SSvM and fit both models
set.seed(123)
n_sim <- 200
true_mu_ssvm <- pi
true_kappa_ssvm <- 5
true_eta_ssvm <- 0.4
sim_data_ssvm <- rssvm(n_sim, true_mu_ssvm, true_kappa_ssvm, true_eta_ssvm)

# Fit SSvM
fit_ssvm_to_ssvm <- mle_ssvm(sim_data_ssvm)
loglik_ssvm_to_ssvm <- loglik_ssvm(fit_ssvm_to_ssvm, sim_data_ssvm)
criteria_ssvm_to_ssvm <- compute_criteria(loglik_ssvm_to_ssvm, 3, n_sim)

# Fit AGvM
fit_agvm_to_ssvm <- mle_agvm(sim_data_ssvm)
loglik_agvm_to_ssvm <- loglik_agvm(fit_agvm_to_ssvm, sim_data_ssvm)
criteria_agvm_to_ssvm <- compute_criteria(loglik_agvm_to_ssvm, 3, n_sim)

# Table for simulation from SSvM
sim_ssvm_table <- data.frame(
  Model = c("SSvM", "AGvM"),
  mu = c(fit_ssvm_to_ssvm[1], fit_agvm_to_ssvm[1]),
  kappa_k1 = c(fit_ssvm_to_ssvm[2], fit_agvm_to_ssvm[2]),
  eta_k2 = c(fit_ssvm_to_ssvm[3], fit_agvm_to_ssvm[3]),
  LogLik = c(loglik_ssvm_to_ssvm, loglik_agvm_to_ssvm),
  AIC = c(criteria_ssvm_to_ssvm$AIC, criteria_agvm_to_ssvm$AIC),
  BIC = c(criteria_ssvm_to_ssvm$BIC, criteria_agvm_to_ssvm$BIC)
)
print("Simulation from SSvM: Fit Comparison")
print(sim_ssvm_table)

# 3.2 Simulate from AGvM and fit both models
true_mu_agvm <- pi
true_k1_agvm <- 0.5
true_k2_agvm <- 0.3
sim_data_agvm <- ragvm(n_sim, true_mu_agvm, true_k1_agvm, true_k2_agvm)

# Fit AGvM
fit_agvm_to_agvm <- mle_agvm(sim_data_agvm)
loglik_agvm_to_agvm <- loglik_agvm(fit_agvm_to_agvm, sim_data_agvm)
criteria_agvm_to_agvm <- compute_criteria(loglik_agvm_to_agvm, 3, n_sim)

# Fit SSvM
fit_ssvm_to_agvm <- mle_ssvm(sim_data_agvm)
loglik_ssvm_to_agvm <- loglik_ssvm(fit_ssvm_to_agvm, sim_data_agvm)
criteria_ssvm_to_agvm <- compute_criteria(loglik_ssvm_to_agvm, 3, n_sim)

# Table for simulation from AGvM
sim_agvm_table <- data.frame(
  Model = c("AGvM", "SSvM"),
  mu = c(fit_agvm_to_agvm[1], fit_ssvm_to_agvm[1]),
  k1_kappa = c(fit_agvm_to_agvm[2], fit_ssvm_to_agvm[2]),
  k2_eta = c(fit_agvm_to_agvm[3], fit_ssvm_to_agvm[3]),
  LogLik = c(loglik_agvm_to_agvm, loglik_ssvm_to_agvm),
  AIC = c(criteria_agvm_to_agvm$AIC, criteria_ssvm_to_agvm$AIC),
  BIC = c(criteria_agvm_to_agvm$BIC, criteria_ssvm_to_agvm$BIC)
)
print("Simulation from AGvM: Fit Comparison")
print(sim_agvm_table)

# Simulation Plots: Density overlays on histograms (linear and polar)
# For SSvM sim data
df_sim_ssvm <- data.frame(theta = sim_data_ssvm)
p_hist_ssvm <- ggplot(df_sim_ssvm, aes(x = theta)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", alpha = 0.5) +
  stat_function(fun = function(x) dssvm(x, fit_ssvm_to_ssvm[1], fit_ssvm_to_ssvm[2], fit_ssvm_to_ssvm[3]), color = "blue", size = 1) +
  stat_function(fun = function(x) dagvm(x, fit_agvm_to_ssvm[1], fit_agvm_to_ssvm[2], fit_agvm_to_ssvm[3]), color = "red", size = 1) +
  labs(title = "Histogram with Densities (SSvM Sim Data)", subtitle = "Blue: SSvM, Red: AGvM") +
  theme_minimal()

p_polar_ssvm <- ggplot() +
  stat_function(fun = function(x) dssvm(x, fit_ssvm_to_ssvm[1], fit_ssvm_to_ssvm[2], fit_ssvm_to_ssvm[3]), geom = "line", color = "blue", xlim = c(0, 2*pi)) +
  stat_function(fun = function(x) dagvm(x, fit_agvm_to_ssvm[1], fit_agvm_to_ssvm[2], fit_agvm_to_ssvm[3]), geom = "line", color = "red", xlim = c(0, 2*pi)) +
  coord_polar() + theme_minimal() + labs(title = "Polar Density Plot (SSvM Sim Data)")

print(p_hist_ssvm)
print(p_polar_ssvm)

# Section 4: Real Data Analysis

# 4.1 mydata from openair (wind direction)
data(mydata)
theta_mydata <- mydata$wd * pi / 180
theta_mydata <- theta_mydata[!is.na(theta_mydata)]

# Fit SSvM
fit_ssvm_mydata <- mle_ssvm(theta_mydata)
loglik_ssvm_mydata <- loglik_ssvm(fit_ssvm_mydata, theta_mydata)
criteria_ssvm_mydata <- compute_criteria(loglik_ssvm_mydata, 3, length(theta_mydata))

# Fit AGvM
fit_agvm_mydata <- mle_agvm(theta_mydata)
loglik_agvm_mydata <- loglik_agvm(fit_agvm_mydata, theta_mydata)
criteria_agvm_mydata <- compute_criteria(loglik_agvm_mydata, 3, length(theta_mydata))

# Table for mydata
mydata_table <- data.frame(
  Model = c("SSvM", "AGvM"),
  mu = c(fit_ssvm_mydata[1], fit_agvm_mydata[1]),
  kappa_k1 = c(fit_ssvm_mydata[2], fit_agvm_mydata[2]),
  eta_k2 = c(fit_ssvm_mydata[3], fit_agvm_mydata[3]),
  LogLik = c(loglik_ssvm_mydata, loglik_agvm_mydata),
  AIC = c(criteria_ssvm_mydata$AIC, criteria_agvm_mydata$AIC),
  BIC = c(criteria_ssvm_mydata$BIC, criteria_agvm_mydata$BIC)
)
print("mydata (openair) Fit Comparison")
print(mydata_table)

# 4.2 fisherB2 from circular
data(fisherB2)
theta_fisherB2 <- as.numeric(fisherB2)  # Ensure numeric

# Fit SSvM
fit_ssvm_fisher <- mle_ssvm(theta_fisherB2)
loglik_ssvm_fisher <- loglik_ssvm(fit_ssvm_fisher, theta_fisherB2)
criteria_ssvm_fisher <- compute_criteria(loglik_ssvm_fisher, 3, length(theta_fisherB2))

# Fit AGvM
fit_agvm_fisher <- mle_agvm(theta_fisherB2)
loglik_agvm_fisher <- loglik_agvm(fit_agvm_fisher, theta_fisherB2)
criteria_agvm_fisher <- compute_criteria(loglik_agvm_fisher, 3, length(theta_fisherB2))

# Table for fisherB2
fisher_table <- data.frame(
  Model = c("SSvM", "AGvM"),
  mu = c(fit_ssvm_fisher[1], fit_agvm_fisher[1]),
  kappa_k1 = c(fit_ssvm_fisher[2], fit_agvm_fisher[2]),
  eta_k2 = c(fit_ssvm_fisher[3], fit_agvm_fisher[3]),
  LogLik = c(loglik_ssvm_fisher, loglik_agvm_fisher),
  AIC = c(criteria_ssvm_fisher$AIC, criteria_agvm_fisher$AIC),
  BIC = c(criteria_ssvm_fisher$BIC, criteria_agvm_fisher$BIC)
)
print("fisherB2 (circular) Fit Comparison")
print(fisher_table)

# Real Data Plots: Similar to simulation (example for mydata)
df_mydata <- data.frame(theta = theta_mydata)
p_hist_mydata <- ggplot(df_mydata, aes(x = theta)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "lightgreen", alpha = 0.5) +
  stat_function(fun = function(x) dssvm(x, fit_ssvm_mydata[1], fit_ssvm_mydata[2], fit_ssvm_mydata[3]), color = "blue", size = 1) +
  stat_function(fun = function(x) dagvm(x, fit_agvm_mydata[1], fit_agvm_mydata[2], fit_agvm_mydata[3]), color = "red", size = 1) +
  labs(title = "Histogram with Densities (mydata)", subtitle = "Blue: SSvM, Red: AGvM") +
  theme_minimal()

p_polar_mydata <- ggplot() +
  stat_function(fun = function(x) dssvm(x, fit_ssvm_mydata[1], fit_ssvm_mydata[2], fit_ssvm_mydata[3]), geom = "line", color = "blue") +
  stat_function(fun = function(x) dagvm(x, fit_agvm_mydata[1], fit_agvm_mydata[2], fit_agvm_mydata[3]), geom = "line", color = "red") +
  coord_polar() + theme_minimal() + labs(title = "Polar Density Plot (mydata)")

print(p_hist_mydata)
print(p_polar_mydata)

# Repeat similar plots for fisherB2
df_fisher <- data.frame(theta = theta_fisherB2)
p_hist_fisher <- ggplot(df_fisher, aes(x = theta)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "lightyellow", alpha = 0.5) +
  stat_function(fun = function(x) dssvm(x, fit_ssvm_fisher[1], fit_ssvm_fisher[2], fit_ssvm_fisher[3]), color = "blue", size = 1) +
  stat_function(fun = function(x) dagvm(x, fit_agvm_fisher[1], fit_agvm_fisher[2], fit_agvm_fisher[3]), color = "red", size = 1) +
  labs(title = "Histogram with Densities (fisherB2)", subtitle = "Blue: SSvM, Red: AGvM") +
  theme_minimal()

p_polar_fisher <- ggplot() +
  stat_function(fun = function(x) dssvm(x, fit_ssvm_fisher[1], fit_ssvm_fisher[2], fit_ssvm_fisher[3]), geom = "line", color = "blue") +
  stat_function(fun = function(x) dagvm(x, fit_agvm_fisher[1], fit_agvm_fisher[2], fit_agvm_fisher[3]), geom = "line", color = "red") +
  coord_polar() + theme_minimal() + labs(title = "Polar Density Plot (fisherB2)")

print(p_hist_fisher)
print(p_polar_fisher)
# Section 1: Plots for Simulated SSvM Data
# Note: The data frame `df_sim_ssvm` and model fits `fit_ssvm_to_ssvm`, `fit_agvm_to_ssvm`
# are assumed to be generated by the preceding code.
p_polar_ssvm <- ggplot() +
  stat_function(fun = function(x) dssvm(x, fit_ssvm_to_ssvm[1], fit_ssvm_to_ssvm[2], fit_ssvm_to_ssvm[3]), geom = "line", color = "blue", xlim = c(0, 2*pi)) +
  stat_function(fun = function(x) dagvm(x, fit_agvm_to_ssvm[1], fit_agvm_to_ssvm[2], fit_agvm_to_ssvm[3]), geom = "line", color = "red", xlim = c(0, 2*pi)) +
  coord_polar() + 
  theme_minimal() + 
  labs(title = "Polar Density Plot (SSvM Sim Data)")
print(p_polar_ssvm)

# Section 2: Plots for Simulated AGvM Data
# Note: The data frame `df_sim_agvm` and model fits `fit_ssvm_to_agvm`, `fit_agvm_to_agvm`
# are assumed to be generated by the preceding code.
df_sim_agvm <- data.frame(theta = sim_data_agvm)
p_polar_agvm <- ggplot() +
  stat_function(fun = function(x) dagvm(x, fit_agvm_to_agvm[1], fit_agvm_to_agvm[2], fit_agvm_to_agvm[3]), geom = "line", color = "red", xlim = c(0, 2*pi)) +
  stat_function(fun = function(x) dssvm(x, fit_ssvm_to_agvm[1], fit_ssvm_to_agvm[2], fit_ssvm_to_agvm[3]), geom = "line", color = "blue", xlim = c(0, 2*pi)) +
  coord_polar() + 
  theme_minimal() + 
  labs(title = "Polar Density Plot (AGvM Sim Data)")
print(p_polar_agvm)

# Section 3: Plots for `mydata`
# Note: The data frame `df_mydata` and model fits `fit_ssvm_mydata`, `fit_agvm_mydata`
# are assumed to be generated by the preceding code.
p_polar_mydata <- ggplot() +
  stat_function(fun = function(x) dssvm(x, fit_ssvm_mydata[1], fit_ssvm_mydata[2], fit_ssvm_mydata[3]), geom = "line", color = "blue") +
  stat_function(fun = function(x) dagvm(x, fit_agvm_mydata[1], fit_agvm_mydata[2], fit_agvm_mydata[3]), geom = "line", color = "red") +
  coord_polar() + 
  theme_minimal() + 
  labs(title = "Polar Density Plot (mydata)")
print(p_polar_mydata)

# Section 4: Plots for `fisherB2` data
# Note: The data frame `df_fisher` and model fits `fit_ssvm_fisher`, `fit_agvm_fisher`
# are assumed to be generated by the preceding code.
p_polar_fisher <- ggplot() +
  stat_function(fun = function(x) dssvm(x, fit_ssvm_fisher[1], fit_ssvm_fisher[2], fit_ssvm_fisher[3]), geom = "line", color = "blue") +
  stat_function(fun = function(x) dagvm(x, fit_agvm_fisher[1], fit_agvm_fisher[2], fit_agvm_fisher[3]), geom = "line", color = "red") +
  coord_polar() + 
  theme_minimal() + 
  labs(title = "Polar Density Plot (fisherB2)")
print(p_polar_fisher)
