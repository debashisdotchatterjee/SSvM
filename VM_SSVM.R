################################################################################
# Fit VM and Sine-Skewed VM (SSVM) to circular wind-direction data
# Produces tables, comparative metrics, and colourful circular plots
#
# Assumes: mydata$wd exists (degrees). If not, replace with your wind vector.
################################################################################

# --- Install / load required packages ----------------------------------------
pkgs <- c("circular","numDeriv","ggplot2","cowplot","dplyr","gridExtra","scales")
for(p in pkgs) if(!requireNamespace(p, quietly=TRUE)) install.packages(p)
library(circular); library(numDeriv); library(ggplot2); library(cowplot)
library(dplyr); library(gridExtra); library(scales)
library(openair)

# --- Data: load your wind directions (degrees) -------------------------------
# Replace this if your data is elsewhere
if(!exists("mydata")) stop("Please load your data.frame 'mydata' with column 'wd' (wind direction in degrees).")
wd_deg <- mydata$wd
wd_deg <- wd_deg[!is.na(wd_deg)]
if(length(wd_deg) < 10) stop("Need at least ~10 observations for stable estimation.")

# convert to radians in (-pi, pi]
theta <- (wd_deg %% 360) * pi/180
theta <- ifelse(theta > pi, theta - 2*pi, theta)
n <- length(theta)

# --- Helper: circular sample mean & resultant -------------------------------
circ_summary <- function(t){
  C <- mean(cos(t)); S <- mean(sin(t))
  mu <- atan2(S, C)
  Rbar <- sqrt(C^2 + S^2)
  list(mu = mu, Rbar = Rbar, C = C, S = S)
}
cs <- circ_summary(theta)

# --- 1) Fit von Mises: moment-based (mu = atan2, solve for kappa) -----------
A1 <- function(k) besselI(k,1)/besselI(k,0)
A1inv_approx <- function(R){
  if(R < 0.53) return(2*R + R^3 + (5*R^5)/6)
  if(R < 0.85) return(-0.4 + 1.39*R + 0.43/(1-R))
  return(1/(R^3 - 4*R^2 + 3*R))
}
mu_vm <- cs$mu
Rbar <- cs$Rbar
kappa_vm <- A1inv_approx(Rbar)
# Newton refine
for(i in 1:10){
  dA <- (A1(kappa_vm + 1e-6) - A1(kappa_vm - 1e-6)) / (2e-6)
  if(dA==0) break
  kappa_vm <- max(1e-8, kappa_vm - (A1(kappa_vm) - Rbar)/dA)
}
logLik_vm <- function(mu, kappa, t){
  n <- length(t)
  n * (-log(2*pi) - log(besselI(kappa,0))) + kappa * sum(cos(t - mu))
}
ll_vm <- logLik_vm(mu_vm, kappa_vm, theta)
aic_vm <- -2*ll_vm + 2*2
bic_vm <- -2*ll_vm + log(n)*2

# --- 2) Fit SSVM by MLE (parameter transforms + penalty) --------------------
# params: mu (real mapped to (-pi,pi]), logkappa (real->kappa>0), eta_raw (real->eta=tanh)
neg_loglik_ssvm_raw <- function(par, t){
  mu_raw <- par[1]; logk <- par[2]; eta_raw <- par[3]
  mu <- ((mu_raw + pi) %% (2*pi)) - pi
  kappa <- exp(logk)
  eta <- tanh(eta_raw)
  vals <- 1 + eta * sin(t - mu)
  # penalize if any vals <= 0 (density must be positive for observed t)
  if(any(vals <= 0)) return(1e12 + sum(pmax(0, -vals))*1e6)
  ll <- length(t) * (-log(2*pi) - log(besselI(kappa,0))) + kappa * sum(cos(t - mu)) + sum(log(vals))
  -ll
}

# start from VM estimates, small eta
start_raw <- c(mu_vm, log(max(1e-8,kappa_vm)), atanh(0.05))
opt <- optim(par = start_raw, fn = neg_loglik_ssvm_raw, t = theta, method = "BFGS",
             control = list(maxit = 5000), hessian = TRUE)
if(opt$convergence != 0) warning("optim warning: convergence code ", opt$convergence)

par_hat_raw <- opt$par
mu_ssvm <- ((par_hat_raw[1] + pi) %% (2*pi)) - pi
kappa_ssvm <- exp(par_hat_raw[2])
eta_ssvm <- tanh(par_hat_raw[3])
ll_ssvm <- -neg_loglik_ssvm_raw(par_hat_raw, theta)
aic_ssvm <- -2*ll_ssvm + 2*3
bic_ssvm <- -2*ll_ssvm + log(n)*3

# Delta method SEs from Hessian (raw params -> transform)
hess <- opt$hessian
cov_raw <- try(solve(hess), silent = TRUE)
if(inherits(cov_raw, "try-error")) {
  se_mu_raw <- se_logk_raw <- se_eta_raw <- NA
  se_mu <- se_kappa <- se_eta <- NA
} else {
  se_raw <- sqrt(diag(cov_raw))
  se_mu_raw <- se_raw[1]
  se_logk_raw <- se_raw[2]
  se_eta_raw <- se_raw[3]
  # delta: mu = identity
  se_mu <- se_mu_raw
  # kappa = exp(logk) -> se(kappa) = kappa * se(logk)
  se_kappa <- abs(kappa_ssvm) * se_logk_raw
  # eta = tanh(eta_raw) -> derivative = 1 - tanh^2
  se_eta <- abs((1 - tanh(par_hat_raw[3])^2) * se_eta_raw)
}

# --- 3) Bootstrap for parameter CIs (nonparametric & parametric) ----------
# parametric bootstrap: simulate from fitted models and re-fit, used for CI and LRT p-value
nboot <- 500   # lower default for speed; increase to 2000 for publication
set.seed(2025)
boot_par_vm <- matrix(NA, nrow = nboot, ncol = 2)      # mu,kappa
boot_par_ssvm <- matrix(NA, nrow = nboot, ncol = 3)    # mu,kappa,eta
lr_stats <- numeric(nboot)

# helper: simulate from VM
rvm <- function(n, mu, kappa){
  # use circular::rvonmises for stable sampling
  circular::rvonmises(n, mu = circular(mu), kappa = kappa)
}

# helper: simulate from SSVM via rejection sampling:
rssvm <- function(n, mu, kappa, eta){
  # target density f(theta) = base_vm * (1 + eta sin(theta-mu))
  # use envelope = c * vm density with c chosen: since |eta|<1, 1+eta sin in [1-|eta|, 1+|eta|]
  # maximum of modifier is 1+|eta| ; use rejection with vm proposals at same mu,kappa
  out <- numeric(0); tries <- 0
  M <- 1 + abs(eta)
  while(length(out) < n){
    prop <- as.numeric(rvm(n*2, mu, kappa))
    u <- runif(length(prop))
    accept <- u < (1 + eta * sin(prop - mu)) / M
    out <- c(out, prop[accept])
    tries <- tries + length(prop)
    if(tries > 1e6) stop("Rejection sampling failed")
  }
  out[1:n]
}

# Boot loop (may take time)
pb <- txtProgressBar(min=0, max=nboot, style=3)
for(b in 1:nboot){
  # simulate under VM to get bootstrap LRT dist (for parametric boot under null),
  # and also re-fit for parametric CIs under each model.
  sim_vm <- as.numeric(rvm(n, mu_vm, kappa_vm))
  # fit vm on sim data
  cs_sim <- circ_summary(sim_vm)
  mu_vm_sim <- cs_sim$mu
  Rbar_sim <- cs_sim$Rbar
  kappa_vm_sim <- A1inv_approx(Rbar_sim)
  for(i in 1:5){
    dA <- (A1(kappa_vm_sim + 1e-6) - A1(kappa_vm_sim - 1e-6)) / (2e-6)
    if(dA==0) break
    kappa_vm_sim <- max(1e-8, kappa_vm_sim - (A1(kappa_vm_sim) - Rbar_sim)/dA)
  }
  boot_par_vm[b,] <- c(mu_vm_sim, kappa_vm_sim)
  
  # re-fit SSVM to sim_vm (may have issues if 1+eta sin <=0; fit starting near 0)
  start_b <- c(mu_vm_sim, log(max(1e-8,kappa_vm_sim)), atanh(0.01))
  opt_b <- try(optim(start_b, neg_loglik_ssvm_raw, t = sim_vm, method = "BFGS", control=list(maxit=2000)), silent=TRUE)
  if(!inherits(opt_b, "try-error") && opt_b$convergence==0){
    parb <- opt_b$par
    mu_b <- ((parb[1]+pi) %% (2*pi)) - pi
    kappa_b <- exp(parb[2])
    eta_b <- tanh(parb[3])
    boot_par_ssvm[b,] <- c(mu_b, kappa_b, eta_b)
    # compute LR stat on this simulated sample
    ll_vm_b <- logLik_vm(mu_vm_sim, kappa_vm_sim, sim_vm)
    ll_ssvm_b <- -neg_loglik_ssvm_raw(parb, sim_vm)
    lr_stats[b] <- 2*(ll_ssvm_b - ll_vm_b)
  } else {
    boot_par_ssvm[b,] <- NA
    lr_stats[b] <- NA
  }
  setTxtProgressBar(pb, b)
}
close(pb)

# Remove failed boot rows
ok_ssvm <- apply(boot_par_ssvm, 1, function(x) all(!is.na(x)))
boot_par_ssvm <- boot_par_ssvm[ok_ssvm, , drop=FALSE]
lr_stats <- lr_stats[ok_ssvm]
boot_par_vm <- boot_par_vm[1:length(lr_stats), , drop=FALSE]  # align lengths

# Parametric bootstrap p-value for LRT: compare observed LR to bootstrap LR distribution
lr_obs <- 2*(ll_ssvm - ll_vm)
pval_lr <- mean(lr_stats >= lr_obs, na.rm = TRUE)

# parametric CIs: percentiles
ci_vm_mu <- quantile(boot_par_vm[,1], c(0.025, 0.975))
ci_vm_kappa <- quantile(boot_par_vm[,2], c(0.025, 0.975))
ci_ssvm_mu <- quantile(boot_par_ssvm[,1], c(0.025, 0.975))
ci_ssvm_kappa <- quantile(boot_par_ssvm[,2], c(0.025, 0.975))
ci_ssvm_eta <- quantile(boot_par_ssvm[,3], c(0.025, 0.975))

# --- 4) Goodness-of-fit tests (Rayleigh for unimodal orientation) ----------
# Rayleigh test for non-uniformity:
ray <- suppression <- try(circular::rayleigh.test(circular(theta)), silent=TRUE)
ray_p <- if(inherits(ray,"try-error")) NA else ray$p.value

# Note: for Kuiper or Watson tests for specific fitted distributions, packages vary.
# We'll compute Kuiper statistic comparing empirical CDF to fitted CDF for VM and SSVM
ecdf_circ <- function(t){
  # returns empirical circular CDF at t in (-pi,pi]
  # use angles sorted and empirical cumulative counts
  t_sorted <- sort(theta)
  function(x){
    # map x to (-pi,pi]
    x2 <- ifelse(x > pi, x - 2*pi, x)
    mean( ( (t_sorted <= x2) ) )
  }
}
# theoretical circular CDF for vm and ssvm (numerical integration)
cdf_vm <- function(x, mu, kappa){
  # integrate density from -pi to x
  sapply(x, function(xx) integrate(function(z) dvm(z, mu, kappa), lower = -pi, upper = xx, rel.tol=1e-8)$value)
}
dvm <- function(t, mu, kappa) (1/(2*pi*besselI(kappa,0))) * exp(kappa * cos(t - mu))
dssvm <- function(t, mu, kappa, eta) dvm(t, mu, kappa) * (1 + eta * sin(t - mu))

# compute Kuiper-like statistic (max diffs)
angles_sorted <- sort(theta)
ecdf_vals <- sapply(angles_sorted, function(x) mean(theta <= x))
# theoretical CDFs at these angles
cdf_vm_vals <- sapply(angles_sorted, function(x) integrate(function(z) dvm(z, mu_vm, kappa_vm), lower = -pi, upper = x)$value)
cdf_ssvm_vals <- sapply(angles_sorted, function(x) integrate(function(z) dssvm(z, mu_ssvm, kappa_ssvm, eta_ssvm), lower = -pi, upper = x)$value)

Kuiper_vm <- (max(ecdf_vals - cdf_vm_vals) + max(cdf_vm_vals - ecdf_vals))
Kuiper_ssvm <- (max(ecdf_vals - cdf_ssvm_vals) + max(cdf_ssvm_vals - ecdf_vals))

# --- 5) Prepare tables ------------------------------------------------------
param_table <- data.frame(
  Model = c("Von Mises", "Sine-skewed VM"),
  mu_rad = c(mu_vm, mu_ssvm),
  mu_deg = c(mu_vm*180/pi, mu_ssvm*180/pi),
  kappa = c(kappa_vm, kappa_ssvm),
  eta = c(NA, eta_ssvm),
  logLik = c(ll_vm, ll_ssvm),
  AIC = c(aic_vm, aic_ssvm),
  BIC = c(bic_vm, bic_ssvm),
  stringsAsFactors = FALSE
)

se_table <- data.frame(
  Param = c("mu_rad","kappa","eta"),
  VM_est = c(mu_vm, kappa_vm, NA),
  VM_SE  = c(NA, NA, NA),
  SSVM_est = c(mu_ssvm, kappa_ssvm, eta_ssvm),
  SSVM_SE  = c(se_mu, se_kappa, se_eta),
  stringsAsFactors = FALSE
)

ci_table <- data.frame(
  Param = c("mu_deg","kappa","eta"),
  VM_L = c(ci_vm_mu[1]*180/pi, ci_vm_kappa[1], NA),
  VM_U = c(ci_vm_mu[2]*180/pi, ci_vm_kappa[2], NA),
  SSVM_L = c(ci_ssvm_mu[1]*180/pi, ci_ssvm_kappa[1], ci_ssvm_eta[1]),
  SSVM_U = c(ci_ssvm_mu[2]*180/pi, ci_ssvm_kappa[2], ci_ssvm_eta[2]),
  stringsAsFactors = FALSE
)

comp_metrics <- data.frame(
  Metric = c("logLik","AIC","BIC","LR_stat","LR_pval_parametric_boot","Kuiper_stat"),
  VM = c(ll_vm, aic_vm, bic_vm, NA, NA, Kuiper_vm),
  SSVM = c(ll_ssvm, aic_ssvm, bic_ssvm, lr_obs, pval_lr, Kuiper_ssvm)
)

# --- 6) Plots: colourful circular plots -------------------------------------
# Prepare data.frame for ggplot (degrees, 0-360)
plot_df <- data.frame(angle_deg = (theta * 180/pi) %% 360)
# circular histogram (rose) with ggplot2
rose <- ggplot(plot_df, aes(x = angle_deg)) +
  geom_histogram(aes(y = ..count..), breaks = seq(0,360,by=10), fill = "#1f77b4", color = "black", alpha = 0.7) +
  coord_polar(start = -pi/2) +
  scale_x_continuous(breaks = seq(0,330, by = 30)) +
  labs(title = "Rose diagram (10° bins)", x = "", y = "Counts") +
  theme_minimal()

plot(rose)
# density overlay: compute fitted densities on fine grid
grid_deg <- seq(0, 360, length.out = 720)
grid_rad <- (grid_deg/180)*pi
vm_dens_grid <- dvm(grid_rad, mu_vm, kappa_vm)
ssvm_dens_grid <- dssvm(grid_rad, mu_ssvm, kappa_ssvm, eta_ssvm)
dens_df <- data.frame(deg = grid_deg, vm = vm_dens_grid, ssvm = ssvm_dens_grid)
# scale densities to match histogram's max for overlay (visual only)
max_hist <- max(ggplot_build(ggplot(plot_df, aes(x=angle_deg)) + geom_histogram(breaks=seq(0,360,by=10)))$data[[1]]$count)
scale_factor <- max_hist / max(c(dens_df$vm, dens_df$ssvm)) * 0.9


polar_overlay <- ggplot() +
  geom_histogram(data=plot_df, aes(x = angle_deg, y = ..count..),
                 breaks = seq(0,360,by=10), fill = alpha("#1f77b4",0.5), color="black") +
  geom_line(data = dens_df, aes(x=deg, y=vm*scale_factor), color = "#ff7f0e", size = 1.2) +
  geom_line(data = dens_df, aes(x=deg, y=ssvm*scale_factor), color = "#2ca02c", size = 1.2, linetype = "dashed") +
  coord_polar(start = -pi/2) +
  labs(title = "Wind directions: histogram + VM (orange) and SSVM (green dashed)",
       x = "", y = "Counts (scaled density)") +
  theme_minimal()

plot(polar_overlay)

# kernel circular density (using circular::density.circular)
dens_circ <- density.circular(circular(theta), bw = 20, kernel = "vonmises")
# convert to data.frame degrees for ggplot
kd_deg <- (dens_circ$x * 180/pi) %% 360
kd_df <- data.frame(deg = kd_deg, y = dens_circ$y)

kernel_plot <- ggplot() +
  geom_line(data = kd_df, aes(x=deg, y=y), size=1.2, color="#1f77b4") +
  geom_line(data = dens_df, aes(x=deg, y=vm), color = "#ff7f0e", size=1.1) +
  geom_line(data = dens_df, aes(x=deg, y=ssvm), color = "#2ca02c", linetype="dashed", size=1.1) +
  labs(title = "Kernel circular density (blue) vs fitted VM (orange) & SSVM (green dashed)",
       x = "Angle (deg)", y = "Density") +
  theme_minimal()
plot(kernel_plot)

# diagnostic: 1 + eta sin(theta-mu) across angle
diag_df <- data.frame(rad = grid_rad,
                      modifier_vm = 1 + 0 * sin(grid_rad - mu_vm),
                      modifier_ssvm = 1 + eta_ssvm * sin(grid_rad - mu_ssvm))
diag_plot <- ggplot(diag_df, aes(x = rad*180/pi)) +
  geom_line(aes(y = modifier_ssvm), color="#2ca02c", size=1.2) +
  geom_hline(yintercept = 0, linetype="dashed", color="red") +
  labs(title = "Modifier (1 + η sin(θ-μ)) for SSVM", x = "Angle (deg)", y = "modifier") +
  theme_minimal()
plot(diag_plot)

# QQ-like: empirical CDF vs fitted CDF
ecdf_df <- data.frame(angle = (angles_sorted*180/pi) %% 360,
                      ecdf = ecdf_vals,
                      cdf_vm = cdf_vm_vals,
                      cdf_ssvm = cdf_ssvm_vals)
qq_plot <- ggplot(ecdf_df, aes(x = angle)) +
  geom_line(aes(y = ecdf, color = "Empirical")) +
  geom_line(aes(y = cdf_vm, color = "VM")) +
  geom_line(aes(y = cdf_ssvm, color = "SSVM")) +
  labs(title = "Empirical CDF vs Fitted CDFs (angles sorted)", x = "Angle (deg)", y = "CDF") +
  scale_color_manual(values = c("Empirical" = "#1f77b4", "VM" = "#ff7f0e", "SSVM" = "#2ca02c")) +
  theme_minimal()
plot(qq_plot)

# Arrange plots
top_row <- plot_grid(rose, polar_overlay, ncol=2, rel_widths = c(1,1))
bottom_row <- plot_grid(kernel_plot, diag_plot, qq_plot, ncol=3, rel_widths = c(1,1,1))
all_plots <- plot_grid(top_row, bottom_row, ncol=1, rel_heights = c(1,0.9))
print(top_row)
print(bottom_row)
print(all_plots)

# --- 7) Print Tables and Summary --------------------------------------------
cat("\n--- Parameter estimates ---\n")
print(param_table)

cat("\n--- SE table (approx) ---\n")
print(se_table)

cat("\n--- Bootstrap CIs (parametric) ---\n")
print(ci_table)

cat("\n--- Comparative metrics ---\n")
print(comp_metrics)
cat(sprintf("\nRayleigh test p-value (non-uniformity): %.4g\n", ray_p))
cat(sprintf("LR obs = %.4f, parametric bootstrap p-value = %.4g\n", lr_obs, pval_lr))
cat(sprintf("Kuiper statistic VM = %.4g, SSVM = %.4g (smaller is better fit to empirical CDF)\n",
            Kuiper_vm, Kuiper_ssvm))


################################################################################
# End of script
################################################################################
