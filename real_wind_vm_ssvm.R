################################################################################
# Fit von Mises (vM) and Sine-Skewed von Mises (SSvM) to circular wind direction
# - Accepts mydata$wd (degrees). If absent, loads openair::mydata.
# - Outputs high-quality colourful plots and CSV tables.
################################################################################

# --- Install / load required packages ----------------------------------------
pkgs <- c(
  "circular","numDeriv","ggplot2","cowplot","dplyr","gridExtra",
  "scales","readr","tibble","viridis","stringr","openair"
)
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
lapply(pkgs, require, character.only = TRUE)
theme_set(theme_minimal(base_size = 12))

# --- Output directories -------------------------------------------------------
root_out <- file.path(getwd(), "paper_outputs", "real_wind")
fig_dir  <- file.path(root_out, "figures")
tab_dir  <- file.path(root_out, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
cat("Outputs:\n  Figures:", fig_dir, "\n  Tables :", tab_dir, "\n\n")

# --- Data: load wind directions (degrees) ------------------------------------
if (!exists("mydata")) {
  message("No 'mydata' in workspace. Loading openair::mydata …")
  data("mydata", package = "openair")
}
stopifnot("wd" %in% names(mydata))
wd_deg <- as.numeric(mydata$wd)
wd_deg <- wd_deg[is.finite(wd_deg)]
if (length(wd_deg) < 30) stop("Need at least ~30 observations for stable estimation.")

# Convert to radians in (-pi, pi]
theta <- (wd_deg %% 360) * pi/180
theta <- ifelse(theta > pi, theta - 2*pi, theta)
n <- length(theta)
cat("Sample size n =", n, " (radians in (-pi,pi])\n\n")

# --- Helpers: circular sample mean & resultant -------------------------------
circ_summary <- function(t){
  C <- mean(cos(t)); S <- mean(sin(t))
  mu <- atan2(S, C)
  Rbar <- sqrt(C^2 + S^2)
  list(mu = mu, Rbar = Rbar, C = C, S = S)
}
cs <- circ_summary(theta)

# --- 1) Fit von Mises: moment-based start + Newton refine --------------------
A1  <- function(k) besselI(k,1)/besselI(k,0)
A1inv_approx <- function(R){
  if (R < 0.53) return(2*R + R^3 + (5*R^5)/6)
  if (R < 0.85) return(-0.4 + 1.39*R + 0.43/(1-R))
  1/(R^3 - 4*R^2 + 3*R)
}
mu_vm   <- cs$mu
Rbar    <- cs$Rbar
kappa_vm <- max(1e-8, A1inv_approx(Rbar))
# Newton refine kappa
for (i in 1:12){
  dA <- (A1(kappa_vm + 1e-6) - A1(kappa_vm - 1e-6)) / (2e-6)
  if (abs(dA) < 1e-10) break
  kappa_vm <- max(1e-8, kappa_vm - (A1(kappa_vm) - Rbar)/dA)
}

# vM log-likelihood
ll_vm_fun <- function(mu, kappa, t){
  n <- length(t)
  n * (-log(2*pi) - log(besselI(kappa,0))) + kappa * sum(cos(t - mu))
}
ll_vm  <- ll_vm_fun(mu_vm, kappa_vm, theta)
aic_vm <- -2*ll_vm + 2*2
bic_vm <- -2*ll_vm + log(n)*2

# --- 2) SSvM MLE (mu, logk, eta_raw with tanh link) --------------------------
# Density helpers
dvm   <- function(t, mu, kappa) (1/(2*pi*besselI(kappa,0))) * exp(kappa * cos(t - mu))
dssvm <- function(t, mu, kappa, eta) {
  tilt <- 1 + eta * sin(t - mu)
  base <- dvm(t, mu, kappa)
  base * pmax(tilt, .Machine$double.eps) # safe floor
}

# Negative log-likelihood on raw scale
nll_ssvm_raw <- function(par, t){
  mu_raw <- par[1]; logk <- par[2]; eta_raw <- par[3]
  mu    <- ((mu_raw + pi) %% (2*pi)) - pi
  kappa <- exp(logk)
  eta   <- tanh(eta_raw)
  tilt  <- 1 + eta * sin(t - mu)
  if (any(tilt <= 0)) return(1e12 + sum(pmax(0, -tilt))*1e6)
  ll <- length(t) * (-log(2*pi) - log(besselI(kappa,0))) +
    kappa * sum(cos(t - mu)) + sum(log(tilt))
  -ll
}

# Start values from vM + small skew
start_raw <- c(mu_vm, log(max(1e-8, kappa_vm)), atanh(0.05))
opt <- optim(par = start_raw, fn = nll_ssvm_raw, t = theta, method = "BFGS",
             control = list(maxit = 6000, reltol = 1e-10), hessian = TRUE)
if (opt$convergence != 0) warning("optim convergence code ", opt$convergence)

par_hat_raw  <- opt$par
mu_ssvm      <- ((par_hat_raw[1] + pi) %% (2*pi)) - pi
kappa_ssvm   <- exp(par_hat_raw[2])
eta_ssvm     <- tanh(par_hat_raw[3])
ll_ssvm      <- -nll_ssvm_raw(par_hat_raw, theta)
aic_ssvm     <- -2*ll_ssvm + 2*3
bic_ssvm     <- -2*ll_ssvm + log(n)*3

# Delta method SEs from Hessian
hess <- opt$hessian
cov_raw <- try(solve(hess), silent = TRUE)
if (inherits(cov_raw,"try-error")) {
  se_mu <- se_kappa <- se_eta <- NA_real_
} else {
  se_mu_raw   <- sqrt(cov_raw[1,1])
  se_logk_raw <- sqrt(cov_raw[2,2])
  se_eta_raw  <- sqrt(cov_raw[3,3])
  se_mu    <- se_mu_raw
  se_kappa <- abs(kappa_ssvm) * se_logk_raw
  se_eta   <- abs((1 - tanh(par_hat_raw[3])^2) * se_eta_raw)
}

# LR stat (asymptotic and parametric bootstrap)
lr_obs <- 2 * (ll_ssvm - ll_vm)

# --- 3) Parametric bootstrap (under vM null) for LR p-value + CIs -----------
set.seed(2025)
nboot <- 400  # bump to 1500–2000 for the manuscript run
rvm_fast <- function(n, mu, kappa) as.numeric(circular::rvonmises(n, mu = circular(mu), kappa = kappa))

boot_vm   <- matrix(NA_real_, nboot, 2)
boot_ssvm <- matrix(NA_real_, nboot, 3)
lr_stats  <- rep(NA_real_, nboot)

pb <- txtProgressBar(min = 0, max = nboot, style = 3)
for (b in 1:nboot) {
  sim <- rvm_fast(n, mu_vm, kappa_vm)
  # vM re-estimate
  csb <- circ_summary(sim); mu_b <- csb$mu; Rb <- csb$Rbar
  kap_b <- max(1e-8, A1inv_approx(Rb))
  for (i in 1:6) {
    dA <- (A1(kap_b + 1e-6) - A1(kap_b - 1e-6))/(2e-6)
    if (abs(dA) < 1e-10) break
    kap_b <- max(1e-8, kap_b - (A1(kap_b) - Rb)/dA)
  }
  boot_vm[b,] <- c(mu_b, kap_b)
  # SSvM fit on same sample
  st  <- c(mu_b, log(max(1e-8, kap_b)), atanh(0.01))
  fit <- try(optim(st, nll_ssvm_raw, t = sim, method = "BFGS",
                   control = list(maxit = 3000, reltol = 1e-10)),
             silent = TRUE)
  if (!inherits(fit,"try-error") && fit$convergence == 0) {
    mu_hat_b    <- ((fit$par[1]+pi) %% (2*pi)) - pi
    kap_hat_b   <- exp(fit$par[2])
    eta_hat_b   <- tanh(fit$par[3])
    boot_ssvm[b,] <- c(mu_hat_b, kap_hat_b, eta_hat_b)
    ll_vm_b    <- ll_vm_fun(mu_b, kap_b, sim)
    ll_ssvm_b  <- -nll_ssvm_raw(fit$par, sim)
    lr_stats[b] <- 2 * (ll_ssvm_b - ll_vm_b)
  }
  setTxtProgressBar(pb, b)
}
close(pb)

ok <- complete.cases(boot_ssvm) & is.finite(lr_stats)
boot_vm   <- boot_vm[ok, , drop = FALSE]
boot_ssvm <- boot_ssvm[ok, , drop = FALSE]
lr_stats  <- lr_stats[ok]

pval_lr_boot <- mean(lr_stats >= lr_obs)

# Percentile parametric CIs
q2 <- function(x) if (length(x)) quantile(x, c(.025,.975), na.rm = TRUE) else c(NA,NA)
ci_vm_mu     <- q2(boot_vm[,1]);  ci_vm_kappa  <- q2(boot_vm[,2])
ci_ssvm_mu   <- q2(boot_ssvm[,1]);ci_ssvm_kappa<- q2(boot_ssvm[,2]); ci_ssvm_eta <- q2(boot_ssvm[,3])

# --- 4) Simple goodness-of-fit tests ----------------------------------------
ray <- try(circular::rayleigh.test(circular(theta)), silent = TRUE)
ray_p <- if (inherits(ray,"try-error")) NA_real_ else as.numeric(ray$p.value)

# Fast model CDFs on a grid (avoid many integrate() calls)
grid_deg <- seq(0, 360, by = 0.5)
grid_rad <- (grid_deg %% 360) * pi/180
grid_rad <- ifelse(grid_rad > pi, grid_rad - 2*pi, grid_rad)
vm_grid  <- dvm(grid_rad,  mu_vm,  kappa_vm)
ss_grid  <- dssvm(grid_rad, mu_ssvm, kappa_ssvm, eta_ssvm)
F_vm_grid <- cumsum(vm_grid); F_vm_grid <- F_vm_grid / max(F_vm_grid)
F_ss_grid <- cumsum(ss_grid); F_ss_grid <- F_ss_grid / max(F_ss_grid)

# Empirical CDF on the same grid (P–P view)
emp_deg <- (theta * 180/pi) %% 360
F_emp_grid <- ecdf(emp_deg)(grid_deg)

# Kuiper-like metric vs each model
Kuiper_vm   <- (max(F_emp_grid - F_vm_grid) + max(F_vm_grid - F_emp_grid))
Kuiper_ssvm <- (max(F_emp_grid - F_ss_grid) + max(F_ss_grid - F_emp_grid))

# --- 5) Tables ---------------------------------------------------------------
param_table <- tibble::tibble(
  Model    = c("Von Mises","Sine-skewed VM"),
  mu_rad   = c(mu_vm, mu_ssvm),
  mu_deg   = c(mu_vm*180/pi, mu_ssvm*180/pi),
  kappa    = c(kappa_vm, kappa_ssvm),
  eta      = c(NA, eta_ssvm),
  logLik   = c(ll_vm, ll_ssvm),
  AIC      = c(aic_vm, aic_ssvm),
  BIC      = c(bic_vm, bic_ssvm)
)

se_table <- tibble::tibble(
  Param    = c("mu_rad","kappa","eta"),
  VM_est   = c(mu_vm, kappa_vm, NA),
  VM_SE    = c(NA, NA, NA),
  SSVM_est = c(mu_ssvm, kappa_ssvm, eta_ssvm),
  SSVM_SE  = c(se_mu, se_kappa, se_eta)
)

ci_table <- tibble::tibble(
  Param = c("mu_deg","kappa","eta"),
  VM_L  = c(ci_vm_mu[1]*180/pi,  ci_vm_kappa[1], NA),
  VM_U  = c(ci_vm_mu[2]*180/pi,  ci_vm_kappa[2], NA),
  SSVM_L= c(ci_ssvm_mu[1]*180/pi,ci_ssvm_kappa[1], ci_ssvm_eta[1]),
  SSVM_U= c(ci_ssvm_mu[2]*180/pi,ci_ssvm_kappa[2], ci_ssvm_eta[2])
)

comp_metrics <- tibble::tibble(
  Metric = c("logLik","AIC","BIC","LR_stat","LR_pval_param_boot","Kuiper_stat"),
  VM     = c(ll_vm, aic_vm, bic_vm, NA, NA, Kuiper_vm),
  SSVM   = c(ll_ssvm, aic_ssvm, bic_ssvm, lr_obs, pval_lr_boot, Kuiper_ssvm)
)

# Save CSVs
readr::write_csv(param_table, file.path(tab_dir, "parameters.csv"))
readr::write_csv(se_table,    file.path(tab_dir, "ses.csv"))
readr::write_csv(ci_table,    file.path(tab_dir, "cis_parametric.csv"))
readr::write_csv(comp_metrics,file.path(tab_dir, "comparison.csv"))

# Print to console (nice)
cat("\n--- Parameter estimates ---\n"); print(param_table)
cat("\n--- SEs (delta method for SSvM) ---\n"); print(se_table)
cat("\n--- Parametric bootstrap CIs ---\n"); print(ci_table)
cat("\n--- Comparative metrics ---\n"); print(comp_metrics)
cat(sprintf("\nRayleigh test p-value (non-uniformity): %.4g\n", ray_p))
cat(sprintf("LR obs = %.4f; parametric bootstrap p-value = %.4g\n", lr_obs, pval_lr_boot))
cat(sprintf("Kuiper VM = %.4g, SSVM = %.4g (smaller is better)\n\n", Kuiper_vm, Kuiper_ssvm))

# --- 6) Plots: colourful circular set ---------------------------------------
# Base histogram (“rose”) + overlay curves scaled to counts
plot_df   <- tibble::tibble(angle_deg = emp_deg)
hist_base <- ggplot(plot_df, aes(x = angle_deg)) +
  geom_histogram(aes(y = after_stat(count)), breaks = seq(0,360,by=10),
                 fill = alpha("#3b82f6", 0.65), color = "grey20", linewidth = 0.2) +
  coord_polar(start = -pi/2) +
  scale_x_continuous(breaks = seq(0,330, by = 30)) +
  labs(title = "Wind directions: rose diagram (10° bins)", x = NULL, y = "Counts")

# Overlay densities scaled to histogram
tmp_counts <- ggplot_build(
  ggplot(plot_df, aes(x = angle_deg)) + geom_histogram(breaks = seq(0,360,by=10))
)$data[[1]]$count
max_hist   <- max(tmp_counts)
dens_df    <- tibble::tibble(deg = grid_deg, vm = vm_grid, ssvm = ss_grid)
scale_fac  <- max_hist / max(dens_df$vm, dens_df$ssvm) * 0.9

polar_overlay <- ggplot() +
  geom_histogram(data = plot_df, aes(x = angle_deg, y = after_stat(count)),
                 breaks = seq(0,360,by=10), fill = alpha("#3b82f6",0.45),
                 color = "grey20", linewidth = 0.2) +
  geom_line(data = dens_df, aes(x = deg, y = vm * scale_fac, color = "vM"), linewidth = 1.2) +
  geom_line(data = dens_df, aes(x = deg, y = ssvm * scale_fac, color = "SSvM"),
            linewidth = 1.2, linetype = "longdash") +
  coord_polar(start = -pi/2) +
  scale_color_manual(values = c("vM" = "#fb923c", "SSvM" = "#22c55e")) +
  labs(title = "Histogram with fitted densities (scaled)", x = NULL, y = "Counts / scaled pdf", color = NULL)

# Kernel circular density vs fitted curves (linear view)
dens_circ <- density.circular(circular(theta), bw = 20*pi/180, kernel = "vonmises")
kd_df <- tibble::tibble(deg = (dens_circ$x*180/pi) %% 360, kde = dens_circ$y)

kernel_plot <- ggplot() +
  geom_line(data = kd_df, aes(x = deg, y = kde, color = "Kernel vm"), linewidth = 1.2) +
  geom_line(data = dens_df, aes(x = deg, y = vm,  color = "vM"), linewidth = 1.1) +
  geom_line(data = dens_df, aes(x = deg, y = ssvm, color = "SSvM"),
            linewidth = 1.1, linetype = "longdash") +
  scale_color_manual(values = c("Kernel vm" = "#3b82f6", "vM" = "#fb923c", "SSvM" = "#22c55e")) +
  labs(title = "Kernel circular density vs fitted vM / SSvM",
       x = "Angle (deg)", y = "Density", color = NULL)

# SSvM modifier: 1 + eta sin(theta - mu)
diag_df <- tibble::tibble(deg = grid_deg,
                          modifier = 1 + eta_ssvm * sin(grid_rad - mu_ssvm))
diag_plot <- ggplot(diag_df, aes(x = deg, y = modifier)) +
  geom_line(linewidth = 1.2, color = "#22c55e") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "SSvM modifier 1 + η sin(θ − μ)", x = "Angle (deg)", y = "Modifier")

# P–P view (empirical vs model CDF)
pp_df <- tibble::tibble(F_emp = F_emp_grid,
                        F_vM  = F_vm_grid,
                        F_SS  = F_ss_grid)
pp_plot <- ggplot(pp_df) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey40") +
  geom_path(aes(x = F_emp, y = F_vM, color = "vM"), linewidth = 1) +
  geom_path(aes(x = F_emp, y = F_SS, color = "SSvM"), linewidth = 1) +
  scale_color_manual(values = c("vM" = "#fb923c", "SSvM" = "#22c55e")) +
  labs(title = "Empirical vs model CDF (P–P)",
       x = "Empirical CDF", y = "Model CDF", color = NULL)

# Arrange and save all figures
top_row    <- plot_grid(hist_base, polar_overlay, ncol = 2)
bottom_row <- plot_grid(kernel_plot, diag_plot, pp_plot, ncol = 3)
all_plots  <- plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1, 1.05))

ggsave(file.path(fig_dir, "rose.png"),        hist_base,   width = 7.5, height = 6,  dpi = 300)
ggsave(file.path(fig_dir, "overlay.png"),     polar_overlay,width = 7.5, height = 6,  dpi = 300)
ggsave(file.path(fig_dir, "kernel.png"),      kernel_plot, width = 7.2, height = 4.6, dpi = 300)
ggsave(file.path(fig_dir, "modifier.png"),    diag_plot,   width = 7.2, height = 4.2, dpi = 300)
ggsave(file.path(fig_dir, "pp.png"),          pp_plot,     width = 6.2, height = 5.2, dpi = 300)
ggsave(file.path(fig_dir, "panel.png"),       all_plots,   width = 12,  height = 10, dpi = 300)

ggsave(file.path(fig_dir, "rose.pdf"),        hist_base,   width = 7.5, height = 6)
ggsave(file.path(fig_dir, "overlay.pdf"),     polar_overlay,width = 7.5, height = 6)
ggsave(file.path(fig_dir, "kernel.pdf"),      kernel_plot, width = 7.2, height = 4.6)
ggsave(file.path(fig_dir, "modifier.pdf"),    diag_plot,   width = 7.2, height = 4.2)
ggsave(file.path(fig_dir, "pp.pdf"),          pp_plot,     width = 6.2, height = 5.2)
ggsave(file.path(fig_dir, "panel.pdf"),       all_plots,   width = 12,  height = 10)

cat("Saved figures to:", fig_dir, "\nSaved tables to:", tab_dir, "\n\n")
################################################################################
# End of script
################################################################################
