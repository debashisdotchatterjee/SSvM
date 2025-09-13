## ============================================================
## SSvM vs VM (robust) on real data: mydata + fisherB2
## - Stable MLE for VM and SSvM
## - Robust parameter transforms and initialisation (no non-finite in optim)
## - Rose plots + fitted density overlays (VM in purple, SSvM in steelblue)
## ============================================================

suppressPackageStartupMessages({
  library(circular)
  library(openair)   # mydata
  library(CircStats) # rvonmises etc.
})

wrap_0_2pi <- function(x) x %% (2*pi)

# -------------------------
# SSvM density (stable form)
# -------------------------
dssvm <- function(theta, mu, kappa, eta) {
  theta <- as.numeric(theta)
  if (!is.finite(mu) || !is.finite(kappa) || !is.finite(eta)) return(rep(NA_real_, length(theta)))
  if (kappa <= 0 || abs(eta) >= 1) return(rep(NA_real_, length(theta)))
  # use scaled Bessel for numerical stability
  I0s <- besselI(kappa, 0, expon.scaled = TRUE)
  base <- exp(kappa * (cos(theta - mu) - 1)) / (2*pi*I0s)
  skew <- 1 + eta * sin(theta - mu)
  base * skew
}

# -------------------------
# Helper: A1 and inverse approx for von Mises
# -------------------------
A1 <- function(k) besselI(k,1)/besselI(k,0)
A1inv_approx <- function(R){
  if (R < 0.53) return(2*R + R^3 + (5*R^5)/6)
  if (R < 0.85) return(-0.4 + 1.39*R + 0.43/(1-R))
  return(1/(R^3 - 4*R^2 + 3*R))
}

# -------------------------
# MLE for VM (stable reparameterisation)
# -------------------------
mle_vm <- function(data) {
  data <- wrap_0_2pi(as.numeric(data))
  Cbar <- mean(cos(data)); Sbar <- mean(sin(data))
  mu0  <- atan2(Sbar, Cbar)
  R    <- sqrt(Cbar^2 + Sbar^2)
  kappa0 <- A1inv_approx(R)
  kappa0 <- max(kappa0, 1e-4)
  
  nll <- function(p) {
    mu <- p[1]
    kappa <- exp(p[2])
    f <- dvonmises(circular(data), mu = circular(mu), kappa = kappa)
    if (any(!is.finite(f) | f <= 0)) return(1e10)
    -sum(log(f))
  }
  
  par0 <- c(mu0, log(kappa0))
  fit <- try(optim(par0, nll, method = "L-BFGS-B",
                   lower = c(-Inf, log(1e-6)), upper = c(Inf, log(1e3)),
                   control = list(maxit = 2000)), silent = TRUE)
  if (inherits(fit, "try-error")) {
    # fallback: return method of moments estimate
    return(list(par = c(mu = mu0, kappa = kappa0), logLik = NA_real_))
  }
  mu_hat <- wrap_0_2pi(fit$par[1])
  kappa_hat <- exp(fit$par[2])
  list(par = c(mu = mu_hat, kappa = kappa_hat), logLik = -fit$value)
}

# -------------------------
# MLE for SSvM (uses multiple starts on eta) - robust
# -------------------------
mle_ssvm <- function(data, n_starts = 9L) {
  data <- wrap_0_2pi(as.numeric(data))
  Cbar <- mean(cos(data)); Sbar <- mean(sin(data))
  mu0  <- atan2(Sbar, Cbar)
  R    <- sqrt(Cbar^2 + Sbar^2)
  kappa0 <- A1inv_approx(R)
  kappa0 <- max(kappa0, 1e-3)
  
  eta_grid <- seq(-0.8, 0.8, length.out = n_starts)
  best <- list(logLik = -Inf, par = c(mu=NA_real_, kappa=NA_real_, eta=NA_real_))
  for (eta0 in eta_grid) {
    par0 <- c(mu = mu0, xi = log(kappa0), zeta = atanh(eta0))
    nll <- function(p) {
      mu <- p[1]; kappa <- exp(p[2]); eta <- tanh(p[3])
      f <- dssvm(data, mu, kappa, eta)
      if (any(!is.finite(f) | f <= 0)) return(1e10)
      -sum(log(f))
    }
    fit <- try(optim(par0, nll, method = "L-BFGS-B",
                     lower = c(-Inf, log(1e-6), -5), upper = c(Inf, log(1e3), 5),
                     control = list(maxit = 3000)), silent = TRUE)
    if (inherits(fit, "try-error")) next
    mu_hat <- wrap_0_2pi(fit$par[1]); kappa_hat <- exp(fit$par[2]); eta_hat <- tanh(fit$par[3])
    ll <- -fit$value
    if (is.finite(ll) && ll > best$logLik) best <- list(logLik = ll, par = c(mu=mu_hat,kappa=kappa_hat,eta=eta_hat))
  }
  best
}

# -------------------------
# Robust rose + overlay plotting (manual polar polygons)
# -------------------------
draw_rose_with_fits <- function(theta,
                                ssvm_par = NULL, vm_par = NULL,
                                binwidth_deg = 10, prop = 1.5, main = "Rose + fits") {
  theta <- wrap_0_2pi(as.numeric(theta))
  bins  <- max(1L, round(360 / binwidth_deg))
  brk <- seq(0, 2*pi, length.out = bins+1)
  idx <- cut(theta, breaks = brk, include.lowest = TRUE, right = FALSE, labels = FALSE)
  tab <- tabulate(idx, nbins = bins)
  rc <- list(counts = tab, breaks = brk, mid = brk[-length(brk)] + (pi/bins), bw = 2*pi/bins)
  n <- length(theta)
  max_ct <- max(rc$counts, na.rm = TRUE)
  if (max_ct <= 0) stop("All counts are zero; check data.")
  
  oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar), add = TRUE)
  par(mar = c(2,2,3,2))
  plot(0,0, type='n', xlim=c(-prop*1.15, prop*1.15), ylim=c(-prop*1.15, prop*1.15),
       xlab="", ylab="", axes = FALSE, asp = 1, main = main)
  rs <- pretty(c(0, prop), n=3)
  for (r in rs) symbols(0,0, circles = r, inches=FALSE, add = TRUE, fg = "grey85")
  txt_pos <- c(pi/2, 0, 3*pi/2, pi)
  text( (prop+0.05)*cos(txt_pos), (prop+0.05)*sin(txt_pos), labels = c("N","E","S","W"))
  
  for (i in seq_len(bins)) {
    mid <- rc$mid[i]
    ct  <- rc$counts[i]
    r_outer <- prop * (ct / max_ct)
    theta1 <- mid - rc$bw/2
    theta2 <- mid + rc$bw/2
    phis <- seq((pi/2)-theta1, (pi/2)-theta2, length.out = 12)
    xs <- c(0, r_outer*cos(phis), 0)
    ys <- c(0, r_outer*sin(phis), 0)
    polygon(xs, ys, col = "grey90", border = "grey60")
  }
  
  add_curve <- function(pdf_fun, par_vec, col, grid = seq(0, 2*pi - 1e-6, length.out = 720L)) {
    dens <- pdf_fun(grid, par_vec)
    y <- n * rc$bw * dens
    r <- prop * (y / max_ct)
    r[!is.finite(r) | r < 0] <- NA
    phi <- (pi/2) - grid
    x <- r * cos(phi); yv <- r * sin(phi)
    lines(x, yv, col = col, lwd = 2)
  }
  
  pdf_ssvm <- function(th, par) dssvm(th, par["mu"], par["kappa"], par["eta"]) 
  pdf_vm   <- function(th, par) as.numeric(dvonmises(circular(th), mu = circular(par["mu"]), kappa = par["kappa"]))
  
  if (!is.null(ssvm_par) && all(is.finite(ssvm_par[c("mu","kappa","eta")])) ) add_curve(pdf_ssvm, ssvm_par, col = "steelblue")
  if (!is.null(vm_par) && all(is.finite(vm_par[c("mu","kappa")])) ) add_curve(pdf_vm, vm_par, col = "purple")
  
  arrow_one <- function(mu, col) {
    r_mu <- prop * 1.05; phi_mu <- (pi/2) - mu
    arrows(0, 0, r_mu*cos(phi_mu), r_mu*sin(phi_mu), length = 0.08, col = col, lwd = 1.5)
  }
  if (!is.null(ssvm_par) && is.finite(ssvm_par["mu"])) arrow_one(ssvm_par["mu"], "steelblue")
  if (!is.null(vm_par) && is.finite(vm_par["mu"])) arrow_one(vm_par["mu"], "purple")
  
  legend("topright", inset = 0.02, bty = "n",
         legend = c("SSvM  n·Δθ·f(θ)", "VM  n·Δθ·f(θ)"),
         col = c("steelblue","purple"), lwd = 2)
}

# ============================================================
# Run: mydata (openair) and fisherB2 (circular)
# ============================================================

# 1) mydata (wind directions in degrees)
data(mydata, package = "openair")
theta_my <- mydata$wd
theta_my <- theta_my[is.finite(theta_my)]
theta_my <- wrap_0_2pi(theta_my * pi/180)

fit_vm_my <- mle_vm(theta_my)
fit_ssvm_my <- mle_ssvm(theta_my)

ssvm_my <- fit_ssvm_my$par; names(ssvm_my) <- c("mu","kappa","eta")
vm_my   <- fit_vm_my$par;   names(vm_my)   <- c("mu","kappa")

ll_ssvm_my <- fit_ssvm_my$logLik
ll_vm_my   <- fit_vm_my$logLik

crit_ssvm_my <- c(AIC = -2*ll_ssvm_my + 2*3, BIC = -2*ll_ssvm_my + log(length(theta_my))*3)
crit_vm_my   <- c(AIC = -2*ll_vm_my + 2*2, BIC = -2*ll_vm_my + log(length(theta_my))*2)

cat("\nopenair::mydata — Fit Comparison\n")
print(data.frame(
  Model = c("SSvM","VM"),
  mu = c(ssvm_my["mu"], vm_my["mu"]),
  kappa = c(ssvm_my["kappa"], vm_my["kappa"]),
  eta = c(ssvm_my["eta"], NA),
  logLik = c(ll_ssvm_my, ll_vm_my),
  AIC = c(crit_ssvm_my["AIC"], crit_vm_my["AIC"]),
  BIC = c(crit_ssvm_my["BIC"], crit_vm_my["BIC"]) ), row.names = FALSE)

# Plot
if (dev.interactive()) dev.new(noRStudioGD = TRUE)
draw_rose_with_fits(theta_my, ssvm_par = ssvm_my, vm_par = vm_my,
                    binwidth_deg = 10, prop = 1.5,
                    main = "Wind directions (mydata): rose + SSvM & VM")

# 2) fisherB2 (already in radians)
data(fisherB2, package = "circular")
theta_f <- as.numeric(fisherB2)
theta_f <- theta_f[is.finite(theta_f)]
theta_f <- wrap_0_2pi(theta_f)

fit_vm_f <- mle_vm(theta_f)
fit_ssvm_f <- mle_ssvm(theta_f)

ssvm_f <- fit_ssvm_f$par; names(ssvm_f) <- c("mu","kappa","eta")
vm_f   <- fit_vm_f$par;   names(vm_f)   <- c("mu","kappa")

ll_ssvm_f <- fit_ssvm_f$logLik
ll_vm_f   <- fit_vm_f$logLik

crit_ssvm_f <- c(AIC = -2*ll_ssvm_f + 2*3, BIC = -2*ll_ssvm_f + log(length(theta_f))*3)
crit_vm_f   <- c(AIC = -2*ll_vm_f + 2*2, BIC = -2*ll_vm_f + log(length(theta_f))*2)

cat("\ncircular::fisherB2 — Fit Comparison\n")
print(data.frame(
  Model = c("SSvM","VM"),
  mu = c(ssvm_f["mu"], vm_f["mu"]),
  kappa = c(ssvm_f["kappa"], vm_f["kappa"]),
  eta = c(ssvm_f["eta"], NA),
  logLik = c(ll_ssvm_f, ll_vm_f),
  AIC = c(crit_ssvm_f["AIC"], crit_vm_f["AIC"]),
  BIC = c(crit_ssvm_f["BIC"], crit_vm_f["BIC"]) ), row.names = FALSE)

if (dev.interactive()) dev.new(noRStudioGD = TRUE)
draw_rose_with_fits(theta_f, ssvm_par = ssvm_f, vm_par = vm_f,
                    binwidth_deg = 10, prop = 1.5,
                    main = "Fisher B2: rose + SSvM & VM")

# End -----------------------------------------------------------------------
