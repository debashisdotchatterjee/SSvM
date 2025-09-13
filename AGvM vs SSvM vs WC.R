## ============================================================
## SSvM vs AGvM vs Wrapped Cauchy on real data
## Rose plots + fitted circular densities (base + circular only)
## Added: Wrapped Cauchy density, MLE and overlay
## ============================================================

suppressPackageStartupMessages({
  library(circular)
  library(openair)   # mydata
})

wrap_0_2pi <- function(x) x %% (2*pi)

## ---------- SSvM (stable closed-form) ----------
dssvm <- function(theta, mu, kappa, eta) {
  theta <- as.numeric(theta)
  if (!is.finite(mu) || !is.finite(kappa) || !is.finite(eta)) return(rep(NA_real_, length(theta)))
  if (kappa <= 0 || abs(eta) >= 1) return(rep(NA_real_, length(theta)))
  I0s <- besselI(kappa, 0, expon.scaled = TRUE)               # I0(k)*exp(-k)
  base <- exp(kappa*(cos(theta - mu) - 1)) / (2*pi*I0s)       # stable
  skew <- pmax(1 + eta*sin(theta - mu), 0)                    # clamp tiny negatives
  base * skew
}
loglik_ssvm <- function(par, data) {
  mu <- par[1]; kappa <- par[2]; eta <- par[3]
  if (!is.finite(mu) || !is.finite(kappa) || !is.finite(eta)) return(-Inf)
  if (kappa <= 0 || abs(eta) >= 1) return(-Inf)
  f <- dssvm(data, mu, kappa, eta)
  if (any(!is.finite(f) | f <= 0)) return(-Inf)
  sum(log(f))
}
mle_ssvm <- function(data, n_starts = 9L) {
  data <- wrap_0_2pi(data)
  Cbar <- mean(cos(data)); Sbar <- mean(sin(data))
  mu0  <- atan2(Sbar, Cbar)
  R    <- sqrt(Cbar^2 + Sbar^2)
  kappa0 <- if (R < 0.53) 2*R + R^3 + (5*R^5)/6 else if (R < 0.85) -0.4 + 1.39*R + 0.43/(1-R) else 1/(R^3 - 4*R^2 + 3*R)
  kappa0 <- max(kappa0, 1e-3)
  eta_grid <- seq(-0.8, 0.8, length.out = n_starts)
  best <- list(logLik = -Inf, par = c(mu=NA_real_, kappa=NA_real_, eta=NA_real_))
  for (eta0 in eta_grid) {
    par0 <- c(mu=mu0, xi=log(kappa0), zeta=atanh(eta0))
    nll <- function(p) {
      mu <- p[1]; kappa <- exp(p[2]); eta <- tanh(p[3])
      val <- -loglik_ssvm(c(mu,kappa,eta), data)
      if (!is.finite(val)) 1e10 else val
    }
    fit <- try(optim(par0, nll, method="L-BFGS-B",
                     lower=c(-Inf, log(1e-4), -5),
                     upper=c( Inf, log(50),    5),
                     control=list(maxit=2000)), silent=TRUE)
    if (inherits(fit,"try-error")) next
    mu_hat    <- wrap_0_2pi(fit$par[1])
    kappa_hat <- exp(fit$par[2])
    eta_hat   <- tanh(fit$par[3])
    ll <- -fit$value
    if (is.finite(ll) && ll > best$logLik) best <- list(logLik = ll, par = c(mu=mu_hat,kappa=kappa_hat,eta=eta_hat))
  }
  best
}

## ---------- AGvM (robust numeric normalizer) ----------
.Zcache <- new.env(parent = emptyenv())
Z_agvm <- function(k1, k2) {
  if (!is.finite(k1) || !is.finite(k2)) return(NA_real_)
  key <- paste0(round(k1,7), "_", round(k2,7))
  if (exists(key, envir=.Zcache, inherits=FALSE)) return(get(key, envir=.Zcache))
  fun <- function(t) exp(k1*cos(t) + k2*sin(2*t))
  val <- try(integrate(fun, 0, 2*pi, subdivisions=3000L, rel.tol=1e-9)$value, silent=TRUE)
  if (inherits(val,"try-error") || !is.finite(val) || val <= 0) val <- NA_real_
  assign(key, val, envir=.Zcache); val
}
dagvm <- function(theta, mu, k1, k2) {
  theta <- as.numeric(theta)
  if (!is.finite(mu) || !is.finite(k1) || !is.finite(k2) || k1 <= 0) return(rep(NA_real_, length(theta)))
  Z <- Z_agvm(k1, k2)
  if (!is.finite(Z) || Z <= 0) return(rep(NA_real_, length(theta)))
  exp(k1*cos(theta - mu) + k2*sin(2*(theta - mu))) / Z
}
loglik_agvm <- function(par, data) {
  mu <- par[1]; k1 <- par[2]; k2 <- par[3]
  if (!is.finite(mu) || !is.finite(k1) || !is.finite(k2) || k1 <= 0) return(-Inf)
  Z <- Z_agvm(k1, k2); if (!is.finite(Z) || Z <= 0) return(-Inf)
  n  <- length(data)
  s1 <- sum(cos(data - mu)); s2 <- sum(sin(2*(data - mu)))
  -n*log(Z) + k1*s1 + k2*s2
}
mle_agvm <- function(data, n_starts = 9L) {
  data <- wrap_0_2pi(data)
  Cbar <- mean(cos(data)); Sbar <- mean(sin(data))
  mu0  <- atan2(Sbar, Cbar)
  k1_grid <- exp(seq(log(0.2), log(5), length.out=n_starts))
  k2_grid <- seq(-2, 2, length.out=n_starts)
  best <- list(logLik = -Inf, par = c(mu=NA_real_, k1=NA_real_, k2=NA_real_))
  for (k1s in k1_grid) for (k2s in k2_grid) {
    par0 <- c(mu=mu0, log_k1=log(k1s), k2=k2s)
    nll <- function(p) {
      mu <- p[1]; k1 <- exp(p[2]); k2 <- p[3]
      val <- -loglik_agvm(c(mu,k1,k2), data)
      if (!is.finite(val)) 1e10 else val
    }
    fit <- try(optim(par0, nll, method="L-BFGS-B",
                     lower=c(-Inf, log(1e-4), -10),
                     upper=c( Inf, log(50),   10),
                     control=list(maxit=2000)), silent=TRUE)
    if (inherits(fit,"try-error")) next
    mu_hat <- wrap_0_2pi(fit$par[1]); k1_hat <- exp(fit$par[2]); k2_hat <- fit$par[3]
    ll <- -fit$value
    if (is.finite(ll) && ll > best$logLik) best <- list(logLik = ll, par = c(mu=mu_hat,k1=k1_hat,k2=k2_hat))
  }
  best
}

## ---------- Wrapped Cauchy (truncated infinite sum) ----------
## density: f_{WC}(theta; mu, gamma) = sum_{n=-Inf}^{Inf} gamma / [ pi (gamma^2 + (theta-mu+2*pi*n)^2) ]
## we truncate the sum at n = -N..N; N=50 is typically ample for reasonable gamma
dwc_trunc <- function(theta, mu, gamma, N = 50) {
  theta <- as.numeric(theta)
  if (!is.finite(mu) || !is.finite(gamma) || gamma <= 0) return(rep(NA_real_, length(theta)))
  # center difference
  d <- outer(theta, (-N:N), function(t, n) t - mu + 2*pi*n)
  num <- gamma
  denom <- pi * (gamma^2 + d^2)
  dens <- rowSums(num / denom)
  dens
}
loglik_wc <- function(par, data, N = 50) {
  mu <- par[1]; gamma <- par[2]
  if (!is.finite(mu) || !is.finite(gamma) || gamma <= 0) return(-Inf)
  f <- dwc_trunc(data, mu, gamma, N = N)
  if (any(!is.finite(f) | f <= 0)) return(-Inf)
  sum(log(f))
}
mle_wc <- function(data, n_starts = 9L, N = 50) {
  data <- wrap_0_2pi(data)
  Cbar <- mean(cos(data)); Sbar <- mean(sin(data))
  mu0  <- atan2(Sbar, Cbar)
  # gamma initial: use a small positive value; try a grid on log-scale
  gamma_grid <- exp(seq(log(0.01), log(2), length.out = n_starts))
  best <- list(logLik = -Inf, par = c(mu=NA_real_, gamma=NA_real_))
  for (g0 in gamma_grid) {
    par0 <- c(mu = mu0, log_g = log(g0))
    nll <- function(p) {
      mu <- p[1]; gamma <- exp(p[2])
      val <- -loglik_wc(c(mu, gamma), data, N = N)
      if (!is.finite(val)) 1e10 else val
    }
    fit <- try(optim(par0, nll, method = "L-BFGS-B",
                     lower = c(-Inf, log(1e-6)),
                     upper = c( Inf, log(10)),
                     control = list(maxit = 2000)), silent = TRUE)
    if (inherits(fit, "try-error")) next
    mu_hat <- wrap_0_2pi(fit$par[1]); gamma_hat <- exp(fit$par[2])
    ll <- -fit$value
    if (is.finite(ll) && ll > best$logLik) best <- list(logLik = ll, par = c(mu = mu_hat, gamma = gamma_hat))
  }
  best
}

## ---------- AIC/BIC ----------
criteria <- function(loglik, p, n) c(AIC = -2*loglik + 2*p, BIC = -2*loglik + log(n)*p)

## ---------- Helper: rose counts (for scaling overlays) ----------
rose_counts <- function(theta, bins) {
  theta <- wrap_0_2pi(theta)
  brk <- seq(0, 2*pi, length.out = bins+1)
  idx <- cut(theta, breaks = brk, include.lowest = TRUE, right = FALSE, labels = FALSE)
  tab <- tabulate(idx, nbins = bins)
  list(counts = tab, breaks = brk, mid = brk[-length(brk)] + (pi/bins), bw = 2*pi/bins)
}

## ---------- Draw rose + density overlays (base + circular only) ----------
## Orientation: zero at North, clockwise increase (like wind roses)
draw_rose_with_fits <- function(theta,
                                ssvm_par = NULL, agvm_par = NULL, wc_par = NULL,
                                binwidth_deg = 10,
                                prop = 1.5,
                                main = "Rose with fitted circular densities") {
  
  theta <- wrap_0_2pi(theta)
  bins  <- max(1L, round(360 / binwidth_deg))
  rc    <- rose_counts(theta, bins = bins)
  n     <- length(theta)
  max_ct <- max(rc$counts, na.rm = TRUE)
  if (max_ct <= 0) stop("All counts are zero; check data.")
  
  ## 1) Base rose histogram (drawn by circular)
  rose.diag(circular(theta, type="angles", units="radians", template="none"),
            bins = bins, col = "grey90", border = "grey50",
            prop = prop,  # max bar length on the plot
            shrink = 1,   # full radius
            axes = TRUE,  # draw circular axes
            zero = pi/2,  # 0° at North
            rotation = "clock", main = main)
  
  ## 2) Build dense grid for overlays (avoid exactly 0 and 2π to keep continuity)
  grid_th <- seq(0, 2*pi - 1e-6, length.out = 720L)
  
  ## Expected counts per angle:  y = n * bw * pdf(theta)
  add_curve <- function(y, col, grid = grid_th) {
    y[!is.finite(y) | y < 0] <- NA
    if (all(is.na(y))) return(invisible(NULL))
    ## Convert counts to the same radial scale used by rose.diag:
    ## bar length r = prop * (count / max_ct). So scale expected counts similarly.
    r <- prop * (y / max_ct)
    
    ## Convert (theta,r) to Cartesian with same orientation as rose.diag:
    ## zero = pi/2, rotation = "clock"  => plotting angle phi = zero - theta
    phi <- (pi/2) - grid
    x <- r * cos(phi)
    y <- r * sin(phi)
    
    lines(x, y, col = col, lwd = 2)
  }
  
  ## SSvM overlay
  if (!is.null(ssvm_par) && all(is.finite(ssvm_par[c("mu","kappa","eta")]))) {
    y_ssvm <- n * rc$bw * dssvm(grid_th, ssvm_par["mu"], ssvm_par["kappa"], ssvm_par["eta"]) 
    add_curve(y_ssvm, col = "steelblue")
    ## μ arrow
    r_mu <- prop * 1.05
    phi_mu <- (pi/2) - ssvm_par["mu"]
    arrows(0, 0, r_mu*cos(phi_mu), r_mu*sin(phi_mu), length = 0.08, col = "steelblue", lwd = 1.5)
  }
  
  ## AGvM overlay
  if (!is.null(agvm_par) && all(is.finite(agvm_par[c("mu","k1","k2")]))) {
    y_agvm <- n * rc$bw * dagvm(grid_th, agvm_par["mu"], agvm_par["k1"], agvm_par["k2"]) 
    add_curve(y_agvm, col = "firebrick")
    ## μ arrow
    r_mu <- prop * 1.05
    phi_mu <- (pi/2) - agvm_par["mu"]
    arrows(0, 0, r_mu*cos(phi_mu), r_mu*sin(phi_mu), length = 0.08, col = "firebrick", lwd = 1.5)
  }
  
  ## Wrapped Cauchy overlay
  if (!is.null(wc_par) && all(is.finite(wc_par[c("mu","gamma")]))) {
    y_wc <- n * rc$bw * dwc_trunc(grid_th, wc_par["mu"], wc_par["gamma"], N = 80)
    add_curve(y_wc, col = "darkgreen")
    ## μ arrow
    r_mu <- prop * 1.05
    phi_mu <- (pi/2) - wc_par["mu"]
    arrows(0, 0, r_mu*cos(phi_mu), r_mu*sin(phi_mu), length = 0.08, col = "darkgreen", lwd = 1.5)
  }
  
  legend("topright", inset = 0.02, bty = "n",
         col = c("steelblue","firebrick","darkgreen"), lwd = 2,
         legend = c("SSvM  n·Δθ·f(θ)", "AGvM  n·Δθ·f(θ)", "Wrapped Cauchy  n·Δθ·f(θ)"))
}

## ============================================================
## REAL DATA 1: openair::mydata (wind direction in degrees)
## ============================================================
data(mydata, package = "openair")
theta_my <- mydata$wd
theta_my <- theta_my[is.finite(theta_my)]
theta_my <- wrap_0_2pi(theta_my * pi/180)

fit_ssvm_my <- mle_ssvm(theta_my)
fit_agvm_my <- mle_agvm(theta_my)
fit_wc_my   <- mle_wc(theta_my)

ssvm_my <- fit_ssvm_my$par; names(ssvm_my) <- c("mu","kappa","eta")
agvm_my <- fit_agvm_my$par; names(agvm_my) <- c("mu","k1","k2")
wc_my   <- fit_wc_my$par;   names(wc_my)   <- c("mu","gamma")

ll_ssvm_my  <- fit_ssvm_my$logLik
ll_agvm_my  <- fit_agvm_my$logLik
ll_wc_my    <- fit_wc_my$logLik
crit_ssvm_my<- criteria(ll_ssvm_my, p = 3, n = length(theta_my))
crit_agvm_my<- criteria(ll_agvm_my, p = 3, n = length(theta_my))
crit_wc_my  <- criteria(ll_wc_my,   p = 2, n = length(theta_my))

cat("\nopenair::mydata — Fit Comparison\n")
print(data.frame(
  Model    = c("SSvM","AGvM","WrappedCauchy"),
  mu       = c(ssvm_my["mu"],  agvm_my["mu"], wc_my["mu"]),
  kappa_k1 = c(ssvm_my["kappa"], agvm_my["k1"], wc_my["gamma"]),
  eta_k2   = c(ssvm_my["eta"],   agvm_my["k2"], NA),
  logLik   = c(ll_ssvm_my, ll_agvm_my, ll_wc_my),
  AIC      = c(crit_ssvm_my["AIC"], crit_agvm_my["AIC"], crit_wc_my["AIC"]),
  BIC      = c(crit_ssvm_my["BIC"], crit_agvm_my["BIC"], crit_wc_my["BIC"]) 
), row.names = FALSE)

## Plot (base + circular only)
if (dev.interactive()) dev.new(noRStudioGD = TRUE)  # open new device if interactive
draw_rose_with_fits(theta_my, ssvm_par = ssvm_my, agvm_par = agvm_my, wc_par = wc_my,
                    binwidth_deg = 10,
                    prop = 1.5,
                    main = "Wind directions (mydata): rose + fitted densities")

## ============================================================
## REAL DATA 2: circular::fisherB2 (radians)
## ============================================================
data(fisherB2, package = "circular")
theta_f <- as.numeric(fisherB2)
theta_f <- theta_f[is.finite(theta_f)]
theta_f <- wrap_0_2pi(theta_f)

fit_ssvm_f <- mle_ssvm(theta_f)
fit_agvm_f <- mle_agvm(theta_f)
fit_wc_f   <- mle_wc(theta_f)

ssvm_f <- fit_ssvm_f$par; names(ssvm_f) <- c("mu","kappa","eta")
agvm_f <- fit_agvm_f$par; names(agvm_f) <- c("mu","k1","k2")
wc_f   <- fit_wc_f$par;   names(wc_f)   <- c("mu","gamma")

ll_ssvm_f  <- fit_ssvm_f$logLik
ll_agvm_f  <- fit_agvm_f$logLik
ll_wc_f    <- fit_wc_f$logLik
crit_ssvm_f<- criteria(ll_ssvm_f, p = 3, n = length(theta_f))
crit_agvm_f<- criteria(ll_agvm_f, p = 3, n = length(theta_f))
crit_wc_f  <- criteria(ll_wc_f,   p = 2, n = length(theta_f))

cat("\ncircular::fisherB2 — Fit Comparison\n")
print(data.frame(
  Model    = c("SSvM","AGvM","WrappedCauchy"),
  mu       = c(ssvm_f["mu"],  agvm_f["mu"], wc_f["mu"]),
  kappa_k1 = c(ssvm_f["kappa"], agvm_f["k1"], wc_f["gamma"]),
  eta_k2   = c(ssvm_f["eta"],   agvm_f["k2"], NA),
  logLik   = c(ll_ssvm_f, ll_agvm_f, ll_wc_f),
  AIC      = c(crit_ssvm_f["AIC"], crit_agvm_f["AIC"], crit_wc_f["AIC"]),
  BIC      = c(crit_ssvm_f["BIC"], crit_agvm_f["BIC"], crit_wc_f["BIC"]) 
), row.names = FALSE)

if (dev.interactive()) dev.new(noRStudioGD = TRUE)
draw_rose_with_fits(theta_f, ssvm_par = ssvm_f, agvm_par = agvm_f, wc_par = wc_f,
                    binwidth_deg = 10,
                    prop = 1.5,
                    main = "Fisher B2: rose + fitted densities")

## End of script
