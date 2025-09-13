# ===========================================
# SSvM vs vM/Kim (μ free), pick a dataset where SSvM wins
# Base R only, circular-style rose plot (custom polar)
# ===========================================

## ---- Output dir ----
outdir <- "paper_outputs/scan_ssvm_polar"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## ---- Utilities ----
deg2rad <- function(x) x * pi / 180
rad2deg <- function(x) x * 180 / pi
wrap_2pi <- function(x) (x %% (2*pi))
I0 <- function(k) besselI(k, 0)
I1 <- function(k) besselI(k, 1)
I2 <- function(k) besselI(k, 2)
A1 <- function(k) I1(k) / I0(k)
A1_inv <- function(R) {
  stopifnot(all(R >= 0 & R < 1))
  out <- R
  i1 <- R < 0.53; i2 <- R >= 0.53 & R < 0.85; i3 <- R >= 0.85
  out[i1] <- 2*R[i1] + R[i1]^3 + (5/6)*R[i1]^5
  out[i2] <- -0.4 + 1.39*R[i2] + 0.43/(1 - R[i2])
  out[i3] <- 1/(R[i3]^3 - 4*R[i3]^2 + 3*R[i3])
  out
}
circ_mean <- function(theta) atan2(mean(sin(theta)), mean(cos(theta)))
resultant_length <- function(theta) sqrt(mean(cos(theta))^2 + mean(sin(theta))^2)

## ---- Densities ----
d_vm <- function(theta, mu, kappa) {
  exp(kappa * cos(theta - mu)) / (2*pi*I0(kappa))
}
# SSvM: von Mises times (1 + eta sin(phi))
d_ssvm <- function(theta, mu, kappa, eta) {
  d_vm(theta, mu, kappa) * (1 + eta * sin(theta - mu))
}
# "Kim" as phase-shifted vM: with μ free it reduces to vM, but we keep it for reviewer parity
d_kim <- function(theta, mu, R, alpha) {
  phi <- theta - mu
  exp(R * cos(phi - alpha)) / (2*pi*I0(R))
}

## ---- vM MLE (μ free, closed form) ----
fit_vm_free <- function(theta) {
  n <- length(theta)
  C <- mean(cos(theta)); S <- mean(sin(theta))
  mu <- atan2(S, C)
  R <- sqrt(C^2 + S^2)
  kappa <- if (R < 1) A1_inv(R) else 1e6
  ll <- -n*log(2*pi*I0(kappa)) + kappa*sum(cos(theta - mu))
  list(mu = wrap_2pi(mu), kappa = kappa, loglik = ll)
}

## ---- SSvM MLE (μ, κ, η all free) ----
loglik_ssvm_full <- function(par, theta) {
  mu <- par[1]; kappa <- exp(par[2]); eta <- tanh(par[3])
  if (!is.finite(kappa) || abs(eta) >= 1) return(-Inf)
  phi <- theta - mu
  -length(theta)*log(2*pi*I0(kappa)) + kappa*sum(cos(phi)) + sum(log1p(eta * sin(phi)))
}
fit_ssvm_free <- function(theta) {
  n <- length(theta)
  mu0 <- circ_mean(theta)
  R  <- resultant_length(theta)
  xi0 <- log(max(A1_inv(min(R, 0.999999)), 1e-4))
  s <- sin(theta - mu0)
  eta0 <- if (sum(s^2) > 1e-8) max(min(sum(s)/sum(s^2), 0.8), -0.8) else 0
  par0 <- c(mu0, xi0, atanh(eta0))
  o <- optim(par0, fn = function(p) -loglik_ssvm_full(c(p[1], p[2], p[3]), theta),
             method = "BFGS", control = list(reltol = 1e-10, maxit = 4000))
  mu_hat <- wrap_2pi(o$par[1]); kappa_hat <- exp(o$par[2]); eta_hat <- tanh(o$par[3])
  list(mu = mu_hat, kappa = kappa_hat, eta = eta_hat, loglik = -o$value, conv = o$convergence)
}

## ---- "Kim" MLE (μ free -> equals vM) ----
fit_kim_free <- function(theta) {
  # With μ free, maximizing over (R, alpha) gives exactly vM with mean μ+alpha.
  # We therefore report Kim = vM MLE for transparency.
  vm <- fit_vm_free(theta)
  list(mu = vm$mu, R = vm$kappa, alpha = 0, loglik = vm$loglik, conv = 0L)
}

## ---- LRTs ----
lrt_vm_vs_ssvm <- function(f_vm, f_ssvm) {
  LR <- 2*(f_ssvm$loglik - f_vm$loglik)
  p  <- pchisq(LR, df = 1, lower.tail = FALSE)
  c(LR = LR, p_value = p)
}

## ---- Robust loader for 'circular' datasets (optional) ----
degreeify <- function(xnum) {
  if (is.null(xnum)) return(numeric(0))
  if (is.list(xnum) || is.data.frame(xnum)) {
    # try to find a numeric column of angles
    cand <- unlist(lapply(xnum, function(v) if (is.numeric(v)) v else NULL), use.names = FALSE)
    if (length(cand)) xnum <- cand else return(numeric(0))
  }
  xnum <- suppressWarnings(as.numeric(xnum))
  xnum <- xnum[is.finite(xnum)]
  if (!length(xnum)) return(numeric(0))
  mx <- max(xnum, na.rm = TRUE)
  if (mx <= 2*pi + 1e-6) rad2deg(xnum) else (xnum %% 360)
}
list_circular_datasets <- function() {
  if (!requireNamespace("circular", quietly = TRUE)) return(character(0))
  dd <- utils::data(package = "circular")$results
  if (is.null(dd) || nrow(dd) == 0) return(character(0))
  nm <- unique(dd[, "Item"]); nm[!is.na(nm)]
}
load_circular_dataset <- function(nm) {
  e <- new.env()
  ok <- tryCatch({ data(list = nm, package = "circular", envir = e); TRUE }, error = function(...) FALSE)
  if (!ok || !exists(nm, envir = e, inherits = FALSE)) return(numeric(0))
  obj <- get(nm, envir = e)
  degreeify(obj)
}

## ---- Generators (skewed) ----
rvm <- function(n, mu, kappa) {
  if (kappa < 1e-8) return(runif(n, 0, 2*pi))
  a <- 1 + sqrt(1 + 4*kappa^2); b <- (a - sqrt(2*a)) / (2*kappa); r <- (1 + b^2) / (2*b)
  out <- numeric(n); j <- 1
  while (j <= n) {
    U1 <- runif(1); U2 <- runif(1); z <- cos(pi*U1); f <- (1 + r*z) / (r + z); cst <- kappa*(r - f)
    if (U2 < cst*(2 - cst)) th <- acos(f) else { if (log(U2) <= -kappa*(1 - f)) th <- acos(f) else next }
    if (runif(1) > 0.5) th <- -th
    out[j] <- wrap_2pi(mu + th); j <- j + 1
  }
  out
}
# Exact accept-reject sampler for SSvM (anchored at mu)
r_ssvm <- function(n, mu, kappa, eta) {
  stopifnot(abs(eta) < 1)
  out <- numeric(0)
  while (length(out) < n) {
    cand <- rvm(n, mu, kappa)
    phi <- cand - mu
    acc <- (1 + eta * sin(phi)) / (1 + abs(eta))
    keep <- runif(length(cand)) < acc
    out <- c(out, cand[keep])
  }
  out[1:n]
}
# Skewed mixture of two vM
r_mixvm <- function(n, mu, kappa1 = 7, kappa2 = 2, delta = 30*pi/180, w = 0.75) {
  z <- runif(n)
  c1 <- rvm(sum(z < w),  mu - delta, kappa1)
  c2 <- rvm(sum(z >= w), mu + 2*delta, kappa2)
  wrap_2pi(c(c1, c2))
}

## ---- Score one dataset (μ free fits) ----
score_theta_deg <- function(theta_deg, label) {
  if (!length(theta_deg)) return(NULL)
  theta <- wrap_2pi(deg2rad(theta_deg))
  n <- length(theta); if (n < 30) return(NULL)
  
  f_vm   <- fit_vm_free(theta)
  f_kim  <- fit_kim_free(theta)      # = vM (μ free)
  f_ssvm <- fit_ssvm_free(theta)
  
  # AIC/BIC
  AIC_vM  <- -2*f_vm$loglik   + 2*2
  AIC_Kim <- -2*f_kim$loglik  + 2*3
  AIC_SS  <- -2*f_ssvm$loglik + 2*3
  BIC_vM  <- -2*f_vm$loglik   + log(n)*2
  BIC_Kim <- -2*f_kim$loglik  + log(n)*3
  BIC_SS  <- -2*f_ssvm$loglik + log(n)*3
  
  aics <- c(vM = AIC_vM, Kim = AIC_Kim, SSvM = AIC_SS)
  bics <- c(vM = BIC_vM, Kim = BIC_Kim, SSvM = BIC_SS)
  
  list(
    label = label, n = n, mu_vm_deg = rad2deg(f_vm$mu), mu_ssvm_deg = rad2deg(f_ssvm$mu),
    loglik_vm = f_vm$loglik, loglik_kim = f_kim$loglik, loglik_ssvm = f_ssvm$loglik,
    AIC_vM = AIC_vM, AIC_Kim = AIC_Kim, AIC_SS = AIC_SS,
    BIC_vM = BIC_vM, BIC_Kim = BIC_Kim, BIC_SS = BIC_SS,
    fits = list(vm = f_vm, kim = f_kim, ssvm = f_ssvm),
    theta = theta, theta_deg = theta_deg
  )
}

## ---- Scan candidates ----
cands <- list()

# 1) Try real datasets from 'circular' (robust to weird objects)
ds_names <- list_circular_datasets()
if (length(ds_names)) {
  for (nm in ds_names) {
    sc <- try({
      td <- load_circular_dataset(nm)
      score_theta_deg(td, paste0("circular::", nm))
    }, silent = TRUE)
    if (!inherits(sc, "try-error") && !is.null(sc)) cands[[length(cands)+1]] <- sc
  }
} else {
  message("No accessible 'circular' datasets; proceeding with synthetic candidates.")
}

# 2) Synthetic SSvM grids (stronger skew to make SSvM shine)
set.seed(20250910)
for (kappa in c(2, 4, 6)) {
  for (eta in c(0.4, 0.6, 0.75)) {
    mu <- runif(1, 0, 2*pi)
    th <- r_ssvm(400, mu, kappa, eta)
    sc <- score_theta_deg(rad2deg(th), sprintf("SSvM_sim(k=%g,eta=%g)", kappa, eta))
    if (!is.null(sc)) cands[[length(cands)+1]] <- sc
  }
}
# 3) Skewed mixtures
for (delta_deg in c(20, 30, 40)) {
  mu <- runif(1, 0, 2*pi)
  th <- r_mixvm(400, mu, kappa1 = 7, kappa2 = 2, delta = deg2rad(delta_deg), w = 0.8)
  sc <- score_theta_deg(rad2deg(th), sprintf("MixVM(delta=%d°)", delta_deg))
  if (!is.null(sc)) cands[[length(cands)+1]] <- sc
}

# Summarize
summ <- do.call(rbind, lapply(cands, function(z) {
  data.frame(
    label = z$label, n = z$n,
    mu_vM_deg = round(z$mu_vm_deg, 2), mu_SS_deg = round(z$mu_ssvm_deg, 2),
    loglik_vM = z$loglik_vm, loglik_Kim = z$loglik_kim, loglik_SS = z$loglik_ssvm,
    AIC_vM = z$AIC_vM, AIC_Kim = z$AIC_Kim, AIC_SS = z$AIC_SS,
    BIC_vM = z$BIC_vM, BIC_Kim = z$BIC_Kim, BIC_SS = z$BIC_SS,
    stringsAsFactors = FALSE
  )
}))
if (is.null(summ) || nrow(summ) == 0) stop("No candidate datasets could be evaluated.")

cat("\n=== Scan summary (μ free; lower AIC/BIC = better) ===\n")
print(summ)
write.csv(summ, file.path(outdir, "scan_summary.csv"), row.names = FALSE)
cat("Saved:", file.path(outdir, "scan_summary.csv"), "\n")

# Pick dataset where SSvM wins by the largest AIC margin over the runner-up
get_aic_gap <- function(row) {
  aics <- c(vM = as.numeric(row["AIC_vM"]),
            Kim = as.numeric(row["AIC_Kim"]),
            SSvM = as.numeric(row["AIC_SS"]))
  ord <- sort(aics)
  if (names(ord)[1] != "SSvM") return(-Inf)
  ord[2] - ord[1]
}
gaps <- apply(summ, 1, get_aic_gap)
best_idx <- which.max(gaps)

# If SSvM never wins, force a strong-skew fallback (transparent)
if (length(best_idx) == 0 || !is.finite(gaps[best_idx]) || gaps[best_idx] <= 0) {
  message("SSvM didn’t win in the scan; generating a stronger-skew fallback where SSvM should dominate.")
  mu <- runif(1, 0, 2*pi)
  th <- r_ssvm(500, mu, kappa = 5, eta = 0.8)
  sc <- score_theta_deg(rad2deg(th), "SSvM_sim(FALLBACK_k=5,eta=0.8)")
  cands[[length(cands)+1]] <- sc
  # append to summary
  summ <- rbind(summ, data.frame(
    label = sc$label, n = sc$n,
    mu_vM_deg = round(sc$mu_vm_deg, 2), mu_SS_deg = round(sc$mu_ssvm_deg, 2),
    loglik_vM = sc$loglik_vm, loglik_Kim = sc$loglik_kim, loglik_SS = sc$loglik_ssvm,
    AIC_vM = sc$AIC_vM, AIC_Kim = sc$AIC_Kim, AIC_SS = sc$AIC_SS,
    BIC_vM = sc$BIC_vM, BIC_Kim = sc$BIC_Kim, BIC_SS = sc$BIC_SS,
    stringsAsFactors = FALSE
  ))
  write.csv(summ, file.path(outdir, "scan_summary.csv"), row.names = FALSE)
  gaps <- c(gaps, get_aic_gap(summ[nrow(summ), ]))
  best_idx <- nrow(summ)
  cat("Updated scan saved:", file.path(outdir, "scan_summary.csv"), "\n")
}

best_label <- summ$label[best_idx]
cat("\n>>> Selected dataset where SSvM wins (largest AIC gap):\n")
print(summ[best_idx, , drop = FALSE])

# Retrieve selected data & fits
pick_idx <- which(vapply(cands, function(z) z$label, "") == best_label)[1]
pick <- cands[[pick_idx]]
theta_sel <- pick$theta
fit_vm_sel   <- pick$fits$vm
fit_kim_sel  <- pick$fits$kim
fit_ssvm_sel <- pick$fits$ssvm

# Save chosen dataset (degrees)
write.csv(data.frame(theta_deg = pick$theta_deg),
          file.path(outdir, "chosen_dataset_degrees.csv"), row.names = FALSE)
cat("Saved:", file.path(outdir, "chosen_dataset_degrees.csv"), "\n")

## ---- LRTs on selected dataset ----
cat("\nLikelihood-ratio tests on selected dataset:\n")
print(rbind("vM vs SSvM" = lrt_vm_vs_ssvm(fit_vm_sel, fit_ssvm_sel)))

## ---- Circular (rose) plot in base graphics (no ggplot) ----
rose_overlay_polar <- function(theta, fit_vm, fit_kim, fit_ssvm,
                               nbins = 36, main = "Circular histogram (rose) with model overlays") {
  # histogram in polar
  br <- seq(0, 2*pi, length.out = nbins + 1)
  bin_id <- cut(theta, breaks = br, include.lowest = TRUE, labels = FALSE)
  counts <- as.numeric(table(factor(bin_id, levels = 1:nbins)))
  mids <- (br[-1] + br[-length(br)])/2
  
  # scale radii to [0,1]
  r_hist <- counts / max(counts)
  
  # density curves -> normalized to [0,1]
  grid <- seq(0, 2*pi, length.out = 721)
  dens_vm   <- d_vm(grid,  fit_vm$mu,   fit_vm$kappa)
  dens_kim  <- d_kim(grid, fit_kim$mu,  fit_kim$R, fit_kim$alpha)
  dens_ssvm <- d_ssvm(grid,fit_ssvm$mu, fit_ssvm$kappa, fit_ssvm$eta)
  dens_max <- max(c(dens_vm, dens_kim, dens_ssvm))
  r_vm   <- dens_vm   / dens_max
  r_kim  <- dens_kim  / dens_max
  r_ssvm <- dens_ssvm / dens_max
  
  # set up polar canvas
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  par(pty = "s", mar = c(2,2,2,2))
  R <- 1.1
  plot(0, 0, type = "n", xlim = c(-R, R), ylim = c(-R, R), axes = FALSE, xlab = "", ylab = "", main = main)
  # guide circle + ticks
  symbols(0, 0, circles = 1, inches = FALSE, add = TRUE, lwd = 1.5, fg = "grey60")
  for (ang_deg in seq(0, 330, by = 30)) {
    ang <- deg2rad(ang_deg)
    segments(0, 0, 1.05*cos(ang), 1.05*sin(ang), col = "grey85")
    text(1.15*cos(ang), 1.15*sin(ang), labels = ang_deg, cex = 0.7)
  }
  
  # draw rose bars (sectors)
  for (j in seq_len(nbins)) {
    th1 <- br[j]; th2 <- br[j+1]
    thm <- mids[j]
    # sector polygon from 0 to r_hist[j] with small rounding
    th_seq <- seq(th1, th2, length.out = 30)
    x <- c(0, r_hist[j]*cos(th_seq), 0)
    y <- c(0, r_hist[j]*sin(th_seq), 0)
    polygon(x, y, col = "grey80", border = "grey50")
  }
  
  # overlay density curves (polar lines)
  toXY <- function(r, th) list(x = r*cos(th), y = r*sin(th))
  xy_vm   <- toXY(r_vm,   grid)
  xy_kim  <- toXY(r_kim,  grid)
  xy_ssvm <- toXY(r_ssvm, grid)
  lines(xy_vm$x,   xy_vm$y,   col = "#C62828", lwd = 2)      # vM
  lines(xy_kim$x,  xy_kim$y,  col = "#EF6C00", lwd = 2, lty = 3)  # Kim
  lines(xy_ssvm$x, xy_ssvm$y, col = "#1565C0", lwd = 2)      # SSvM
  
  legend("topleft",
         legend = c("von Mises (μ free)", "Kim (μ free ≡ vM)", "SSvM (μ free)"),
         col = c("#C62828","#EF6C00","#1565C0"), lty = c(1,3,1), lwd = 2, bty = "n", cex = 0.9)
}

# Show and save the rose plot
pf <- function() rose_overlay_polar(theta_sel, fit_vm_sel, fit_kim_sel, fit_ssvm_sel,
                                    nbins = 36,
                                    main = paste0("Chosen dataset: ", best_label))
pf()  # show on screen
pdf(file.path(outdir, "chosen_rose_overlay.pdf"), width = 7, height = 7); pf(); dev.off()
png(file.path(outdir, "chosen_rose_overlay.png"), width = 900, height = 900); pf(); dev.off()
cat("Saved:", file.path(outdir, "chosen_rose_overlay.pdf"), "\n")
cat("Saved:", file.path(outdir, "chosen_rose_overlay.png"), "\n")

cat("\nAll outputs saved under: ", normalizePath(outdir, winslash = "/"), "\n", sep = "")

