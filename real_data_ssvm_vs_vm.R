# ===========================
# Real-data SSvM vs von Mises
# ===========================

# ---- Packages ----
req <- c(
  "circular","CircStats","Directional","EnvStats","openair",
  "dplyr","tidyr","purrr","stringr","tibble","readr","ggplot2",
  "ggthemes","viridis","scales","patchwork","numDeriv","gridExtra","cowplot","kableExtra"
)
miss <- setdiff(req, rownames(installed.packages()))
if (length(miss)) install.packages(miss, quiet = TRUE)
invisible(lapply(req, library, character.only = TRUE))
theme_set(ggplot2::theme_minimal(base_size = 12))

# ---- Output dirs ----
root_out <- file.path(getwd(), "paper_outputs", "realdata")
dir.create(root_out, recursive = TRUE, showWarnings = FALSE)
fig_dir <- file.path(root_out, "figures"); dir.create(fig_dir, showWarnings = FALSE)
tab_dir <- file.path(root_out, "tables");  dir.create(tab_dir, showWarnings = FALSE)

cat("Output directories:\n",
    "  Figures ->", fig_dir, "\n",
    "  Tables  ->", tab_dir, "\n\n")

# ---------------------------
# Helper: tidy printing to console
# ---------------------------
print_tbl <- function(df, caption = NULL) {
  knitr::kable(df, format = "simple", caption = caption, digits = 4)
}

save_csv_and_print <- function(df, path, caption = NULL) {
  readr::write_csv(df, path)
  cat("Saved table:", path, "\n")
  print(print_tbl(df, caption))
  cat("\n")
}

# ---------------------------
# 1) Dataset auto-discovery
# ---------------------------
candidate_pkgs <- c("circular","CircStats","Directional","EnvStats","openair","amt","moveHMM","circularGLM")

find_circular_dataset <- function() {
  patt <- "(wind|direction|bearing|angle|orientation)"
  for (pk in candidate_pkgs) {
    suppressWarnings(ds <- try(utils::data(package = pk)$results, silent = TRUE))
    if (inherits(ds, "try-error") || is.null(ds) || !NROW(ds)) next
    hits <- subset(as.data.frame(ds, stringsAsFactors = FALSE),
                   grepl(patt, Item, ignore.case = TRUE) | grepl(patt, Title, ignore.case = TRUE))
    if (NROW(hits)) {
      # Try each hit in order and attempt to extract an angle vector
      for (i in seq_len(NROW(hits))) {
        dsname <- hits$Item[i]
        # Snapshot objects before/after load to detect what's created
        before <- ls(envir = .GlobalEnv, all.names = TRUE)
        ok <- try(utils::data(list = dsname, package = pk, envir = .GlobalEnv), silent = TRUE)
        after <- ls(envir = .GlobalEnv, all.names = TRUE)
        if (inherits(ok, "try-error")) next
        new_objs <- setdiff(after, before)
        # Try to locate angle data inside any new object
        for (obj in new_objs) {
          x <- get(obj, envir = .GlobalEnv)
          # Direct numeric vector?
          if (is.numeric(x) && length(x) >= 60 && all(is.finite(x))) {
            # Heuristics: is it degrees, radians, etc. Keep as-is; conversion later.
            return(list(pkg = pk, dataset = dsname, obj = obj, data = x, source = "numeric-vector"))
          }
          # Data frame / tibble?
          if (is.data.frame(x)) {
            cn <- names(x)
            angle_cols <- cn[grepl("(dir|bearing|angle|orient|theta|phi)", cn, ignore.case = TRUE)]
            # Prefer single angle-like column
            if (length(angle_cols)) {
              # Pick the first that looks numeric with many finite values
              for (cc in angle_cols) {
                v <- suppressWarnings(as.numeric(x[[cc]]))
                if (sum(is.finite(v)) >= 60) {
                  return(list(pkg = pk, dataset = dsname, obj = obj, data = v, source = paste0("data.frame$", cc)))
                }
              }
            }
          }
          # List with obvious angle component?
          if (is.list(x)) {
            nm <- names(x)
            angle_elt <- nm[grepl("(dir|bearing|angle|orient|theta|phi)", nm, ignore.case = TRUE)]
            if (length(angle_elt)) {
              for (kk in angle_elt) {
                v <- suppressWarnings(as.numeric(x[[kk]]))
                if (sum(is.finite(v)) >= 60) {
                  return(list(pkg = pk, dataset = dsname, obj = obj, data = v, source = paste0("list$", kk)))
                }
              }
            }
          }
        }
        # If we get here, this hit didn't yield usable angles. Clean & try next hit.
        rm(list = new_objs, envir = .GlobalEnv)
      }
    }
  }
  stop("No suitable real dataset with angle-like measurements was found among installed packages.\n",
       "Try installing `circular`, `CircStats`, `Directional`, or `openair` and re-run.")
}

found <- find_circular_dataset()

cat("DATASET FOUND\n",
    "  Package : ", found$pkg, "\n",
    "  Dataset : ", found$dataset, "\n",
    "  Object  : ", found$obj, "\n",
    "  Source  : ", found$source, "\n\n")

raw_angles <- as.numeric(found$data)
raw_angles <- raw_angles[is.finite(raw_angles)]

# ---------------------------
# 2) Coerce to radians in [0, 2π)
# ---------------------------
to_radians <- function(v) {
  # Heuristic: if many values > 2*pi, assume degrees
  if (stats::median(v, na.rm = TRUE) > 6.5) {
    v <- (v %% 360) * pi / 180
  } else {
    # Already radians? wrap to [0, 2π)
    v <- (v %% (2*pi))
  }
  v
}
theta <- to_radians(raw_angles)
theta <- theta[is.finite(theta)]
n <- length(theta)

cat("Sample size used (finite angles):", n, "\n\n")
stopifnot(n >= 60)

# ---------------------------
# 3) Utilities: von Mises & SSvM
# ---------------------------
A1_inv <- function(R) {
  # Best-pieces approximation for inverse A1
  if (R < 0.53) {
    2*R + R^3 + (5*R^5)/6
  } else if (R < 0.85) {
    -0.4 + 1.39*R + 0.43/(1 - R)
  } else {
    1/(R^3 - 4*R^2 + 3*R)
  }
}

# First and second sample trigonometric moments
Cbar <- mean(cos(theta)); Sbar <- mean(sin(theta))
Rbar <- sqrt(Cbar^2 + Sbar^2)
mu_init <- atan2(Sbar, Cbar) %% (2*pi)
kappa_init <- max(1e-8, A1_inv(Rbar))

# SSvM one-step Newton start for eta at mu_init
S <- sin(theta - mu_init)
eta_init <- sum(S) / sum(S^2)
eta_init <- max(min(eta_init, 0.5), -0.5)

# Density helpers
d_vm <- function(th, mu, kappa, log = FALSE) {
  z <- kappa * cos(th - mu)
  const <- -log(2*pi) - log(besselI(kappa, nu = 0, expon.scaled = FALSE))
  if (log) const + z else exp(const + z)
}

d_ssvm <- function(th, mu, kappa, eta, log = FALSE) {
  # Ensure positivity 1 + eta sin(.) > 0
  tilt <- 1 + eta * sin(th - mu)
  if (any(tilt <= 0)) {
    if (log) return(rep(-Inf, length(th))) else return(rep(0, length(th)))
  }
  log_base <- -log(2*pi) - log(besselI(kappa, nu = 0)) + kappa * cos(th - mu)
  out <- log_base + log(tilt)
  if (log) out else exp(out)
}

# Log-likelihoods (sum)
ll_vm <- function(mu, kappa) sum(d_vm(theta, mu, kappa, log = TRUE))
ll_ssvm_core <- function(mu, kappa, eta) sum(d_ssvm(theta, mu, kappa, eta, log = TRUE))

# Safe transforms for optimization
par_vm_pack   <- function(mu, kappa) c(mu = (mu %% (2*pi)), xi = log(kappa))
par_vm_unpack <- function(par) list(mu = (par[1] %% (2*pi)), kappa = exp(par[2]))

par_ss_pack   <- function(mu, kappa, eta) c(mu = (mu %% (2*pi)), xi = log(kappa), zeta = atanh(eta))
par_ss_unpack <- function(par) list(mu = (par[1] %% (2*pi)), kappa = exp(par[2]), eta = tanh(par[3]))

# Objective for optim
nll_vm <- function(par) {
  p <- par_vm_unpack(par)
  -ll_vm(p$mu, p$kappa)
}
nll_ss <- function(par) {
  p <- par_ss_unpack(par)
  # penalize slight violations to keep optimizer stable near boundary
  if (abs(p$eta) >= 0.999) return(1e6)
  base <- d_ssvm(theta, p$mu, p$kappa, p$eta, log = TRUE)
  -sum(base)
}

# ---------------------------
# 4) Fit models (MLE)
# ---------------------------
set.seed(1)
fit_vm <- optim(
  par = par_vm_pack(mu_init, max(kappa_init, 1e-6)),
  fn  = nll_vm,
  method = "BFGS",
  control = list(reltol = 1e-10, maxit = 2000)
)

stopifnot(fit_vm$convergence == 0)
p_vm <- par_vm_unpack(fit_vm$par)
ll_vm_hat <- -fit_vm$value

# SSvM: robust multi-start in zeta
starts_eta <- unique(sign(c(-1,0,1)) * c(0, 0.25, 0.5, abs(eta_init)))
fits <- list(); vals <- numeric()
for (z0 in atanh(pmax(pmin(starts_eta, 0.95), -0.95))) {
  st <- c(mu = mu_init, xi = log(max(kappa_init, 1e-6)), zeta = z0)
  of <- try(optim(
    par = st, fn = nll_ss, method = "BFGS",
    control = list(reltol = 1e-10, maxit = 4000)
  ), silent = TRUE)
  if (!inherits(of, "try-error") && of$convergence == 0) {
    fits[[length(fits)+1]] <- of
    vals[length(vals)+1] <- of$value
  }
}
stopifnot(length(fits) >= 1)
best <- fits[[ which.min(vals) ]]
p_ss <- par_ss_unpack(best$par)
ll_ssvm_hat <- -best$value

# ---------------------------
# 5) SEs via numDeriv on transformed scale -> delta method
# ---------------------------
vm_hess <- numDeriv::hessian(nll_vm, fit_vm$par)
vm_cov  <- tryCatch(solve(vm_hess), error = function(e) { matrix(NA, 2, 2) })
# Jacobian diag(1, kappa) for (mu, xi) -> (mu, kappa)
J_vm <- diag(c(1, p_vm$kappa))
vm_cov_orig <- J_vm %*% vm_cov %*% J_vm
vm_se_mu <- sqrt(vm_cov_orig[1,1])
vm_se_k  <- sqrt(vm_cov_orig[2,2])

ss_hess <- numDeriv::hessian(nll_ss, best$par)
ss_cov  <- tryCatch(solve(ss_hess), error = function(e) { matrix(NA, 3, 3) })
# Jacobian diag(1, kappa, 1-eta^2) for (mu, xi, zeta) -> (mu, kappa, eta)
J_ss <- diag(c(1, p_ss$kappa, 1 - p_ss$eta^2))
ss_cov_orig <- J_ss %*% ss_cov %*% J_ss
ss_se_mu <- sqrt(ss_cov_orig[1,1])
ss_se_k  <- sqrt(ss_cov_orig[2,2])
ss_se_eta<- sqrt(ss_cov_orig[3,3])

# ---------------------------
# 6) Model comparison (AIC/BIC, LR test)
# ---------------------------
aic_vm   <- 2*2  - 2*ll_vm_hat
aic_ssvm <- 2*3  - 2*ll_ssvm_hat
bic_vm   <- log(n)*2 - 2*ll_vm_hat
bic_ssvm <- log(n)*3 - 2*ll_ssvm_hat

lr        <- 2*(ll_ssvm_hat - ll_vm_hat)
p_lr      <- 1 - pchisq(lr, df = 1)

est_table <- tibble::tibble(
  model = c("von Mises", "SSvM"),
  mu_hat = c(p_vm$mu, p_ss$mu),
  se_mu  = c(vm_se_mu, ss_se_mu),
  kappa_hat = c(p_vm$kappa, p_ss$kappa),
  se_kappa  = c(vm_se_k,  ss_se_k),
  eta_hat   = c(NA, p_ss$eta),
  se_eta    = c(NA, ss_se_eta),
  logLik    = c(ll_vm_hat, ll_ssvm_hat),
  AIC       = c(aic_vm, aic_ssvm),
  BIC       = c(bic_vm, bic_ssvm)
)

comp_table <- tibble::tibble(
  test = "LR: SSvM vs von Mises",
  df   = 1L,
  stat = lr,
  pval = p_lr
)

save_csv_and_print(est_table, file.path(tab_dir, "estimates_comparison.csv"),
                   "Parameter estimates and information criteria")
save_csv_and_print(comp_table, file.path(tab_dir, "lr_test.csv"),
                   "Likelihood-ratio test (df=1)")

# ---------------------------
# 7) Plots: rose + overlaid density curves; density-on-circle
# ---------------------------
# Create a tidy angle data frame (degrees for nicer axis breaks)
deg <- (theta * 180/pi) %% 360
df_angles <- tibble::tibble(angle_deg = deg)

# Compute model curves at fine grid (degrees)
grid_deg <- seq(0, 359, by = 1)
grid_rad <- grid_deg * pi/180
vm_pdf   <- d_vm(grid_rad, p_vm$mu, p_vm$kappa, log = FALSE)
ss_pdf   <- d_ssvm(grid_rad, p_ss$mu, p_ss$kappa, p_ss$eta, log = FALSE)

# For rose overlay: scale density to histogram height
# Make a reference histogram to extract max count
h <- hist(deg, breaks = seq(0, 360, by = 10), plot = FALSE, include.lowest = TRUE)
max_count <- max(h$counts); max_pdf <- max(c(vm_pdf, ss_pdf))
scale_fac <- max_count / max_pdf

df_curves <- tibble::tibble(
  angle_deg = grid_deg,
  vm_scaled = vm_pdf * scale_fac,
  ss_scaled = ss_pdf * scale_fac
)

# Rose plot with overlay (10-degree bins)
p_rose <- ggplot() +
  geom_histogram(data = df_angles, aes(x = angle_deg, y = ..count..),
                 breaks = seq(0, 360, by = 10), alpha = 0.6) +
  geom_line(data = df_curves, aes(x = angle_deg, y = vm_scaled),
            linewidth = 0.9) +
  geom_line(data = df_curves, aes(x = angle_deg, y = ss_scaled),
            linewidth = 0.9, linetype = "dashed") +
  coord_polar(start = -pi/2, direction = -1) + # 0° at North; clockwise
  labs(title = "Rose diagram with overlaid fitted densities",
       subtitle = sprintf("Dataset: %s::%s (%s)\nSolid=vM; Dashed=SSvM",
                          found$pkg, found$dataset, found$source),
       x = "Direction (deg)", y = "Count") +
  theme_minimal(base_size = 12)

ggsave(file.path(fig_dir, "rose_overlay_vm_ssvm.png"), p_rose, width = 7.5, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "rose_overlay_vm_ssvm.pdf"), p_rose, width = 7.5, height = 6)
cat("Saved rose overlay plots:\n",
    "  -", file.path(fig_dir, "rose_overlay_vm_ssvm.png"), "\n",
    "  -", file.path(fig_dir, "rose_overlay_vm_ssvm.pdf"), "\n\n")

# Density on circle (radius = scaled density)
df_ring <- df_curves |>
  tidyr::pivot_longer(cols = c(vm_scaled, ss_scaled), names_to = "model", values_to = "rad") |>
  dplyr::mutate(model = dplyr::recode(model, vm_scaled = "vM", ss_scaled = "SSvM"))

p_ring <- ggplot(df_ring, aes(x = angle_deg, y = rad, linetype = model)) +
  geom_path(linewidth = 0.9) +
  coord_polar(start = -pi/2, direction = -1) +
  labs(title = "Fitted densities (scaled) traced on the circle",
       subtitle = sprintf("Dataset: %s::%s", found$pkg, found$dataset),
       x = "Direction (deg)", y = "Scaled density") +
  theme_minimal(base_size = 12)

ggsave(file.path(fig_dir, "density_on_circle.png"), p_ring, width = 7.5, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "density_on_circle.pdf"), p_ring, width = 7.5, height = 6)
cat("Saved circle-density plots:\n",
    "  -", file.path(fig_dir, "density_on_circle.png"), "\n",
    "  -", file.path(fig_dir, "density_on_circle.pdf"), "\n\n")

# Diagnostics: QQ-style angle CDF vs model (optional, simple)
cdf_emp <- ecdf(deg)
cdf_grid <- tibble::tibble(
  angle_deg = grid_deg,
  F_emp = cdf_emp(grid_deg),
  F_vm  = (cumsum(d_vm(grid_rad, p_vm$mu, p_vm$kappa, log = FALSE))/sum(d_vm(grid_rad, p_vm$mu, p_vm$kappa, log = FALSE))),
  F_ss  = (cumsum(d_ssvm(grid_rad, p_ss$mu, p_ss$kappa, p_ss$eta, log = FALSE))/sum(d_ssvm(grid_rad, p_ss$mu, p_ss$kappa, p_ss$eta, log = FALSE)))
)
cdf_long <- cdf_grid |>
  tidyr::pivot_longer(cols = c(F_vm, F_ss), names_to = "model", values_to = "F_mod") |>
  dplyr::mutate(model = recode(model, F_vm="vM", F_ss="SSvM"))

p_cdf <- ggplot(cdf_long, aes(x = F_emp, y = F_mod, color = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  geom_path(linewidth = 0.9) +
  labs(title = "Empirical CDF vs model CDF",
       x = "Empirical CDF", y = "Model CDF") +
  theme_minimal(base_size = 12)

ggsave(file.path(fig_dir, "cdf_qq.png"), p_cdf, width = 6.8, height = 5.2, dpi = 300)
ggsave(file.path(fig_dir, "cdf_qq.pdf"), p_cdf, width = 6.8, height = 5.2)
cat("Saved CDF diagnostics:\n",
    "  -", file.path(fig_dir, "cdf_qq.png"), "\n",
    "  -", file.path(fig_dir, "cdf_qq.pdf"), "\n\n")

# ---------------------------
# 8) Console recap
# ---------------------------
recap <- tibble::tibble(
  package = found$pkg,
  dataset = found$dataset,
  object  = found$obj,
  source  = found$source,
  n       = n,
  vm_logLik = ll_vm_hat,
  ssvm_logLik = ll_ssvm_hat,
  LR = lr, df = 1, p_value = p_lr,
  AIC_vm = aic_vm, AIC_ssvm = aic_ssvm,
  BIC_vm = bic_vm, BIC_ssvm = bic_ssvm
)
save_csv_and_print(recap, file.path(tab_dir, "recap.csv"), "Model comparison recap")

cat("Analysis complete.\n",
    "Use the printed tables above for manuscript copy.\n",
    "Figures saved under:", fig_dir, "\n",
    "Tables  saved under:", tab_dir, "\n")
