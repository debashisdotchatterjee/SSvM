# =========================================================
# Polar density rings for vM (multi-kappa) and SSvM (multi-eta)
# Base R only. Draws to screen and saves two PDFs without overlapping legends.
# =========================================================

## --------- setup & helpers ----------
I0 <- function(x) besselI(x, nu = 0)

vm_pdf <- function(theta, mu, kappa) {
  exp(kappa * cos(theta - mu)) / (2*pi*I0(kappa))
}

ssvm_pdf <- function(theta, mu, kappa, eta) {
  vm_pdf(theta, mu, kappa) * (1 + eta * sin(theta - mu))
}

# Angle grid (fine)
n_seg <- 960
theta <- seq(0, 2*pi, length.out = n_seg + 1)[-1]  # drop duplicate 2π

# Convert "math" theta to screen angle so that:
#   0° is at North and angles increase clockwise
to_screen_angle <- function(theta_rad) (pi/2) - theta_rad

# Color map (cool -> warm), dependency-free
make_col_fun <- function() {
  function(v) {
    v <- pmax(0, pmin(1, v))
    hsv(h = (2/3) * (1 - v), s = 1, v = 1)
  }
}
col_fun <- make_col_fun()

# Draw one colorful density ring as a set of thin quadrilaterals
draw_density_ring <- function(theta, dens, r_inner, band, border = NA) {
  dmax    <- max(dens, na.rm = TRUE)
  r_norm  <- if (dmax > 0) dens / dmax else dens
  r_outer <- r_inner + band * r_norm
  cols    <- col_fun(r_norm)
  
  th1 <- theta
  th2 <- c(theta[-1], theta[1])  # wrap
  a1  <- to_screen_angle(th1)
  a2  <- to_screen_angle(th2)
  
  x1i <- r_inner * cos(a1); y1i <- r_inner * sin(a1)
  x2i <- r_inner * cos(a2); y2i <- r_inner * sin(a2)
  x1o <- r_outer * cos(a1); y1o <- r_outer * sin(a1)
  x2o <- r_outer * cos(a2); y2o <- r_outer * sin(a2)
  
  for (i in seq_along(theta)) {
    polygon(
      x = c(x1i[i], x2i[i], x2o[i], x1o[i]),
      y = c(y1i[i], y2i[i], y2o[i], y1o[i]),
      col = cols[i], border = border
    )
  }
}

# Polar canvas (unit circle), ticks each 30°
polar_canvas <- function(title = NULL, subtitle = NULL,
                         xlim = c(-1.55, 1.65), ylim = c(-1.1, 1.1)) {
  plot.new()
  plot.window(xlim = xlim, ylim = ylim, asp = 1)
  
  # Outer circle
  symbols(0, 0, circles = 1, inches = FALSE, add = TRUE,
          bg = NA, fg = "grey60", lwd = 1)
  
  # Angular spokes every 30°
  for (d in seq(0, 330, by = 30)) {
    th <- d * pi/180
    a  <- to_screen_angle(th)
    segments(0, 0, cos(a), sin(a), col = "grey85", lwd = 1)
    # Tick labels just outside the circle
    tx <- 1.06 * cos(a); ty <- 1.06 * sin(a)
    text(tx, ty, labels = paste0(d, "°"), cex = 0.8, col = "grey25")
  }
  
  if (!is.null(title))    mtext(title,    side = 3, line = 1.1, cex = 1.1, font = 2)
  if (!is.null(subtitle)) mtext(subtitle, side = 3, line = 0.2, cex = 0.9, font = 1)
}

# Arrow at mu
draw_mu_arrow <- function(mu, len = 1.05, lwd = 2, col = "black") {
  a <- to_screen_angle(mu)
  arrows(0, 0, len * cos(a), len * sin(a), length = 0.12, lwd = lwd, col = col)
}

# Vertical colorbar (0..1), placed OUTSIDE the circle on the far right.
draw_colorbar <- function(xleft, ybottom, width = 0.12, height = 0.90, n = 220,
                          title = "Normalized density") {
  y <- seq(0, 1, length.out = n)
  cols <- col_fun(y)
  # white background pad so it never "touches" the plot
  rect(xleft - 0.03, ybottom - 0.05, xleft + width + 0.15, ybottom + height + 0.10,
       col = "white", border = NA)
  for (i in seq_len(n)) {
    rect(
      xleft, ybottom + (i-1) * height/n,
      xleft + width, ybottom + i * height/n,
      col = cols[i], border = NA
    )
  }
  rect(xleft, ybottom, xleft + width, ybottom + height, border = "grey30")
  # ticks
  axis_ticks <- seq(0, 1, by = 0.25)
  for (t in axis_ticks) {
    yy <- ybottom + t * height
    segments(xleft + width, yy, xleft + width + 0.02, yy, col = "grey30")
    text(xleft + width + 0.05, yy, labels = sprintf("%.2f", t), adj = 0, cex = 0.85)
  }
  text(xleft + width/2, ybottom + height + 0.08, title, cex = 0.85, font = 2)
}

## --------- parameters ----------
mu     <- pi/2
kappas <- c(2, 5, 10)      # von Mises panel
kappaR <- 5                # SSvM fixed kappa
etas   <- c(-0.4, 0, 0.4)  # SSvM panel

# Precompute densities
dens_vm_list <- lapply(kappas, function(k) vm_pdf(theta, mu, k))
dens_ss_list <- lapply(etas,   function(e) ssvm_pdf(theta, mu, kappaR, e))

# Ring placements (concentric; non-overlapping)
rinners <- c(0.20, 0.48, 0.76)
bands   <- c(0.24, 0.24, 0.20)

# --------- console tables (quick reference) ----------
vm_tbl  <- data.frame(kappa = kappas,
                      max_density = sapply(dens_vm_list, max))
ss_tbl  <- data.frame(eta = etas,
                      kappa_fixed = kappaR,
                      max_density = sapply(dens_ss_list, max))
cat("\nvon Mises settings & max densities:\n"); print(vm_tbl, row.names = FALSE)
cat("\nSSvM settings (kappa fixed) & max densities:\n"); print(ss_tbl, row.names = FALSE)

# --------- draw: von Mises panel (screen) ----------
par(mar = c(1.6, 1.6, 4.0, 6.5), xpd = NA)  # big right margin; draw outside allowed
polar_canvas(
  title    = "von Mises:  \u03bc = \u03c0/2  (multiple \u03ba)",
  subtitle = "\u03ba \u2208 {2, 5, 10}",
  xlim = c(-1.55, 1.75)   # extra width to the right for colorbar
)
for (i in seq_along(kappas)) {
  draw_density_ring(theta, dens_vm_list[[i]], r_inner = rinners[i], band = bands[i])
  # Ring labels OUTSIDE, top-left gutter (won't overlap)
  lbl <- bquote(kappa == .(kappas[i]))
  text(-1.45, 1.00 - 0.15*(i-1), labels = lbl, adj = 0, cex = 0.95)
}
draw_mu_arrow(mu, len = 1.05)
draw_colorbar(xleft = 1.44, ybottom = -0.45, width = 0.12, height = 0.90,
              title = "Norm. density")

# --------- draw: SSvM panel (screen) ----------
par(mar = c(1.6, 1.6, 4.0, 6.5), xpd = NA)
polar_canvas(
  title    = "SSvM:  \u03bc = \u03c0/2,  \u03ba = 5  (vary \u03b7)",
  subtitle = "\u03b7 \u2208 {-0.4, 0, 0.4}",
  xlim = c(-1.55, 1.75)
)
for (i in seq_along(etas)) {
  draw_density_ring(theta, dens_ss_list[[i]], r_inner = rinners[i], band = bands[i])
  # Ring labels OUTSIDE, top-left gutter (won't overlap)
  lbl <- bquote(eta == .(etas[i]))
  text(-1.45, 1.00 - 0.15*(i-1), labels = lbl, adj = 0, cex = 0.95)
}
draw_mu_arrow(mu, len = 1.05)
draw_colorbar(xleft = 1.44, ybottom = -0.45, width = 0.12, height = 0.90,
              title = "Norm. density")

# --------- save both panels to PDFs ----------
out_dir <- "paper_outputs/simulations"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

file_vm   <- file.path(out_dir, "curve_vm_multi.pdf")
file_ssvm <- file.path(out_dir, "curve_ssvm_multi.pdf")

# Save von Mises panel
pdf(file_vm, width = 6.25, height = 6.25)  # square page; roomy margins
par(mar = c(1.6, 1.6, 4.0, 6.5), xpd = NA)
polar_canvas(
  title    = "von Mises:  \u03bc = \u03c0/2  (multiple \u03ba)",
  subtitle = "\u03ba \u2208 {2, 5, 10}",
  xlim = c(-1.55, 1.75)
)
for (i in seq_along(kappas)) {
  draw_density_ring(theta, dens_vm_list[[i]], r_inner = rinners[i], band = bands[i])
  lbl <- bquote(kappa == .(kappas[i]))
  text(-1.45, 1.00 - 0.15*(i-1), labels = lbl, adj = 0, cex = 0.95)
}
draw_mu_arrow(mu, len = 1.05)
draw_colorbar(xleft = 1.44, ybottom = -0.45, width = 0.12, height = 0.90,
              title = "Norm. density")
dev.off()

# Save SSvM panel
pdf(file_ssvm, width = 6.25, height = 6.25)
par(mar = c(1.6, 1.6, 4.0, 6.5), xpd = NA)
polar_canvas(
  title    = "SSvM:  \u03bc = \u03c0/2,  \u03ba = 5  (vary \u03b7)",
  subtitle = "\u03b7 \u2208 {-0.4, 0, 0.4}",
  xlim = c(-1.55, 1.75)
)
for (i in seq_along(etas)) {
  draw_density_ring(theta, dens_ss_list[[i]], r_inner = rinners[i], band = bands[i])
  lbl <- bquote(eta == .(etas[i]))
  text(-1.45, 1.00 - 0.15*(i-1), labels = lbl, adj = 0, cex = 0.95)
}
draw_mu_arrow(mu, len = 1.05)
draw_colorbar(xleft = 1.44, ybottom = -0.45, width = 0.12, height = 0.90,
              title = "Norm. density")
dev.off()

cat("\nSaved PDFs:\n  - ", file_vm, "\n  - ", file_ssvm, "\n", sep = "")

