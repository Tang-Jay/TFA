#============================================================
# SEL local asymptotic potential function
#
# Local parametrization:
#   mu* = mu0 + tau * sigma / sqrt(n), tau in R
#   tau <= 0 is the local null region; tau > 0 is the local alternative.
#
# Notation:
#   z_alpha2      = qnorm(1 - alpha2)
#   z_alpha1_over2 = qnorm(1 - alpha1 / 2)
#   chi_1_alpha1^2 = z_alpha1_over2^2 = qchisq(1 - alpha1, df = 1)
#   alpha          = nominal target size used to optionally calibrate alpha1
#
# SEL local asymptotic rejection probability / potential:
#   pi_s(tau) =
#   P{ U + V^2 1(V > z_alpha2) > chi_1_alpha1^2 },
#   U ~ chisq_1, V ~ N(tau, 1), U independent of V.
#
# ASEL local asymptotic rejection probability / potential:
#   pi_a(tau) =
#     pi_s(tau),                                      if mu0 = 0,
#     1 - Phi(z_alpha2 - tau) F_chisq1(chi_1_alpha1^2), if mu0 != 0.
#   Since chi_1_alpha1^2 is the upper alpha1 quantile,
#     F_chisq1(chi_1_alpha1^2) = 1 - alpha1.
#
# ASEL2 uses the same ASEL formula with a separate corrected alpha2.
#
# CEL local asymptotic rejection probability:
#   T_c -> V^2 1(V > 0), V ~ N(tau, 1).
#   With the mixture null critical value chi^2_{1,2 alpha},
#   pi_c(tau) = P{V > sqrt(chi^2_{1,2 alpha})}.
#
# CEL2 local asymptotic rejection probability:
#   T_cel2 = ell(0) 1(hat_mu > 0), i.e. CEL with fixed mu0 = 0.
#   With z_{2 alpha}(1) = chi^2_{1, 2 alpha},
#   pi_cel2(tau) = Phi(tau - sqrt(z_{2 alpha}(1))).
#============================================================

rm(list = ls())

#-----------------------------#
# Basic notation
#-----------------------------#
sel_constants <- function(alpha1 = 0.05, alpha2 = 0.08) {
  stopifnot(length(alpha1) == 1, length(alpha2) == 1)
  if (alpha1 <= 0 || alpha1 >= 1) stop("alpha1 must be in (0, 1).")
  if (alpha2 <= 0 || alpha2 >= 1) stop("alpha2 must be in (0, 1).")

  z_alpha2 <- qnorm(1 - alpha2)
  z_alpha1_over2 <- qnorm(1 - alpha1 / 2)
  chi_1_alpha1_sq <- qchisq(1 - alpha1, df = 1)

  list(
    alpha1 = alpha1,
    alpha2 = alpha2,
    z_alpha2 = z_alpha2,
    z_alpha1_over2 = z_alpha1_over2,
    chi_1_alpha1_sq = chi_1_alpha1_sq
  )
}

chisq1_upper_tail <- function(x) {
  pchisq(x, df = 1, lower.tail = FALSE)
}

tau_to_mu <- function(tau, mu0 = 0, sigma = 1, n = 1000) {
  if (sigma <= 0) stop("sigma must be positive.")
  if (n <= 0) stop("n must be positive.")
  mu0 + tau * sigma / sqrt(n)
}

mu_to_tau <- function(mu, mu0 = 0, sigma = 1, n = 1000) {
  if (sigma <= 0) stop("sigma must be positive.")
  if (n <= 0) stop("n must be positive.")
  sqrt(n) * (mu - mu0) / sigma
}

#-----------------------------#
# Direct integral formula
#-----------------------------#
sel_power_integral_one <- function(tau, alpha1 = 0.05, alpha2 = 0.08,
                                   rel.tol = 1e-10) {
  cc <- sel_constants(alpha1, alpha2)
  q <- cc$chi_1_alpha1_sq
  z2 <- cc$z_alpha2

  integrand <- function(v) {
    chisq1_upper_tail(q - v^2) * dnorm(v, mean = tau, sd = 1)
  }

  value <- alpha1 * pnorm(z2 - tau) +
    integrate(
      integrand,
      lower = z2,
      upper = Inf,
      rel.tol = rel.tol,
      subdivisions = 1000
    )$value

  min(max(value, 0), 1)
}

sel_power_integral <- function(tau, alpha1 = 0.05, alpha2 = 0.08,
                               rel.tol = 1e-10) {
  vapply(
    tau,
    sel_power_integral_one,
    numeric(1),
    alpha1 = alpha1,
    alpha2 = alpha2,
    rel.tol = rel.tol
  )
}

#-----------------------------#
# Piecewise formula in manuscript
#-----------------------------#
sel_power_piecewise_one <- function(tau, alpha1 = 0.05, alpha2 = 0.08,
                                    rel.tol = 1e-10) {
  cc <- sel_constants(alpha1, alpha2)
  z2 <- cc$z_alpha2
  zhalf <- cc$z_alpha1_over2
  q <- cc$chi_1_alpha1_sq

  if (alpha1 >= 2 * alpha2) {
    return(1 - (1 - alpha1) * pnorm(z2 - tau))
  }

  integrand <- function(v) {
    chisq1_upper_tail(q - v^2) * dnorm(v, mean = tau, sd = 1)
  }

  middle <- integrate(
    integrand,
    lower = z2,
    upper = zhalf,
    rel.tol = rel.tol,
    subdivisions = 1000
  )$value

  value <- alpha1 * pnorm(z2 - tau) +
    middle +
    1 - pnorm(zhalf - tau)

  min(max(value, 0), 1)
}

sel_power_piecewise <- function(tau, alpha1 = 0.05, alpha2 = 0.08,
                                rel.tol = 1e-10) {
  vapply(
    tau,
    sel_power_piecewise_one,
    numeric(1),
    alpha1 = alpha1,
    alpha2 = alpha2,
    rel.tol = rel.tol
  )
}

# Main wrapper.
# method = "piecewise" uses the simplified formula in the note.
# method = "integral" uses the original integral over [z_alpha2, Inf).
sel_power <- function(tau, alpha1 = 0.05, alpha2 = 0.08,
                      method = c("piecewise", "integral"),
                      rel.tol = 1e-10) {
  method <- match.arg(method)

  if (method == "piecewise") {
    sel_power_piecewise(tau, alpha1, alpha2, rel.tol)
  } else {
    sel_power_integral(tau, alpha1, alpha2, rel.tol)
  }
}

sel_boundary_power <- function(alpha1 = 0.05, alpha2 = 0.08,
                               method = c("piecewise", "integral")) {
  method <- match.arg(method)
  sel_power(0, alpha1 = alpha1, alpha2 = alpha2, method = method)
}

# Calibrate the chi-square critical tail alpha1 so that
#   pi_s(0; alpha1, alpha2) = target_alpha.
# SEL uses alpha1 = target_alpha.
# SEL2 uses the calibrated alpha1 below.
calibrate_sel_alpha1 <- function(target_alpha = 0.05, alpha2 = 0.08,
                                 rel.tol = 1e-10) {
  if (target_alpha <= 0 || target_alpha >= 1) {
    stop("target_alpha must be in (0, 1).")
  }
  if (alpha2 <= 0 || alpha2 >= 0.5) {
    stop("alpha2 must be in (0, 0.5) for this calibration formula.")
  }

  threshold <- 3 * alpha2 - 2 * alpha2^2

  if (target_alpha >= threshold) {
    return((target_alpha - alpha2) / (1 - alpha2))
  }

  f <- function(alpha1) {
    sel_power(0, alpha1 = alpha1, alpha2 = alpha2, method = "piecewise",
              rel.tol = rel.tol) - target_alpha
  }

  upper <- 2 * alpha2 * (1 - sqrt(.Machine$double.eps))
  uniroot(f, lower = .Machine$double.eps, upper = upper, tol = rel.tol)$root
}

# T2 alone:
#   T2 -> chi_1^2(tau^2)
# If the rejection rule is T2 > chi_1_alpha^2, then
#   power_T2(tau) = P{chi_1^2(tau^2) > chi_1_alpha^2}.
t2_power <- function(tau, alpha = 0.05) {
  if (alpha <= 0 || alpha >= 1) stop("alpha must be in (0, 1).")
  chi_1_alpha_sq <- qchisq(1 - alpha, df = 1)
  pchisq(
    chi_1_alpha_sq,
    df = 1,
    ncp = tau^2,
    lower.tail = FALSE
  )
}

# ASEL:
#   if mu0 = 0, ASEL coincides with SEL;
#   if mu0 != 0, the second empirical likelihood term diverges after screening.
asel_power <- function(tau, alpha1 = 0.05, alpha2 = 0.08, mu0 = 1,
                       sel_method = c("piecewise", "integral"),
                       rel.tol = 1e-10) {
  cc <- sel_constants(alpha1, alpha2)
  sel_method <- match.arg(sel_method)

  if (isTRUE(all.equal(mu0, 0))) {
    return(sel_power(tau, alpha1, alpha2, method = sel_method, rel.tol = rel.tol))
  }

  value <- 1 - pnorm(cc$z_alpha2 - tau) * (1 - alpha1)
  pmin(pmax(value, 0), 1)
}

# CEL:
#   T_c = ell(mu0) 1(hat_mu > mu0)
# Under mu* = mu0 + tau sigma / sqrt(n),
#   T_c -> V^2 1(V > 0), V ~ N(tau, 1).
# The least-favorable null limit is 0.5 chi^2_0 + 0.5 chi^2_1,
# so a size-alpha rejection rule uses chi^2_{1, 2 alpha}.
cel_power <- function(tau, alpha = 0.05) {
  if (alpha <= 0 || alpha >= 0.5) stop("alpha must be in (0, 0.5).")
  crit <- qchisq(1 - 2 * alpha, df = 1)
  1 - pnorm(sqrt(crit) - tau)
}

# CEL2:
#   T_cel2 = ell(0) 1(hat_mu > 0), fixed mu0 = 0.
# Under mu* = tau sigma / sqrt(n), tau > 0,
#   pi_cel2(tau) = Phi(tau - sqrt(z_{2 alpha}(1))),
# where z_{2 alpha}(1) is the upper 2 alpha quantile of chi^2_1.
cel2_power <- function(tau, alpha = 0.05) {
  if (alpha <= 0 || alpha >= 0.5) stop("alpha must be in (0, 0.5).")
  crit <- qchisq(1 - 2 * alpha, df = 1)
  value <- pnorm(tau - sqrt(crit))
  pmin(pmax(value, 0), 1)
}

#-----------------------------#
# Plot helper
#-----------------------------#
plot_sel_power <- function(mu_grid, power_sel, power_t2 = NULL,
                           power_sel2 = NULL, power_asel = NULL,
                           power_asel2 = NULL,
                           power_cel = NULL,
                           power_cel2 = NULL,
                           alpha1_sel2 = NULL,
                           alpha = 0.05, alpha1 = 0.05, alpha2 = 0.08,
                           mu0 = 0, n = 1000, sigma = 1,
                           show_boundary_labels = TRUE,
                           main = NULL) {
  if (is.null(main)) {
    main <- sprintf(
      "Power  alpha = %.3f, alpha1 = %.3f, alpha2 = %.3f",
      alpha,
      alpha1,
      alpha2
    )
  }

  old_par <- par(no.readonly = TRUE)  # 保存当前图形参数，作图结束后恢复
  on.exit(par(old_par), add = TRUE)   # 函数退出时自动还原 par，避免污染交互会话
  par(
    oma = c(0, 0, 0, 0),              # 外层边距：下、左、上、右；全 0 去掉多余留白
    mar = c(3.2, 3.5, 1.5, 0.5),      # 绘图区边距：下、左、上、右；比默认更紧凑
    mgp = c(2.2, 0.6, 0),             # 轴标题/刻度/次刻度与绘图区的距离
    las = 1,                          # y 轴刻度水平显示，便于阅读
    cex.main = 1.05,                  # 主标题字号
    cex.lab = 1,                      # 轴标题字号
    cex.axis = 0.95                   # 轴刻度字号
  )

  plot(
    mu_grid,
    power_sel,
    type = "l",                       # 折线图
    lwd = 2,                          # 线宽
    col = "firebrick3",               # SEL 曲线颜色
    xlim = c(0, max(mu_grid)),         # 横轴从 0 开始显示
    ylim = c(0, 1),                   # 拒绝概率范围
    xaxs = "i",                       # x 轴紧贴数据范围，去掉默认 4% 内边距
    yaxs = "i",                       # y 轴紧贴 [0, 1]，减少上下留白
    yaxt = "n",                       # 关闭默认 y 轴，下面手动指定 0.1 间隔刻度
    xlab = expression(mu),            # x 轴：真实均值 mu
    ylab = "rejection probability",   # y 轴：拒绝概率
    main = main                       # 标题：名义水平与 alpha1、alpha2
  )
  axis(2, at = seq(0, 1, by = 0.1), las = 1)
  if (!is.null(power_t2)) {
    lines(mu_grid, power_t2, lwd = 2, lty = 2, col = "gray45")
  }
  if (!is.null(power_sel2)) {
    lines(mu_grid, power_sel2, lwd = 2, lty = 1, col = "forestgreen")
  }
  if (!is.null(power_asel)) {
    lines(mu_grid, power_asel, lwd = 2, lty = 3, col = "purple4")
  }
  if (!is.null(power_asel2)) {
    lines(mu_grid, power_asel2, lwd = 2, lty = 5, col = "dodgerblue3")
  }
  if (!is.null(power_cel)) {
    lines(mu_grid, power_cel, lwd = 2, lty = 4, col = "darkorange3")
  }
  if (!is.null(power_cel2)) {
    lines(mu_grid, power_cel2, lwd = 2, lty = 6, col = "#E6B800")
  }
  grid(col = "gray85")
  abline(h = alpha, lty = 2, col = "gray40")
  abline(v = mu0, lty = 3, col = "gray30")
  boundary_idx <- which.min(abs(mu_grid - mu0))
  boundary_pt_cex <- 0.6                 # 临界值点字号
  points(mu0, power_sel[boundary_idx], pch = 19, col = "firebrick3", cex = boundary_pt_cex)
  if (!is.null(power_t2)) {
    points(mu0, power_t2[boundary_idx], pch = 17, col = "gray45", cex = boundary_pt_cex)
  }
  if (!is.null(power_sel2)) {
    points(mu0, power_sel2[boundary_idx], pch = 15, col = "forestgreen", cex = boundary_pt_cex)
  }
  if (!is.null(power_asel)) {
    points(mu0, power_asel[boundary_idx], pch = 18, col = "purple4", cex = boundary_pt_cex)
  }
  if (!is.null(power_asel2)) {
    points(mu0, power_asel2[boundary_idx], pch = 8, col = "dodgerblue3", cex = boundary_pt_cex)
  }
  if (!is.null(power_cel)) {
    points(mu0, power_cel[boundary_idx], pch = 16, col = "darkorange3", cex = boundary_pt_cex)
  }
  if (!is.null(power_cel2)) {
    points(mu0, power_cel2[boundary_idx], pch = 4, col = "#E6B800", cex = boundary_pt_cex)
  }

  plot_items <- data.frame(
    label = "SEL",
    value = power_sel[boundary_idx],
    col = "firebrick3",
    lty = 1,
    pch = 19,
    stringsAsFactors = FALSE
  )
  if (!is.null(power_t2)) {
    plot_items <- rbind(
      plot_items,
      data.frame(label = "EL", value = power_t2[boundary_idx],
                 col = "gray45", lty = 2, pch = 17,
                 stringsAsFactors = FALSE)
    )
  }
  if (!is.null(power_sel2)) {
    plot_items <- rbind(
      plot_items,
      data.frame(label = "SEL2", value = power_sel2[boundary_idx],
                 col = "forestgreen", lty = 1, pch = 15,
                 stringsAsFactors = FALSE)
    )
  }
  if (!is.null(power_asel)) {
    plot_items <- rbind(
      plot_items,
      data.frame(label = "ASEL", value = power_asel[boundary_idx],
                 col = "purple4", lty = 3, pch = 18,
                 stringsAsFactors = FALSE)
    )
  }
  if (!is.null(power_asel2)) {
    plot_items <- rbind(
      plot_items,
      data.frame(label = "ASEL2", value = power_asel2[boundary_idx],
                 col = "dodgerblue3", lty = 5, pch = 8,
                 stringsAsFactors = FALSE)
    )
  }
  if (!is.null(power_cel)) {
    plot_items <- rbind(
      plot_items,
      data.frame(label = "CEL", value = power_cel[boundary_idx],
                 col = "darkorange3", lty = 4, pch = 16,
                 stringsAsFactors = FALSE)
    )
  }
  if (!is.null(power_cel2)) {
    plot_items <- rbind(
      plot_items,
      data.frame(label = "CEL2", value = power_cel2[boundary_idx],
                 col = "#E6B800", lty = 6, pch = 4,
                 stringsAsFactors = FALSE)
    )
  }
  plot_items <- plot_items[order(plot_items$value), ]

  if (show_boundary_labels) {
    label_items <- plot_items[rev(seq_len(nrow(plot_items))), ]
    label_x <- mu0
    label_y <- seq(0.92, by = -0.065, length.out = nrow(label_items))

    text(
      x = label_x,
      y = label_y,
      labels = sprintf("%s %.3f", label_items$label, label_items$value),
      adj = c(0, 0.5),
      col = label_items$col,
      cex = 0.8,                     # 边界数值标注字号
      xpd = NA
    )
  }

  legend_items <- plot_items[rev(seq_len(nrow(plot_items))), ]

  legend(
    "bottomright",
    legend = legend_items$label,
    col = legend_items$col,
    lty = legend_items$lty,
    lwd = rep(2, nrow(legend_items)),
    pch = legend_items$pch,
    inset = c(0.02, 0.06),
    cex = 0.9,
    bty = "o"
  )
}

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  idx <- grep(file_arg, args, fixed = TRUE)
  if (length(idx) == 0) return(getwd())
  path <- sub(file_arg, "", args[idx[1]], fixed = TRUE)
  if (!file.exists(path)) path <- gsub("~\\+~", " ", path)
  if (!file.exists(path)) return(getwd())
  dirname(normalizePath(path))
}

draw_or_save <- function(plot_fun, save_plot = FALSE, file_path = NULL,
                         width = 7, height = 5, res = 300) {
  if (save_plot) {
    if (is.null(file_path)) stop("file_path is required when save_plot = TRUE.")
    ext <- tolower(tools::file_ext(file_path))
    if (ext == "eps") {
      postscript(
        file_path,
        width = width,
        height = height,
        horizontal = FALSE,
        onefile = FALSE,
        paper = "special",
        family = "Helvetica"
      )
    } else if (ext == "png") {
      png(file_path, width = width, height = height, units = "in", res = res)
    } else {
      stop("Unsupported file extension: ", ext, ". Use .eps or .png.")
    }
    on.exit(dev.off(), add = TRUE)
    plot_fun()
    return(invisible(file_path))
  }

  plot_fun()
  invisible(NULL)
}

#-----------------------------#
# Example run
#-----------------------------#
alpha <- 0.05
alpha2 <- 0.08
alpha2_asel2 <- 0.049 # 小于 0.05
alpha1 <- alpha
alpha1_asel2 <- calibrate_sel_alpha1(alpha, alpha2_asel2)
mu0 <- 0.1
sigma <- 1
n <- 1000
tau_grid <- seq(-5, 5, by = 0.02)
mu_grid <- tau_to_mu(tau_grid, mu0 = mu0, sigma = sigma, n = n)
save_plot <- F

power_piecewise <- sel_power(tau_grid, alpha1, alpha2, method = "piecewise")
power_integral <- sel_power(tau_grid, alpha1, alpha2, method = "integral")
alpha1_sel2 <- calibrate_sel_alpha1(alpha, alpha2)
power_sel2 <- sel_power(tau_grid, alpha1_sel2, alpha2, method = "piecewise")
power_t2 <- t2_power(tau_grid, alpha)
power_asel <- asel_power(tau_grid, alpha1, alpha2, mu0 = mu0)
power_asel2 <- asel_power(tau_grid, alpha1_asel2, alpha2_asel2, mu0 = mu0)
power_cel <- cel_power(tau_grid, alpha)
tau_grid_cel2 <- mu_to_tau(mu_grid, mu0 = 0, sigma = sigma, n = n)
power_cel2 <- cel2_power(tau_grid_cel2, alpha)

check <- max(abs(power_piecewise - power_integral))
cat("max |piecewise - integral| =", signif(check, 6), "\n")
cat("mu0                         =", signif(mu0, 8), "\n")
cat("sigma                       =", signif(sigma, 8), "\n")
cat("n                           =", signif(n, 8), "\n")
cat("nominal alpha               =", signif(alpha, 8), "\n")
cat("SEL/ASEL alpha1             =", signif(alpha1, 8), "\n")
cat("boundary rejection pi_s(0)  =", signif(power_piecewise[which.min(abs(tau_grid))], 6), "\n")
cat("boundary rejection pi_a(0)  =", signif(power_asel[which.min(abs(tau_grid))], 6), "\n")
cat("boundary rejection pi_a2(0) =", signif(power_asel2[which.min(abs(tau_grid))], 6), "\n")
cat("ASEL mu0 case               =", if (isTRUE(all.equal(mu0, 0))) "mu0 = 0" else "mu0 != 0", "\n")
cat("ASEL2 alpha2                =", signif(alpha2_asel2, 8), "\n")
cat("ASEL2 alpha1                =", signif(alpha1_asel2, 8), "\n")
cat("SEL2 calibrated alpha1      =", signif(alpha1_sel2, 8), "\n")
cat("SEL2 critical value         =", signif(qchisq(1 - alpha1_sel2, 1), 8), "\n")
cat("boundary rejection SEL2     =", signif(power_sel2[which.min(abs(tau_grid))], 6), "\n")
cat("boundary rejection CEL      =", signif(power_cel[which.min(abs(tau_grid))], 6), "\n")
cat("CEL2 at plotted mu0         =", signif(power_cel2[which.min(abs(tau_grid))], 6), "\n")
cat("CEL2 fixed mu0              = 0\n")
cat("boundary rejection EL       =", signif(power_t2[which.min(abs(tau_grid))], 6), "\n")

out_dir <- script_dir()
out_eps <- file.path(
  out_dir,
  sprintf("SEL-Theory-Reject-Ratio-alpha1-%.3f-alpha2-%.3f.eps",
          alpha1, alpha2)
)

if (!save_plot && !interactive()) {
  cat("save_plot = FALSE; plots are drawn on the active graphics device if available.\n")
  cat("Set save_plot <- TRUE to save EPS files when running by Rscript.\n")
}

draw_or_save(
  function() {
    plot_sel_power(
      mu_grid = mu_grid,
      power_sel = power_piecewise,
      power_t2 = power_t2,
      power_sel2 = power_sel2,
      power_asel = power_asel,
      power_asel2 = power_asel2,
      power_cel = power_cel,
      power_cel2 = power_cel2,
      alpha1_sel2 = alpha1_sel2,
      alpha = alpha,
      alpha1 = alpha1,
      alpha2 = alpha2,
      mu0 = mu0,
      n = n,
      sigma = sigma
    )
  },
  save_plot = save_plot,
  file_path = out_eps
)

if (save_plot) {
  cat("potential plot saved:", out_eps, "\n")
}
