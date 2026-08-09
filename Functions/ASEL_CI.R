#===============================================
#                    Compas  ASEL              #
#===============================================
# epsilonG 的 ASEL 置信区间：g_i = Y_i - theta_P - epsilonG
# 分拆经验似然 T_a(mu)；临界值 z_{alpha}(1)（SEL.tex §3.2 / §3.3）
# gate 激活时使用 tilde T_2 = ell_2(0)（与 SEL 的 ell_2(2 mu) 不同）
#
# 三种检验（alternative 别名 test1 / test2 / test3）：
#   test-1 / greater   : H0 0<=mu<=mu0 vs H1 mu>mu0        -> CI [lb, +infty)
#   test-2 / less      : H0 -mu0<=mu<=0 vs H1 mu<-mu0      -> CI (-infty, ub]
#   test-3 / two.sided : H0 -mu1<=mu<=mu2 vs H1 outside   -> CI [lb, ub]
#===============================================

source("Functions/GlambdaChen.R")

# ---------------------------------------------------------------------------
# alternative 规范化
# ---------------------------------------------------------------------------

.normalize_asel_alternative <- function(alternative) {
  alt <- tolower(trimws(as.character(alternative)))
  switch(alt,
    "greater" = "greater",
    ">" = "greater",
    "right" = "greater",
    "test1" = "greater",
    "大于" = "greater",
    "less" = "less",
    "<" = "less",
    "left" = "less",
    "test2" = "less",
    "小于" = "less",
    "two.sided" = "two.sided",
    "twosided" = "two.sided",
    "two-sided" = "two.sided",
    "interval" = "two.sided",
    "test3" = "two.sided",
    "or" = "two.sided",
    "双边" = "two.sided",
    stop(
      "alternative 须为 greater / less / two.sided 之一，当前为: ",
      alternative
    )
  )
}

# ---------------------------------------------------------------------------
# 基础：单样本经验似然
# ---------------------------------------------------------------------------

el_mean_ratio <- function(x, mu) {
  x <- as.numeric(x)
  if (length(x) < 2L) {
    stop("样本量至少为 2。")
  }

  z <- matrix(x - mu, ncol = 1L)
  lam <- c(lambdaChen(z))
  aa <- 1 + t(lam) %*% t(z)
  if (any(aa <= 0)) {
    return(0)
  }
  npi <- 1 / (1 + lam * z)
  if (any(npi <= 0)) {
    return(0)
  }
  prod(npi)
}

el_mean_stat <- function(x, mu) {
  ratio <- el_mean_ratio(x, mu)
  if (ratio <= 0) {
    return(Inf)
  }
  -2 * log(ratio)
}

# ---------------------------------------------------------------------------
# 分拆样本 d1, d2、gate 常数 c_n 与 tilde T_2
# ---------------------------------------------------------------------------

.split_deltas <- function(x) {
  x <- as.numeric(x)
  n <- length(x)
  if (n %% 2L != 0L) {
    x <- x[seq_len(n - 1L)]
    n <- length(x)
  }
  n1 <- n %/% 2L
  x1 <- x[seq_len(n1)]
  x2 <- x[(n1 + 1L):n]
  m <- min(length(x1), length(x2))
  list(
    d1 = x1[seq_len(m)] - x2[seq_len(m)],
    d2 = x1[seq_len(m)] + x2[seq_len(m)],
    m = m
  )
}

.asel_cn <- function(d2, alpha2 = 0.08) {
  sigma2_hat <- sd(d2)
  if (is.na(sigma2_hat) || sigma2_hat <= 0) {
    sigma2_hat <- max(abs(d2), 1e-8)
  }
  qnorm(1 - alpha2) * sigma2_hat
}

.asel_t2 <- function(d2) {
  el_mean_stat(d2, 0)
}

# ---------------------------------------------------------------------------
# test-1 / greater：T_a = T_1 + tilde T_2 * 1(T_3 > c_n)
# ---------------------------------------------------------------------------

el_asel_stat_greater <- function(d1, d2, mu, cn, T1 = NULL, T2 = NULL) {
  mu <- as.numeric(mu)
  m <- length(d1)
  if (is.null(T1)) {
    T1 <- el_mean_stat(d1, 0)
  }
  if (is.null(T2)) {
    T2 <- .asel_t2(d2)
  }
  T3 <- sqrt(m) * (mean(d2) - 2 * mu)
  if (T3 > cn) {
    T1 + T2
  } else {
    T1
  }
}

# ---------------------------------------------------------------------------
# test-2 / less：T_a = T_1 + tilde T_2 * 1(T_3 < -c_n)
# ---------------------------------------------------------------------------

el_asel_stat_less <- function(d1, d2, mu, cn, T1 = NULL, T2 = NULL) {
  mu <- as.numeric(mu)
  m <- length(d1)
  if (is.null(T1)) {
    T1 <- el_mean_stat(d1, 0)
  }
  if (is.null(T2)) {
    T2 <- .asel_t2(d2)
  }
  T3 <- sqrt(m) * (mean(d2) - 2 * mu)
  if (T3 < -cn) {
    T1 + T2
  } else {
    T1
  }
}

# ---------------------------------------------------------------------------
# test-3 / two.sided：T_a = T_1 + tilde T_2 * 1(|T_3| > c_n)
# ---------------------------------------------------------------------------

el_asel_stat_interval <- function(d1, d2, mu, cn, eps_lo, eps_hi,
                                 T1 = NULL, T2 = NULL) {
  mu <- as.numeric(mu)
  m <- length(d1)
  if (is.null(T1)) {
    T1 <- el_mean_stat(d1, 0)
  }
  if (is.null(T2)) {
    T2 <- .asel_t2(d2)
  }
  T3 <- sqrt(m) * (mean(d2) - 2 * mu)
  if (T3 < -cn || T3 > cn) {
    T1 + T2
  } else {
    T1
  }
}

# ---------------------------------------------------------------------------
# test-1 / greater：置信区间下界
# ---------------------------------------------------------------------------

.asel_boundary_greater <- function(d1, d2, m, cn, T1, T2, crit) {
  Ta <- T1 + T2
  if (T1 > crit) {
    return(NA_real_)
  }
  if (Ta <= crit) {
    return(-Inf)
  }
  (mean(d2) - cn / sqrt(m)) / 2
}

# ---------------------------------------------------------------------------
# test-2 / less：置信区间上界
# ---------------------------------------------------------------------------

.asel_boundary_less <- function(d1, d2, m, cn, T1, T2, crit) {
  Ta <- T1 + T2
  if (T1 > crit) {
    return(NA_real_)
  }
  if (Ta <= crit) {
    return(Inf)
  }
  (mean(d2) + cn / sqrt(m)) / 2
}

# ---------------------------------------------------------------------------
# test-3 / two.sided：置信区间 [lb, ub]
# ---------------------------------------------------------------------------

.asel_interval_mu_range <- function(mu_hat, eps_lo, eps_hi, d1, d2,
                                    mu_lo = -1, mu_hi = 1, n_grid = 800L) {
  list(
    lo = mu_lo,
    hi = mu_hi,
    mu_grid = seq(mu_lo, mu_hi, length.out = n_grid)
  )
}

.asel_boundary_interval <- function(d1, d2, m, cn, crit, mu_hat, eps_lo, eps_hi,
                                    T1, T2) {
  # T1 > crit → NA；T2 < crit → (-Inf, Inf)；否则在 [-1, 1] 网格上求 lb/ub
  if (T1 > crit) {
    return(list(
      lb = NA_real_, ub = NA_real_,
      anchor_lo = eps_lo, anchor_hi = eps_hi
    ))
  }
  if (T2 < crit) {
    return(list(
      lb = -Inf, ub = Inf,
      anchor_lo = eps_lo, anchor_hi = eps_hi
    ))
  }

  rg <- .asel_interval_mu_range(mu_hat, eps_lo, eps_hi, d1, d2)
  mu_grid <- rg$mu_grid
  Ts <- vapply(
    mu_grid, el_asel_stat_interval, numeric(1L),
    d1 = d1, d2 = d2, cn = cn,
    eps_lo = eps_lo, eps_hi = eps_hi,
    T1 = T1, T2 = T2
  )
  ok <- is.finite(Ts) & Ts <= crit
  if (!any(ok)) {
    return(list(
      lb = NA_real_, ub = NA_real_,
      anchor_lo = eps_lo, anchor_hi = eps_hi
    ))
  }

  lb <- min(mu_grid[ok])
  ub <- max(mu_grid[ok])
  if (ok[1L]) {
    lb <- -Inf
  }
  if (ok[length(ok)]) {
    ub <- Inf
  }

  list(lb = lb, ub = ub, anchor_lo = eps_lo, anchor_hi = eps_hi)
}

# ---------------------------------------------------------------------------
# 三种检验：统一求置信区间边界
# ---------------------------------------------------------------------------

.asel_ci_bounds <- function(d1, d2, m, cn, T1, T2, crit, mu_hat,
                            alternative, mu0, eps_lo, eps_hi) {
  alternative <- .normalize_asel_alternative(alternative)

  if (alternative == "greater") {
    lb <- .asel_boundary_greater(d1, d2, m, cn, T1, T2, crit)
    if (is.na(lb)) {
      return(list(
        lb = NA_real_, ub = NA_real_,
        anchor_lo = mu0, anchor_hi = NA_real_
      ))
    }
    return(list(lb = lb, ub = Inf, anchor_lo = mu0, anchor_hi = NA_real_))
  }

  if (alternative == "less") {
    ub <- .asel_boundary_less(d1, d2, m, cn, T1, T2, crit)
    if (is.na(ub)) {
      return(list(
        lb = NA_real_, ub = NA_real_,
        anchor_lo = NA_real_, anchor_hi = mu0
      ))
    }
    return(list(lb = -Inf, ub = ub, anchor_lo = NA_real_, anchor_hi = mu0))
  }

  if (eps_lo >= eps_hi) {
    stop("two.sided 要求 epsilon_lo < epsilon_hi（即 -mu_1 < mu_2）。")
  }
  .asel_boundary_interval(
    d1, d2, m, cn, crit, mu_hat, eps_lo, eps_hi, T1, T2
  )
}

.asel_al_length <- function(lb, ub) {
  if (is.na(lb) || is.na(ub)) {
    return(NA_real_)
  }
  if (!is.finite(lb) || !is.finite(ub)) {
    return(Inf)
  }
  ub - lb
}

# ---------------------------------------------------------------------------
# 作图
# ---------------------------------------------------------------------------

.finite_plot_ylim <- function(values, crit, default_max = 1) {
  vals <- c(as.numeric(values), crit)
  vals <- vals[is.finite(vals)]
  ymax <- if (length(vals)) max(vals) else max(crit, default_max, na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0) {
    ymax <- max(crit, default_max, na.rm = TRUE)
  }
  ymax * 1.05
}

.asel_plot_annotate_bound <- function(bound, y_bound, label, pos = 2) {
  if (!is.finite(bound)) {
    return(invisible(NULL))
  }
  abline(v = bound, col = "blue", lty = 1, lwd = 1.5)
  points(bound, y_bound, pch = 19, col = "blue", cex = 0.9)
  text(bound, y_bound, labels = label, pos = pos, col = "blue", cex = 0.9)
}

.asel_plot <- function(d1, d2, m, cn, crit, bounds, alternative,
                       mu_hat, T1, T2, sex_name, age_cat_name,
                       alpha, alpha2, plot_label = NULL) {
  alternative <- .normalize_asel_alternative(alternative)
  default_main <- paste0(
    "ASEL(", alternative, ") alpha=", alpha, " alpha2=", alpha2, " ",
    sex_name, "  ", age_cat_name
  )
  plot_main <- if (!is.null(plot_label) && nzchar(trimws(plot_label))) {
    trimws(plot_label)
  } else {
    default_main
  }
  Ta <- T1 + T2
  pad <- max(0.05, 0.15 * diff(range(c(d1, d2))), na.rm = TRUE)

  if (alternative == "greater") {
    lb <- bounds$lb
    mu_switch <- (mean(d2) - cn / sqrt(m)) / 2
    if (is.finite(lb)) {
      pad <- max(pad, abs(mu_hat - lb), na.rm = TRUE)
    }
    x_left <- if (is.finite(lb)) lb - pad * 0.5 else mu_hat - pad
    x_right <- mu_hat + pad
    y_max <- .finite_plot_ylim(c(Ta, T1), crit)

    par(mfrow = c(1, 1))
    plot(
      c(x_left, mu_switch, x_right), c(Ta, T1, T1),
      type = "n", xlim = c(x_left, x_right), ylim = c(0, y_max),
      main = plot_main,
      xlab = "epsilonG", ylab = expression(T[a])
    )
    segments(x_left, Ta, mu_switch, Ta, lwd = 1.5)
    segments(mu_switch, T1, x_right, T1, lwd = 1.5)
    segments(mu_switch, T1, mu_switch, Ta, lty = 3, col = "gray50")
  } else if (alternative == "less") {
    ub <- bounds$ub
    mu_switch <- (mean(d2) + cn / sqrt(m)) / 2
    if (is.finite(ub)) {
      pad <- max(pad, abs(mu_hat - ub), na.rm = TRUE)
    }
    x_left <- mu_hat - pad
    x_right <- if (is.finite(ub)) ub + pad * 0.5 else mu_hat + pad
    x_lo <- min(x_left, mu_switch, if (is.finite(ub)) ub else mu_hat) - pad * 0.25
    x_hi <- max(x_right, mu_switch, mu_hat) + pad * 0.25
    epsilonGs <- seq(x_lo, x_hi, length.out = 200L)
    Ta_vals <- vapply(
      epsilonGs, el_asel_stat_less, numeric(1L),
      d1 = d1, d2 = d2, cn = cn, T1 = T1, T2 = T2
    )
    y_max <- .finite_plot_ylim(Ta_vals, crit)

    par(mfrow = c(1, 1))
    plot(
      epsilonGs, Ta_vals,
      type = "l", ylim = c(0, y_max),
      main = plot_main,
      xlab = "epsilonG", ylab = expression(T[a])
    )
    if (is.finite(mu_switch) && mu_switch >= x_lo && mu_switch <= x_hi) {
      abline(v = mu_switch, col = "gray50", lty = 3)
    }
  } else {
    rg <- .asel_interval_mu_range(
      mu_hat, bounds$anchor_lo, bounds$anchor_hi, d1, d2,
      n_grid = 400L
    )
    epsilonGs <- rg$mu_grid
    Ta_vals <- vapply(
      epsilonGs, el_asel_stat_interval, numeric(1L),
      d1 = d1, d2 = d2, cn = cn,
      eps_lo = bounds$anchor_lo, eps_hi = bounds$anchor_hi,
      T1 = T1, T2 = T2
    )
    y_max <- .finite_plot_ylim(Ta_vals, crit)

    par(mfrow = c(1, 1))
    plot(
      epsilonGs, Ta_vals,
      type = "l", ylim = c(0, y_max),
      main = plot_main,
      xlab = "epsilonG", ylab = expression(T[a])
    )
  }

  abline(h = crit, col = "red", lty = 2)
  x_usr <- par("usr")
  text(
    x_usr[1] + 0.06 * diff(x_usr[1:2]), crit,
    labels = paste0("crit=", round(crit, 4)),
    pos = 3, col = "red", cex = 0.9
  )

  if (alternative == "greater" && is.finite(bounds$lb)) {
    .asel_plot_annotate_bound(
      bounds$lb, T1,
      paste0("lb=", round(bounds$lb, 4)),
      pos = 2
    )
  }
  if (alternative == "less" && is.finite(bounds$ub)) {
    y_ub <- el_asel_stat_less(d1, d2, bounds$ub, cn, T1, T2)
    .asel_plot_annotate_bound(
      bounds$ub, y_ub,
      paste0("ub=", round(bounds$ub, 4)),
      pos = 4
    )
  }
  if (alternative == "two.sided") {
    if (is.finite(bounds$lb)) {
      y_lb <- el_asel_stat_interval(
        d1, d2, bounds$lb, cn, bounds$anchor_lo, bounds$anchor_hi, T1, T2
      )
      .asel_plot_annotate_bound(
        bounds$lb, y_lb,
        paste0("lb=", round(bounds$lb, 4)),
        pos = 2
      )
    }
    if (is.finite(bounds$ub)) {
      y_ub <- el_asel_stat_interval(
        d1, d2, bounds$ub, cn, bounds$anchor_lo, bounds$anchor_hi, T1, T2
      )
      .asel_plot_annotate_bound(
        bounds$ub, y_ub,
        paste0("ub=", round(bounds$ub, 4)),
        pos = 4
      )
    }
  }

  abline(v = mu_hat, col = "darkgreen", lty = 2)
  text(
    mu_hat, y_max * 0.85,
    labels = paste0("ep_h=", round(mu_hat, 4)),
    pos = 4, col = "darkgreen", cex = 0.9
  )
}

# ---------------------------------------------------------------------------
# 主入口
# ---------------------------------------------------------------------------

ASEL_CI <- function(group_specific_data,
                    target = theta_P,
                    alpha = 0.95,
                    alpha2 = 0.08,
                    alternative = "greater",
                    mu0 = 0,
                    epsilon_lo = NULL,
                    epsilon_hi = NULL,
                    plot = TRUE,
                    label = NULL,
                    plot_label = NULL,
                    verbose = TRUE) {
  sex_name <- unique(group_specific_data$sex)
  age_cat_name <- unique(group_specific_data$age_cat)

  if (length(sex_name) != 1) {
    sex_name <- ""
  }
  if (length(age_cat_name) != 1) {
    age_cat_name <- ""
  }

  y <- as.numeric(group_specific_data$Y)
  n <- length(y)
  if (n < 4L) {
    stop("ASEL 至少需要 4 个观测（n = 2m, m >= 2）。")
  }
  if (alpha <= 0 || alpha >= 1) {
    stop("alpha 必须介于 (0, 1) 之间。")
  }

  alt <- .normalize_asel_alternative(alternative)
  x <- y - target
  split <- .split_deltas(x)
  d1 <- split$d1
  d2 <- split$d2
  m <- split$m
  mu_hat <- mean(x)
  cn <- .asel_cn(d2, alpha2)
  crit <- qchisq(alpha, 1L)
  T1 <- el_mean_stat(d1, 0)
  T2 <- .asel_t2(d2)

  if (alt == "two.sided") {
    if (is.null(epsilon_lo) || is.null(epsilon_hi)) {
      stop("two.sided 须同时提供 epsilon_lo 与 epsilon_hi（容忍区间端点）。")
    }
    eps_lo <- as.numeric(epsilon_lo)
    eps_hi <- as.numeric(epsilon_hi)
  } else {
    eps_lo <- NA_real_
    eps_hi <- NA_real_
    mu0 <- as.numeric(mu0)
  }

  bounds <- .asel_ci_bounds(
    d1, d2, m, cn, T1, T2, crit, mu_hat,
    alt, mu0, eps_lo, eps_hi
  )
  lb <- bounds$lb
  ub <- bounds$ub
  ASELAL <- .asel_al_length(lb, ub)

  if (plot) {
    .asel_plot(
      d1, d2, m, cn, crit, bounds, alt,
      mu_hat, T1, T2, sex_name, age_cat_name, alpha, alpha2,
      plot_label = plot_label
    )
  }

  if (verbose) {
    if (!is.null(label) && nzchar(label)) {
      cat('"', label, '"\n', sep = "")
    }
    cat("alternative ", "n ", " ep_h ", " lb ", " ub ", " ASELAL", "\n")
    cat(alt, n, mu_hat, lb, ub, ASELAL, "\n")
    if (!is.null(label) && nzchar(label)) {
      cat("\n")
    }
    if (is.null(label) || !nzchar(label)) {
      cat("样本epsilonG", mu_hat, "\n")
    }
  }

  invisible(list(
    n = n,
    epsilonG = mu_hat,
    lb = lb,
    ub = ub,
    ASELAL = ASELAL,
    crit = crit,
    cn = cn,
    T1 = T1,
    T2 = T2,
    alpha = alpha,
    alpha2 = alpha2,
    alternative = alt,
    mu0 = if (alt != "two.sided") mu0 else NA_real_,
    epsilon_lo = if (alt == "two.sided") eps_lo else NA_real_,
    epsilon_hi = if (alt == "two.sided") eps_hi else NA_real_
  ))
}
