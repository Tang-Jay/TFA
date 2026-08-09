#===============================================
#                    Compas  CEL               #
#===============================================
# epsilonG 的 CEL 置信区间：g_i = Y_i - theta_P - epsilonG
# 临界值 z_{2a}(1)（混合 chi^2 上 2a 分位，ELFA Theorem 5/7/8）
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

.normalize_cel_alternative <- function(alternative) {
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
# test-1 / greater：H0 mu <= mu0 vs H1 mu > mu0（ELFA Sec 4.2）
# ---------------------------------------------------------------------------

el_cel_stat_greater <- function(x, mu0) {
  mu0 <- as.numeric(mu0)
  xbar <- mean(x)
  if (xbar <= mu0) {
    return(0)
  }
  el_mean_stat(x, mu0)
}

# ---------------------------------------------------------------------------
# test-2 / less：H0 mu >= mu0 vs H1 mu < mu0
# ---------------------------------------------------------------------------

el_cel_stat_less <- function(x, mu0) {
  mu0 <- as.numeric(mu0)
  xbar <- mean(x)
  if (xbar >= mu0) {
    return(0)
  }
  el_mean_stat(x, mu0)
}

# ---------------------------------------------------------------------------
# test-3 / two.sided：H0 eps_lo <= mu <= eps_hi vs H1 outside
# ---------------------------------------------------------------------------

el_cel_stat_interval <- function(x, eps_lo, eps_hi) {
  eps_lo <- as.numeric(eps_lo)
  eps_hi <- as.numeric(eps_hi)
  xbar <- mean(x)
  if (xbar >= eps_lo && xbar <= eps_hi) {
    return(0)
  }
  if (xbar < eps_lo) {
    return(el_mean_stat(x, eps_lo))
  }
  el_mean_stat(x, eps_hi)
}

# ---------------------------------------------------------------------------
# 数值求根：EL 剖面 {mu : ell(mu) >= cut}
# ---------------------------------------------------------------------------

.expand_bracket <- function(f, left, right, max_expand = 20L) {
  fl <- f(left)
  fr <- f(right)
  k <- 0L
  while (fl * fr > 0 && k < max_expand) {
    width <- right - left
    if (width <= 0) {
      width <- max(abs(left), abs(right), 1)
    }
    left <- left - width
    right <- right + width
    fl <- f(left)
    fr <- f(right)
    k <- k + 1L
  }
  if (fl * fr > 0) {
    return(NULL)
  }
  c(left, right)
}

.boundary_root <- function(x, cut, start, direction = c("left", "right")) {
  direction <- match.arg(direction)
  f <- function(mu) el_mean_ratio(x, mu) - cut

  if (direction == "left") {
    bracket <- .expand_bracket(f, start - diff(range(x)) / 2, start)
  } else {
    bracket <- .expand_bracket(f, start, start + diff(range(x)) / 2)
  }
  if (is.null(bracket)) {
    return(NA_real_)
  }

  uniroot(f, bracket)$root
}

# ---------------------------------------------------------------------------
# test-1 / greater：置信区间下界
# ---------------------------------------------------------------------------

.cel_boundary_greater <- function(x, mu_hat, cut) {
  lb <- .boundary_root(x, cut, mu_hat, "left")
  if (is.na(lb)) {
    lb <- mu_hat
  }
  lb
}

# ---------------------------------------------------------------------------
# test-2 / less：置信区间上界
# ---------------------------------------------------------------------------

.cel_boundary_less <- function(x, mu_hat, cut) {
  ub <- .boundary_root(x, cut, mu_hat, "right")
  if (is.na(ub)) {
    ub <- mu_hat
  }
  ub
}

# ---------------------------------------------------------------------------
# test-3 / two.sided：置信区间 [lb, ub]
# ---------------------------------------------------------------------------

.cel_boundary_interval <- function(x, mu_hat, cut) {
  lb <- .boundary_root(x, cut, mu_hat, "left")
  ub <- .boundary_root(x, cut, mu_hat, "right")
  if (is.na(lb)) {
    lb <- mu_hat
  }
  if (is.na(ub)) {
    ub <- mu_hat
  }
  list(lb = lb, ub = ub)
}

# ---------------------------------------------------------------------------
# 三种检验：统一求置信区间边界
# ---------------------------------------------------------------------------

.cel_ci_bounds <- function(x, mu_hat, crit, cut, alternative, mu0, eps_lo, eps_hi) {
  alternative <- .normalize_cel_alternative(alternative)

  if (alternative == "greater") {
    lb <- .cel_boundary_greater(x, mu_hat, cut)
    return(list(lb = lb, ub = Inf, anchor_lo = mu0, anchor_hi = NA_real_))
  }

  if (alternative == "less") {
    ub <- .cel_boundary_less(x, mu_hat, cut)
    return(list(lb = -Inf, ub = ub, anchor_lo = NA_real_, anchor_hi = mu0))
  }

  if (eps_lo >= eps_hi) {
    stop("two.sided 要求 epsilon_lo < epsilon_hi（即 -mu_1 < mu_2）。")
  }
  bounds <- .cel_boundary_interval(x, mu_hat, cut)
  list(
    lb = bounds$lb,
    ub = bounds$ub,
    anchor_lo = eps_lo,
    anchor_hi = eps_hi
  )
}

.cel_al_length <- function(lb, ub) {
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

.cel_plot_title <- function(alternative, alpha, sex_name, age_cat_name, plot_label = NULL) {
  if (!is.null(plot_label) && nzchar(trimws(plot_label))) {
    return(trimws(plot_label))
  }
  paste0(
    "CEL(", alternative, ") alpha=", alpha, " ",
    sex_name, "  ", age_cat_name
  )
}

.cel_plot_annotate_bound <- function(bound, crit, y_max, label, pos = 4) {
  if (!is.finite(bound)) {
    return(invisible(NULL))
  }
  abline(v = bound, col = "blue", lty = 1, lwd = 1.5)
  points(bound, crit, pch = 19, col = "blue", cex = 0.9)
  text(
    bound, crit + 0.05 * y_max,
    labels = label, pos = pos, col = "blue", cex = 0.9
  )
}

.cel_plot <- function(x, mu_hat, crit, cut, bounds, alternative, sex_name, age_cat_name, alpha,
                      plot_label = NULL) {
  alternative <- .normalize_cel_alternative(alternative)
  pad <- max(0.05, 0.15 * diff(range(x)), na.rm = TRUE)

  if (alternative == "greater") {
    lb <- bounds$lb
    pad <- max(pad, abs(mu_hat - lb), na.rm = TRUE)
    epsilonGs <- seq(lb - pad * 0.5, mu_hat + pad, length.out = 200L)
    Tc <- vapply(epsilonGs, el_cel_stat_greater, numeric(1L), x = x)
    ylab <- expression(T[c])
  } else if (alternative == "less") {
    ub <- bounds$ub
    pad <- max(pad, abs(mu_hat - ub), na.rm = TRUE)
    epsilonGs <- seq(mu_hat - pad, ub + pad * 0.5, length.out = 200L)
    Tc <- vapply(epsilonGs, el_cel_stat_less, numeric(1L), x = x)
    ylab <- expression(T[c])
  } else {
    lb <- bounds$lb
    ub <- bounds$ub
    pad <- max(pad, abs(mu_hat - lb), abs(mu_hat - ub), na.rm = TRUE)
    lo_plot <- if (is.finite(lb)) lb - pad * 0.5 else mu_hat - 2 * pad
    hi_plot <- if (is.finite(ub)) ub + pad * 0.5 else mu_hat + 2 * pad
    epsilonGs <- seq(lo_plot, hi_plot, length.out = 200L)
    Tc <- vapply(epsilonGs, el_mean_stat, numeric(1L), x = x)
    ylab <- expression(ell(mu))
  }

  par(mfrow = c(1, 1))
  y_max <- .finite_plot_ylim(Tc, crit)
  plot(
    epsilonGs, Tc,
    type = "l", ylim = c(0, y_max),
    main = .cel_plot_title(alternative, alpha, sex_name, age_cat_name, plot_label),
    xlab = "epsilonG", ylab = ylab
  )
  abline(h = crit, col = "red", lty = 2)
  text(
    par("usr")[1] + 0.03, crit,
    labels = paste0("crit=", round(crit, 4)),
    pos = 3, col = "red", cex = 0.9
  )

  if (alternative == "greater") {
    .cel_plot_annotate_bound(
      bounds$lb, crit, y_max,
      paste0("lb=", round(bounds$lb, 4)),
      pos = 4
    )
  } else if (alternative == "less") {
    .cel_plot_annotate_bound(
      bounds$ub, crit, y_max,
      paste0("ub=", round(bounds$ub, 4)),
      pos = 4
    )
  } else {
    if (is.finite(bounds$lb)) {
      .cel_plot_annotate_bound(
        bounds$lb, crit, y_max,
        paste0("lb=", round(bounds$lb, 4)),
        pos = 4
      )
    }
    if (is.finite(bounds$ub)) {
      .cel_plot_annotate_bound(
        bounds$ub, crit, y_max,
        paste0("ub=", round(bounds$ub, 4)),
        pos = 2
      )
    }
    abline(v = bounds$anchor_lo, col = "purple", lty = 3)
    abline(v = bounds$anchor_hi, col = "purple", lty = 3)
    text(
      bounds$anchor_lo, y_max * 0.65,
      labels = paste0("eps_lo=", round(bounds$anchor_lo, 4)),
      pos = 4, col = "purple", cex = 0.85
    )
    text(
      bounds$anchor_hi, y_max * 0.55,
      labels = paste0("eps_hi=", round(bounds$anchor_hi, 4)),
      pos = 2, col = "purple", cex = 0.85
    )
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

CEL_CI <- function(group_specific_data,
                   target = theta_P,
                   alpha = 0.95,
                   alternative = "greater",
                   mu0 = 0,
                   epsilon_lo = NULL,
                   epsilon_hi = NULL,
                   plot = TRUE,
                   plot_label = NULL) {
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
  if (n < 2L) {
    stop("分组样本量至少为 2。")
  }
  if (alpha <= 0 || alpha >= 1) {
    stop("alpha 必须介于 (0, 1) 之间。")
  }

  alt <- .normalize_cel_alternative(alternative)
  x <- y - target
  mu_hat <- mean(x)
  sig_level <- 1 - alpha
  crit <- qchisq(1 - 2 * sig_level, 1L)
  cut <- exp(-crit / 2)

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

  bounds <- .cel_ci_bounds(
    x, mu_hat, crit, cut, alt, mu0, eps_lo, eps_hi
  )
  lb <- bounds$lb
  ub <- bounds$ub
  CELAL <- .cel_al_length(lb, ub)

  if (plot) {
    .cel_plot(
      x, mu_hat, crit, cut, bounds, alt,
      sex_name, age_cat_name, alpha, plot_label = plot_label
    )
  }

  cat("alternative ", "n ", " ep_h ", " lb ", " ub ", " CELAL", "\n")
  cat(alt, n, mu_hat, lb, ub, CELAL, "\n")
  cat("样本epsilonG", mu_hat, "\n")

  invisible(list(
    n = n,
    epsilonG = mu_hat,
    lb = lb,
    ub = ub,
    CELAL = CELAL,
    crit = crit,
    cut = cut,
    alpha = alpha,
    alternative = alt,
    mu0 = if (alt != "two.sided") mu0 else NA_real_,
    epsilon_lo = if (alt == "two.sided") eps_lo else NA_real_,
    epsilon_hi = if (alt == "two.sided") eps_hi else NA_real_
  ))
}
