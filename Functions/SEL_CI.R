#===============================================
#                    Compas  SEL               #
#===============================================
# epsilonG 的 SEL 置信区间：g_i = Y_i - theta_P - epsilonG
# 分拆经验似然 T_s(mu)；临界值 z_{alpha}(1)（SEL.tex Theorem 4 / §3.3）
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

.normalize_sel_alternative <- function(alternative) {
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
# 分拆样本 d1, d2 与 gate 常数 c_n
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

.sel_cn <- function(d2, alpha2 = 0.08) {
  sigma2_hat <- sd(d2)
  if (is.na(sigma2_hat) || sigma2_hat <= 0) {
    sigma2_hat <- max(abs(d2), 1e-8)
  }
  qnorm(1 - alpha2) * sigma2_hat
}

# ---------------------------------------------------------------------------
# test-1 / greater：T_s = T_1 + T_2 * 1(T_3 > c_n)，T_3 = sqrt(m)(bar d2 - 2 mu)
# ---------------------------------------------------------------------------

el_sel_stat_greater <- function(d1, d2, mu, cn) {
  mu <- as.numeric(mu)
  m <- length(d1)
  T1 <- el_mean_stat(d1, 0)
  T2 <- el_mean_stat(d2, 2 * mu)
  T3 <- sqrt(m) * (mean(d2) - 2 * mu)
  if (T3 > cn) {
    T1 + T2
  } else {
    T1
  }
}

# ---------------------------------------------------------------------------
# test-2 / less：T_s = T_1 + T_2 * 1(T_3 < -c_n)
# ---------------------------------------------------------------------------

el_sel_stat_less <- function(d1, d2, mu, cn, mu0 = NULL) {
  if (!is.null(mu0)) {
    mu <- mu0
  }
  mu <- as.numeric(mu)
  m <- length(d1)
  T1 <- el_mean_stat(d1, 0)
  T2 <- el_mean_stat(d2, 2 * mu)
  T3 <- sqrt(m) * (mean(d2) - 2 * mu)
  if (T3 < -cn) {
    T1 + T2
  } else {
    T1
  }
}

# ---------------------------------------------------------------------------
# test-3 / two.sided (or)：T_s = T_1 + T_2 * 1(|T_3| > c_n)
# ---------------------------------------------------------------------------

el_sel_stat_interval <- function(d1, d2, mu, cn) {
  mu <- as.numeric(mu)
  m <- length(d1)
  T1 <- el_mean_stat(d1, 0)
  T2 <- el_mean_stat(d2, 2 * mu)
  T3 <- sqrt(m) * (mean(d2) - 2 * mu)
  if (T3 < -cn || T3 > cn) {
    T1 + T2
  } else {
    T1
  }
}

# ---------------------------------------------------------------------------
# 数值求根
# ---------------------------------------------------------------------------

.expand_bracket <- function(f, left, right, max_expand = 20L) {
  fl <- f(left)
  fr <- f(right)
  k <- 0L
  while (is.finite(fl) && is.finite(fr) && fl * fr > 0 && k < max_expand) {
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
  if (!is.finite(fl) || !is.finite(fr) || fl * fr > 0) {
    return(NULL)
  }
  c(left, right)
}

.sel_find_bracket <- function(f, lo, hi, n = 200L) {
  if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
    return(NULL)
  }
  xs <- seq(lo, hi, length.out = n)
  ys <- vapply(xs, f, numeric(1))
  ok <- is.finite(ys)
  if (sum(ok) < 2L) {
    return(NULL)
  }
  xs <- xs[ok]
  ys <- ys[ok]
  for (i in seq_len(length(ys) - 1L)) {
    if (ys[i] * ys[i + 1L] <= 0) {
      return(c(xs[i], xs[i + 1L]))
    }
  }
  NULL
}

.sel_f_active <- function(d2, T1, crit) {
  function(mu0) {
    mu0 <- as.numeric(mu0)
    ratio <- el_mean_ratio(d2, 2 * mu0)
    if (ratio <= 0) {
      return(NA_real_)
    }
    T1 - crit - 2 * log(ratio)
  }
}

.sel_root_on_bracket <- function(f, lo, hi) {
  bracket <- .sel_find_bracket(f, lo, hi)
  if (is.null(bracket)) {
    bracket <- .expand_bracket(f, lo, hi)
  }
  if (is.null(bracket)) {
    return(NA_real_)
  }
  tryCatch(
    uniroot(f, bracket, extendInt = "no")$root,
    error = function(e) NA_real_
  )
}

# ---------------------------------------------------------------------------
# test-1 / greater：置信区间下界
# ---------------------------------------------------------------------------

.sel_boundary_greater <- function(d1, d2, m, cn, T1, crit, mu_hat) {
  if (T1 > crit) {
    return(NA_real_)
  }

  mu_switch <- (mean(d2) - cn / sqrt(m)) / 2
  f <- .sel_f_active(d2, T1, crit)
  hi <- mu_switch - 1e-8
  span <- max(diff(range(d2)), abs(mu_hat - mu_switch), 1, na.rm = TRUE)
  lo <- max(min(d2) / 2 + 1e-6 * span, hi - 5 * span)

  if (hi <= lo) {
    return(if (T1 <= crit) -Inf else NA_real_)
  }

  mu_star <- .sel_root_on_bracket(f, lo, hi)
  if (!is.finite(mu_star)) {
    f_hi <- f(hi)
    if (is.finite(f_hi) && f_hi > 0) {
      return(mu_switch)
    }
    return(-Inf)
  }
  mu_star
}

# ---------------------------------------------------------------------------
# test-2 / less：置信区间上界
# ---------------------------------------------------------------------------

.sel_boundary_less <- function(d1, d2, m, cn, T1, crit, mu_hat) {
  if (T1 > crit) {
    return(NA_real_)
  }

  mu_switch <- (mean(d2) + cn / sqrt(m)) / 2
  span <- max(diff(range(d2)), abs(mu_hat - mu_switch), 1, na.rm = TRUE)
  g <- function(mu) el_sel_stat_less(d1, d2, mu, cn) - crit

  lo <- min(mu_hat, mu_switch) - span
  hi <- max(mu_hat + span, mu_switch + span, 0.25)

  if (g(lo) > 0) {
    lo <- min(mu_hat, mu_switch) - 5 * span
    if (g(lo) > 0) {
      return(Inf)
    }
  }

  if (g(hi) <= 0) {
    for (k in seq_len(15L)) {
      hi <- hi + span
      if (g(hi) > 0) {
        break
      }
    }
    if (g(hi) <= 0) {
      return(Inf)
    }
  }

  if (hi <= lo) {
    return(Inf)
  }

  bracket <- .sel_find_bracket(g, lo, hi)
  if (is.null(bracket)) {
    bracket <- .expand_bracket(g, lo, hi)
  }

  mu_star <- NA_real_
  if (!is.null(bracket)) {
    mu_star <- tryCatch(
      uniroot(g, bracket, extendInt = "no")$root,
      error = function(e) NA_real_
    )
  }

  if (!is.finite(mu_star)) {
    return(Inf)
  }

  mu_star
}

# ---------------------------------------------------------------------------
# test-3 / two.sided：置信区间 [lb, ub]
# ---------------------------------------------------------------------------

.sel_boundary_interval <- function(d1, d2, m, cn, T1, crit, mu_hat, eps_lo, eps_hi) {
  if (T1 > crit) {
    return(list(lb = NA_real_, ub = NA_real_))
  }

  span <- max(
    abs(mu_hat - eps_lo), abs(mu_hat - eps_hi),
    abs(eps_hi - eps_lo), diff(range(d2)) / 4, 0.25,
    na.rm = TRUE
  )
  lo <- min(mu_hat, eps_lo) - 5 * span
  hi <- max(mu_hat, eps_hi) + 5 * span
  mu_grid <- seq(lo, hi, length.out = 600L)
  Ts <- vapply(
    mu_grid, el_sel_stat_interval, numeric(1L),
    d1 = d1, d2 = d2, cn = cn
  )
  ok <- is.finite(Ts) & Ts <= crit
  if (!any(ok)) {
    return(list(lb = NA_real_, ub = NA_real_))
  }

  lb <- min(mu_grid[ok])
  ub <- max(mu_grid[ok])
  if (ok[1L]) {
    lb <- -Inf
  }
  if (ok[length(ok)]) {
    ub <- Inf
  }

  list(lb = lb, ub = ub)
}

# ---------------------------------------------------------------------------
# 三种检验：统一求置信区间边界
# ---------------------------------------------------------------------------

.sel_ci_bounds <- function(d1, d2, m, cn, T1, crit, mu_hat,
                           alternative, mu0, eps_lo, eps_hi) {
  alternative <- .normalize_sel_alternative(alternative)

  if (alternative == "greater") {
    lb <- .sel_boundary_greater(d1, d2, m, cn, T1, crit, mu_hat)
    if (is.na(lb)) {
      return(list(
        lb = NA_real_, ub = NA_real_,
        anchor_lo = mu0, anchor_hi = NA_real_
      ))
    }
    return(list(lb = lb, ub = Inf, anchor_lo = mu0, anchor_hi = NA_real_))
  }

  if (alternative == "less") {
    ub <- .sel_boundary_less(d1, d2, m, cn, T1, crit, mu_hat)
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
  bounds <- .sel_boundary_interval(
    d1, d2, m, cn, T1, crit, mu_hat, eps_lo, eps_hi
  )
  list(
    lb = bounds$lb,
    ub = bounds$ub,
    anchor_lo = eps_lo,
    anchor_hi = eps_hi
  )
}

.sel_al_length <- function(lb, ub) {
  if (is.na(lb) || is.na(ub)) {
    return(NA_real_)
  }
  if (!is.finite(lb) || !is.finite(ub)) {
    return(Inf)
  }
  ub - lb
}

.sel_Ts_at_bound <- function(d1, d2, bound, cn,
                             side = c("greater", "less", "interval", "or")) {
  side <- match.arg(side)
  if (side == "or") {
    side <- "interval"
  }
  if (!is.finite(bound)) {
    return(NA_real_)
  }
  switch(
    side,
    greater = el_sel_stat_greater(d1, d2, bound, cn),
    less = el_sel_stat_less(d1, d2, bound, cn),
    interval = el_sel_stat_interval(d1, d2, bound, cn)
  )
}

.sel_bound_y <- function(Ts_bound, crit) {
  if (!is.finite(Ts_bound)) {
    return(crit)
  }
  if (abs(Ts_bound - crit) < 0.05) {
    crit
  } else {
    Ts_bound
  }
}

# ---------------------------------------------------------------------------
# 剖面网格（作图用）
# ---------------------------------------------------------------------------

.sel_greater_profile <- function(d1, d2, cn, mu_hat, lb,
                                 lb_plot = NULL, pad = NULL, n_grid = 200L) {
  lb_use <- if (is.finite(lb)) lb else if (!is.null(lb_plot)) lb_plot else mu_hat
  if (is.null(pad)) {
    pad <- max(0.5, abs(mu_hat - lb_use), na.rm = TRUE)
  }
  mu_grid <- seq(lb_use - pad * 0.5, mu_hat + pad, length.out = n_grid)
  list(
    mu_grid = mu_grid,
    Ts = vapply(
      mu_grid, el_sel_stat_greater, numeric(1L),
      d1 = d1, d2 = d2, cn = cn
    )
  )
}

.sel_less_profile <- function(d1, d2, cn, mu_hat, ub, mu0_ref = 0,
                              ub_plot = NULL, pad = NULL, n_grid = 200L) {
  mu0_ref <- as.numeric(mu0_ref)
  ub_use <- if (is.finite(ub)) ub else if (!is.null(ub_plot)) ub_plot else mu_hat
  if (is.null(pad)) {
    pad <- max(0.5, 0.15 * diff(range(d2)), abs(mu_hat - ub_use), mu0_ref, na.rm = TRUE)
    if (is.finite(ub_use) && ub_use > 0) {
      pad <- max(pad, ub_use, na.rm = TRUE)
    }
  }
  mu_lo <- min(c(-ub_use, ub_use, mu_hat, 0), na.rm = TRUE) - pad
  mu_hi <- max(c(ub_use, mu_hat, 0), na.rm = TRUE) + pad * 0.25
  mu_grid <- seq(mu_lo, mu_hi, length.out = n_grid)
  list(
    mu0_grid = mu_grid,
    mu_grid = mu_grid,
    Ts = vapply(
      mu_grid, el_sel_stat_less, numeric(1L),
      d1 = d1, d2 = d2, cn = cn
    )
  )
}

.sel_interval_profile <- function(d1, d2, cn, mu_hat, lb, ub, eps_lo, eps_hi,
                                  lb_plot = NULL, ub_plot = NULL, pad = NULL,
                                  x_span = NULL, n_grid = 200L) {
  lb_use <- if (is.finite(lb)) {
    lb
  } else if (!is.null(lb_plot)) {
    lb_plot
  } else {
    min(eps_lo, mu_hat)
  }
  ub_use <- if (is.finite(ub)) {
    ub
  } else if (!is.null(ub_plot)) {
    ub_plot
  } else {
    max(eps_hi, mu_hat)
  }
  if (is.null(pad)) {
    pad <- max(
      0.5,
      if (!is.null(x_span)) 0.15 * x_span else 0,
      abs(mu_hat - lb_use), abs(mu_hat - ub_use),
      abs(eps_hi - eps_lo), na.rm = TRUE
    )
  }
  mu_lo <- min(c(lb_use, ub_use, mu_hat, eps_lo), na.rm = TRUE) - pad
  mu_hi <- max(c(lb_use, ub_use, mu_hat, eps_hi), na.rm = TRUE) + pad
  mu_grid <- seq(mu_lo, mu_hi, length.out = n_grid)
  list(
    mu_grid = mu_grid,
    Ts = vapply(
      mu_grid, el_sel_stat_interval, numeric(1L),
      d1 = d1, d2 = d2, cn = cn
    )
  )
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

.sel_plot_title <- function(alpha, alpha2, sex_name, age_cat_name, plot_label = NULL) {
  if (!is.null(plot_label) && nzchar(trimws(plot_label))) {
    return(trimws(plot_label))
  }
  suffix <- trimws(paste(sex_name, age_cat_name))
  paste0(
    "SEL alpha =", alpha, " alpha2 =", alpha2,
    if (nzchar(suffix)) paste0(" ", suffix) else ""
  )
}

.sel_plot_pad <- function(x, mu_hat) {
  x_sd <- sd(x)
  if (!is.finite(x_sd) || x_sd <= 0) {
    x_sd <- max(abs(x - mu_hat), 0.05, na.rm = TRUE)
  }
  x_span <- diff(range(x))
  if (!is.finite(x_span) || x_span <= 0) {
    x_span <- max(2 * x_sd, 0.05, na.rm = TRUE)
  }
  list(x_sd = x_sd, x_span = x_span)
}

.sel_plot_annotate_bound <- function(bound, Ts_bound, crit, y_max, label, pos = 4,
                                     align_y = FALSE) {
  if (!is.finite(bound)) {
    return(invisible(NULL))
  }
  y_bound <- .sel_bound_y(Ts_bound, crit)
  abline(v = bound, col = "blue", lty = 1, lwd = 1.5)
  points(bound, y_bound, pch = 19, col = "blue", cex = 0.9)
  y_text <- if (align_y) y_bound else y_bound + 0.05 * y_max
  text(
    bound, y_text,
    labels = label, pos = pos, col = "blue", cex = 0.9
  )
}

.sel_plot <- function(x, d1, d2, cn, crit, bounds, alternative,
                      mu_hat, sex_name, age_cat_name, alpha, alpha2,
                      plot_label = NULL, mu0 = 0) {
  alternative <- .normalize_sel_alternative(alternative)
  x <- as.numeric(x)
  pad_info <- .sel_plot_pad(x, mu_hat)

  if (alternative == "greater") {
    lb_plot <- if (is.finite(bounds$lb)) {
      bounds$lb
    } else {
      mu_hat - 2 * pad_info$x_sd / sqrt(length(x))
    }
    pad <- max(0.05, 0.15 * pad_info$x_span, abs(mu_hat - lb_plot), na.rm = TRUE)
    prof <- .sel_greater_profile(
      d1, d2, cn, mu_hat, bounds$lb,
      lb_plot = lb_plot, pad = pad
    )
    epsilonGs <- prof$mu_grid
    Ts <- prof$Ts
  } else if (alternative == "less") {
    ub_plot <- if (is.finite(bounds$ub)) {
      bounds$ub
    } else {
      mu_hat + 2 * pad_info$x_sd / sqrt(length(x))
    }
    pad <- max(
      0.05, 0.15 * pad_info$x_span,
      abs(mu_hat - ub_plot), mu0, na.rm = TRUE
    )
    if (is.finite(ub_plot) && ub_plot > 0) {
      pad <- max(pad, ub_plot, na.rm = TRUE)
    }
    prof <- .sel_less_profile(
      d1, d2, cn, mu_hat, bounds$ub,
      mu0_ref = mu0, ub_plot = ub_plot, pad = pad
    )
    epsilonGs <- prof$mu_grid
    Ts <- prof$Ts
  } else {
    lb_plot <- if (is.finite(bounds$lb)) {
      bounds$lb
    } else {
      min(bounds$anchor_lo, mu_hat) - 2 * pad_info$x_sd / sqrt(length(x))
    }
    ub_plot <- if (is.finite(bounds$ub)) {
      bounds$ub
    } else {
      max(bounds$anchor_hi, mu_hat) + 2 * pad_info$x_sd / sqrt(length(x))
    }
    prof <- .sel_interval_profile(
      d1, d2, cn, mu_hat, bounds$lb, bounds$ub,
      bounds$anchor_lo, bounds$anchor_hi,
      lb_plot = lb_plot, ub_plot = ub_plot,
      x_span = pad_info$x_span
    )
    epsilonGs <- prof$mu_grid
    Ts <- prof$Ts
  }

  par(mfrow = c(1, 1))
  Ts_plot <- ifelse(is.finite(Ts), Ts, NA_real_)
  y_max <- .finite_plot_ylim(Ts_plot, crit)
  plot(
    epsilonGs, Ts_plot,
    type = "l", ylim = c(0, y_max),
    main = .sel_plot_title(alpha, alpha2, sex_name, age_cat_name, plot_label),
    xlab = "epsilonG", ylab = expression(T[s])
  )
  abline(h = crit, col = "red", lty = 2)
  x_usr <- par("usr")
  crit_x <- if (alternative == "two.sided") {
    x_usr[1] + 0.06 * diff(x_usr[1:2])
  } else {
    x_usr[1] + 0.03
  }
  text(
    crit_x, crit,
    labels = paste0("crit=", round(crit, 4)),
    pos = 3, col = "red", cex = 0.9
  )
  if (alternative == "less") {
    abline(v = 0, col = "gray60", lty = 3)
  }

  if (alternative == "greater") {
    .sel_plot_annotate_bound(
      bounds$lb,
      .sel_Ts_at_bound(d1, d2, bounds$lb, cn, "greater"),
      crit, y_max,
      paste0("lb=", round(bounds$lb, 4)),
      pos = 4
    )
  } else if (alternative == "less") {
    if (is.finite(bounds$ub)) {
      .sel_plot_annotate_bound(
        bounds$ub,
        .sel_Ts_at_bound(d1, d2, bounds$ub, cn, "less"),
        crit, y_max,
        paste0("ub=", round(bounds$ub, 4)),
        pos = 2
      )
    } else {
      text(mu_hat, y_max * 0.55, labels = "ub=Inf", pos = 4, col = "blue", cex = 0.9)
    }
  } else {
    if (is.finite(bounds$lb)) {
      .sel_plot_annotate_bound(
        bounds$lb,
        .sel_Ts_at_bound(d1, d2, bounds$lb, cn, "interval"),
        crit, y_max,
        paste0("lb=", round(bounds$lb, 4)),
        pos = 2,
        align_y = TRUE
      )
    }
    if (is.finite(bounds$ub)) {
      .sel_plot_annotate_bound(
        bounds$ub,
        .sel_Ts_at_bound(d1, d2, bounds$ub, cn, "interval"),
        crit, y_max,
        paste0("ub=", round(bounds$ub, 4)),
        pos = 4,
        align_y = TRUE
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

SEL_CI <- function(group_specific_data,
                   target = theta_P,
                   alpha = 0.95,
                   alpha2 = 0.08,
                   alternative = "greater",
                   mu0 = 0,
                   mu1 = NULL,
                   mu2 = NULL,
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
  if (n < 4L) {
    stop("SEL 至少需要 4 个观测（n = 2m, m >= 2）。")
  }
  if (alpha <= 0 || alpha >= 1) {
    stop("alpha 必须介于 (0, 1) 之间。")
  }

  alt <- .normalize_sel_alternative(alternative)
  x <- y - target
  split <- .split_deltas(x)
  d1 <- split$d1
  d2 <- split$d2
  m <- split$m
  mu_hat <- mean(x)
  cn <- .sel_cn(d2, alpha2)
  crit <- qchisq(alpha, 1L)
  T1 <- el_mean_stat(d1, 0)

  if (alt == "two.sided") {
    if (is.null(epsilon_lo) && !is.null(mu1)) {
      epsilon_lo <- -as.numeric(mu1)
    }
    if (is.null(epsilon_hi) && !is.null(mu2)) {
      epsilon_hi <- as.numeric(mu2)
    }
    if (is.null(epsilon_lo) || is.null(epsilon_hi)) {
      stop("two.sided 须同时提供 epsilon_lo 与 epsilon_hi（或 mu1 与 mu2）。")
    }
    eps_lo <- as.numeric(epsilon_lo)
    eps_hi <- as.numeric(epsilon_hi)
  } else {
    eps_lo <- NA_real_
    eps_hi <- NA_real_
    mu0 <- as.numeric(mu0)
  }

  bounds <- .sel_ci_bounds(
    d1, d2, m, cn, T1, crit, mu_hat,
    alt, mu0, eps_lo, eps_hi
  )
  lb <- bounds$lb
  ub <- bounds$ub
  SELAL <- .sel_al_length(lb, ub)

  Ts_lb <- if (alt == "two.sided") {
    .sel_Ts_at_bound(d1, d2, lb, cn, "interval")
  } else {
    NA_real_
  }
  Ts_ub <- if (alt == "less") {
    .sel_Ts_at_bound(d1, d2, ub, cn, "less")
  } else if (alt == "two.sided") {
    .sel_Ts_at_bound(d1, d2, ub, cn, "interval")
  } else {
    NA_real_
  }

  if (plot) {
    .sel_plot(
      x, d1, d2, cn, crit, bounds, alt,
      mu_hat, sex_name, age_cat_name, alpha, alpha2, plot_label, mu0
    )
  }

  cat("alternative ", "n ", " ep_h ", " lb ", " ub ", " SELAL", "\n")
  cat(alt, n, mu_hat, lb, ub, SELAL, "\n")
  cat("样本epsilonG", mu_hat, "\n")

  invisible(list(
    n = n,
    epsilonG = mu_hat,
    lb = lb,
    ub = ub,
    SELAL = SELAL,
    crit = crit,
    cn = cn,
    T1 = T1,
    Ts_lb = Ts_lb,
    Ts_ub = Ts_ub,
    alpha = alpha,
    alpha2 = alpha2,
    alternative = alt,
    mu0 = if (alt != "two.sided") mu0 else NA_real_,
    mu1 = if (alt == "two.sided") -eps_lo else NA_real_,
    mu2 = if (alt == "two.sided") eps_hi else NA_real_,
    epsilon_lo = if (alt == "two.sided") eps_lo else NA_real_,
    epsilon_hi = if (alt == "two.sided") eps_hi else NA_real_
  ))
}
