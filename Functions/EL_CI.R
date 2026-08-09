#===============================================
#                    Compas  AL                #
#===============================================
# epsilonG 置信区间：估计方程 g_i = Y_i - theta_P - epsilonG
# 等价于对 (Y_i - theta_P) 求均值的经验似然（uniroot 求根）
#===============================================

source("Functions/GlambdaChen.R")

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

EL_CI <- function(group_specific_data, target = theta_P, alpha = 0.95,
                   plot = TRUE, plot_label = NULL) {
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

  x <- y - target
  cut <- exp(-qchisq(alpha, 1L) / 2)
  elrMax_epsilonG <- mean(x)
  lb <- .boundary_root(x, cut, elrMax_epsilonG, "left")
  ub <- .boundary_root(x, cut, elrMax_epsilonG, "right")
  ELAL <- ub - lb

  if (plot) {
    pad <- max(0.05, 0.1 * ELAL, na.rm = TRUE)
    a <- lb - pad
    b <- ub + pad
    if (is.na(a) || is.na(b)) {
      a <- elrMax_epsilonG - 0.2
      b <- elrMax_epsilonG + 0.2
    }
    epsilonGs <- seq(a, b, length.out = 200L)
    elRatio <- vapply(epsilonGs, el_mean_ratio, numeric(1L), x = x)

    par(mfrow = c(1, 1))
    plot_title <- if (!is.null(plot_label) && nzchar(trimws(plot_label))) {
      trimws(plot_label)
    } else {
      paste0("alpha =", alpha, " ", sex_name, "  ", age_cat_name)
    }
    plot(
      epsilonGs, elRatio,
      type = "l", xlim = c(a, b), ylim = c(0, 1),
      main = plot_title,
      xlab = "epsilonG", ylab = "elr"
    )
    abline(h = cut, col = "red")
    abline(v = elrMax_epsilonG, col = "blue1")
    abline(v = lb, col = "blue1")
    abline(v = ub, col = "blue1")
    text(elrMax_epsilonG + 0.03 * (b - a), 0.95, labels = paste0("el_h=", round(elrMax_epsilonG, 4)), col = "blue1")
    text(ub + 0.02 * (b - a), 0.2, labels = paste0("ub=", round(ub, 4)), col = "blue1")
    text(lb - 0.02 * (b - a), 0.2, labels = paste0("lb=", round(lb, 4)), col = "blue1")
  }

  cat("n ", " ep_h1", " lb ", " ub ", " ELAL", "\n")
  cat(n, elrMax_epsilonG, lb, ub, ELAL, "\n")
  cat("样本epsilonG", mean(y) - target, "\n")

  invisible(list(
    n = n,
    epsilonG = elrMax_epsilonG,
    lb = lb,
    ub = ub,
    ELAL = ELAL,
    alpha = alpha,
    cut = cut
  ))
}
