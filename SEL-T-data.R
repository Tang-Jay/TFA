#==================================================#
#  SEL-T-data：重尾 DGP 功效模拟（SEL.tex §4, test-1）#
#  X_i = mu + Z_i, Z_i ~ t(3)                       #
#==================================================#
rm(list = ls())
library(emplik)

normalize_test <- function(test) {
  key <- tolower(trimws(as.character(test)))
  switch(
    key,
    "test1" = "test1", "1" = "test1", "greater" = "test1", "大于" = "test1",
    "test2" = "test2", "2" = "test2", "less" = "test2", "小于" = "test2",
    "test3" = "test3", "3" = "test3", "or" = "test3", "双侧" = "test3",
    stop("test 须为 test1 / test2 / test3 之一，当前为: ", test)
  )
}

test_config <- function(test, mu0 = 0.25, mu1 = NULL, mu2 = NULL, mu_grid = NULL) {
  test <- normalize_test(test)
  mu1 <- if (is.null(mu1)) mu0 else mu1
  mu2 <- if (is.null(mu2)) mu0 else mu2
  cfg <- switch(
    test,
    test1 = list(
      test = test, tex_ref = "§4, test-1, t(3)+mu",
      mu0 = mu0, mu1 = mu0, mu2 = mu0,
      mu_l = 0.05, mu_u = 0.45, mu_step = 0.01,
      vline = mu0
    ),
    test2 = list(
      test = test, tex_ref = "§3.3, test-2",
      mu0 = mu0, mu1 = mu0, mu2 = mu0,
      mu_l = -0.5, mu_u = 0, mu_step = 0.01,
      vline = -mu0
    ),
    test3 = list(
      test = test, tex_ref = "§3.3, test-3",
      mu0 = mu0, mu1 = mu1, mu2 = mu2,
      mu_l = -0.5, mu_u = 0.5, mu_step = 0.01,
      vline = c(-mu1, mu2)
    )
  )
  if (!is.null(mu_grid)) {
    cfg$mu_l <- mu_grid$mu_l
    cfg$mu_u <- mu_grid$mu_u
    cfg$mu_step <- mu_grid$mu_step
  }
  cfg
}

resolve_an <- function(n, an_method, a2, a) {
  switch(
    an_method,
    "qnorm" = qnorm(1 - a2),
    "theory" = qnorm((1 + sqrt(9 - 8 * a)) / 4)^2,
    "logn" = log(n),
    "logn/2" = log(n) / 6,
    "logn/n" = log(n) / n,
    "loglogn" = log(log(n)),
    stop("an_method 必须是 'qnorm'、'theory'、'logn'、'logn/2'、'logn/n' 或 'loglogn'")
  )
}

make_dgp_rng <- function(dgp) {
  switch(
    dgp,
    t3 = function(n, mu) rt(n, df = 3) + mu,
    stop("未知 dgp: ", dgp)
  )
}

simulate_reject <- function(test, n, mu, mu0, mu1, mu2, an, a, cut_T1, cut_T2, rng = NULL) {
  test <- normalize_test(test)
  if (is.null(rng)) rng <- make_dgp_rng("t3")
  n1 <- floor(n / 2)
  X <- rng(n, mu)
  X1 <- X[seq_len(n1)]
  X2 <- X[(n1 + 1):n]
  m <- min(length(X1), length(X2))
  D1 <- X1[seq_len(m)] - X2[seq_len(m)]
  D2 <- X1[seq_len(m)] + X2[seq_len(m)]
  s2 <- sd(D2)
  xbar <- mean(X)
  T1 <- emplik::el.test(x = as.vector(D1), mu = 0)$`-2LLR`

  if (test == "test1") {
    cn <- an * s2
    T2 <- emplik::el.test(x = as.vector(D2), mu = 2 * mu0)$`-2LLR`
    T2b <- emplik::el.test(x = as.vector(D2), mu = 0)$`-2LLR`
    T3 <- sqrt(m) * (mean(D2) - 2 * mu0)
    T_SEL <- T1 + T2 * (T3 > cn)
    T_ASEL <- T1 + T2b * (T3 > cn)
    T_EL <- emplik::el.test(x = as.vector(X), mu = mu0)$`-2LLR`
    T_CEL <- T_EL * (xbar > mu0)
    rej_t <- t.test(X, mu = mu0, alternative = "greater")$p.value < a
  } else if (test == "test2") {
    cn <- -an * s2
    T2 <- emplik::el.test(x = as.vector(D2), mu = -2 * mu0)$`-2LLR`
    T2b <- emplik::el.test(x = as.vector(D2), mu = 0)$`-2LLR`
    T3 <- sqrt(m) * (mean(D2) + 2 * mu0)
    T_SEL <- T1 + T2 * (T3 < cn)
    T_ASEL <- T1 + T2b * (T3 < cn)
    T_EL <- emplik::el.test(x = as.vector(X), mu = -mu0)$`-2LLR`
    T_CEL <- T_EL * (xbar < -mu0)
    rej_t <- t.test(X, mu = -mu0, alternative = "less")$p.value < a
  } else {
    cn <- an * s2
    T2p <- emplik::el.test(x = as.vector(D2), mu = 2 * mu2)$`-2LLR`
    T2m <- emplik::el.test(x = as.vector(D2), mu = -2 * mu1)$`-2LLR`
    T2b <- emplik::el.test(x = as.vector(D2), mu = 0)$`-2LLR`
    T3p <- sqrt(m) * (mean(D2) - 2 * mu2)
    T3m <- sqrt(m) * (mean(D2) + 2 * mu1)
    T_SEL <- T1 + T2p * (T3p > cn) + T2m * (T3m < -cn)
    T_ASEL <- T1 + T2b * ((T3p > cn) + (T3m < -cn))
    if (xbar > mu2) {
      T_EL <- emplik::el.test(x = as.vector(X), mu = mu2)$`-2LLR`
      T_CEL <- T_EL
    } else if (xbar < -mu1) {
      T_EL <- emplik::el.test(x = as.vector(X), mu = -mu1)$`-2LLR`
      T_CEL <- T_EL
    } else {
      T_EL <- 0
      T_CEL <- 0
    }
    rej_t <- (xbar > mu2 &&
                t.test(X, mu = mu2, alternative = "greater")$p.value < a) ||
      (xbar < -mu1 &&
         t.test(X, mu = -mu1, alternative = "less")$p.value < a)
  }

  c(sel = T_SEL > cut_T1, asel = T_ASEL > cut_T1, cel = T_CEL > cut_T2,
    el = T_EL > cut_T1, t = rej_t)
}

plot_power_curve <- function(mus, Power, cfg, n, alpha,
                             Power_ASEL = NULL, Power_CEL = NULL,
                             Power_EL = NULL, Power_t = NULL) {
  par(oma = c(0, 0, 0, 0), mar = c(2, 3, 1.5, 1))
  plot(mus, Power, xaxs = "i", yaxs = "i", type = "n", xaxt = "n", yaxt = "n",
       ylim = c(0, 1), xlab = "", ylab = "", col = "white")
  axis(1, pretty(c(cfg$mu_l, cfg$mu_u)), mgp = c(4, 0, 0), tcl = 0.5,
       font = 2, lwd = 2, las = 1)
  axis(2, seq(0, 1, 0.1), mgp = c(4, 0.2, 0), tcl = 0.5, font = 2, lwd = 2, las = 1)
  axis(3, pretty(c(cfg$mu_l, cfg$mu_u)), tck = 0.01, labels = FALSE, tcl = 0, lwd = 2)
  axis(4, seq(0, 1, 0.1), labels = FALSE, tcl = 0, lwd = 2)
  main_title <- switch(
    cfg$test,
    test1 = bquote(H[0]~":"~0 <= ~mu <=  mu[0]~" vs "~H[1]~":"~mu > mu[0]~"  n = "~.(n)),
    test2 = bquote(H[0]~":"-mu[0] <= ~mu <= ~0~" vs "~H[1]~":"~mu < -mu[0]~"  n = "~.(n)),
    test3 = bquote(H[0]~":"~-mu[1] <= ~ mu <= ~ mu[2]~" vs "~H[1]~":"~mu < -mu[1]~" or "~mu > mu[2]~"  n = "~.(n))
  )
  title(main = main_title, cex.main = 1.1, font.main = 2)
  title(xlab = "true mean (mu)", line = 1, col.lab = 1, font.lab = 2, cex.lab = 1)
  title(ylab = "power", line = 2, col.lab = 1, font.lab = 2, cex.lab = 1)
  abline(h = alpha, lwd = 1.5, col = "gray")
  abline(v = cfg$vline, lty = 2, col = "gray")
  lines(mus, Power, col = "#ff8c00", lty = 2, lwd = 2)
  points(mus, Power, pch = 17, cex = 0.6, col = "#ff8c00")
  if (!is.null(Power_ASEL)) {
    lines(mus, Power_ASEL, col = "#cc0000", lty = 1, lwd = 2)
    points(mus, Power_ASEL, pch = 19, cex = 0.55, col = "#cc0000")
  }
  if (!is.null(Power_CEL)) {
    lines(mus, Power_CEL, col = "#8844cc", lty = 5, lwd = 2)
    points(mus, Power_CEL, pch = 2, cex = 0.55, col = "#8844cc")
  }
  if (!is.null(Power_EL)) {
    lines(mus, Power_EL, col = "gray40", lty = 3, lwd = 2)
    points(mus, Power_EL, pch = 0, cex = 0.55, col = "gray40")
  }
  if (!is.null(Power_t)) {
    lines(mus, Power_t, col = "#22aa66", lty = 6, lwd = 2)
    points(mus, Power_t, pch = 3, cex = 0.55, col = "#22aa66")
  }
  leg <- col <- lty <- pch <- NULL
  if (!is.null(Power_ASEL)) {
    leg <- c(leg, "ASEL"); col <- c(col, "#cc0000"); lty <- c(lty, 1); pch <- c(pch, 19)
  }
  leg <- c(leg, "SEL"); col <- c(col, "#ff8c00"); lty <- c(lty, 2); pch <- c(pch, 17)
  if (!is.null(Power_CEL)) {
    leg <- c(leg, "CEL"); col <- c(col, "#8844cc"); lty <- c(lty, 5); pch <- c(pch, 2)
  }
  if (!is.null(Power_EL)) {
    leg <- c(leg, "EL"); col <- c(col, "gray40"); lty <- c(lty, 3); pch <- c(pch, 0)
  }
  if (!is.null(Power_t)) {
    leg <- c(leg, "t"); col <- c(col, "#22aa66"); lty <- c(lty, 6); pch <- c(pch, 3)
  }
  legend("topleft", legend = leg, col = col, lty = lty, lwd = 2,
         pch = pch, bty = "o", cex = 0.85)
}

draw_power_plot <- function(mus, Power_SEL, cfg, n, alpha, Power_ASEL, Power_CEL, Power_EL, Power_t) {
  par(mfrow = c(1, 1), mar = c(3, 3, 2, 1))
  plot_power_curve(
    mus, Power_SEL, cfg, n, alpha = alpha,
    Power_ASEL = Power_ASEL, Power_CEL = Power_CEL,
    Power_EL = Power_EL, Power_t = Power_t
  )
}

sim_file_stem <- function(sim_prefix, test) {
  paste0(sim_prefix, "-", normalize_test(test))
}

run_power_sim <- function(test, ns, nsim, a, a2, mu0, mu1, mu2,
                          an_method, save_data, show_plot, save_plot, power_dir,
                          sim_prefix = "SEL-T",
                          mu_grid = NULL,
                          rng = NULL) {
  test <- normalize_test(test)
  cfg <- test_config(test, mu0 = mu0, mu1 = mu1, mu2 = mu2, mu_grid = mu_grid)
  mus <- seq(cfg$mu_l, cfg$mu_u, cfg$mu_step)
  l <- length(mus)
  cut_T1 <- qchisq(1 - a, 1)
  cut_T2 <- qchisq(1 - 2 * a, 1)
  stem <- sim_file_stem(sim_prefix, test)

  Power_SEL <- Power_ASEL <- Power_CEL <- Power_EL <- Power_t <- rep(0, l)
  power <- matrix(0, nrow = l, ncol = 5)
  colnames(power) <- c("SEL", "ASEL", "CEL", "EL", "t")

  cat("\n==========", sim_prefix, test, "（SEL.tex）==========\n")
  for (n in ns) {
    an <- resolve_an(n, an_method, a2, a)
    cat("a2=", a2, " an =", an, " c =", cut_T1, " n =", n, "\n")
    for (j in seq_along(mus)) {
      mu <- mus[j]
      rej <- replicate(
        nsim,
        simulate_reject(test, n, mu, cfg$mu0, cfg$mu1, cfg$mu2, an, a, cut_T1, cut_T2, rng = rng),
        simplify = "array"
      )
      Power_SEL[j] <- mean(rej["sel", ])
      Power_ASEL[j] <- mean(rej["asel", ])
      Power_CEL[j] <- mean(rej["cel", ])
      Power_EL[j] <- mean(rej["el", ])
      Power_t[j] <- mean(rej["t", ])
      power[j, ] <- c(Power_SEL[j], Power_ASEL[j], Power_CEL[j], Power_EL[j], Power_t[j])
      cat("mu =", round(mu, 3),
          " SEL =", round(Power_SEL[j], 4), " ASEL =", round(Power_ASEL[j], 4),
          " CEL =", round(Power_CEL[j], 4),
          " EL =", round(Power_EL[j], 4),
          " t =", round(Power_t[j], 4), "\n")
    }
    cat("完成", stem, " n =", n, "\n\n")

    a2_file <- formatC(a2, format = "f", digits = 4)
    if (save_data) {
      data_path <- file.path(power_dir, paste0(stem, "-nsim-", nsim, "-a2-", a2_file, "-n-", n, ".RData"))
      mu_l <- cfg$mu_l
      mu_u <- cfg$mu_u
      save(
        mus, power, mu_l, mu_u, mu0, mu1, mu2,
        test, n, a, a2, nsim, an_method, an,
        file = data_path
      )
      cat("数据已保存:", data_path, "\n\n")
    }
    if (show_plot) {
      draw_power_plot(mus, Power_SEL, cfg, n, alpha = a,
                      Power_ASEL = Power_ASEL, Power_CEL = Power_CEL,
                      Power_EL = Power_EL, Power_t = Power_t)
      cat("已绘制功效曲线:", stem, " n =", n, "\n\n")
    }
    if (save_plot) {
      file_path <- file.path(power_dir, paste0(stem, "-nsim-", nsim, "-a2-", a2_file, "-n-", n, ".png"))
      png(file_path, width = 7, height = 5, units = "in", res = 300)
      draw_power_plot(mus, Power_SEL, cfg, n, alpha = a,
                      Power_ASEL = Power_ASEL, Power_CEL = Power_CEL,
                      Power_EL = Power_EL, Power_t = Power_t)
      dev.off()
      cat("图片已保存:", file_path, "\n\n")
    }
  }
  invisible(list(cfg = cfg, mus = mus, power = power))
}

#================================================#
#  全局参数
#================================================#
set.seed(2)
sim_prefix <- "SEL-T"
dgp <- "t3"
ns <- c(30,100,300,1000)
a2 <- 0.08
mu0 <- 0.25
mu1 <- 0.25
mu2 <- 0.25
an_method <- "qnorm"
nsim <- 2000
a <- 0.05
save_data <- T
show_plot <- TRUE
save_plot <- FALSE
power_dir <- "."
test <- "test1"
run_all_tests <- T
mu_grid_test1 <- list(mu_l = 0, mu_u = 0.5, mu_step = 0.01)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  test <- args[1]
  run_all_tests <- FALSE
}
if (length(args) >= 2) {
  run_all_tests <- as.logical(args[2])
}

tests <- if (run_all_tests) c("test1", "test2", "test3") else normalize_test(test)
rng <- make_dgp_rng(dgp)

for (test_cur in tests) {
  mg <- if (test_cur == "test1") mu_grid_test1 else NULL
  run_power_sim(
    test = test_cur, ns = ns, nsim = nsim, a = a, a2 = a2,
    mu0 = mu0, mu1 = mu1, mu2 = mu2, an_method = an_method,
    save_data = save_data, show_plot = show_plot, save_plot = save_plot,
    power_dir = power_dir, sim_prefix = sim_prefix,
    mu_grid = mg, rng = rng
  )
}
