#==================================================#
#  读取 *-data.R 生成的 .RData 并绘制功效图          #
#  先运行 SEL-E/N/T/C-data.R（文件名含 test）      #
#  图片保存至当前目录，格式 eps                          #
#==================================================#
rm(list = ls())
#================================================#
#  读图参数（与 SEL-*-data.R 保存文件名一致）
#  数据：data/SEL-N-test1-nsim-2000-a2-0.0800-n-100.RData
#  出图：SEL-N-test1-nsim-2000-a2-0.0800-n-100.eps
#================================================#

test <- "test1" # "test1","test2","test3"
sim_prefixes <- c("SEL-E", "SEL-N", "SEL-T", "SEL-C") # "SEL-E", "SEL-N", "SEL-T", "SEL-C"  
nsims <- c(2000)
ns <- c(30,60,100, 300, 1000) # ns <- c(30,60,100, 300, 1000)
a2 <- 0.08
power_dir <- "data"
fig_dir <- "."
fig_ext <- "eps"
save_plot <- F

power_rdata_path <- function(sim_prefix, test, nsim, a2, n, dir = power_dir) {
  a2_file <- formatC(a2, format = "f", digits = 4)
  file.path(dir, paste0(sim_prefix, "-", test, "-nsim-", nsim, "-a2-", a2_file, "-n-", n, ".RData"))
}

power_fig_path <- function(sim_prefix, test, nsim, a2, n,
                           dir = fig_dir, ext = fig_ext) {
  a2_file <- formatC(a2, format = "f", digits = 4)
  file.path(
    dir,
    paste0(sim_prefix, "-", test, "-nsim-", nsim, "-a2-", a2_file, "-n-", n, ".", ext)
  )
}

.open_fig_device <- function(fig_file, width = 7, height = 5) {
  postscript(
    fig_file,
    width = width,
    height = height,
    horizontal = FALSE,
    onefile = FALSE,
    paper = "special",
    family = "Helvetica"
  )
}

#=================================================#
# 功效图（与 *-data.R 一致：名义水平横线 + 临界值垂线）
#=================================================#
normalize_test_id <- function(test) {
  switch(
    tolower(trimws(as.character(test))),
    test1 = "test1", "1" = "test1", greater = "test1",
    test2 = "test2", "2" = "test2", less = "test2",
    test3 = "test3", "3" = "test3", or = "test3",
    tolower(trimws(as.character(test)))
  )
}

power_test_vline <- function(test, mu0, mu1, mu2) {
  switch(
    normalize_test_id(test),
    test1 = mu0,
    test2 = -mu0,
    test3 = c(-mu1, mu2),
    mu0
  )
}

plot_power_curve <- function(mus, Power, mu_l, mu_u, mu0, n, alpha,
                             test = "test1",
                             mu1 = mu0, mu2 = mu0,
                             vline = NULL,
                             Power_ASEL = NULL,
                             Power_CEL = NULL,
                             Power_EL = NULL,
                             Power_t = NULL) {
  if (is.null(vline)) {
    vline <- power_test_vline(test, mu0, mu1, mu2)
  }
  par(oma = c(0, 0, 0, 0), mar = c(2, 3, 2, 1))
  plot(mus, Power, xaxs = "i", yaxs = "i", type = "n", xaxt = "n", yaxt = "n",
       ylim = c(0, 1), xlab = "", ylab = "", col = "white")
  axis(1, pretty(c(mu_l, mu_u)), mgp = c(4, 0, 0), tcl = 0.5, font = 2, lwd = 2, las = 1)
  axis(2, seq(0, 1, 0.1), mgp = c(4, 0.2, 0), tcl = 0.5, font = 2, lwd = 2, las = 1)
  axis(3, pretty(c(mu_l, mu_u)), tck = 0.01, labels = FALSE, tcl = 0, lwd = 2)
  axis(4, seq(0, 1, 0.1), labels = FALSE, tcl = 0, lwd = 2)
  main_title <- switch(
    normalize_test_id(test),
    test1 = bquote(H[0]~":"~0 <= ~mu <=  mu[0]~" vs "~H[1]~":"~mu > mu[0]~"  n = "~.(n)),
    test2 = bquote(H[0]~":"-mu[0] <= ~mu <= ~0~" vs "~H[1]~":"~mu < -mu[0]~"  n = "~.(n)),
    test3 = bquote(H[0]~":"~-mu[1] <= ~ mu <= ~ mu[2]~" vs "~H[1]~":"~mu < -mu[1]~" or "~mu > mu[2]~"  n = "~.(n)),
    bquote(H[0]~":"~(0 <= mu) <= mu[0]~" vs "~H[1]~":"~mu > mu[0]~"  n = "~.(n))
  )
  title(main = main_title, cex.main = 1.25, font.main = 2)
  title(xlab = "true mean (mu)", line = 1, col.lab = 1, font.lab = 2, cex.lab = 1)
  title(ylab = "rejection frequency", line = 2, col.lab = 1, font.lab = 2, cex.lab = 1)
  abline(h = alpha, lwd = 1.5, col = "gray")
  abline(v = vline, lty = 2, col = "gray")
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

#================================================#
#  循环读取并作图
#================================================#
if (save_plot && !dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

for (sim_prefix in sim_prefixes) {
  for (nsim in nsims) {
    for (n in ns) {
      power_file <- power_rdata_path(sim_prefix, test, nsim, a2, n, power_dir)
      if (!file.exists(power_file)) {
        warning(
          "未找到 ", power_file, "，跳过 ", sim_prefix, " ", test, " n = ", n,
          "（请先运行 ", sim_prefix, "-data.R，test = ", test, "，nsim = ", nsim, "）。"
        )
        next
      }
      load(power_file)
      cat("已读取:", power_file, "  power 维度:", paste(dim(power), collapse = " x "),
          "  test =", test, "  nsim =", nsim, "\n")

      Power_SEL  <- power[, "SEL"]
      Power_ASEL <- power[, "ASEL"]
      Power_CEL  <- power[, "CEL"]
      Power_EL   <- power[, "EL"]
      Power_t    <- power[, "t"]

      if (save_plot) {
        fig_file <- power_fig_path(sim_prefix, test, nsim, a2, n, fig_dir, fig_ext)
        .open_fig_device(fig_file)
      }
      par(mfrow = c(1, 1), mar = c(3, 3, 2, 1))
      mu1_plot <- if (exists("mu1")) mu1 else mu0
      mu2_plot <- if (exists("mu2")) mu2 else mu0
      test_plot <- if (exists("test")) test else "test1"
      plot_power_curve(
        mus, Power_SEL, mu_l, mu_u, mu0, n, alpha = a,
        test = test_plot, mu1 = mu1_plot, mu2 = mu2_plot,
        Power_ASEL = Power_ASEL,
        Power_CEL = Power_CEL, Power_EL = Power_EL,
        Power_t = Power_t
      )
      if (save_plot) {
        dev.off()
        cat("图片已保存:", fig_file, "\n\n")
      }
    }
  }
}
