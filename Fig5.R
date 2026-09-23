#==============================================#
#       读取 SEL 模拟结果并绘制 FFR–Power 图     #
#==============================================#
rm(list = ls())
# 独立运行，仅使用 R 自带包；在 RStudio 修改末尾参数后点击 Source。

#==============================================#
#                 数据读取                     #
#==============================================#
read_simulation_data <- function(data_dir, dgp, n, k, methods,
                                 effect_sizes, null_epsilons) {
  codes <- c(Normal = "N", Exponential = "E")
  if (!dgp %in% names(codes)) stop("dgp 请选择 Normal 或 Exponential。")
  if (!length(methods) || anyNA(methods) || anyDuplicated(methods) ||
      any(!methods %in% c("CEL", "SEL", "ASEL", "t"))) stop("methods 设置无效。")
  tables <- list()
  required <- c("method", "dgp", "N", "K", "target_pi1", "null_epsilon",
                "effect_size", "fdr", "power", "m1", "q", "epsilon0", "theta0",
                "a2", "nsim", "seed", "code_version")
  for (method in methods) for (pi1 in c(0, 0.5, 1)) {
    file_path <- file.path(data_dir, sprintf("%s-%s-n%s-K%s-pi%s.csv",
                                           method, codes[[dgp]], n, k, pi1))
    if (!file.exists(file_path)) stop("找不到数据文件：", file_path)
    d <- utils::read.csv(file_path, stringsAsFactors = FALSE)
    missing <- setdiff(required, names(d))
    if (length(missing)) stop(basename(file_path), " 缺少列：", paste(missing, collapse = ", "))
    if (!nrow(d) || anyNA(d[setdiff(required, "power")])) stop("数据为空或含缺失值：", file_path)
    if (any(d$method != paste0(method, "-BH") | d$dgp != dgp |
            d$N != n | d$K != k | d$target_pi1 != pi1 | d$m1 != k * pi1)) {
      stop("文件内容与方法或场景设置不符：", file_path)
    }
    numeric_columns <- c("null_epsilon", "effect_size", "fdr", "q", "epsilon0", "theta0", "a2")
    if (any(!vapply(d[numeric_columns], function(x) is.numeric(x) && all(is.finite(x)), logical(1))) ||
        any(d$fdr < 0 | d$fdr > 1) ||
        (pi1 > 0 && (anyNA(d$power) || any(!is.finite(d$power) | d$power < 0 | d$power > 1)))) {
      stop("FFR、Power 或场景参数包含无效数值：", file_path)
    }
    if (pi1 == 0) {
      if (any(d$effect_size != 0 | d$null_epsilon < 0 | d$null_epsilon > d$epsilon0)) {
        stop("全原假设参数不在容忍区间内：", file_path)
      }
    } else if (any(d$effect_size <= 0 | d$null_epsilon != d$epsilon0)) {
      stop("备择配置必须以容忍边界为基准且 delta > 0：", file_path)
    }
    x_name <- if (pi1 == 0) "null_epsilon" else "effect_size"
    expected <- if (pi1 == 0) null_epsilons else effect_sizes
    if (pi1 == 0 && dgp == "Exponential") expected <- expected[expected > 0]
    x <- round(d[[x_name]], 12)
    if (!length(expected) || anyDuplicated(x) || !all(round(expected, 12) %in% x)) {
      stop("横轴网格缺失或重复：", file_path)
    }
    d <- d[x %in% round(expected, 12), required, drop = FALSE]
    tables[[length(tables) + 1L]] <- d[order(d[[x_name]]), , drop = FALSE]
  }
  results <- do.call(rbind, tables)
  for (column in c("nsim", "seed", "q", "epsilon0", "theta0", "a2", "code_version")) {
    if (length(unique(results[[column]])) != 1L) stop("数据批次参数不一致：", column)
  }
  rownames(results) <- NULL
  message(dgp, "，K = ", k, "：读取 ", length(tables), " 张表，共 ", nrow(results), " 行。")
  results
}

#==============================================#
#                 结果作图                     #
#==============================================#
draw_ffr_power <- function(results, methods) {
  # 固定顺序：pi1=0 的 FFR；pi1=0.5 的 FFR/Power；pi1=1 的 Power。
  # CSV 的 fdr 列作为 FFR 纵轴，不重新计算指标。
  panel_specs <- data.frame(pi1 = c(0, 0.5, 0.5, 1),
                            metric = c("fdr", "fdr", "power", "power"))
  methods_plot <- intersect(c("CEL", "SEL", "ASEL", "t"), methods)
  colors <- c(CEL = "red", SEL = "orange", ASEL = "blue", t = "green3")
  points <- c(CEL = 2, SEL = 17, ASEL = 19, t = 3)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(1, 4), oma = c(0, 0, 1.8, 0),
                mar = c(2.8, 3.2, 1.7, 0.75), mgp = c(1.8, 0.6, 0), las = 1,
                cex = 0.95, cex.main = 1.05, cex.lab = 1, cex.axis = 0.9)
  for (i in seq_len(nrow(panel_specs))) {
    pi1 <- panel_specs$pi1[i]
    metric <- panel_specs$metric[i]
    panel <- results[results$target_pi1 == pi1, , drop = FALSE]
    x_name <- if (pi1 == 0) "null_epsilon" else "effect_size"
    x_minor <- sort(unique(panel[[x_name]]))
    x_range <- c(0, if (pi1 == 0) results$epsilon0[1] else max(x_minor))
    x_major <- if (pi1 == 0) x_minor else seq(0, x_range[2], by = 0.1)
    graphics::plot(NA, xlim = x_range, ylim = c(0, 1), xaxs = "i", yaxs = "i",
                   xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                   main = sprintf("target pi1 = %.2f", pi1))
    graphics::abline(v = x_major, h = seq(0, 1, 0.2), col = "gray88", lty = 1)
    graphics::axis(1, at = x_minor, labels = FALSE, lwd = 1.5, tcl = 0.35)
    graphics::axis(1, at = x_major, labels = sprintf(if (pi1 == 0) "%.2f" else "%.1f", x_major),
                   font = 2, lwd = 1.5, tcl = 0.35)
    graphics::axis(2, at = seq(0, 1, 0.2), labels = sprintf("%.1f", seq(0, 1, 0.2)),
                   font = 2, lwd = 1.5, tcl = 0.35)
    graphics::axis(3, labels = FALSE, lwd = 1.5, tcl = 0)
    graphics::axis(4, labels = FALSE, lwd = 1.5, tcl = 0)
    graphics::title(xlab = if (pi1 == 0) expression(epsilon) else expression(delta), line = 1.6, font.lab = 2)
    graphics::title(ylab = if (metric == "fdr") "FFR" else "Power", line = 2.1, font.lab = 2)
    if (metric == "fdr") graphics::abline(h = results$q[1], col = "gray35", lty = 2, lwd = 1.5)
    for (method in methods_plot) {
      curve <- panel[panel$method == paste0(method, "-BH"), , drop = FALSE]
      curve <- curve[order(curve[[x_name]]), , drop = FALSE]
      graphics::lines(curve[[x_name]], curve[[metric]], col = colors[[method]], lty = 1, lwd = 2)
      graphics::points(curve[[x_name]], curve[[metric]], col = colors[[method]], pch = points[[method]], cex = 0.65)
    }
    if (i == 1L) {
      graphics::legend("topright", legend = paste0(methods_plot, "-BH"), col = unname(colors[methods_plot]),
                       lty = 1, pch = unname(points[methods_plot]), lwd = 1.8,
                       bg = "white", cex = 0.70, inset = 0.005, x.intersp = 0.7, y.intersp = 0.9)
    }
  }
  title <- paste0(paste(methods_plot, collapse = "/"), " + BH")
  graphics::mtext(sprintf("%s (%s, K = %d, N = %s)", title, results$dgp[1], results$K[1], results$N[1]),
                  side = 3, outer = TRUE, line = 0.3, font = 2, cex = 1.35)
  invisible(results)
}

save_ffr_power_plot <- function(results, methods, file_path, width, height) {
  grDevices::postscript(file_path, width = width, height = height, paper = "special",
                       horizontal = FALSE, onefile = FALSE, family = "Helvetica")
  on.exit(grDevices::dev.off(), add = TRUE)
  draw_ffr_power(results, methods)
  invisible(file_path)
}

#==============================================#
#                 参数与运行                   #
#==============================================#
process_results <- function(results, methods, output_dir, save_plot, plot_width, plot_height) {
  if (any(!is.finite(c(plot_width, plot_height))) || any(c(plot_width, plot_height) <= 0)) {
    stop("图片宽高必须为正数，单位为英寸。")
  }
  codes <- c(Normal = "N", Exponential = "E")
  file_path <- file.path(output_dir, sprintf("all-%s-n%s-K%s.eps",
                                           codes[[results$dgp[1]]], results$N[1], results$K[1]))
  if (save_plot) {
    if (!dir.exists(output_dir)) stop("输出目录不存在：", output_dir)
    save_ffr_power_plot(results, methods, file_path, width = plot_width, height = plot_height)
    message("图片已保存：", file_path)
  }
  if (interactive()) draw_ffr_power(results, methods)
  invisible(list(table = results, plot_path = if (save_plot) file_path else character()))
}

#==============================================#
#                 参数设置                     #
#==============================================#
project_dir <- path.expand("~/zotero-obsidian/05 归档/SEL/3 定稿/v2/FDR")
data_dir <- file.path(project_dir, "data")
output_dir <- data_dir
n <- 1000                           # 每张图对应一个总样本量
ks <- c(26)
dgps <- c("Normal", "Exponential")
methods <- c("CEL", "SEL", "ASEL", "t") # 可选单个方法，如 "CEL"
effect_sizes <- seq(0.05, 0.60, 0.05)
null_epsilons <- seq(0, 0.25, 0.05)   # 指数设计自动排除 epsilon=0
save_plot <- TRUE                   # FALSE：不保存；仍可在 RStudio 中显示
plot_width <- 16                 # 整张 EPS 的宽度，英寸
plot_height <- 5.5                   # 整张 EPS 的高度，英寸

for (k in ks) for (dgp in dgps) {
  results <- read_simulation_data(data_dir, dgp, n, k, methods, effect_sizes, null_epsilons)
  process_results(results, methods, output_dir, save_plot, plot_width, plot_height)
}
