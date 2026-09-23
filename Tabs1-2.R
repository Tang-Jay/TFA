#==============================================#
#       读取已有 CSV，在控制台打印表 1 和表 2    #
#==============================================#
# 仅使用 base R；在 RStudio 点击 Source 即可运行。
# fdr 是逐次 FDP 的均值，对应手稿中的 FFR；不重新运行模拟。

# ---- 数据读取 ----
read_table_data <- function(data_dir, n, k, methods) {
  codes <- c(Normal = "N", Exponential = "E")
  batch_columns <- c("nsim", "seed", "q", "epsilon0", "theta0", "a2",
                     "code_version", "sampling_design", "theta_mode",
                     "replicate_start", "replicate_end", "replicate_policy")
  required <- c("method", "dgp", "N", "K", "target_pi1", "null_epsilon",
                "effect_size", "fdr", "power", "m1", batch_columns)
  tables <- list()
  for (dgp in names(codes)) for (pi1 in c(0, 1)) for (method in methods) {
    file_path <- file.path(data_dir, sprintf("%s-%s-n%s-K%s-pi%s.csv",
                                           method, codes[[dgp]], n, k, pi1))
    if (!file.exists(file_path)) stop("找不到文件：", file_path)
    d <- read.csv(file_path, stringsAsFactors = FALSE)
    missing <- setdiff(required, names(d))
    if (length(missing)) stop(file_path, " 缺少列：", paste(missing, collapse = ", "))
    if (!nrow(d) || anyNA(d[setdiff(required, "power")])) {
      stop("数据为空或必要列含缺失值：", file_path)
    }
    if (any(d$method != paste0(method, "-BH") | d$dgp != dgp |
            d$N != n | d$K != k | d$target_pi1 != pi1 | d$m1 != k * pi1)) {
      stop("文件内容与方法或场景不符：", file_path)
    }
    metric <- if (pi1 == 0) "fdr" else "power"
    x_name <- if (pi1 == 0) "null_epsilon" else "effect_size"
    if (!is.numeric(d[[metric]]) || any(!is.finite(d[[metric]])) ||
        any(d[[metric]] < 0 | d[[metric]] > 1) ||
        any(!is.finite(d[[x_name]])) || anyDuplicated(round(d[[x_name]], 12))) {
      stop("指标无效或横轴网格重复：", file_path)
    }
    if (pi1 == 0) {
      if (any(!is.na(d$power)) || any(d$effect_size != 0 |
          d$null_epsilon < 0 | d$null_epsilon > d$epsilon0)) {
        stop("全原假设配置无效，Power 应为 NA：", file_path)
      }
    } else if (any(d$effect_size <= 0 | d$null_epsilon != d$epsilon0)) {
      stop("备择配置应满足 delta > 0 且 null_epsilon = epsilon0：", file_path)
    }
    tables[[length(tables) + 1L]] <- d[required]
  }
  results <- do.call(rbind, tables)
  for (column in batch_columns) {
    if (length(unique(results[[column]])) != 1L) stop("数据批次不一致：", column)
  }
  rownames(results) <- NULL
  results
}

# ---- 按手稿指定网格提取数值 ----
build_panel <- function(results, dgp, pi1, grid, methods) {
  x_name <- if (pi1 == 0) "null_epsilon" else "effect_size"
  metric <- if (pi1 == 0) "fdr" else "power"
  panel <- data.frame(grid)
  names(panel) <- if (pi1 == 0) "epsilon" else "delta"
  for (method in methods) {
    d <- results[results$dgp == dgp & results$target_pi1 == pi1 &
                   results$method == paste0(method, "-BH"), , drop = FALSE]
    index <- match(round(grid, 12), round(d[[x_name]], 12))
    # 指数分布 epsilon=0 未设计实验，只以 NA / "-" 占位。
    absent_by_design <- dgp == "Exponential" & pi1 == 0 & grid == 0
    if (any(is.na(index) & !absent_by_design)) {
      stop(dgp, " / ", method, " 缺少所需网格点。")
    }
    values <- d[[metric]][index]
    values[absent_by_design] <- NA_real_
    panel[[paste0(method, "-BH")]] <- values
  }
  panel
}

# ---- 控制台输出：两种设计并排，保留手稿的小数位数 ----
print_table <- function(normal, exponential, title, digits) {
  format_panel <- function(panel) {
    panel[[1]] <- sprintf("%.2f", panel[[1]])
    for (j in seq.int(2L, ncol(panel))) {
      panel[[j]] <- ifelse(is.na(panel[[j]]), "-",
                           sprintf(paste0("%.", digits, "f"), panel[[j]]))
    }
    panel
  }
  left <- capture.output(print(format_panel(normal), row.names = FALSE, right = TRUE))
  right <- capture.output(print(format_panel(exponential), row.names = FALSE, right = TRUE))
  width <- max(nchar(left), nchar("N(epsilon, 1)"))
  cat("\n", title, "\n", sep = "")
  cat(sprintf("%-*s    %s\n", width, "N(epsilon, 1)", "Exp(1/epsilon)"))
  cat(sprintf("%-*s    %s", width, left, right), sep = "\n")
  cat("\n")
  invisible(list(Normal = normal, Exponential = exponential))
}

# ---- 参数设置与运行 ----
project_dir <- path.expand("~/zotero-obsidian/05 归档/SEL/3 定稿/v2/FDR")
data_dir <- file.path(project_dir, "data")
n <- 1000
k <- 26
methods <- c("CEL", "SEL", "ASEL", "t")
null_epsilons <- seq(0, 0.25, 0.05)
normal_deltas <- c(0.05, 0.20, 0.25, 0.35, 0.50, 0.60)
exponential_deltas <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.60)

results <- read_table_data(data_dir, n, k, methods)
cat(sprintf("n = %d, K = %d, B = %d, q = %.2f, epsilon0 = %.2f, alpha2 = %.2f\n",
            n, k, results$nsim[1], results$q[1], results$epsilon0[1], results$a2[1]))
table1 <- print_table(
  build_panel(results, "Normal", 0, null_epsilons, methods),
  build_panel(results, "Exponential", 0, null_epsilons, methods),
  title = "Table 1: Empirical FFR under the global null (pi1 = 0)", digits = 2
)
cat("Tolerance boundary: epsilon = ", results$epsilon0[1],
    "; '-' = configuration not simulated.\n", sep = "")
table2 <- print_table(
  build_panel(results, "Normal", 1, normal_deltas, methods),
  build_panel(results, "Exponential", 1, exponential_deltas, methods),
  title = "Table 2: Empirical power for selected configurations (pi1 = 1)", digits = 4
)
# table1 / table2 保留未四舍五入的数值，供后续检查；本脚本不写结果文件。
