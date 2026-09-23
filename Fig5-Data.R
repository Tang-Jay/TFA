#==============================================#
#                 SEL 多组 FDR 与 Power 模拟实验
#==============================================#
rm(list = ls())

# 单文件实现：按算法顺序阅读；在 RStudio 修改末尾参数后点击 Source。
# 经验似然使用 emplik，不加载外部函数文件。
if (!requireNamespace("emplik", quietly = TRUE)) {
  stop("需要安装 R 包 'emplik'。", call. = FALSE)
}

#==============================================#
#                 辅助函数
#==============================================#
validate_scalar <- function(x, name, lower = -Inf, upper = Inf, integer = FALSE) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
      x < lower || x > upper || (integer && x != floor(x))) {
    stop(name, " 参数不合法。", call. = FALSE)
  }
  invisible(x)
}

# 与 Generate-Data.R 一致：取 emplik::el.test() 的 -2LLR；报错或非有限值视为失败。
compute_el_statistic <- function(x, mu = 0) {
  result <- function(statistic, status) list(statistic = statistic, status = status)
  if (length(x) < 2L) return(result(NA_real_, "insufficient_n"))
  if (length(mu) != 1L || !is.finite(mu) || any(!is.finite(x))) {
    return(result(NA_real_, "nonfinite_input"))
  }
  fit <- tryCatch(suppressWarnings(emplik::el.test(x, mu = mu)),
                  error = function(e) NULL)
  if (is.null(fit) || !is.finite(fit$`-2LLR`)) {
    return(result(NA_real_, "solver_fail"))
  }
  # 凸包状态只作诊断；与 A 一样接受包返回的有限统计量，不另行替换为 Inf。
  hull_infeasible <- (min(x) >= mu || max(x) <= mu) && !all(x == mu)
  result(max(0, as.numeric(fit$`-2LLR`)),
         if (hull_infeasible) "hull_infeasible" else "ok")
}

# SEL 边界尾概率。仅作为明确命名的可选 SEL-BC 方法。
compute_sel_bc_pvalue <- function(statistic, a2) {
  validate_scalar(a2, "a2", 0, 0.499999)
  if (is.na(statistic) || statistic < 0) return(NA_real_)
  if (is.infinite(statistic)) return(0)
  if (statistic == 0) return(1)
  base <- stats::pchisq(statistic, 1, lower.tail = FALSE)
  if (a2 == 0) return(base)
  gate_cut <- stats::qnorm(1 - a2)
  root_t <- sqrt(statistic)
  # 在 sqrt(t) 处分段，尾部解析积分，避免极小 p 值的无穷积分误差。
  if (root_t <= gate_cut) return((1 - a2) * base + a2)
  integral <- tryCatch(stats::integrate(function(z) {
    stats::pchisq(pmax(0, statistic - z^2), 1, lower.tail = FALSE) * stats::dnorm(z)
  }, lower = gate_cut, upper = root_t, rel.tol = 1e-9, abs.tol = 1e-13,
  subdivisions = 500L, stop.on.error = TRUE)$value, error = function(e) NA_real_)
  min(1, max(0, (1 - a2) * base + integral + stats::pnorm(root_t, lower.tail = FALSE)))
}

resolve_methods <- function(methods) {
  allowed <- c("CEL", "SEL", "ASEL", "t", "SEL-BC")
  if (!is.character(methods) || !length(methods) || anyNA(methods) ||
      anyDuplicated(methods) || any(!methods %in% allowed)) {
    stop("methods 使用 CEL、SEL、ASEL、t、SEL-BC 中不重复的名称。")
  }
  methods
}

# 同一组观测只计算所选方法的 BH 输入。
# CEL 在均值不超边界时取 1，否则取半个尾概率；t 使用单侧 t 尾概率。
# SEL/ASEL 的原始筛查分数沿用现有定义，其校准由模拟结果评估。
compute_group_scores <- function(x, epsilon0 = 0.25, theta0 = 0, a2 = 0.08,
                                 methods = c("CEL", "SEL", "ASEL", "t")) {
  validate_scalar(epsilon0, "epsilon0", 0)
  validate_scalar(theta0, "theta0")
  validate_scalar(a2, "a2", 0, 0.499999)
  methods <- resolve_methods(methods)
  out <- data.frame(method = methods, score = NA_real_, status = "ok", failed = FALSE,
                    hull_infeasible = FALSE, solver_fail = FALSE,
                    gate = NA, n = length(x), n_pairs = floor(length(x) / 2),
                    aux_zero_hull = NA, stringsAsFactors = FALSE)
  invalid <- if (length(x) < 2L) "insufficient_n" else if (any(!is.finite(x))) {
    "nonfinite_input"
  } else if (!is.finite(stats::sd(x)) || stats::sd(x) == 0) "zero_variance" else NULL
  if (!is.null(invalid)) {
    out$status <- invalid
    out$failed <- TRUE
    return(out)
  }
  boundary <- theta0 + epsilon0
  # 失败记为 NA；主循环丢弃整次模拟，包括其他方法已算出的结果。
  set_score <- function(method, value, fits = list(), failure = NULL) {
    if (!method %in% methods) return(invisible(NULL))
    index <- match(method, out$method)
    states <- vapply(fits, function(fit) fit$status, character(1))
    bad <- states[!states %in% c("ok", "hull_infeasible")]
    if (length(bad)) failure <- bad[1]
    if (is.null(failure) && (!is.finite(value) || value < 0 || value > 1)) {
      failure <- "nonfinite_score"
    }
    out$hull_infeasible[index] <<- any(states == "hull_infeasible")
    out$solver_fail[index] <<- any(states == "solver_fail")
    out$failed[index] <<- !is.null(failure)
    out$status[index] <<- if (!is.null(failure)) failure else if (any(states == "hull_infeasible")) {
      "hull_infeasible"
    } else "ok"
    out$score[index] <<- if (is.null(failure)) value else NA_real_
  }
  tail_chi <- function(t) stats::pchisq(t, df = 1, lower.tail = FALSE)
  if ("CEL" %in% methods) {
    if (mean(x) <= boundary) set_score("CEL", 1) else {
      raw <- compute_el_statistic(x, boundary)
      set_score("CEL", 0.5 * tail_chi(raw$statistic), list(raw))
    }
  }
  t_statistic <- (mean(x) - boundary) / (stats::sd(x) / sqrt(length(x)))
  set_score("t", stats::pt(t_statistic, df = length(x) - 1, lower.tail = FALSE))
  h <- floor(length(x) / 2)
  screening <- intersect(c("SEL", "ASEL", "SEL-BC"), methods)
  if (!length(screening)) return(out)
  if (h < 2L) {
    for (method in screening) set_score(method, NA_real_, failure = "insufficient_pairs")
    return(out)
  }
  # x 顺序在 DGP 内已预定；奇数时仅舍弃最后一个观测。
  # 配对差 D1 消去共同均值；配对和 D2 保留相对外部参考值的偏离。
  d1 <- x[seq_len(h)] - x[h + seq_len(h)]
  d2 <- x[seq_len(h)] + x[h + seq_len(h)] - 2 * theta0
  s2 <- stats::sd(d2)
  if (!is.finite(s2) || s2 == 0) {
    for (method in screening) set_score(method, NA_real_, failure = "zero_pair_variance")
    return(out)
  }
  # gate: sqrt(h)*(mean(D2)-2*epsilon0) > z_(1-a2)*sd(D2)。
  # SEL = EL(D1;0) + gate*EL(D2;2*epsilon0)；ASEL 将第二个约束改为 0。
  gate <- sqrt(h) * (mean(d2) - 2 * epsilon0) > stats::qnorm(1 - a2) * s2
  out$gate[out$method %in% screening] <- gate
  out$aux_zero_hull[out$method %in% screening] <-
    (min(d2) >= 0 || max(d2) <= 0) && !all(d2 == 0)
  baseline <- compute_el_statistic(d1, 0)
  fits_sel <- fits_asel <- list(baseline)
  ts <- ta <- baseline$statistic
  if (gate) {
    if (any(methods %in% c("SEL", "SEL-BC"))) {
      boundary_fit <- compute_el_statistic(d2, 2 * epsilon0)
      fits_sel <- c(fits_sel, list(boundary_fit))
      ts <- ts + boundary_fit$statistic
    }
    if ("ASEL" %in% methods) {
      zero_fit <- compute_el_statistic(d2, 0)
      fits_asel <- c(fits_asel, list(zero_fit))
      ta <- ta + zero_fit$statistic
    }
  }
  set_score("SEL", tail_chi(ts), fits_sel)
  set_score("ASEL", tail_chi(ta), fits_asel)
  if ("SEL-BC" %in% methods) set_score("SEL-BC", compute_sel_bc_pvalue(ts, a2), fits_sel)
  out
}

# 与 A 一样手写 BH 排序阈值；补充 BY 使用调低后的 q。
bh_reject <- function(p, q = 0.05, adjustment = "BH") {
  validate_scalar(q, "q", .Machine$double.eps, 1)
  if (!length(p) || any(!is.finite(p)) || any(p < 0 | p > 1)) {
    stop("BH 输入必须为非空的 [0,1] 有限向量。", call. = FALSE)
  }
  if (!adjustment %in% c("BH", "BY")) stop("adjustment 必须为 BH 或 BY。")
  m <- length(p)
  if (adjustment == "BY") q <- q / sum(1 / seq_len(m))
  ord <- order(p)
  k <- which(p[ord] <= seq_len(m) * q / m)
  if (length(k) == 0L) return(rep(FALSE, m))
  p <= p[ord][max(k)]
}

#==============================================================#
# 场景网格与真实参数；场景编号只作标签，不管理随机数。
#==============================================================#
build_figure_grid <- function(n, k, dgp, effect_sizes, epsilon0 = 0.25,
                              theta0 = 0, null_epsilons = seq(0, epsilon0, 0.05)) {
  validate_scalar(n, "n", 2, integer = TRUE)
  validate_scalar(k, "k", 2, integer = TRUE)
  if (k %% 2 != 0) stop("五面板的 pi1=0.5 要求 k 为偶数。")
  if (length(dgp) != 1L || is.na(dgp) ||
      !dgp %in% c("Normal", "t3", "Exponential")) {
    stop("dgp 必须是单个已知分布名称。")
  }
  if (anyDuplicated(effect_sizes)) stop("effect_sizes 不能重复。")
  for (delta in effect_sizes) validate_scalar(delta, "effect_size", 0)
  positive_effects <- sort(effect_sizes[effect_sizes > 0])
  if (!length(positive_effects)) stop("五面板需要至少一个正效应量。")

  validate_scalar(epsilon0, "epsilon0", 0)
  validate_scalar(theta0, "theta0")
  if (!length(null_epsilons) || anyDuplicated(null_epsilons)) stop("null_epsilons 必须非空且不重复。")
  for (epsilon in null_epsilons) validate_scalar(epsilon, "null_epsilon", 0, epsilon0)
  # 全原假设下所有组共享 epsilon；指数设计排除零值。
  if (dgp == "Exponential") null_epsilons <- null_epsilons[null_epsilons > 0]
  if (!length(null_epsilons)) stop("当前分布没有可用的全原假设 epsilon 网格。")
  if (dgp == "Exponential" && any(theta0 + c(null_epsilons, epsilon0) <= 0)) {
    stop("指数分布的条件均值必须为正。")
  }
  alternatives <- expand.grid(target_pi1 = c(0.5, 1), effect_size = positive_effects)
  alternatives$null_epsilon <- epsilon0
  alternatives$null_config <- "boundary"
  effects <- rbind(data.frame(target_pi1 = 0, effect_size = 0,
                              null_epsilon = sort(null_epsilons), null_config = "common"),
                   alternatives)
  grid <- data.frame(dgp = dgp, K = k, n_group_target = n / k,
                     effects, N = n,
                     block = "figure", sampling_design = "random",
                     scenario_index = seq_len(nrow(effects)))
  grid$scenario_id <- sprintf("K%s-figure-%03d", k, grid$scenario_index)
  # 固定执行顺序；新增场景会消耗随机数，因此后续样本与旧版本不同。
  grid <- grid[order(grid$target_pi1, grid$null_epsilon, grid$effect_size), , drop = FALSE]
  rownames(grid) <- NULL
  grid
}

build_group_truth <- function(scenario, epsilon0 = 0.25, theta0 = 0) {
  K <- scenario$K
  validate_scalar(K, "K", 1, integer = TRUE)
  validate_scalar(epsilon0, "epsilon0", 0)
  validate_scalar(scenario$target_pi1, "target_pi1", 0, 1)
  validate_scalar(scenario$effect_size, "effect_size", 0)
  if (scenario$sampling_design == "overlap") {
    if (K != 19L || scenario$dgp != "Normal" ||
        !scenario$target_pi1 %in% c(0, 0.5) || scenario$null_config != "boundary") {
      stop("overlap 配置固定为 10 个正态原子组和 9 个相邻并集。")
    }
    atoms <- rep(epsilon0, 10)
    if (scenario$target_pi1 > 0) atoms[6:10] <- epsilon0 + scenario$effect_size
    epsilon <- c(atoms, (atoms[1:9] + atoms[2:10]) / 2)
    group_names <- c(paste0("A", 1:10), paste0("A", 1:9, "+A", 2:10))
    probability <- c(rep(0.1, 10), rep(0.2, 9))
  } else {
    planned_m1 <- round(K * scenario$target_pi1)
    if (abs(planned_m1 - K * scenario$target_pi1) > 1e-9) {
      stop("K * target_pi1 必须是整数，以免静默改变设计。")
    }
    planned_m0 <- K - planned_m1
    if (scenario$null_config == "common") {
      validate_scalar(scenario$null_epsilon, "null_epsilon", 0, epsilon0)
      if (scenario$target_pi1 != 0 || scenario$effect_size != 0) {
        stop("common 配置只用于 pi1=0、delta=0 的全原假设扫描。")
      }
    }
    null_values <- switch(scenario$null_config, common = scenario$null_epsilon,
                          boundary = epsilon0, interior = 0.10,
                          mixed_null = c(0.05, 0.15, epsilon0), parity = 0,
                          stop("未知 null_config。"))
    if (any(null_values < 0 | null_values > epsilon0)) stop("null_config 超出容忍区间。")
    epsilon <- c(rep(null_values, length.out = planned_m0),
                 rep(epsilon0 + scenario$effect_size, planned_m1))
    group_names <- paste0("G", seq_len(K))
    probability <- rep(1 / K, K)
  }
  data.frame(group_id = seq_len(K), group_name = group_names, epsilon = epsilon,
             mean = theta0 + epsilon, null_true = epsilon <= epsilon0,
             probability = probability, expected_n = scenario$N * probability,
             stringsAsFactors = FALSE)
}

generate_group_data <- function(scenario, truth) {
  validate_scalar(scenario$N, "N", 1, integer = TRUE)
  K <- nrow(truth)
  if (!scenario$sampling_design %in% c("fixed", "random", "overlap")) stop("未知抽样设计。")
  if (scenario$sampling_design == "fixed") {
    validate_scalar(scenario$n_group_target, "n_group_target", 2, integer = TRUE)
    if (scenario$N != K * scenario$n_group_target) stop("fixed 需要 N=K*r。")
    atom <- rep(seq_len(K), each = scenario$n_group_target)
  } else {
    atom_count <- if (scenario$sampling_design == "overlap") 10L else K
    atom <- pmin(floor(stats::runif(scenario$N) * atom_count) + 1L, atom_count)
  }
  # ID / 分组顺序在生成 M 之前已确定，奇数组最后一个 ID 与观测大小无关。
  mu <- truth$mean[atom]
  if (scenario$dgp == "Exponential" && any(mu <= 0)) {
    stop("指数分布的条件均值必须为正。")
  }
  x <- switch(scenario$dgp,
              Normal = stats::rnorm(length(atom), mean = mu),
              t3 = stats::rt(length(atom), df = 3) + mu,
              Exponential = stats::rexp(length(atom), rate = 1 / mu),
              stop("未知 dgp。"))
  ids <- lapply(seq_len(if (scenario$sampling_design == "overlap") 10L else K),
                function(g) which(atom == g))
  if (scenario$sampling_design == "overlap") {
    ids <- c(ids, lapply(1:9, function(g) which(atom %in% c(g, g + 1L))))
  }
  list(observations = data.frame(id = seq_along(x), atom = atom, M = x),
       group_ids = ids, groups = lapply(ids, function(index) x[index]),
       group_n = lengths(ids))
}

#==============================================================#
# family 指标、MCSE 与汇总
#==============================================================#
compute_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
compute_mcse <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2L) NA_real_ else stats::sd(x) / sqrt(length(x))
}

# 每次重复在同一份数据上比较各方法：V 为错拒数，S 为正确拒绝数，R=V+S。
# FDP=V/max(R,1)；TPR=S/m1；没有真实备择时 TPR 为 NA。
compute_family_results <- function(scores, truth, group_n, replicate, qs, overlap = FALSE) {
  specs <- data.frame(base = unique(scores$method), adjustment = "BH", stringsAsFactors = FALSE)
  by_methods <- intersect(c("CEL", "t"), specs$base)
  if (overlap && length(by_methods)) {
    specs <- rbind(specs, data.frame(base = by_methods, adjustment = "BY"))
  }
  m0 <- sum(truth$null_true)
  m1 <- sum(!truth$null_true)
  K <- nrow(truth)
  rows <- vector("list", nrow(specs) * length(qs))
  index <- 0L
  for (i in seq_len(nrow(specs))) {
    d <- scores[scores$method == specs$base[i], , drop = FALSE]
    d <- d[match(truth$group_id, d$group_id), , drop = FALSE]
    for (q in qs) {
      rejected <- bh_reject(d$score, q, specs$adjustment[i])
      V <- sum(rejected & truth$null_true)
      S <- sum(rejected & !truth$null_true)
      R <- V + S
      index <- index + 1L
      rows[[index]] <- data.frame(
        replicate = replicate, method = paste(specs$base[i], specs$adjustment[i], sep = "-"), q = q,
        R = R, V = V, S = S, fdp = V / max(R, 1), tpr = if (m1 > 0) S / m1 else NA_real_,
        gate_rate_null = compute_mean(d$gate[truth$null_true]),
        gate_rate_alt = compute_mean(d$gate[!truth$null_true]),
        gate_missing_rate = if (all(is.na(d$gate)) && !specs$base[i] %in% c("SEL", "ASEL", "SEL-BC")) {
          NA_real_
        } else mean(is.na(d$gate)),
        hull_infeasible_rate = mean(d$hull_infeasible),
        aux_zero_hull_rate = compute_mean(d$aux_zero_hull),
        solver_fail_rate = mean(d$solver_fail), failed_group_rate = mean(d$failed),
        affected_family = any(d$failed), min_group_n = min(group_n),
        median_group_n = stats::median(group_n), max_group_n = max(group_n),
        min_pair_n = min(floor(group_n / 2)),
        null_tail_q_over_K = if (m0) mean(d$score[truth$null_true] <= q / K) else NA_real_,
        null_tail_q_over_2 = if (m0) mean(d$score[truth$null_true] <= q / 2) else NA_real_,
        null_tail_q = if (m0) mean(d$score[truth$null_true] <= q) else NA_real_)
      for (status in c("insufficient_n", "nonfinite_input", "zero_variance", "insufficient_pairs",
                       "zero_pair_variance", "nonfinite_score")) {
        rows[[index]][[paste0(status, "_rate")]] <- mean(d$status == status)
      }
    }
  }
  do.call(rbind, rows)
}

# 仅汇总有效重复：FDR=mean(FDP)，Power=mean(TPR)。
summarize_fdr_results <- function(family) {
  keys <- interaction(family$method, family$q, drop = TRUE)
  rows <- lapply(split(family, keys), function(d) {
    fdr <- mean(d$fdp)
    power <- compute_mean(d$tpr)
    out <- data.frame(method = d$method[1], q = d$q[1], nsim = nrow(d),
                      fdr = fdr, power = power,
                      mean_rejections = mean(d$R), mean_false = mean(d$V), mean_true = mean(d$S),
                      prob_any_rejection = mean(d$R > 0),
                      affected_family_rate = mean(d$affected_family),
                      min_group_n = min(d$min_group_n), min_pair_n = min(d$min_pair_n),
                      mean_median_group_n = mean(d$median_group_n))
    means <- c("gate_rate_null", "gate_rate_alt", "gate_missing_rate", "hull_infeasible_rate",
               "aux_zero_hull_rate", "solver_fail_rate", "failed_group_rate",
               "null_tail_q_over_K", "null_tail_q_over_2", "null_tail_q",
               "insufficient_n_rate", "nonfinite_input_rate", "zero_variance_rate",
               "insufficient_pairs_rate", "zero_pair_variance_rate", "nonfinite_score_rate")
    for (column in means) {
      out[[column]] <- compute_mean(d[[column]])
      out[[paste0(column, "_mcse")]] <- compute_mcse(d[[column]])
    }
    out
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#==============================================#
#                 数据保存
#==============================================#
# 方法 × 分布 × n × K × pi 决定文件名；N 是 Normal 的缩写，n 表示总样本量。
# BH 文件省略 -BH 后缀；补充 BY 方法保留 -BY，避免与 BH 文件重名。
build_method_summary_filename <- function(method, dgp, k, pi1, n) {
  codes <- c(Normal = "N", t3 = "T", Exponential = "E")
  allowed <- c("CEL", "SEL", "ASEL", "t", "SEL-BC", "CEL-BY", "t-BY")
  method <- sub("-BH$", "", method)
  if (length(method) != 1L || is.na(method) || !method %in% allowed) stop("未知方法名称。")
  if (length(dgp) != 1L || is.na(dgp) || !dgp %in% names(codes)) stop("未知分布。")
  validate_scalar(n, "n", 2, integer = TRUE)
  validate_scalar(k, "K", 1, integer = TRUE)
  validate_scalar(pi1, "pi1", 0, 1)
  paste0(method, "-", codes[[dgp]], "-n", format(n, scientific = FALSE, trim = TRUE),
         "-K", format(k, scientific = FALSE, trim = TRUE),
         "-pi", format(pi1, digits = 12, scientific = FALSE, trim = TRUE), ".csv")
}

# 每张表合并该方法、分布、n、K、pi 下的所有网格点；保留运行元数据。
# null_epsilon 是原假设组的偏离，effect_size 是备择组超过阈值的 delta。
# 同名文件直接覆盖，不另外保存逐效应量子表、总大表或逐次数据。
save_simulation_data <- function(results, output_dir = ".") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- character()
  groups <- unique(results[, c("method", "dgp", "N", "K", "target_pi1"), drop = FALSE])
  for (i in seq_len(nrow(groups))) {
    g <- groups[i, ]
    table_data <- results[results$method == g$method & results$dgp == g$dgp &
                            results$N == g$N & results$K == g$K &
                            results$target_pi1 == g$target_pi1, , drop = FALSE]
    table_data <- table_data[order(table_data$N, table_data$null_epsilon,
                                  table_data$effect_size, table_data$q), , drop = FALSE]
    file_path <- file.path(output_dir, build_method_summary_filename(g$method, g$dgp, g$K, g$target_pi1, n = g$N))
    utils::write.csv(table_data, file_path, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
    paths <- c(paths, file_path)
    message("数据已保存：", file_path)
  }
  invisible(paths)
}

#==============================================#
#                 结果作图
#==============================================#
get_method_styles <- function() {
  # 此处的行顺序同时决定曲线、标题和图例的顺序。
  # CEL 红色，SEL 橙色，ASEL 蓝色，t 绿色。
  data.frame(
    method = c("CEL-BH", "SEL-BH", "ASEL-BH", "t-BH",
               "SEL-BC-BH", "CEL-BY", "t-BY"),
    color = c("red", "orange", "blue", "green3",
              "#0072B2", "red", "green3"),
    lty = c(5, 2, 1, 6, 4, 3, 3),
    pch = c(2, 17, 19, 3, 15, 1, 4),
    stringsAsFactors = FALSE
  )
}

split_plot_pages <- function(summary) {
  required <- c("block", "dgp", "K", "sampling_design", "null_config", "q", "epsilon0", "theta0", "a2")
  missing <- setdiff(c(required, "method", "fdr", "power", "m1", "effect_size", "n_group_target"), names(summary))
  if (length(missing)) stop("缺少绘图字段：", paste(missing, collapse = ", "))
  key <- paste(apply(summary[, required, drop = FALSE], 1, paste, collapse = "|"),
               summary$m1 == 0, sep = "|")
  split(summary, factor(key, levels = unique(key)))
}

draw_diagnostic_fdr_power <- function(data, certification_only = FALSE) {
  styles <- get_method_styles()
  keep <- if (certification_only) c("CEL-BH", "t-BH", "CEL-BY", "t-BY") else styles$method
  data <- data[data$method %in% keep, , drop = FALSE]
  if (!nrow(data)) return(invisible(NULL))
  styles <- styles[styles$method %in% data$method, , drop = FALSE]
  null_page <- all(data$m1 == 0)
  budget_page <- all(data$sampling_design %in% c("random", "overlap"))
  panels <- if (null_page || budget_page) data.frame(target_pi1 = unique(data$target_pi1)) else {
    unique(data[, c("target_pi1", "n_group_target"), drop = FALSE])
  }
  panels <- panels[order(panels$target_pi1), , drop = FALSE]
  metrics <- if (null_page) {
    if (certification_only) "fdr" else c("fdr", "gate_rate_null")
  } else c("fdr", "power")
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  legend_lines <- nrow(styles) + 1
  graphics::par(mfrow = c(nrow(panels), length(metrics)), oma = c(legend_lines, 0.3, 4.5, 0.3),
                mar = c(3.5, 4.1, 2.2, 0.9), mgp = c(2.5, 0.7, 0), las = 1)
  for (i in seq_len(nrow(panels))) {
    panel <- data[data$target_pi1 == panels$target_pi1[i], , drop = FALSE]
    if (!null_page && !budget_page) panel <- panel[panel$n_group_target == panels$n_group_target[i], , drop = FALSE]
    x_name <- if (budget_page) "N" else if (!null_page) "effect_size" else "n_group_target"
    x_label <- if (x_name == "N") "Total sample size (N)" else if (x_name == "effect_size") "Effect above tolerance (delta)" else "Sample size per group (r)"
    for (metric in metrics) {
      pd <- panel
      if (metric == "gate_rate_null") pd <- pd[pd$method %in% c("SEL-BH", "ASEL-BH", "SEL-BC-BH"), , drop = FALSE]
      xlim <- range(pd[[x_name]])
      if (diff(xlim) == 0) xlim <- xlim + c(-1, 1) * max(0.01, abs(xlim[1]) * 0.05)
      ylim <- c(0, 1)
      if (certification_only && metric == "fdr") {
        ymax <- max(c(0.10, 2 * pd$q, pd$fdr), na.rm = TRUE)
        ylim <- c(0, min(1, ymax * 1.05))
      }
      if (metric == "gate_rate_null") {
        ymax <- max(c(0.15, 2 * pd$a2, pd$gate_rate_null + 1.96 * pd$gate_rate_null_mcse), na.rm = TRUE)
        ylim <- c(0, min(1, ymax * 1.05))
      }
      title <- if (null_page) "All hypotheses are null" else if (budget_page) {
        sprintf("Realized pi1 = %.3f; delta = %g", panel$realized_pi1[1], panel$effect_size[1])
      } else sprintf("pi1 = %g; expected r = %g", panels$target_pi1[i], panels$n_group_target[i])
      y_label <- switch(metric, fdr = "False discovery rate", power = "Power", gate_rate_null = "Null gate activation")
      graphics::plot(NA, xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i", xlab = x_label,
                     ylab = y_label, main = title, cex.main = 0.95)
      graphics::grid(col = "gray90")
      graphics::abline(h = if (metric == "gate_rate_null") pd$a2[1] else if (metric == "fdr") pd$q[1] else NA_real_,
                       col = "gray40", lty = 2)
      for (j in seq_len(nrow(styles))) {
        curve <- pd[pd$method == styles$method[j], , drop = FALSE]
        if (!nrow(curve)) next
        curve <- curve[order(curve[[x_name]]), , drop = FALSE]
        if (anyDuplicated(curve[[x_name]])) stop("同面板同方法出现重复 x；请先筛选批次或合并 family 后重新汇总。")
        y <- curve[[metric]]
        if (all(is.na(y))) next
        if (metric == "gate_rate_null") {
          lower <- pmax(0, y - 1.96 * curve$gate_rate_null_mcse)
          upper <- pmin(1, y + 1.96 * curve$gate_rate_null_mcse)
          ok <- is.finite(lower) & is.finite(upper) & upper > lower
          if (any(ok)) graphics::arrows(curve[[x_name]][ok], lower[ok], curve[[x_name]][ok], upper[ok],
                                        angle = 90, code = 3, length = 0.025, col = styles$color[j])
        }
        graphics::lines(curve[[x_name]], y, col = styles$color[j], lty = styles$lty[j], lwd = 1.7)
        graphics::points(curve[[x_name]], y, col = styles$color[j], pch = styles$pch[j], cex = 0.7)
      }
      if (null_page && !certification_only && metric == "fdr" &&
          data$null_config[1] == "boundary" && data$sampling_design[1] == "fixed" && "ASEL-BH" %in% data$method) {
        lower_bound <- 1 - (1 - data$a2[1])^data$K[1]
        asel_color <- styles$color[match("ASEL-BH", styles$method)]
        graphics::abline(h = lower_bound, col = asel_color, lty = 3)
        graphics::mtext(sprintf("ASEL asymptotic lower bound: %.3f", lower_bound), side = 3, line = 0.1, cex = 0.62, col = asel_color)
      }
    }
  }
  header <- sprintf("%s | %s | K=%d | %s | q=%g, a2=%g", data$dgp[1], data$sampling_design[1],
                    data$K[1], data$null_config[1], data$q[1], data$a2[1])
  method_title <- paste(sub("-BH$", "", styles$method), collapse = "/")
  graphics::mtext(method_title, side = 3, outer = TRUE, line = 3.2, font = 2, cex = 1)
  graphics::mtext(header, side = 3, outer = TRUE, line = 2, font = 2, cex = 0.95)
  graphics::mtext(sprintf("Tolerance = %g; MC repetitions = %s%s", data$epsilon0[1],
                          paste(sort(unique(data$nsim)), collapse = "/"),
                          if (certification_only) "; certification detail" else "; SEL/ASEL are screening scores"),
                  side = 3, outer = TRUE, line = 0.8, cex = 0.75)
  # 全页图例独占底部区域，避免遮挡曲线或越出最后一个面板。
  legend_height <- (0.2 * legend_lines) / graphics::par("din")[2]
  graphics::par(oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), cex = 1,
                fig = c(0, 1, 0, legend_height), new = TRUE)
  graphics::plot.new()
  graphics::legend("center", styles$method, col = styles$color, lty = styles$lty, pch = styles$pch,
                   ncol = 1L, bty = "n", cex = 0.8, xpd = NA)
  invisible(NULL)
}

save_diagnostic_fdr_power <- function(data, file_path, certification_only = FALSE) {
  panels <- if (all(data$m1 == 0) || all(data$sampling_design %in% c("random", "overlap"))) {
    length(unique(data$target_pi1))
  } else nrow(unique(data[, c("target_pi1", "n_group_target")]))
  grDevices::postscript(file_path, width = 11, height = max(5, 2.8 * panels + 1.6),
                       paper = "special", horizontal = FALSE, onefile = FALSE, family = "Helvetica")
  on.exit(grDevices::dev.off(), add = TRUE)
  draw_diagnostic_fdr_power(data, certification_only)
  invisible(file_path)
}

run_diagnostic_fdr_plots <- function(summary, output_stem, save_plot = TRUE) {
  pages <- split_plot_pages(summary)
  paths <- character()
  for (i in seq_along(pages)) {
    path <- paste0(output_stem, "-page", sprintf("%02d", i), "-all.eps")
    if (save_plot) {
      save_diagnostic_fdr_power(pages[[i]], path)
      paths <- c(paths, path)
    }
    if (interactive()) draw_diagnostic_fdr_power(pages[[i]])
  }
  invisible(paths)
}

#==============================================================#
# 主实验布局：一行五列 (0:FDR, 0.5:FDR/Power, 1:FDR/Power)。
# FDR 面板的虚线标出 BH 水平，颜色区分方法，线型区分总样本量。
# 全原假设面板横轴为 epsilon；其余面板横轴为正效应 delta。
#==============================================================#
prepare_reference_plot_data <- function(data) {
  required <- c("method", "N", "K", "dgp", "target_pi1", "effect_size",
                "null_epsilon", "epsilon0", "theta0", "fdr", "power", "m1")
  if (any(!required %in% names(data))) stop("五面板缺少字段：", paste(setdiff(required, names(data)), collapse = ", "))
  if (any(!data$target_pi1 %in% c(0, 0.5, 1))) stop("五面板只接受 target_pi1=0,0.5,1。")
  # 每个展示点保留其实际模拟场景编号，不复制或补造估计值。
  data$source_scenario_id <- if ("scenario_id" %in% names(data)) data$scenario_id else NA_character_
  data$plot_reused <- FALSE
  rows <- list()
  for (method in unique(data$method)) for (n in sort(unique(data$N))) {
    d <- data[data$method == method & data$N == n, , drop = FALSE]
    null_curve <- d[d$target_pi1 == 0, , drop = FALSE]
    if (!nrow(null_curve) || any(null_curve$m1 != 0) ||
        any(!is.finite(null_curve$null_epsilon)) ||
        any(null_curve$null_epsilon < 0 | null_curve$null_epsilon > null_curve$epsilon0) ||
        any(null_curve$effect_size != 0) || anyDuplicated(null_curve$null_epsilon)) {
      stop("每种 method/N 的 pi1=0 各 epsilon 需要唯一的全原假设模拟结果；请使用新设计重新生成数据。")
    }
    if (any(null_curve$dgp == "Exponential" &
            (null_curve$null_epsilon <= 0 | null_curve$theta0 + null_curve$null_epsilon <= 0))) {
      stop("指数设计的全原假设 epsilon 和条件均值必须为正。")
    }
    positive <- d[d$target_pi1 > 0 & d$effect_size > 0, , drop = FALSE]
    if (!all(c(0.5, 1) %in% positive$target_pi1)) stop("五面板缺少 pi1=0.5 或 1 的正效应结果；请运行完整 figure 配置。")
    if (any(positive$m1 <= 0)) stop("正效应的预设备择配置没有真备择。")
    if (anyDuplicated(paste(positive$target_pi1, positive$effect_size))) stop("重复的五面板数据，请先合并 family 后汇总。")
    null_curve <- null_curve[order(null_curve$null_epsilon), , drop = FALSE]
    null_curve$power <- NA_real_
    # 全原假设各点使用自己的 FDR；备择面板仅使用正效应结果。
    rows[[length(rows) + 1L]] <- rbind(null_curve, positive)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

split_reference_plot_pages <- function(summary) {
  # 同一张图包含 common 全原假设扫描和 boundary 备择配置，不能按 null_config 拆页。
  columns <- c("block", "dgp", "K", "sampling_design", "q", "epsilon0", "theta0", "a2")
  if (any(!columns %in% names(summary))) stop("缺少 SEL 汇总元数据；请使用 SEL-FDR 生成的 summary.csv。")
  key <- apply(summary[, columns, drop = FALSE], 1, paste, collapse = "|")
  split(summary, factor(key, levels = unique(key)))
}

draw_reference_fdr_power <- function(data, certification_only = FALSE) {
  styles <- get_method_styles()
  if (certification_only) {
    styles <- styles[styles$method %in% c("CEL-BH", "t-BH"), , drop = FALSE]
    styles$pch <- c(16, 17)
  }
  styles <- styles[styles$method %in% data$method, , drop = FALSE]
  data <- prepare_reference_plot_data(data[data$method %in% styles$method, , drop = FALSE])
  panel_specs <- data.frame(target_pi1 = c(0, 0.5, 0.5, 1, 1),
                            metric = c("fdr", "fdr", "power", "fdr", "power"))
  ns <- sort(unique(data$N))
  curves <- expand.grid(method = styles$method, N = ns, stringsAsFactors = FALSE)
  # 有多个总样本量时仍按方法分组，保证图例从上到下按 CEL/SEL/ASEL/t 排列。
  curves <- curves[order(match(curves$method, styles$method), curves$N), , drop = FALSE]
  curves$color <- styles$color[match(curves$method, styles$method)]
  curves$pch <- styles$pch[match(curves$method, styles$method)]
  curves$lty <- rep(1:6, length.out = length(ns))[match(curves$N, ns)]
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  graphics::par(mfrow = c(1, 5), oma = c(0, 0, 1.8, 0),
                mar = c(2.8, 3.2, 1.7, 0.75), mgp = c(1.8, 0.6, 0), las = 1,
                cex = 0.95, cex.main = 1.05, cex.lab = 1, cex.axis = 0.9)
  for (i in seq_len(nrow(panel_specs))) {
    pi1 <- panel_specs$target_pi1[i]
    metric <- panel_specs$metric[i]
    panel <- data[data$target_pi1 == pi1, , drop = FALSE]
    x_name <- if (pi1 == 0) "null_epsilon" else "effect_size"
    x_minor <- sort(unique(panel[[x_name]]))
    x_range <- if (pi1 == 0) c(0, data$epsilon0[1]) else c(0, max(x_minor))
    if (diff(x_range) == 0) x_range <- x_range + c(-0.01, 0.01)
    x_major <- if (pi1 == 0) x_minor else seq(0, x_range[2], by = 0.1)
    graphics::plot(NA, xlim = x_range, ylim = c(0, 1), xaxs = "i", yaxs = "i",
                   xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                   main = sprintf("target pi1 = %.2f", pi1))
    graphics::abline(v = x_major, h = seq(0, 1, 0.2), col = "gray88", lty = 1)
    graphics::axis(1, at = x_minor, labels = FALSE, lwd = 1.5, tcl = 0.35)
    graphics::axis(1, at = x_major, labels = sprintf(if (pi1 == 0) "%.2f" else "%.1f", x_major),
                   font = 2, lwd = 1.5, tcl = 0.35)
    graphics::axis(2, at = seq(0, 1, 0.2), labels = sprintf("%.1f", seq(0, 1, 0.2)), font = 2, lwd = 1.5, tcl = 0.35)
    graphics::axis(3, labels = FALSE, lwd = 1.5, tcl = 0)
    graphics::axis(4, labels = FALSE, lwd = 1.5, tcl = 0)
    graphics::title(xlab = if (pi1 == 0) expression(epsilon) else expression(delta), line = 1.6, font.lab = 2)
    graphics::title(ylab = if (metric == "fdr") "FDR" else "Power", line = 2.1, font.lab = 2)
    if (metric == "fdr") graphics::abline(h = data$q[1], col = "gray35", lty = 2, lwd = 1.5)
    for (j in seq_len(nrow(curves))) {
      curve <- panel[panel$method == curves$method[j] & panel$N == curves$N[j], , drop = FALSE]
      curve <- curve[order(curve[[x_name]]), , drop = FALSE]
      graphics::lines(curve[[x_name]], curve[[metric]], col = curves$color[j], lty = curves$lty[j], lwd = 2)
      graphics::points(curve[[x_name]], curve[[metric]], col = curves$color[j], pch = curves$pch[j], cex = 0.65)
    }
    if (i == 1L) {
      graphics::legend("topright", legend = curves$method, col = curves$color,
                       lty = curves$lty, pch = curves$pch, lwd = 1.8,
                       bg = "white", cex = if (nrow(curves) > 6) 0.60 else 0.70,
                       inset = 0.005, x.intersp = 0.7, y.intersp = 0.9)
    }
  }
  title <- if (certification_only) "CEL-BH vs t-BH" else {
    paste0(paste(sub("-BH$", "", styles$method), collapse = "/"), " + BH")
  }
  graphics::mtext(sprintf("%s (%s, K = %d, N = %s)", title, data$dgp[1], data$K[1],
                          paste(ns, collapse = "/")),
                  side = 3, outer = TRUE, line = 0.3, font = 2, cex = 1.35)
  invisible(data)
}

resolve_plot_layout <- function(data, layout = "auto") {
  layout <- match.arg(layout, c("auto", "reference", "diagnostic"))
  if (layout != "auto") return(layout)
  if ("block" %in% names(data) && all(startsWith(data$block, "figure"))) "reference" else "diagnostic"
}

draw_fdr_power <- function(data, certification_only = FALSE, layout = "auto") {
  if (resolve_plot_layout(data, layout) == "reference") {
    draw_reference_fdr_power(data, certification_only)
  } else draw_diagnostic_fdr_power(data, certification_only)
}

save_fdr_power <- function(data, file_path, certification_only = FALSE, layout = "auto") {
  if (resolve_plot_layout(data, layout) == "diagnostic") {
    return(save_diagnostic_fdr_power(data, file_path, certification_only))
  }
  prepare_reference_plot_data(data)
  grDevices::postscript(file_path, width = 18, height = 6, paper = "special",
                       horizontal = FALSE, onefile = FALSE, family = "Helvetica")
  on.exit(grDevices::dev.off(), add = TRUE)
  draw_reference_fdr_power(data, certification_only)
  invisible(file_path)
}

run_fdr_plots <- function(summary, output_stem, save_plot = TRUE, layout = "auto") {
  if (!save_plot && !interactive()) return(invisible(character()))
  if (resolve_plot_layout(summary, layout) == "diagnostic") {
    return(run_diagnostic_fdr_plots(summary, output_stem, save_plot))
  }
  pages <- split_reference_plot_pages(summary)
  # 部分场景运行时只保存统计结果，不把缺失的 pi1 面板画成 0。
  complete <- vapply(pages, function(d) {
    all(c(0, 0.5, 1) %in% d$target_pi1) && any(d$effect_size > 0)
  }, logical(1))
  if (any(!complete)) message("跳过 ", sum(!complete), " 个缺少 pi1=0/0.5/1 的五面板页面；结果 CSV 仍保留。")
  paths <- character()
  # 每个分布与组数组合输出一张图，包含本次所选方法。
  for (i in which(complete)) {
    d <- pages[[i]]
    codes <- c(Normal = "N", t3 = "T", Exponential = "E", Gamma = "G")
    sample_sizes <- paste(format(sort(unique(d$N)), scientific = FALSE, trim = TRUE), collapse = "-")
    path <- file.path(dirname(output_stem),
                      paste0("all-", codes[[d$dgp[1]]], "-n", sample_sizes, "-K", d$K[1], ".eps"))
    if (save_plot) {
      save_fdr_power(d, path, layout = "reference")
      paths <- c(paths, path)
    }
    if (interactive()) draw_reference_fdr_power(d)
  }
  invisible(paths)
}

#==============================================#
#                 模拟主函数
#==============================================#
# 单个场景：生成真实标签 → 重复生成数据 → 各组检验 → BH → 返回场景汇总表。
run_fdr_scenario <- function(scenario, nsim = 20L, qs = 0.05, a2 = 0.08,
                             epsilon0 = 0.25, theta0 = 0,
                             progress = TRUE, methods = c("CEL", "SEL", "ASEL", "t")) {
  methods <- resolve_methods(methods)
  validate_scalar(nsim, "nsim", 2, integer = TRUE)
  validate_scalar(a2, "a2", 0, 0.499999)
  if (!length(qs) || anyDuplicated(qs)) stop("qs 必须为非空且无重复。")
  for (q in qs) validate_scalar(q, "q", .Machine$double.eps, 1)
  truth <- build_group_truth(scenario, epsilon0, theta0)
  started <- proc.time()[["elapsed"]]
  rows <- vector("list", nsim)
  progress_every <- max(1L, floor(nsim / 10L))
  completed <- 0L
  attempted <- 0L
  # 与 A 一致：任意组或任一所选方法失败，整次重抽，直到获得 nsim 次有效重复。
  # 随机数直接接着外部 set.seed() 的序列消耗，不在场景内重设种子。
  while (completed < nsim) {
    attempted <- attempted + 1L
    data <- generate_group_data(scenario, truth)
    group_scores <- vector("list", nrow(truth))
    ok <- TRUE
    for (j in seq_len(nrow(truth))) {
      d <- compute_group_scores(data$groups[[j]], epsilon0, theta0, a2, methods = methods)
      if (any(d$failed) || any(!is.finite(d$score))) {
        ok <- FALSE
        break
      }
      d$group_id <- j
      group_scores[[j]] <- d
    }
    if (!ok) next

    completed <- completed + 1L
    scores <- do.call(rbind, group_scores)
    rows[[completed]] <- compute_family_results(scores, truth, data$group_n,
      completed, qs, scenario$sampling_design == "overlap")
    if (progress && completed %% progress_every == 0L) {
      message("  ", completed, "/", nsim, " valid families; discarded=", attempted - completed)
    }
  }
  family <- do.call(rbind, rows)
  m1 <- sum(!truth$null_true)
  summary <- summarize_fdr_results(family)
  meta <- scenario[, setdiff(names(scenario), c("q", "method")), drop = FALSE]
  summary <- cbind(meta[rep(1L, nrow(summary)), , drop = FALSE], summary)
  summary$m0 <- nrow(truth) - m1
  summary$m1 <- m1
  summary$realized_pi1 <- m1 / nrow(truth)
  summary$theta_mode <- "known_external"
  summary$theta0 <- theta0
  summary$epsilon0 <- epsilon0
  summary$a2 <- a2
  summary$replicate_start <- 1L
  summary$replicate_end <- nsim
  summary$attempted_families <- attempted
  summary$discarded_families <- attempted - completed
  summary$elapsed_seconds <- proc.time()[["elapsed"]] - started
  rownames(summary) <- NULL
  summary
}

# 一次固定一个 n、k、dgp，仅遍历效应量及五面板所需的原假设/备择场景。
# 计算完成后统一保存、打印和作图。
# seed 仅写入结果元数据；真正的 set.seed() 只在末尾执行区调用一次。
run_simulation <- function(n, k, dgp, effect_sizes, nsim, seed,
                           qs = 0.05, a2 = 0.08,
                           epsilon0 = 0.25, theta0 = 0,
                           save_data = FALSE, save_plot = FALSE,
                           plot_dir = ".", progress = TRUE,
                           methods = c("CEL", "SEL", "ASEL", "t"),
                           null_epsilons = seq(0, epsilon0, 0.05)) {
  methods <- resolve_methods(methods)
  validate_scalar(seed, "seed", 0, .Machine$integer.max, integer = TRUE)
  validate_scalar(epsilon0, "epsilon0", 0)
  validate_scalar(theta0, "theta0")

  # ---- 固定 n、k、dgp，构建五面板场景 ----
  grid <- build_figure_grid(n, k, dgp, effect_sizes, epsilon0, theta0, null_epsilons)
  if (save_data || save_plot) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(plot_dir)) stop("无法创建输出目录：", plot_dir)
  }

  # ---- 所选方法共用有效重复；改变方法可能改变重抽次数及后续随机样本 ----
  all_rows <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    s <- grid[i, , drop = FALSE]
    axis_name <- if (s$target_pi1 == 0) "epsilon" else "delta"
    axis_value <- if (s$target_pi1 == 0) s$null_epsilon else s$effect_size
    message(sprintf("[%d/%d] %s, n=%d, k=%d, pi1=%.2f, %s=%.2f",
                    i, nrow(grid), s$dgp, s$N, s$K, s$target_pi1, axis_name, axis_value))
    all_rows[[i]] <- run_fdr_scenario(
      scenario = s, nsim = nsim, qs = qs, a2 = a2,
      epsilon0 = epsilon0, theta0 = theta0,
      progress = progress, methods = methods
    )
  }
  result <- do.call(rbind, all_rows)
  rownames(result) <- NULL
  result$seed <- seed
  result$code_version <- "1.8.0-null-epsilon-sweep"
  result$run_label <- "figure"
  result$include_bc <- "SEL-BC" %in% methods
  result$R_version <- R.version.string
  result$rng <- paste(c(RNGkind(), "sequential; one external set.seed"), collapse = "/")
  result$el_engine <- "emplik::el.test"
  result$emplik_version <- as.character(utils::packageVersion("emplik"))
  result$replicate_policy <- "resample_whole_family_until_nsim_valid"

  # ---- 保存、打印并作图 ----
  if (save_data) save_simulation_data(result, output_dir = plot_dir)
  cat(sprintf("\n k = %d, epsilon0 = %g, theta0 = %g, a2 = %g, nsim = %d\n",
              k, epsilon0, theta0, a2, nsim))
  # 控制台区分全原假设 epsilon 和备择 delta；全原假设的 delta 固定为 0。
  for (method in unique(result$method)) {
    table_data <- result[result$method == method,
                         c("method", "target_pi1", "null_epsilon", "effect_size", "fdr", "power"), drop = FALSE]
    table_data <- table_data[order(table_data$target_pi1, table_data$null_epsilon,
                                  table_data$effect_size), , drop = FALSE]
    table_data$method <- sub("-BH$", "", table_data$method)
    cat("\n", table_data$method[1], "\n", sep = "")
    print(table_data, row.names = FALSE)
  }
  if (save_plot || interactive()) {
    output_stem <- file.path(plot_dir, paste0("SEL-FDR-K", k))
    paths <- run_fdr_plots(result, output_stem = output_stem, save_plot = save_plot)
    if (save_plot) message("图片已保存：", paste(paths, collapse = ", "))
  }
  result
}

#==============================================#
#                 参数设置
#==============================================#
nsim <- 100                          # 每个场景的有效重复次数；失败整次重抽
n <- 1000                           # 单个总样本量；随机分组后每组人数约为 n/k
ks <- c(26)                      # 最外层依次运行的组数
dgps <- c("Normal","Exponential")   # 单个分布，可改为"Normal"  "t3"、"Exponential"
effect_sizes <- seq(0, 0.6, 0.05)     # 备择均值为 theta0 + epsilon0 + delta
qs <- 0.05                          # BH 水平
a2 <- 0.08                          # SEL/ASEL 门控参数，与 BH 水平分开
epsilon0 <- 0.25                     # 容忍阈值
null_epsilons <- seq(0, epsilon0, 0.05) # pi1=0 时所有组共享 epsilon；指数设计自动排除 0
theta0 <- 0                          # 已知外部参考值
methods <- c("CEL", "SEL", "ASEL", "t") # 可设为 "CEL"，或加入 "SEL-BC"
save_data <- T                    # TRUE：按 方法-分布-n-K-pi 保存，合并所有效应量
save_plot <- F                    # TRUE：额外保存 EPS；交互环境仍正常显示
plot_dir <- "."                     # 当前工作目录 getwd()，不是自动定位脚本目录
seed <- 2

#==============================================#
#                 执行模拟
#==============================================#
# 只在这里设置一次种子；随后依次消耗 k、场景和重复所需的随机数。
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(seed)
for (dgp in dgps) {
  for (k in ks) {
    run_simulation(
      n = n, k = k, dgp = dgp, effect_sizes = effect_sizes, nsim = nsim, seed = seed,
      qs = qs, a2 = a2, epsilon0 = epsilon0, theta0 = theta0,
      null_epsilons = null_epsilons,
      methods = methods, save_data = save_data,
      save_plot = save_plot, plot_dir = plot_dir
    )
  }
}

# ---- 小规模检查（保持注释，不参与正式执行）----
# 临时副本中设置 nsim=2、effect_sizes=c(0,0.1)、plot_dir=tempdir() 后 Source。
# 例如 CEL-N-n2000-K30-pi1.csv 包含 CEL/Normal/n2000/K30/pi1=1 的全部效应量。
# 汇总不能恢复逐次模拟记录；methods="CEL" 时只计算、保存和绘制 CEL。
