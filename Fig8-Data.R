# =============================================== #
#   Compas Fig7 Audit Data (Female, EL/CEL/SEL/ASEL) #
# =============================================== #
rm(list = ls())
source("Functions/EL_CI.R")
source("Functions/CEL_CI.R")
source("Functions/SEL_CI.R")
source("Functions/ASEL_CI.R")
source("Functions/audit_run_ci.R")

library(dplyr)
library(readr)

test <- "test2"
ci_methods <- c("ELCI","CELCI","SELCI","ASELCI") # "ELCI","CELCI","SELCI","ASELCI"
alpha <- 0.9
alpha2 <- 0.08
dataname <- "Female-Race-Age"
save_data <- T  # FALSE 时不写 CSV

epsilon0 <- 0
epsilon_tol <- 0.05
epsilon_lo <- -epsilon_tol
epsilon_hi <- epsilon_tol
test <- normalize_audit_test(test)

# ---------------------获得数据------------------- #
compas_data <- read_csv(file("data/compas_data.csv"), show_col_types = FALSE)
table(compas_data$sex)
table(compas_data$age_cat)
table(compas_data$race)

# 一、数据处理
identical(compas_data$decile_score...12, compas_data$decile_score...40)
identical(compas_data$priors_count...15, compas_data$priors_count...49)
compas_data <- compas_data %>%
  rename(decile_score = decile_score...12) %>%
  select(-decile_score...40)  # 删除重复列
compas_data <- compas_data %>%
  rename(priors_count = priors_count...15) %>%
  select(-priors_count...49)  # 删除重复列

# 二、数据筛选
clean_data <- compas_data %>%
  select(sex, age_cat, race, decile_score, two_year_recid) %>%
  filter(
    !is.na(decile_score),
    !is.na(two_year_recid)
  ) %>%
  rename(Y = two_year_recid) %>%
  mutate(Y_h = ifelse(decile_score >= 5, 1, 0)) %>%
  select(-decile_score)

# 生成 L 列（Y == Y_h 时 L=0，否则 L=1）
clean_data <- clean_data %>%
  mutate(L = ifelse(Y == Y_h, 0, 1))
table(clean_data$L)

# 选出目标-女性
Female_data <- clean_data %>%
  filter(sex == "Female")
postive_target_data <- Female_data %>%
  filter(Y_h == 1)

# ---------------------基准线-------------------- #
# 全种族平均值基准线（阳性预测子集）
postive_data <- clean_data %>%
  filter(Y_h == 1)
table(postive_data$Y)
theta_P <- mean(postive_data$Y)
cat("theta_P =", theta_P, "\n")

fig6_min_n <- 8L

# 计算前预检：Fig6 要求 n >= 8，并与 run_ci 其余跳过规则一致
audit_group_check <- function(data, ci_method, target = theta_P) {
  y <- as.numeric(data$Y)
  n <- length(y)

  if (n < fig6_min_n) {
    return(list(
      ok = FALSE, n = n,
      reason = sprintf("n=%d < %d（最小样本量）", n, fig6_min_n)
    ))
  }

  min_n <- ci_min_n(ci_method)
  if (n < min_n) {
    return(list(
      ok = FALSE, n = n,
      reason = sprintf("n=%d < %d（%s 最小样本量）", n, min_n, ci_method)
    ))
  }

  if (ci_method %in% c("ELCI", "CELCI")) {
    x <- y - target
    if (!is.finite(sd(x)) || sd(x) == 0) {
      y_vals <- paste(sort(unique(y)), collapse = ",")
      return(list(
        ok = FALSE, n = n,
        reason = sprintf("n=%d 时 disparity 无变异（Y=%s）", n, y_vals)
      ))
    }
  }

  list(ok = TRUE, n = n, reason = NA_character_)
}

filter_audit_rows <- function(rows, ci_method, target = theta_P) {
  kept <- list()
  skipped <- list()
  for (row in rows) {
    chk <- audit_group_check(row[[2]], ci_method, target)
    if (isTRUE(chk$ok)) {
      kept <- c(kept, list(row))
    } else {
      skipped <- c(skipped, list(list(
        group = row[[1]], n = chk$n, reason = chk$reason
      )))
    }
  }
  list(kept = kept, skipped = skipped)
}

print_skipped_audit_groups <- function(skipped, ci_method) {
  if (length(skipped) == 0L) {
    return(invisible(NULL))
  }
  cat("\n----- 跳过（不参与", ci_method, "计算）-----\n", sep = "")
  for (s in skipped) {
    cat(sprintf("  [%s] %s\n", s$group, s$reason))
  }
}

audit_rows <- list(
  list("F", postive_target_data),
  list("< 25", filter(postive_target_data, age_cat == "Less than 25")),
  list("25 - 45", filter(postive_target_data, age_cat == "25 - 45")),
  list("> 45", filter(postive_target_data, age_cat == "Greater than 45")),
  list("< 25, Caucasian", filter(postive_target_data, age_cat == "Less than 25", race == "Caucasian")),
  list("< 25, African-American", filter(postive_target_data, age_cat == "Less than 25", race == "African-American")),
  # list("< 25, Hispanic", filter(postive_target_data, age_cat == "Less than 25", race == "Hispanic")),  # n=7
  # list("< 25, Other", filter(postive_target_data, age_cat == "Less than 25", race == "Other")),  # n=6
  list("25 - 45, Caucasian", filter(postive_target_data, age_cat == "25 - 45", race == "Caucasian")),
  list("25 - 45, African-American", filter(postive_target_data, age_cat == "25 - 45", race == "African-American")),
  # list("25 - 45, Hispanic", filter(postive_target_data, age_cat == "25 - 45", race == "Hispanic")),  # n=7
  # list("25 - 45, Native American", filter(postive_target_data, age_cat == "25 - 45", race == "Native American")),  # n=2
  list("> 45, Caucasian", filter(postive_target_data, age_cat == "Greater than 45", race == "Caucasian")),
  list("> 45, African-American", filter(postive_target_data, age_cat == "Greater than 45", race == "African-American"))
  # list("> 45, Hispanic", filter(postive_target_data, age_cat == "Greater than 45", race == "Hispanic")),
  # list("> 45, Native American", filter(postive_target_data, age_cat == "Greater than 45", race == "Native American"))
)

audit_plot_label <- function(alpha, group_label, n) {
  sprintf("alpha=%g [%s] n=%d", alpha, group_label, n)
}

# 三、数据计算（默认逐组出图；第四节 CSV 仍 plot = FALSE）
# 汇总森林图见 Fig4-6.R（test2 + 与 alpha 一致）
for (ci_method in ci_methods) {
  fl <- filter_audit_rows(audit_rows, ci_method, theta_P)

  cat("\n==========", ci_method,
      sprintf("（计算 %d 组, with plot）==========\n", length(fl$kept)))

  for (row in fl$kept) {
    group_label <- row[[1]]
    group_specific_data <- row[[2]]
    cat("\n-------------------------------\n")
    cat("[", group_label, "] 分组样本数:", nrow(group_specific_data), "\n", sep = "")
    run_ci(group_specific_data, ci_method, theta_P, alpha, alpha2,
           test = test, mu0 = epsilon0,
           epsilon_lo = epsilon_lo, epsilon_hi = epsilon_hi,
           group = group_label, plot = TRUE,
           plot_label = audit_plot_label(
             alpha, group_label, nrow(group_specific_data)
           ))
  }
}

# 四、数据保存（save_data = TRUE 时写入 data/*.csv）
# Inf/-Inf 存为 "-Inf" / "Inf" 字符串，避免 Excel 将 -Inf 当公式
fmt_audit_csv <- function(v) {
  if (length(v) != 1L || is.na(v)) {
    return(NA_character_)
  }
  if (is.infinite(v)) {
    return(if (v > 0) "Inf" else "-Inf")
  }
  format(round(v, 4), trim = TRUE, scientific = FALSE)
}

if (isTRUE(save_data)) {
  for (ci_method in ci_methods) {
    fl <- filter_audit_rows(audit_rows, ci_method, theta_P)
    al_name <- switch(ci_method,
      "ELCI" = "ELAL",
      "CELCI" = "CELAL",
      "SELCI" = "SELAL",
      "ASELCI" = "ASELAL"
    )
    audit_df <- bind_rows(lapply(fl$kept, function(row) {
      r <- run_ci(row[[2]], ci_method, theta_P, alpha, alpha2,
                  test = test, mu0 = epsilon0,
                  epsilon_lo = epsilon_lo, epsilon_hi = epsilon_hi,
                  group = row[[1]], plot = FALSE)
      tibble(
        group = row[[1]],
        lb = fmt_audit_csv(r$lb),
        ub = fmt_audit_csv(r$ub),
        epsilonG = fmt_audit_csv(r$epsilonG),
        ELAL = fmt_audit_csv(r[[al_name]]),
        n = as.integer(r$n)
      )
    }))
    out_path <- if (ci_method == "ELCI") {
      sprintf("data/%s-Fairness-Audit-%.2f-%s.csv", dataname, alpha, ci_method)
    } else {
      sprintf(
        "data/%s-Fairness-Audit-%.2f-%s-%s.csv",
        dataname, alpha, test, ci_method
      )
    }
    write_csv(audit_df, file(out_path), quote = "all")
    cat("已保存:", out_path, "\n")
  }
} else {
  cat("\n未保存 CSV（save_data = FALSE）\n")
}

# 五、跳过分组汇总（全部结果输出后提醒）
for (ci_method in ci_methods) {
  fl <- filter_audit_rows(audit_rows, ci_method, theta_P)
  print_skipped_audit_groups(fl$skipped, ci_method)
}
