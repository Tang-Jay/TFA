# =============================================== #
#        Compas Fig6 Audit Data (EL/CEL/SEL/ASEL) #
# =============================================== #
rm(list = ls())
source("Functions/EL_CI.R")
source("Functions/CEL_CI.R")
source("Functions/SEL_CI.R")
source("Functions/ASEL_CI.R")
source("Functions/audit_run_ci.R")

library(dplyr)
library(readr)

test <- "test1"
ci_methods <- c("ELCI", "CELCI", "SELCI", "ASELCI")
alpha <- 0.95
alpha2 <- 0.08
dataname <- "Sex-Age"
save_data <- TRUE  # FALSE 时不写 CSV
epsilon0 <- 0
epsilon_tol <- 0.05
epsilon_lo <- -epsilon_tol
epsilon_hi <- epsilon_tol
test <- normalize_audit_test(test)
cat("test =", test,
    "| CEL alternative =", audit_test_alternative(test), "\n")

age_cat_label <- function(x) {
  switch(x,
    "Less than 25" = "< 25",
    "Greater than 45" = "> 45",
    x
  )
}

sex_label <- function(x) {
  switch(x, "Male" = "M", "Female" = "F", x)
}

# ---------------------获得数据------------------- #
compas_data <- read_csv("data/compas_data.csv")
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

# 选出阳预测
postive_clean_data <- clean_data %>%
  filter(Y_h == 1)
table(postive_clean_data$L)

# ---------------------基准线-------------------- #
# 选出白人族
caucasian_data <- clean_data %>%
  filter(race == "Caucasian")
# 选出阳预测白人族
postive_caucasian_data <- caucasian_data %>%
  filter(Y_h == 1)
# 查看数据
table(postive_caucasian_data$Y)
table(postive_caucasian_data$L)
# 基准线
theta_P = mean(postive_caucasian_data$Y)
cat('theta_P =', theta_P, '\n')

# 三、数据计算
# ---------------------开始计算------------------- #
sexs <- unique(postive_clean_data$sex)
age_cats <- unique(postive_clean_data$age_cat)

for (ci_method in ci_methods) {
  cat("\n==========", ci_method, "==========\n")

  # ALL
  group_specific_data <- postive_clean_data
  cat("\n-------------------------------\n")
  cat("ALL,")
  cat("分组样本数:", nrow(group_specific_data), "\n")
  run_ci(group_specific_data, ci_method, theta_P, alpha, alpha2,
         test = test, mu0 = epsilon0,
         epsilon_lo = epsilon_lo, epsilon_hi = epsilon_hi,
         group = "All")

  # Sex
  for (sex_name in sexs) {
    group_specific_data <- postive_clean_data %>%
      filter(sex == sex_name)
    cat("\n-------------------------------\n")
    cat(sex_name, ",")
    cat("分组样本数:", nrow(group_specific_data), "\n")
    run_ci(group_specific_data, ci_method, theta_P, alpha, alpha2,
           test = test, mu0 = epsilon0,
           epsilon_lo = epsilon_lo, epsilon_hi = epsilon_hi,
           group = sex_label(sex_name))
  }

  # age_cat
  for (age_cat_name in age_cats) {
    group_specific_data <- postive_clean_data %>%
      filter(age_cat == age_cat_name)
    cat("\n-------------------------------\n")
    cat(age_cat_name, ",")
    cat("分组样本数:", nrow(group_specific_data), "\n")
    run_ci(group_specific_data, ci_method, theta_P, alpha, alpha2,
           test = test, mu0 = epsilon0,
           epsilon_lo = epsilon_lo, epsilon_hi = epsilon_hi,
           group = age_cat_label(age_cat_name))
  }

  # sex * age_cat
  for (sex_name in sexs) {
    for (age_cat_name in age_cats) {
      group_specific_data <- postive_clean_data %>%
        filter(sex == sex_name, age_cat == age_cat_name)
      cat("\n-------------------------------\n")
      cat(sex_name, ",", age_cat_name, ",")
      cat("分组样本数:", nrow(group_specific_data), "\n")
      run_ci(group_specific_data, ci_method, theta_P, alpha, alpha2,
             test = test, mu0 = epsilon0,
             epsilon_lo = epsilon_lo, epsilon_hi = epsilon_hi,
             group = paste0(sex_label(sex_name), " (", age_cat_label(age_cat_name), ")"))
    }
  }
}

# 四、数据保存（save_data = TRUE 时写入 data/*.csv）
audit_rows <- list(
  list("All", postive_clean_data),
  list("< 25", filter(postive_clean_data, age_cat == "Less than 25")),
  list("25 - 45", filter(postive_clean_data, age_cat == "25 - 45")),
  list("> 45", filter(postive_clean_data, age_cat == "Greater than 45")),
  list("M", filter(postive_clean_data, sex == "Male")),
  list("F", filter(postive_clean_data, sex == "Female")),
  list("M (< 25)", filter(postive_clean_data, sex == "Male", age_cat == "Less than 25")),
  list("F (< 25)", filter(postive_clean_data, sex == "Female", age_cat == "Less than 25")),
  list("M (25 - 45)", filter(postive_clean_data, sex == "Male", age_cat == "25 - 45")),
  list("F (25 - 45)", filter(postive_clean_data, sex == "Female", age_cat == "25 - 45")),
  list("M (> 45)", filter(postive_clean_data, sex == "Male", age_cat == "Greater than 45")),
  list("F (> 45)", filter(postive_clean_data, sex == "Female", age_cat == "Greater than 45"))
)
if (isTRUE(save_data)) {
  for (ci_method in ci_methods) {
    al_name <- switch(ci_method,
      "ELCI" = "ELAL",
      "CELCI" = "CELAL",
      "SELCI" = "SELAL",
      "ASELCI" = "ASELAL"
    )
    audit_df <- bind_rows(lapply(audit_rows, function(row) {
      r <- run_ci(row[[2]], ci_method, theta_P, alpha, alpha2,
                  test = test, mu0 = epsilon0,
                  epsilon_lo = epsilon_lo, epsilon_hi = epsilon_hi,
                  group = row[[1]], plot = FALSE)
      tibble(
        group = row[[1]],
        lb = round(r$lb, 4),
        ub = round(r$ub, 4),
        epsilonG = round(r$epsilonG, 4),
        ELAL = round(r[[al_name]], 4),
        n = r$n
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
    write_csv(audit_df, out_path)
    cat("已保存:", out_path, "\n")
  }
} else {
  cat("\n未保存 CSV（save_data = FALSE）\n")
}
