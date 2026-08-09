# =============================================== #
#   Table 1: African-American All-group LB      #
# =============================================== #
rm(list = ls())

dataname <- "African-American"
test <- "test1" # test1/greater、test2/less、test3/or（CEL/SEL/ASEL 输出含 test）
alphas <- c(0.90, 0.95)
methods <- c(EL = "ELCI", CEL = "CELCI", SEL = "SELCI", ASEL = "ASELCI")

# ELCI 无 test 参数，文件名不含 test；CELCI/SELCI/ASELCI 含 test
audit_csv_path <- function(dataname, alpha, test, ci_method) {
  if (ci_method == "ELCI") {
    sprintf("data/%s-Fairness-Audit-%.2f-%s.csv", dataname, alpha, ci_method)
  } else {
    sprintf(
      "data/%s-Fairness-Audit-%.2f-%s-%s.csv",
      dataname, alpha, test, ci_method
    )
  }
}

read_all_lb <- function(alpha) {
  vapply(methods, function(m) {
    path <- audit_csv_path(dataname, alpha, test, m)
    df <- read.csv(path)
    df$lb[df$group == "All"]
  }, numeric(1))
}

mat <- t(sapply(alphas, read_all_lb))
tab1 <- data.frame(
  `Confidence level` = paste0(as.integer(alphas * 100), "%"),
  mat,
  check.names = FALSE
)
colnames(tab1)[-1] <- names(methods)

print(tab1, digits = 4, row.names = FALSE)


