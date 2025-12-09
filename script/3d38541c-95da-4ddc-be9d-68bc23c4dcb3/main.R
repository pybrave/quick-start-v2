library(tidyverse)
params <- jsonlite::fromJSON("params.json")

df <- read_tsv(params$logistic$content)

independent_variable <- params$logistic$independent_variable$columns_name |>paste0(collapse = " + ")
outcome <- params$logistic$outcome$columns_name
# y ~ x1 + x2
logistic_formula <- as.formula(str_glue("{outcome} ~ {independent_variable}"))
model <- glm(logistic_formula, data = df, family = binomial)

sink(file = str_glue("output/summary.txt"))
summary(model)
sink()
sink(file = str_glue("output/prompt.ai"))
summary(model)
sink()

OR <- exp(coef(model))
OR
CI <- exp(confint.default(model))
CI
p_values <- summary(model)$coefficients[,4]

# 整理成表格
result <- data.frame(
  Variable = names(OR),
  OR = OR,
  Lower95CI = CI[,1],
  Upper95CI = CI[,2],
  Pvalue = p_values
)

write_tsv(result, file = str_glue("output/logistic_res.tsv"))


