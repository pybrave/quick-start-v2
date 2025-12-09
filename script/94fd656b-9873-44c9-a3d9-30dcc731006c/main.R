library(tidyverse)
library(glmnet)

set.seed(2025)

params <- jsonlite::fromJSON("params.json")

# method: "lasso" / "ridge" / "elasticnet"
# alpha_en: Elastic Net alpha，只有 method="elasticnet" 时生效
# standardize: 是否标准化
# family: binomial / gaussian

method  <- params$method
alpha_en <-  params$alpha_en  #0.5
standardize = params$standardize
family <-  params$family #"gaussian"


df <- read_tsv(params$regression$content)
independent_variable <- params$regression$independent_variable$columns_name 
outcome <- params$regression$outcome$columns_name

X <- as.matrix(df[, independent_variable])
y <- df[[outcome]]

alpha <- switch(method,
                "lasso" = 1,
                "ridge" = 0,
                "elasticnet" = alpha_en)

cv_model <- cv.glmnet(X, y, alpha = alpha, family = family, standardize = standardize)

pdf(file = str_glue("output/cv_model.pdf") )
plot(cv_model)
dev.off()


sink(file = str_glue("output/cv_model.txt"))
print(cv_model)
sink()


best_lambda <- cv_model$lambda.min

coef_df <- as.data.frame(as.matrix(coef(cv_model, s = "lambda.min")))
names(coef_df) <- "Coefficient"
coef_df$Variable <- rownames(coef_df)
rownames(coef_df) <- NULL
if (family == "binomial") {
  coef_df$OR <- exp(coef_df$Coefficient)
}

write_tsv(coef_df,file = str_glue("output/coef_df.tsv"))






