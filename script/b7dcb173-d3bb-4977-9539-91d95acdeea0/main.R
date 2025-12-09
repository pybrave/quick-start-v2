# install.packages("semTools")
library(tidyverse)
library(lavaan)
library(MASS)
library(semTools)
library(lavaanPlot)


params <- jsonlite::fromJSON("params.json")

df <- read_tsv(params$mediation_file$content)

# 验证性因素分析（CFA, Confirmatory Factor Analysis
model_cfa <- '
  LS =~ LS1 + LS2 + LS3
  OT =~ OT1 + OT2 + OT3
  JS =~ JS1 + JS2 + JS3
'

fit_cfa <- cfa(model_cfa, data=df)
summary(fit_cfa, fit.measures=TRUE, standardized=TRUE)

reliability(fit_cfa)

model_sem <- '
  # Measurement model
  LS =~ LS1 + LS2 + LS3
  OT =~ OT1 + OT2 + OT3
  JS =~ JS1 + JS2 + JS3
  
  # Structural model
  OT ~ a*LS
  JS ~ b*OT + c*LS
  
  # Indirect effect
  indirect := a*b
  
  # Total effect
  total := c + (a*b)
'

fit_sem <- sem(model_sem, data=df, se="bootstrap", bootstrap=5000)
sink(file = str_glue("output/sem_summary.txt"))
summary(fit_sem, fit.measures=TRUE, standardized=TRUE, ci=TRUE)
sink()

standardizedSolution(fit_sem)  |>
  write_tsv(file = str_glue("output/standardized_solution.tsv"))

# install.packages("lavaanPlot")
lavaanPlot(model = fit_sem)

fit_df<-fitMeasures(fit_sem, c("chisq","df","chisq.scaled","cfi","tli","rmsea","srmr")) 
sink(file = str_glue("output/model_fit_indicators.txt"))
print(fit_df)
sink()

