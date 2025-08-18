## ----fig.width=6, fig.height=4------------------------------------------------
library(SVEMnet)

# Example data
data <- iris
svem_model <- SVEMnet(Sepal.Length ~ ., data = data, nBoot = 300)
coef(svem_model)

## ----fig.width=6, fig.height=4------------------------------------------------
plot(svem_model)

## -----------------------------------------------------------------------------
predictions <- predict(svem_model, data)
print(predictions)

## ----echo=FALSE, out.width='100%', fig.cap="Whole model test result"----------
knitr::include_graphics("figures/whole_model_test.png")

## ----echo=FALSE, out.width='100%', fig.cap="Whole Model Test Results for Example 2"----
knitr::include_graphics("figures/whole_model_2.png")

## ----show-file4, echo=FALSE, results='asis'-----------------------------------
raw_file3 <- readLines("debias_and_objective.R")
cat("```r\n")
cat(raw_file3, sep="\n")
cat("\n```\n")

## ----echo=FALSE, out.width='100%', fig.cap="Simulation-blocked Residual LRMSE for different objective functions and debiasing"----
knitr::include_graphics("figures/14AUG25.png")

## ----show-file5, echo=FALSE, results='asis'-----------------------------------
raw_file4 <- readLines("svem_vs_lassoCV_par_v2.R")
cat("```r\n")
cat(raw_file4, sep="\n")
cat("\n```\n")

## ----show-file6, echo=FALSE, results='asis'-----------------------------------
raw_file6 <- readLines("debias_aicc_cvglm.R")
cat("```r\n")
cat(raw_file6, sep="\n")
cat("\n```\n")

