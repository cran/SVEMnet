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

## ----show-file6, echo=FALSE, results='asis'-----------------------------------
raw_file6 <- readLines("bic_sim.R")
cat("```r\n")
cat(raw_file6, sep="\n")
cat("\n```\n")

