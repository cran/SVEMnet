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

## ----echo=FALSE, out.width='100%', fig.cap="LRMSE for {debias}x{objective}"----
knitr::include_graphics("figures/lrmse.png")

## ----echo=FALSE, out.width='100%', fig.cap="Paired LRMSE for {debias}x{objective}"----
knitr::include_graphics("figures/paired_lrmse.png")

## ----show-file, echo=FALSE, results='asis'------------------------------------
raw_file <- readLines("svemnet_sim.R")
cat("```r\n")
cat(raw_file, sep="\n")
cat("\n```\n")

## ----echo=FALSE, out.width='100%', fig.cap="LRMSE for different weight_scheme"----
knitr::include_graphics("figures/lrmse_identity.png")

## ----echo=FALSE, out.width='100%', fig.cap="Residual LRMSE for different weight_scheme"----
knitr::include_graphics("figures/paired_lmrse_identity.png")

## ----show-file2, echo=FALSE, results='asis'-----------------------------------
raw_file2 <- readLines("svemnet_sim_identity.R")
cat("```r\n")
cat(raw_file2, sep="\n")
cat("\n```\n")

## ----echo=FALSE, out.width='100%', fig.cap="Residual LRMSE for different modeling approaches"----
knitr::include_graphics("figures/mixture_v2.png")

## ----show-file3, echo=FALSE, results='asis'-----------------------------------
raw_file2 <- readLines("svemnet_sim_mixture_v2.R")
cat("```r\n")
cat(raw_file2, sep="\n")
cat("\n```\n")

