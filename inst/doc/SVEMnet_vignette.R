## ----fig.width=6, fig.height=4------------------------------------------------
library(SVEMnet)

# Example data
data <- iris
svem_model <- SVEMnet(Sepal.Length ~ ., data = data, relaxed=FALSE,glmnet_alpha=c(1),nBoot = 50)
coef(svem_model)

## ----fig.width=6, fig.height=4------------------------------------------------
plot(svem_model)

## ----eval=TRUE----------------------------------------------------------------
predictions <- predict(svem_model, data)
print(predictions)

## ----echo=FALSE, out.width='100%', fig.cap="Whole Model Test Results for Example 2"----
knitr::include_graphics("figures/whole_model_2.png")

## ----eval=FALSE---------------------------------------------------------------
# # Simulate simple mixed-type data
# set.seed(1)
# n  <- 120
# X1 <- runif(n)
# X2 <- runif(n)
# F  <- factor(sample(c("lo","hi"), n, replace = TRUE))
# y1 <- 1 + 1.5*X1 - 0.8*X2 + 0.4*(F=="hi") + rnorm(n, 0, 0.2)
# y2 <- 0.7 + 0.4*X1 + 0.4*X2 - 0.2*(F=="hi") + rnorm(n, 0, 0.2)
# dat <- data.frame(y1, y2, X1, X2, F)
# 
# # Fit two SVEM models (keep defaults modest in vignettes)
# m1 <- SVEMnet(y1 ~ X1 + X2 + F, dat, nBoot = 30)
# m2 <- SVEMnet(y2 ~ X1 + X2 + F, dat, nBoot = 30)
# 
# # Predict
# head(predict(m1, newdata = dat))
# head(predict(m2, newdata = dat))

## ----eval=FALSE---------------------------------------------------------------
# res_serial_y1 <- svem_significance_test(
#   y1 ~ X1 + X2 + F, dat,
#   nPoint = 2000, nSVEM = 10, nPerm = 150,
#   nBoot = 80, glmnet_alpha = 1, relaxed = FALSE,
#   verbose = TRUE
# )
# res_serial_y1$p_value

## ----eval=FALSE---------------------------------------------------------------
# objs  <- list(y1 = m1, y2 = m2)
# goals <- list(
#   y1 = list(goal = "max",    weight = 0.6),
#   y2 = list(goal = "target", weight = 0.4, target = 0.9)
# )
# 
# opt_out <- svem_optimize_random(
#   objects      = objs,
#   goals        = goals,
#   n            = 3000,
#   agg          = "mean",
#   debias       = FALSE,
#   ci           = TRUE,
#   level        = 0.95,
#   k_candidates = 5,
#   top_frac     = 0.02,
#   verbose      = TRUE
# )
# 
# opt_out$best_x
# opt_out$best_pred
# head(opt_out$candidates)

