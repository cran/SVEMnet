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

## ----echo=FALSE, out.width='100%', fig.cap="Whole model test result"----------
knitr::include_graphics("figures/whole_model_test.png")

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

## ----eval=FALSE---------------------------------------------------------------
# data(lipid_screen)
# str(lipid_screen)
# 
# spec <- bigexp_terms(
#   Potency ~ PEG + Helper + Ionizable + Cholesterol +
#     Ionizable_Lipid_Type + N_P_ratio + flow_rate,
#   data               = lipid_screen,
#   factorial_order    = 3,
#   include_pure_cubic = TRUE,
#   include_pc_3way    = FALSE
# )
# 
# form_pot <- bigexp_formula(spec, "Potency")
# form_siz <- bigexp_formula(spec, "Size")
# form_pdi <- bigexp_formula(spec, "PDI")

## ----eval=FALSE---------------------------------------------------------------
# set.seed(1)
# fit_pot <- SVEMnet(form_pot, lipid_screen, nBoot = 40)
# fit_siz <- SVEMnet(form_siz, lipid_screen, nBoot = 40)
# fit_pdi <- SVEMnet(form_pdi, lipid_screen, nBoot = 40)
# 
# objs <- list(Potency = fit_pot, Size = fit_siz, PDI = fit_pdi)

## ----eval=FALSE---------------------------------------------------------------
# test_pot <- svem_significance_test(form_pot, lipid_screen, nPoint=2000, nSVEM=10, nPerm=150,
#                                    nBoot=80, glmnet_alpha=1, relaxed=FALSE, verbose=TRUE)
# test_siz <- svem_significance_test(form_siz, lipid_screen, nPoint=2000, nSVEM=10, nPerm=150,
#                                    nBoot=80, glmnet_alpha=1, relaxed=FALSE, verbose=TRUE)
# test_pdi <- svem_significance_test(form_pdi, lipid_screen, nPoint=2000, nSVEM=10, nPerm=150,
#                                    nBoot=80, glmnet_alpha=1, relaxed=FALSE, verbose=TRUE)
# 
# c(Potency = test_pot$p_value, Size = test_siz$p_value, PDI = test_pdi$p_value)

## ----eval=FALSE---------------------------------------------------------------
# # Parallel runs (example; adjust nCore to your machine)
# par_pot <- svem_significance_test_parallel(form_pot, lipid_screen,
#                                            nPoint=3000, nSVEM=10, nPerm=150,
#                                            nCore = max(1L, parallel::detectCores()-1L),
#                                            seed = 123, verbose=TRUE)
# par_siz <- svem_significance_test_parallel(form_siz, lipid_screen,
#                                            nPoint=3000, nSVEM=10, nPerm=150,
#                                            nCore = max(1L, parallel::detectCores()-1L),
#                                            seed = 123, verbose=TRUE)
# par_pdi <- svem_significance_test_parallel(form_pdi, lipid_screen,
#                                            nPoint=3000, nSVEM=10, nPerm=150,
#                                            nCore = max(1L, parallel::detectCores()-1L),
#                                            seed = 123, verbose=TRUE)
# 
# # Plot all three together
# plot(par_pot, par_siz, par_pdi, labels = c("Potency","Size","PDI"))

## ----eval=FALSE---------------------------------------------------------------
# goals <- list(
#   Potency = list(goal = "max", weight = 0.7),
#   Size    = list(goal = "min", weight = 0.2),
#   PDI     = list(goal = "min", weight = 0.1)
# )
# 
# mixL <- list(list(
#   vars  = c("Cholesterol","PEG","Ionizable","Helper"),
#   lower = c(0.10, 0.01, 0.10, 0.10),
#   upper = c(0.60, 0.05, 0.60, 0.60),
#   total = 1
# ))
# 
# opt_out <- svem_optimize_random(
#   objects        = objs,
#   goals          = goals,
#   n              = 10000,
#   mixture_groups = mixL,
#   agg            = "mean",
#   debias         = FALSE,
#   ci             = TRUE,
#   level          = 0.95,
#   k_candidates   = 5,
#   top_frac       = 0.01,
#   verbose        = TRUE
# )
# 
# opt_out$best_x
# opt_out$best_pred
# opt_out$best_ci
# opt_out$candidates

