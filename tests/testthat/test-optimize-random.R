skip_on_cran()
skip_if_not_installed("SVEMnet")
skip_if_not_installed("cluster")

test_that("svem_optimize_random returns best row and diverse candidates", {
  set.seed(42)
  n  <- 150
  X1 <- runif(n)
  X2 <- runif(n)
  F  <- factor(sample(c("lo","hi"), n, replace = TRUE))

  y <- 1 + 1.2*X1 - 0.4*X2 + 0.5*(F=="hi") + rnorm(n, 0, 0.25)
  z <- 0.4 + 0.8*X1 + 0.2*X2 - 0.2*(F=="hi") + rnorm(n, 0, 0.25)
  d <- data.frame(y, z, X1, X2, F)

  fit1 <- SVEMnet::SVEMnet(y ~ X1 + X2 + F, d, nBoot = 30, glmnet_alpha = 1)
  fit2 <- SVEMnet::SVEMnet(z ~ X1 + X2 + F, d, nBoot = 30, glmnet_alpha = 1)

  goals <- list(
    y = list(goal = "max", weight = 0.6),
    z = list(goal = "target", weight = 0.4, target = 0.8)
  )

  out <- SVEMnet::svem_optimize_random(
    objects      = list(y = fit1, z = fit2),
    goals        = goals,
    n            = 600,
    level        = 0.95,
    k_candidates = 3,
    top_frac     = 0.05,
    verbose      = FALSE
  )

  # core structure
  expect_true(is.list(out))
  expect_true(is.data.frame(out$best_x))
  expect_true(is.numeric(out$best_pred["y"]))
  expect_true(is.numeric(out$best_pred["z"]))
  expect_true(all(c("best_idx","best_x","best_pred","score_table","goals") %in% names(out)))

  # new weight fields (old `weights` was removed)
  expect_true(all(c("weights_original","weights_final") %in% names(out)))
  expect_true(is.null(out$wmt_p_values))           # since reweight_by_wmt = FALSE
  expect_true(is.null(out$wmt_multipliers))        # since reweight_by_wmt = FALSE

  # score table sanity
  expect_true(is.data.frame(out$score_table))
  expect_true(all(c("score","uncertainty_measure") %in% colnames(out$score_table)))

  # candidates (if requested/returned)
  if (!is.null(out$candidates)) {
    expect_true(nrow(out$candidates) <= 3)
    expect_true(all(c("y","z") %in% colnames(out$candidates)))
    expect_true(all(c("score","uncertainty_measure") %in% colnames(out$candidates)))
  }

  # exploration artifacts may be present depending on defaults; just sanity check if they exist
  if (!is.null(out$exploration_candidates)) {
    expect_true(is.data.frame(out$exploration_candidates))
    expect_true(all(c("score","uncertainty_measure") %in% colnames(out$exploration_candidates)))
  }
})
