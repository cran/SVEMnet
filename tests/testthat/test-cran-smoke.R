# Keep these public-API checks small enough to run on CRAN. The larger
# bootstrap, cross-validation, plotting and permutation suites skip CRAN.

test_that("small elastic-net ensembles predict, score and select candidates", {
  d <- expand.grid(x = seq(-1, 1, length.out = 8),
                   z = c(-1, 1), g = factor(c("a", "b")))
  d$y <- 2 + 1.5 * d$x - .7 * d$z + (d$g == "b") +
    .1 * sin(seq_len(nrow(d)))
  set.seed(9101)
  fit <- SVEMnet(y ~ x + z + g, d, nBoot = 3,
                 glmnet_alpha = c(.5, 1), relaxed = TRUE)
  expect_s3_class(fit, "svem_model")
  expect_equal(fit$nBoot, 3L)
  expect_true(all(is.finite(fit$coef_matrix)))

  nd <- d[c(1, 7, 18), ]
  X <- model.matrix(y ~ x + z + g, nd)
  member_predictions <- X %*% t(fit$coef_matrix[, colnames(X), drop = FALSE])
  pred <- predict(fit, nd, debias = FALSE,
                  se.fit = TRUE, interval = TRUE, level = .8)
  expect_equal(as.numeric(pred$fit), as.numeric(rowMeans(member_predictions)))
  expect_equal(as.numeric(pred$se.fit), as.numeric(apply(member_predictions, 1, sd)))
  expect_equal(as.numeric(pred$lwr), as.numeric(apply(member_predictions, 1, quantile, .1)))
  expect_equal(as.numeric(pred$upr), as.numeric(apply(member_predictions, 1, quantile, .9)))
  restored <- unserialize(serialize(fit, NULL))
  expect_identical(predict(restored, nd), predict(fit, nd))

  scored <- svem_score_random(
    list(y = fit),
    goals = list(y = list(goal = "max", weight = 1,
                         lower_acceptable = -2, upper_acceptable = 7)),
    n = 24, combine = "mean", numeric_sampler = "uniform", verbose = FALSE
  )
  tab <- scored$score_table
  expect_equal(nrow(tab), 24L)
  expect_true(all(is.finite(tab$y_pred)))
  expect_true(all(tab$x >= -1 & tab$x <= 1 & tab$z >= -1 & tab$z <= 1))
  expect_true(all(as.character(tab$g) %in% levels(d$g)))
  expect_equal(tab$y_pred, as.numeric(predict(fit, tab, debias = FALSE)))
  expected_desirability <- pmin(1, pmax(0, (tab$y_pred + 2) / 9))
  expect_equal(tab$y_des, expected_desirability)
  expect_equal(tab$score, expected_desirability)
  selected <- svem_select_from_score_table(tab, target = "score",
                                            direction = "max", k = 2,
                                            top_type = "n", top = 6)
  expect_equal(selected$best$score, max(tab$score))
  expect_equal(nrow(selected$candidates), 2L)
  expect_true(all(selected$candidates$score %in% tab$score))

  # Reuse the fitted ensemble and small candidate table for the experimental
  # proposer; no additional bootstrap fit or random search is needed.
  batch <- svem_thompson_batch(
    list(y = fit), goals = list(y = list(goal = "max", weight = 1)),
    batch_size = 2, candidates = tab[c("x", "z", "g")], verbose = FALSE
  )
  expect_s3_class(batch, "svem_thompson_batch")
  expect_length(batch$candidate_index, 2L)
  expect_identical(anyDuplicated(batch$candidate_index), 0L)
  chosen <- tab[batch$candidate_index, ]
  X_chosen <- model.matrix(y ~ x + z + g, transform(chosen, y = 0))
  expected_draw <- vapply(seq_len(2), function(i) {
    sum(X_chosen[i, ] * fit$coef_matrix[batch$members[i, "y"], colnames(X_chosen)])
  }, numeric(1L))
  expect_equal(batch$proposals$y_draw, expected_draw)
})

test_that("small binomial ensembles average members on the requested scale", {
  d <- data.frame(x = rep(seq(-1, 1, length.out = 10), 4),
                  z = rep(c(-1, 1), each = 20),
                  y = rep(c(0, 1, 0, 0, 1, 0, 1, 1, 0, 1), 4))
  set.seed(9102)
  fit <- SVEMnet(y ~ x + z, d, nBoot = 3, glmnet_alpha = 1,
                 relaxed = FALSE, family = "binomial")
  nd <- d[1:5, ]
  X <- model.matrix(y ~ x + z, nd)
  member_links <- X %*% t(fit$coef_matrix[, colnames(X), drop = FALSE])
  link <- predict(fit, nd, type = "link", debias = FALSE)
  response <- predict(fit, nd, type = "response", debias = FALSE)
  expect_length(response, 5L)
  expect_true(all(is.finite(link)))
  expect_equal(as.numeric(link), as.numeric(rowMeans(member_links)))
  expect_equal(as.numeric(response), as.numeric(rowMeans(plogis(member_links))))
  expect_true(all(response >= 0 & response <= 1))
  expect_identical(predict(fit, nd, type = "class"), as.integer(response >= .5))
})

test_that("small forward fits agree with an independent least-squares fit", {
  d <- data.frame(x = seq(-1, 1, length.out = 24))
  d$y <- 2 + 3 * d$x + .01 * sin(seq_len(nrow(d)))
  reference <- lm(y ~ x, d)
  fit <- forward_aicc(y ~ x, d)
  ensemble <- svem_forward(y ~ x, d, nBoot = 1, weight_scheme = "Identity")
  expect_identical(fit$selected_terms, "x")
  expect_equal(unname(fit$parms), unname(coef(reference)), tolerance = 1e-8)
  expect_equal(as.numeric(predict(fit, d)), as.numeric(fitted(reference)),
               tolerance = 1e-8)
  expect_equal(unname(ensemble$parms), unname(coef(reference)), tolerance = 1e-8)
})

test_that("a small cross-validated fit predicts from its selected coefficients", {
  d <- expand.grid(x = seq(-1, 1, length.out = 12), z = c(-1, 1))
  d$y <- 2 + 1.5 * d$x - .5 * d$z + .1 * cos(seq_len(nrow(d)))
  fit <- glmnet_with_cv(y ~ x + z, d, nfolds = 3, repeats = 1,
                        glmnet_alpha = 1, seed = 9103)
  nd <- d[c(1, 8, 16), ]
  X <- model.matrix(y ~ x + z, nd)
  expect_true(all(is.finite(fit$parms)))
  expect_equal(as.numeric(predict(fit, nd)),
               as.numeric(X %*% fit$parms[colnames(X)]))
})
