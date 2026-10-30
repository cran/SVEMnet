# Independent glmnet-path selection oracle.

## ---------------------------------------------------------------------------
## Oracle for the published selector (paper Fig. 1, App. A). Should pass now.
## ---------------------------------------------------------------------------

.oracle_member <- function(X, y, wt, wv, objective, family = "gaussian",
                           relaxed = FALSE, gammas = c(0.2, 0.6, 1)) {
  n <- length(y)
  args <- list(x=X,y=y,weights=wt,alpha=1,nlambda=500,
               standardize=TRUE,intercept=TRUE,family=family,relax=relaxed)
  if ("control" %in% names(formals(glmnet::glmnet))) args$control <- list(maxit=1e6)
  else args$maxit <- 1e6
  g <- do.call(glmnet::glmnet,args)
  neff <- sum(wv)^2 / sum(wv^2)
  nadm <- max(2, min(n, neff))
  best <- list(score = Inf)
  for (gam in if (relaxed) gammas else 1) {
    B <- if (relaxed) as.matrix(coef(g, s = g$lambda, gamma = gam))
         else as.matrix(coef(g, s = g$lambda))
    eta <- cbind(1, X) %*% B   # n x L, exact path columns (no interpolation)
    if (family == "gaussian") {
      loss <- colSums(wv * (y - eta)^2)
      fit_term <- n * log(pmax(loss, .Machine$double.xmin) / n)
    } else {
      mu <- pmin(pmax(1 / (1 + exp(-eta)), 1e-8), 1 - 1e-8)
      loss <- pmax(-colSums(wv * (y * log(mu) + (1 - y) * log(1 - mu))),
                   .Machine$double.eps)
      fit_term <- 2 * loss
    }
    k <- 1 + colSums(B[-1, , drop = FALSE] != 0)
    adm <- (k - 1) < nadm
    crit <- switch(objective,
      wSSE = loss,
      wAIC = ifelse(adm, fit_term + 2 * k, Inf),
      wBIC = ifelse(adm, fit_term + log(nadm) * k, Inf))
    j <- which.min(crit)
    if (crit[j] < best$score) {
      best <- list(score = crit[j], beta = B[, j], k = k[j], gamma = gam)
    }
  }
  best
}

test_that("Gaussian SVEMnet members match an independent wAIC/wBIC/wSSE oracle", {
  set.seed(5101)
  n <- 30
  d <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1),
                  x3 = runif(n, -1, 1))
  d$y <- 1 + 2 * d$x1 - d$x2 + 0.5 * d$x1 * d$x2 + rnorm(n, sd = 0.3)
  f <- y ~ (x1 + x2 + x3)^2
  X <- model.matrix(f, d)[, -1, drop = FALSE]

  for (relaxed in c(FALSE, TRUE)) {
    for (obj in c("wAIC", "wBIC", "wSSE")) {
      set.seed(77)
      fit <- SVEMnet(f, d, nBoot = 3, glmnet_alpha = 1, relaxed = relaxed,
                     objective = obj, store_member_weights = TRUE)
      expect_identical(fit$diagnostics$fallback_rate, 0)
      for (b in seq_len(fit$nBoot)) {
        o <- .oracle_member(X, d$y, fit$member_weights$training[b, ],
                            fit$member_weights$validation[b, ], obj,
                            relaxed = relaxed)
        expect_equal(unname(fit$coef_matrix[b, ]), unname(o$beta),
                     tolerance = 1e-10)
        expect_equal(fit$member_diagnostics$objective_value[b], unname(o$score),
                     tolerance = 1e-8)
        expect_identical(as.numeric(fit$member_diagnostics$support_size[b]),
                         as.numeric(o$k))
        if (relaxed) expect_equal(fit$best_relax_gammas[b], o$gamma)
      }
    }
  }
})

test_that("binomial SVEMnet members match the 2*NLL + g*k oracle", {
  set.seed(5103)
  n <- 60
  d <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  d$y <- rbinom(n, 1, plogis(-0.2 + 2 * d$x1 - 1.5 * d$x2))
  f <- y ~ x1 * x2
  X <- model.matrix(f, d)[, -1, drop = FALSE]
  set.seed(78)
  fit <- SVEMnet(f, d, nBoot = 3, glmnet_alpha = 1, family = "binomial",
                 store_member_weights = TRUE)
  expect_identical(fit$objective_used, "wBIC")
  expect_false(fit$relaxed)
  for (b in seq_len(fit$nBoot)) {
    o <- .oracle_member(X, d$y, fit$member_weights$training[b, ],
                        fit$member_weights$validation[b, ], "wBIC",
                        family = "binomial")
    expect_equal(unname(fit$coef_matrix[b, ]), unname(o$beta), tolerance = 1e-8)
  }
})

test_that("auto rules resolve by family as published", {
  set.seed(5104)
  d <- data.frame(x = runif(30)); d$y <- d$x + rnorm(30, sd = 0.1)
  d$yb <- rbinom(30, 1, plogis(3 * (d$x - 0.5)))  # noisy: avoid separation
  g <- SVEMnet(y ~ x, d, nBoot = 2)
  b <- suppressWarnings(SVEMnet(yb ~ x, d, nBoot = 2, family = "binomial"))
  expect_identical(c(g$objective_used, b$objective_used), c("wAIC", "wBIC"))
  expect_identical(c(g$relaxed, b$relaxed), c(TRUE, FALSE))
})
