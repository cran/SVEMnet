review_mixed_grid <- function(n = 150L) {
  set.seed(9021)
  g <- as.data.frame(matrix(runif(n * 5L), n, 5L))
  names(g) <- paste0("x", 1:5)
  g$a <- factor(rep(c("p", "q", "r"), length.out = n))
  g$b <- factor(rep(c("low", "high"), each = 2L, length.out = n))
  g
}

test_that("every budget branch keeps layout, design, batch and serialization consistent", {
  skip_if_not_installed("mgcv")
  g <- review_mixed_grid()
  d <- g[seq_len(40L), ]
  d$y <- with(d, sin(3 * x1) + x2 * (a == "q") + .3 * (b == "high"))
  built <- 0L
  notes <- character()
  # Together these budgets should reach full, pair-slope, curvature, level-slope,
  # cell and marginal reductions. rank = 6 is the d + 1 linear limit.
  settings <- c(lapply(c(200, 110, 90, 70, 40, 30, 20), function(m) list(max_columns = m)),
                list(list(rank = 6L)))
  for (control in settings) {
    s <- do.call(spline_terms, c(list(y ~ x1 + x2 + x3 + x4 + x5 + a + b, g), control))
    built <- built + 1L
    notes <- c(notes, s$basis_notes)
    x <- .svem_spline_surface(s, g)
    expect_identical(ncol(x), length(s$penalty))
    expect_identical(length(s$columns), length(s$penalty))
    expect_identical(length(s$column_families), length(s$penalty))
    expect_identical(colnames(x), s$columns)
    expect_lte(length(s$penalty), s$max_columns)
    expect_true(all(is.finite(x)))
    set.seed(11)
    fit <- svem_spline(s, d, nBoot = 2L, lambda = c(Inf, 1, .1))
    whole <- predict(fit, g)
    pieces <- numeric(nrow(g))
    for (i in split(seq_len(nrow(g)), rep(1:3, length.out = nrow(g)))) pieces[i] <- predict(fit, g[i, ])
    expect_equal(whole, pieces, tolerance = 1e-10)
    expect_identical(predict(unserialize(serialize(fit, NULL)), g), whole)
    edf <- fit$member_diagnostics[grep("^edf_", names(fit$member_diagnostics))]
    expect_equal(rowSums(edf), fit$member_diagnostics$slope_edf, tolerance = 1e-10)
  }
  expect_identical(built, length(settings))
  # Confirm the settings really reach the branches they are meant to cover.
  for (pattern in c("slopes omitted", "curvature rank reduced", "Marginal curvature", "linear directions only")) {
    expect_true(any(grepl(pattern, notes, ignore.case = TRUE)), info = pattern)
  }
})

test_that("an unseen nominal combination has no fitted cell-specific effect", {
  skip_if_not_installed("mgcv")
  g <- expand.grid(u = 0:6 / 6, a = factor(c("A", "B")), b = factor(c("x", "y")))
  d <- g[!(g$a == "B" & g$b == "y"), ]
  d$y <- with(d, u + .4 * (a == "B") + .2 * (b == "y")) + .01 * sin(seq_len(nrow(d)))
  s <- spline_terms(y ~ u + a + b, g)
  set.seed(3)
  f <- svem_spline(s, d, nBoot = 5L)
  # Level order A, B and x, y place combination (B, y) in cell 4.
  unseen <- grep("^cell_a:b_4(_mode[0-9]+)?$", colnames(f$coef_matrix))
  expect_gt(length(unseen), 0L)
  expect_lt(max(abs(f$coef_matrix[, unseen])), 1e-8)
  target <- g[g$a == "B" & g$b == "y", ]
  expect_true(all(is.finite(predict(f, target))))
  expect_true(all(.svem_spline_general_covered(s, f$nominal_levels, target)))
})

test_that("an unmeasured level keeps zero coefficients and valid normal equations", {
  skip_if_not_installed("mgcv")
  g <- expand.grid(u = 0:6 / 6, kind = factor(c("a", "b", "c")))
  d <- g[g$kind != "c", ]
  d$y <- with(d, u^2 + .3 * (kind == "b")) + .01 * cos(seq_len(nrow(d)))
  s <- spline_terms(y ~ u + kind, g)
  set.seed(8)
  f <- svem_spline(s, d, nBoot = 3L, store_member_weights = TRUE, lambda = c(Inf, 1, .1))
  absent <- grep("^nominal_kind_c", colnames(f$coef_matrix))
  expect_gt(length(absent), 0L)
  expect_lt(max(abs(f$coef_matrix[, absent])), 1e-8)
  x <- cbind(1, .svem_spline_design(s, d, f$blocks))
  z <- (d$y - f$response_center) / f$response_scale
  for (i in 1:3) {
    w <- f$member_weights$training[i, ]
    w <- w / sum(w)
    lam <- f$member_diagnostics$lambda[i]
    if (is.infinite(lam)) next
    lhs <- crossprod(x, x * w) + lam * diag(c(0, f$penalty))
    expect_equal(as.numeric(lhs %*% f$coef_matrix[i, ]), as.numeric(crossprod(x, w * z)), tolerance = 1e-8)
  }
})

test_that("limit specifications fit and predict with no saved numeric smooth", {
  # Guards the partial-match failure in which spec$smooth returned spec$smoother.
  line <- expand.grid(u = 0:7 / 7, kind = factor(c("a", "b")))
  nominal <- expand.grid(a = factor(c("p", "q")), b = factor(c("r", "s")))
  fixed <- data.frame(x = .15)
  cases <- list(
    list(spline_terms(y ~ u + kind, line, rank = 2L), transform(line, y = u + .1 * (kind == "b")), "linear limit"),
    list(spline_terms(y ~ a + b, nominal), transform(nominal[rep(1:4, 3), ], y = seq_len(12) %% 5), "categorical shrinkage"),
    list(spline_terms(y ~ x, fixed), data.frame(x = rep(.15, 6), y = 1:6), "intercept limit"))
  for (case in cases) {
    s <- case[[1L]]
    expect_identical(s$smoother, case[[3L]])
    expect_null(s[["smooth"]])
    set.seed(5)
    f <- svem_spline(s, case[[2L]], nBoot = 2L)
    expect_true(all(is.finite(predict(f, case[[2L]]))))
  }
})

# Independent review regression
test_that("rounding-only response variation reports no_signal", {
  skip_if_not_installed("mgcv")
  g <- expand.grid(u = 0:5 / 5, kind = factor(c("a", "b")))
  d <- g
  d$y <- .3
  d$y[seq(1L, nrow(d), by = 2L)] <- .1 + .2
  expect_false(all(d$y == d$y[1L]))
  set.seed(14)
  f <- svem_spline(spline_terms(y ~ u + kind, g), d, nBoot = 5L)
  expect_true(all(f$mean_coefficients[-1L] == 0))
  expect_true(all(f$member_diagnostics$intercept_only))
  expect_equal(predict(f, g), rep(.3, nrow(g)), tolerance = 1e-12)
  des <- function(x, p) pmin(1, pmax(0, x))
  out <- .svem_spline_scores(data.frame(id = seq_len(nrow(g))), g, list(y = f), list(y = list()),
                             des, c(y = 1), "geom", 1e-12)
  expect_identical(out$status$status, "no_signal")
  # A real signal at a large offset must stay nonconstant.
  d$z <- 1e6 + d$u
  set.seed(14)
  h <- svem_spline(spline_terms(y ~ u + kind, g), d, response = "z", nBoot = 5L)
  expect_false(all(h$mean_coefficients[-1L] == 0))
})

# Independent review regression
test_that("rounding noise in a fixed coordinate creates no extra sites", {
  skip_if_not_installed("mgcv")
  clean <- data.frame(u = rep(0:4 / 4, each = 8L), fixed = .1)
  noisy <- clean
  noisy$fixed <- .1 + sin(seq_len(nrow(noisy))) * .Machine$double.eps / 8
  expect_gt(nrow(unique(noisy)), nrow(unique(clean)))
  a <- spline_terms(y ~ u + fixed, clean)
  b <- spline_terms(y ~ u + fixed, noisy)
  expect_identical(b$fixed, "fixed")
  expect_identical(b$rank, a$rank)
  xa <- .svem_spline_surface(a, clean)
  xb <- .svem_spline_surface(b, clean)
  expect_equal(xa %*% (t(xa) / a$penalty), xb %*% (t(xb) / b$penalty), tolerance = 1e-7)
})

# Independent review regression
test_that("a numeric-only specification makes no nominal modelling claim", {
  skip_if_not_installed("mgcv")
  set.seed(4)
  g <- as.data.frame(matrix(runif(400), 100L, 4L))
  names(g) <- letters[1:4]
  out <- capture.output(print(spline_terms(y ~ a + b + c + d, g)))
  expect_false(any(grepl("level shifts|nominal cells", out)))
})

# Independent review regression
test_that("weak support, reduced bases and training extrapolation are disclosed without changing status", {
  skip_if_not_installed("mgcv")
  des <- function(x, p) pmin(1, pmax(0, x / 5))
  score <- function(f, pool) .svem_spline_scores(data.frame(id = seq_len(nrow(pool))), pool, list(y = f),
    list(y = list()), des, c(y = 1), "geom", 1e-12)$status
  g <- data.frame(u = seq(.4, .6, length.out = 21))
  wide <- data.frame(u = seq(0, 1, length.out = 16))
  wide$y <- sin(3 * wide$u)
  set.seed(6)
  f <- svem_spline(spline_terms(y ~ u, g, max_columns = 50L, bounds = rbind(c(u = 0), c(u = 1))), wide, nBoot = 3L)
  status <- score(f, g)
  expect_identical(status$status, "ok")
  expect_match(status$message, "outside the reference grid")

  line <- expand.grid(u = 0:7 / 7, kind = factor(c("a", "b")))
  d <- transform(line, y = u + .1 * (kind == "b") + .01 * sin(seq_along(u)))
  set.seed(6)
  limit <- svem_spline(spline_terms(y ~ u + kind, line, rank = 2L), d, nBoot = 3L)
  status <- score(limit, line)
  expect_identical(status$status, "ok")
  expect_match(status$message, "linear directions only")

  # Hand-built diagnostics exercise the threshold deterministically.
  weak <- limit
  weak$member_diagnostics$intercept_only <- c(TRUE, TRUE, FALSE)
  expect_match(.svem_spline_general_notes(weak), "67% of members")
  weak$member_diagnostics$intercept_only <- c(TRUE, FALSE, FALSE)
  expect_no_match(.svem_spline_general_notes(weak), "of members")
})

test_that("general controls and construction metadata are saved consistently", {
  defaults <- .svem_spline_general_defaults()
  for (fun in list(spline_terms,.svem_spline_general_terms,.svem_spline_spec_from_bigexp))
    for (v in names(defaults)) expect_identical(as.numeric(eval(formals(fun)[[v]])),defaults[[v]])
  g <- expand.grid(u=0:3/3,kind=factor(c("a","b")))
  s <- spline_terms(y~u+kind,g,rank=2,marginal_modes=0)
  expect_identical(s$controls$marginal_modes,0L)
  expect_equal(s$legal_tol,1e-8)
  expect_true(is.na(s$mgcv_version))
  expect_null(s$smooth)
  saved <- s;saved$legal_tol <- .01
  outside <- g[1,,drop=FALSE];outside$u <- -.005
  expect_true(.svem_spline_legal(saved,outside))
  expect_false(.svem_spline_legal(s,outside))
})
