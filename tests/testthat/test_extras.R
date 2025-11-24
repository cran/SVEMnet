test_that("SVEMnet predictions are reproducible given set.seed()", {
  skip_on_cran()

  set.seed(123)
  n  <- 40
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y  <- 1 + 2 * x1 - 0.5 * x2 + rnorm(n)
  dat <- data.frame(y, x1, x2)

  set.seed(2024)
  fit1 <- SVEMnet(y ~ x1 + x2, dat, nBoot = 20)

  set.seed(2024)
  fit2 <- SVEMnet(y ~ x1 + x2, dat, nBoot = 20)

  p1 <- predict(fit1, newdata = dat)
  p2 <- predict(fit2, newdata = dat)

  expect_equal(p1, p2, tolerance = 1e-10)
})

test_that("svem_select_from_score_table handles min and top_type = 'n'", {
  set.seed(4)
  st <- data.frame(
    x = rnorm(50),
    score = rnorm(50)
  )

  sel <- svem_select_from_score_table(
    score_table = st,
    target      = "score",
    direction   = "min",
    k           = 3,
    top_type    = "n",
    top         = 10
  )

  best <- sel$best
  cand <- sel$candidates

  expect_equal(nrow(best), 1L)
  expect_equal(best$score, min(st$score, na.rm = TRUE))


  expect_lte(nrow(cand), 3L)
})


test_that("svem_export_candidates_csv combines selection objects correctly", {
  skip_on_cran()

  set.seed(5)
  n  <- 30
  x1 <- rnorm(n)
  y  <- x1 + rnorm(n)
  dat <- data.frame(y, x1)

  fit <- SVEMnet(y ~ x1, dat, nBoot = 10)
  objs  <- list(y = fit)
  goals <- list(y = list(goal = "max", weight = 1))

  scored <- svem_score_random(
    objects = objs,
    goals   = goals,
    n       = 200,
    verbose = FALSE
  )

  sel1 <- svem_select_from_score_table(
    score_table = scored$score_table,
    target      = "score",
    direction   = "max",
    k           = 3,
    top_type    = "frac",
    top         = 0.2,
    label       = "round1"
  )

  sel2 <- svem_select_from_score_table(
    score_table = scored$score_table,
    target      = "uncertainty_measure",
    direction   = "max",
    k           = 2,
    top_type    = "frac",
    top         = 0.2,
    label       = "round1_explore"
  )

  # As separate args
  out1 <- svem_export_candidates_csv(sel1, sel2, write_file = FALSE)
  # As a list
  out2 <- svem_export_candidates_csv(list(sel1, sel2), write_file = FALSE)

  expect_s3_class(out1, "data.frame")
  expect_true(all(c("candidate_type", "selection_label") %in% names(out1)))
  expect_true(all(out1$candidate_type %in% c("best", "medoid")))
  expect_equal(out1, out2)
})


test_that("svem_wmt_multi returns named p_values and multipliers", {
  skip_on_cran()

  data(lipid_screen)

  set.seed(1)
  spec <- bigexp_terms(
    Potency ~ PEG + Helper + Ionizable + Cholesterol +
      Ionizable_Lipid_Type + N_P_ratio + flow_rate,
    data             = lipid_screen,
    factorial_order  = 1,
    polynomial_order = 1
  )

  form_pot <- bigexp_formula(spec, "Potency")
  form_siz <- bigexp_formula(spec, "Size")

  mix <- list(list(
    vars  = c("PEG", "Helper", "Ionizable", "Cholesterol"),
    lower = c(0.01, 0.10, 0.10, 0.10),
    upper = c(0.05, 0.60, 0.60, 0.60),
    total = 1.0
  ))

  res <- svem_wmt_multi(
    formulas       = list(Potency = form_pot, Size = form_siz),
    data           = lipid_screen,
    mixture_groups = mix,
    plot           = FALSE,
    verbose        = FALSE
  )

  expect_true(all(c("Potency", "Size") %in% names(res$p_values)))
  expect_true(all(c("Potency", "Size") %in% names(res$multipliers)))
  expect_true(all(res$multipliers >= 0))
})

