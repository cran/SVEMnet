portfolio_fixture <- function() {
  tab <- data.frame(x = 0:19, day = factor(rep("A", 20)), recipe_id = letters[1:20],
    y_pred = c(100, 20, rep(0, 18)), y_gaupro_pred = c(0, 0, 30, rep(0, 17)),
    y_extra_trees_pred = c(rep(0, 3), 40, rep(0, 16)),
    y_forward_pred = c(rep(0, 4), 50, rep(0, 15)),
    score = c(1, .2, rep(0, 18)), gaupro_score = c(0, 0, .3, rep(0, 17)),
    extra_trees_score = c(rep(0, 3), .4, rep(0, 16)),
    forward_score = c(rep(0, 4), .5, rep(0, 15)),
    y_lwr = 0, y_upr = 1:20, y_ciw_w = (1:20) / 20,
    uncertainty_measure = (1:20) / 20,
    y_gp_latent_sd = c(rep(0, 5), 100, rep(0, 14)))
  list(tab = tab, history = data.frame(x = 0, day = "B"),
       keys = .svem_recipe_keys(tab, "x", c(x = 1)), increments = c(x = 1),
       schema = list(predictor_vars = c("x", "day"), blocking = "day"),
       engines = c("gaupro", "extra_trees", "forward"),
       goals = list(y = list(goal = "max")), weights = c(y = 1),
       uncertainty_weights = c(y = 1), anchors = list(y = c(1, 20)),
       gp_endpoint = "y")
}

portfolio_run <- function(fixture = portfolio_fixture(), budget = 5L,
                          incumbent = 1L, ...) {
  do.call(.svem_worklist, c(fixture, list(budget = budget, allocation = "portfolio",
                                        incumbent = incumbent), list(...)))
}

test_that("portfolio budgets 1 to 5 have explicit repeat and nested model roles", {
  for (budget in 1:5) {
    result <- portfolio_run(budget = budget)
    expect_identical(result$worklist$recipe_id, letters[seq_len(budget)])
    expect_identical(result$worklist$slot, seq_len(budget))
    expect_identical(result$worklist$new_recipe, c(FALSE, rep(TRUE, budget - 1L)))
    expect_identical(result$worklist$role,
                     c("incumbent_repeat", rep("optimum", budget - 1L)))
    expect_identical(result$worklist$nominated_by[1], "observed + svemnet")
    expect_equal(attr(result$worklist, "budget"), budget)
    expect_equal(attr(result$worklist, "shortfall"), 0)
    expect_equal(result$shortfall, 0)
    expect_equal(result$allocation_details$selected_count, 0L)
    expect_equal(result$nominations$slot, c(seq_len(budget), rep(NA_integer_, 5L - budget)))
    expect_identical(result$nominations$endorses_incumbent,
                     c(FALSE, TRUE, FALSE, FALSE, FALSE))
  }
  result <- portfolio_run(budget = 6L)
  expect_identical(result$worklist$role[6], "high_spread_medoid")
  expect_false(result$worklist$recipe_id[6] == "f")
  expect_equal(result$allocation_details$eligible_count, 15)
  expect_equal(result$allocation_details$shortlist_count, 2)
  expect_equal(result$allocation_details$effective_fraction, 2 / 15)
  expect_equal(result$worklist$nearest_distance[1], 0)
})

test_that("shared optima retain endorsements even after capacity is full", {
  fixture <- portfolio_fixture()
  fixture$tab$y_gaupro_pred <- fixture$tab$y_extra_trees_pred <-
    fixture$tab$y_forward_pred <- fixture$tab$y_pred
  tight <- portfolio_run(fixture, budget = 2L)
  expect_identical(tight$worklist$recipe_id, c("a", "b"))
  expect_identical(tight$worklist$nominated_by,
    c("observed + svemnet + gaupro + extra_trees + forward",
      "svemnet + gaupro + extra_trees + forward"))
  expect_identical(tight$nominations$status, rep("selected", 5))
  expect_identical(tight$nominations$slot, c(1L, rep(2L, 4)))
  filled <- portfolio_run(fixture, budget = 10L, exploration_fraction = .01)
  expect_equal(nrow(filled$worklist), 10)
  expect_equal(filled$allocation_details$eligible_count, 18)
  expect_equal(filled$allocation_details$shortlist_count, 8)
  expect_equal(filled$allocation_details$selected_count, 8)
  expect_equal(filled$allocation_details$effective_fraction, 8 / 18)
  expect_setequal(filled$worklist$recipe_id[-(1:2)], letters[13:20])
  expect_false(anyDuplicated(filled$worklist$recipe_id) > 0)
})

test_that("portfolio ordering, optional incumbent, and exclusions are inspectable", {
  changed <- portfolio_run(budget = 2, priority = c("incumbent", "forward", "gaupro",
                                                   "extra_trees", "svemnet"))
  expect_identical(changed$worklist$recipe_id, c("a", "e"))
  omitted <- portfolio_run(budget = 1, incumbent = NULL)
  expect_identical(omitted$worklist$recipe_id, "b")
  expect_true(omitted$worklist$new_recipe)
  expect_identical(omitted$nominations$status[1], "not_supplied")
  subset <- portfolio_run(budget = 2, priority = c("incumbent", "forward"))
  expect_identical(subset$worklist$recipe_id, c("a", "e"))
  expect_identical(subset$nominations$status[2:4], rep("not_prioritized", 3))
  implicit <- portfolio_run(budget = NULL)
  expect_equal(nrow(implicit$worklist), 5L)

  fixture <- portfolio_fixture()
  fixture$tab <- fixture$tab[20:1, ]
  fixture$keys <- .svem_recipe_keys(fixture$tab, "x", c(x = 1))
  shuffled <- portfolio_run(fixture, budget = 7, incumbent = 20L)
  expect_identical(shuffled$worklist$recipe_id, portfolio_run(budget = 7)$worklist$recipe_id)
})

test_that("forward uses the common raw or joint scientific objective", {
  fixture <- portfolio_fixture()
  fixture$goals$y$goal <- "min"
  for (column in c("y_pred", "y_gaupro_pred", "y_extra_trees_pred", "y_forward_pred")) {
    fixture$tab[[column]] <- -fixture$tab[[column]]
  }
  expect_identical(portfolio_run(fixture)$worklist$recipe_id, letters[1:5])
  fixture <- portfolio_fixture()
  fixture$goals$z <- list(goal = "min")
  fixture$weights <- c(y = .6, z = .4)
  fixture$tab$z_pred <- 1:20
  fixture$tab$z_lwr <- fixture$tab$z_upr <- fixture$tab$z_ciw_w <- 0
  fixture$tab$forward_score <- c(rep(0, 5), 1, rep(0, 14))
  joint <- portfolio_run(fixture)
  expect_identical(joint$worklist$recipe_id, c("a", "b", "c", "d", "f"))
})

test_that("unavailable scores and exploration produce honest shortfalls", {
  fixture <- portfolio_fixture()
  for (column in c("y_pred", "y_gaupro_pred", "y_extra_trees_pred", "y_forward_pred")) {
    fixture$tab[[column]] <- NA_real_
  }
  fixture$anchors <- list(y = NULL)
  expect_warning(result <- portfolio_run(fixture, budget = 3), "short by 2")
  expect_identical(result$worklist$recipe_id, "a")
  expect_equal(result$shortfall, 2)
  expect_equal(attr(result$worklist, "shortfall"), 2)
  expect_true(all(result$nominations$status[-1] == "no_finite_score"))

  fixture <- portfolio_fixture()
  fixture$history <- fixture$tab[c("x", "day")]
  expect_warning(result <- portfolio_run(fixture, budget = 3), "short by 2")
  expect_identical(result$worklist$recipe_id, "a")
  expect_true(all(result$nominations$status[-1] == "no_unprepared_recipe"))
  fixture <- portfolio_fixture()
  fixture$tab$uncertainty_measure <- NA_real_
  expect_warning(result <- portfolio_run(fixture, budget = 6), "short by 1")
  expect_equal(result$allocation_details$eligible_count, 0)
  expect_equal(result$allocation_details$selected_count, 0)
})

test_that("portfolio cutoff ties and invalid controls are explicit", {
  fixture <- portfolio_fixture()
  fixture$tab$uncertainty_measure[6:20] <- 1
  result <- portfolio_run(fixture, budget = 7, exploration_fraction = .25)
  expect_equal(result$allocation_details$shortlist_count, 4)
  expect_equal(result$allocation_details$cutoff_ties, 15)
  expect_equal(result$allocation_details$cutoff_ties_retained, 4)
  expect_error(portfolio_run(priority = c("svemnet", "svemnet")), "distinct supported")
  expect_error(portfolio_run(priority = c("incumbent", "unknown")), "distinct supported")
  expect_error(portfolio_run(exploration_fraction = 0), "exploration_fraction")
  expect_error(portfolio_run(exploration_fraction = 1.1), "exploration_fraction")
  expect_error(portfolio_run(incumbent = 1.5), "incumbent")
  expect_error(portfolio_run(incumbent = 21), "incumbent")
})

test_that("the legacy worklist path keeps its return shape and priority", {
  fixture <- portfolio_fixture()
  fixture$engines <- c("gaupro", "extra_trees")
  implicit <- do.call(.svem_worklist, c(fixture, list(budget = 5L)))
  explicit <- do.call(.svem_worklist, c(fixture, list(budget = 5L, allocation = "legacy")))
  expect_identical(implicit, explicit)
  expect_identical(implicit$worklist$recipe_id, c("b", "c", "d", "t", "f"))
  expect_identical(implicit$worklist$role,
                   c(rep("optimum", 3), "maximum_spread", "latent_sd"))
  expect_false("new_recipe" %in% names(implicit$worklist))
  expect_null(attr(implicit$worklist, "budget"))
  expect_null(attr(implicit$worklist, "shortfall"))
  expect_false("allocation_details" %in% names(implicit))
  expect_false("endorses_incumbent" %in% names(implicit$nominations))
})

test_that("portfolio uncertainty requires finite increasing anchors", {
  fixture <- portfolio_fixture()
  for (anchors in list(c(1, 1), c(1, Inf), c(1, NA_real_), c(2, 1))) {
    fixture$anchors$y <- anchors
    expect_warning(result <- portfolio_run(fixture, budget = 6), "short by 1")
    expect_equal(result$shortfall, 1)
    expect_equal(result$allocation_details$eligible_count, 0)
  }
  fixture$anchors$y <- c(0, 1e-200)
  fixture$tab$y_upr <- (1:20) * 1e-200
  tiny <- portfolio_run(fixture, budget = 6)
  expect_identical(tiny$worklist$recipe_id, portfolio_run(budget = 6)$worklist$recipe_id)
  expect_equal(tiny$shortfall, 0)
})

test_that("portfolio-only metadata cannot overwrite an input predictor", {
  fixture <- portfolio_fixture()
  fixture$schema$predictor_vars <- c("x", "day", "new_recipe")
  fixture$tab$new_recipe <- seq_len(20)
  fixture$history$new_recipe <- 0
  expect_error(portfolio_run(fixture), "Reserved portfolio worklist name: new_recipe")
  legacy <- suppressWarnings(do.call(.svem_worklist, c(fixture, list(budget = 2))))
  expect_true("new_recipe" %in% names(legacy$worklist))
  fixture <- portfolio_fixture()
  fixture["history"] <- list(NULL)
  result <- portfolio_run(fixture, budget = 2)
  expect_true(result$nominations$best_prepared[result$nominations$engine == "svemnet"])
  expect_identical(result$worklist$new_recipe, c(FALSE, TRUE))
})
