test_that("worklist slots use raw optima, shared endorsements, and explicit tie policies", {
  tab <- data.frame(x = 0:6, recipe_id = letters[1:7], y_pred = c(100, 98, 96, 95, 91, 93, 90),
                    y_gaupro_pred = c(101, 99, 95, 92, 90, 93, 81),
                    y_extra_trees_pred = c(80, 81, 99, 95, 92, 93, 90),
                    score = c(1, 1, 1, .99, .9, .95, .8),
                    y_lwr = 0, y_upr = c(2, 3, 4, 6, 14, 9, 2),
                    uncertainty_measure = c(.2, .4, .5, 1, 1, 1, .1),
                    y_gp_latent_sd = c(0, 1, 2, 10, 1, 9, 4))
  run <- function(tab, budget, history = data.frame(x = 0)) {
    .svem_worklist(tab, history, .svem_recipe_keys(tab, "x", c(x = 1)), c(x = 1),
                   list(predictor_vars = "x"), c("gaupro", "extra_trees"),
                   list(y = list(goal = "max")), c(y = 1), c(y = 1),
                   list(y = c(1, 6)), "y", budget)
  }
  out <- run(tab, 4)
  expect_equal(out$worklist$recipe_id, c("b", "c", "d", "f"))
  expect_equal(out$worklist$nominated_by[1], "svemnet + gaupro")
  expect_equal(out$worklist$role, c("optimum", "optimum", "best_predicted_tie", "latent_sd"))
  expect_equal(out$worklist$nearest_distance[1], 1/6)
  expect_equal(out$nominations$best_prepared, c(TRUE, TRUE, FALSE))
  tight <- run(tab, 1)
  expect_equal(tight$nominations$status, c("selected", "selected", "dropped_by_budget"))
  expect_equal(tight$nominations$slot, c(1L, 1L, NA_integer_))
  unique_max <- tab
  unique_max$uncertainty_measure[c(5, 6)] <- .99
  expect_equal(run(unique_max, 3)$worklist$role[3], "maximum_spread")
  expect_equal(out$exploration$recipe_id, c("d", "e"))
  expect_equal(out$exploration$tie_count, c(3L, 3L))
  expect_equal(out$exploration$exact_ties, c(3L, 3L))
  expect_equal(out$tie_endpoint_ranges$max_width, 14)
  expect_equal(out$exploration$min_desirability, rep(.9, 2L))
  expect_equal(run(tab[nrow(tab):1, ], 4)$worklist$recipe_id, out$worklist$recipe_id)
  for (n in 1:5) expect_equal(nrow(run(tab, n)$worklist), n)
  csv <- svem_export_candidates_csv(out$selections, write_file = FALSE)
  expect_equal(csv$recipe_id, out$worklist$recipe_id)
  expect_warning(no_new <- run(tab, 1, tab["x"]), "short by 1")
  expect_equal(nrow(no_new$worklist), 0L)
  expect_true(all(no_new$nominations$status == "no_unprepared_recipe"))
  # A previously prepared recipe is still the same physical preparation
  # when the operator/plate changes; a cross-block repeat needs a dedicated slot.
  tab$operator <- "A"
  prior <- data.frame(x = 0, operator = "B")
  index <- .svem_pool_index(tab[c("x", "operator")], tab, c(x = 1), "operator")
  blocked <- .svem_worklist(tab, prior, index$keys, c(x = 1),
                            list(predictor_vars = c("operator", "x"), blocking = "operator"),
                            character(), list(y = list(goal = "max")), c(y = 1), c(y = 1),
                            list(y = c(1, 6)), NULL, 1)
  expect_true(blocked$nominations$best_prepared)
  expect_identical(blocked$worklist$recipe_id, "b")
})

test_that("constrained supplied pools are checked and support physical IDs", {
  dat <- data.frame(x = rep(0:4 / 4, 3), a = rep(c(.25, .5, .75), each = 5))
  dat$b <- 1 - dat$a
  dat$y <- dat$x + 3 * dat$a + .01 * cos(seq_len(nrow(dat)))
  set.seed(55)
  fit <- SVEMnet(y ~ a + b + x, dat, nBoot = 4, glmnet_alpha = 1)
  groups <- list(list(vars = c("a", "b"), lower = c(0, 0), upper = c(1, 1), total = 1))
  pool <- data.frame(a = c(0, .25, .50, 1, .25), x = c(.25, .5, .75, .5, .5))
  pool$b <- 1 - pool$a
  pool$recipe_id <- c("A", "B", "C", "D", "E")
  fixed <- list(y = list(goal = "max", weight = 1, lower_acceptable = -10, upper_acceptable = -9))
  out <- svem_score_random(fit, fixed, candidates = pool, mixture_groups = groups,
                           increments = c(a = .25, b = .25), budget = 1L, verbose = FALSE)
  expect_equal(out$score_table$recipe_id, LETTERS[1:4])
  expect_equal(out$score_table$score, rep(1, 4))
  expect_equal(out$worklist$recipe_id, out$score_table$recipe_id[which.max(out$score_table$y_pred)])
  flipped <- svem_score_random(fit, fixed, candidates = pool[5:1, ], mixture_groups = groups,
                               budget = 1, verbose = FALSE)
  expect_identical(flipped$worklist$recipe_id, out$worklist$recipe_id)
  expect_identical(out$history_source, "no_preparation_history")
  expect_true(is.na(out$nominations$best_prepared))
  expect_error(svem_score_random(fit, fixed, candidates = transform(pool, b = .5),
                                 mixture_groups = groups, verbose = FALSE), "mixture bounds or total")
  expect_error(svem_score_random(fit, fixed, candidates = transform(pool, x = 100),
                                 mixture_groups = groups, verbose = FALSE), "outside its fitted")
  expect_error(svem_score_random(fit, fixed, candidates = pool, mixture_groups = groups,
                                 increments = c(a = .2), verbose = FALSE), "increments")
})

test_that("paired engines score the same raw recipes, anchors and weights", {
  skip_if_not_installed("GauPro")
  skip_if_not_installed("ranger")
  dat <- expand.grid(a = seq(.2, .8, length.out = 5), x = c(-1, 0, 1))
  dat$b <- 1 - dat$a
  dat$y <- 3 * dat$a - dat$x^2 + sin(seq_len(nrow(dat))) / 20
  dat$z <- dat$x - dat$a^2 + cos(seq_len(nrow(dat))) / 10
  groups <- list(list(vars = c("a", "b"), lower = c(.1, .1), upper = c(.9, .9), total = 1))
  set.seed(123)
  fits <- lapply(c("y", "z"), function(resp) {
    SVEMnet(stats::reformulate(c("a", "b", "x"), response = resp), dat,
            nBoot = 4, glmnet_alpha = 1, companions = TRUE, mixture_groups = groups)
  })
  goals <- list(y = list(goal = "max", weight = 2, lower_acceptable = -2, upper_acceptable = 3),
                z = list(goal = "min", weight = 1, lower_acceptable = -2, upper_acceptable = 2))
  x <- data.frame(a = seq(.25, .75, length.out = 8), x = seq(-.6, .7, length.out = 8))
  x$b <- 1 - x$a
  out <- suppressWarnings(svem_score_random(fits, goals, candidates = x, gp_endpoint = "y",
                                            budget = 3, verbose = FALSE))
  tab <- out$score_table
  for (engine in c("gaupro", "extra_trees")) {
    mu_y <- predict(fits[[1L]], x, engine = engine)
    mu_z <- predict(fits[[2L]], x, engine = engine)
    expect_equal(tab[[paste0("y_", engine, "_pred")]], as.numeric(mu_y))
    expect_equal(tab[[paste0("z_", engine, "_pred")]], as.numeric(mu_z))
    pos <- 1e-6 + (1-1e-6) * pmin(1, pmax(0, (mu_y + 2)/5))
    neg <- 1e-6 + (1-1e-6) * pmin(1, pmax(0, (2 - mu_z)/4))
    expect_equal(tab[[paste0(engine, "_score")]], as.numeric((pos^2 * neg)^(1/3)))
  }
  expect_equal(tab$y_gp_latent_sd, as.numeric(predict(fits[[1]], x,
                                     engine = "gaupro", se.fit = TRUE)$se.fit))
  prepared_pool <- dat[c("a", "b", "x")]
  prepared_pool$recipe_id <- paste0("old", seq_len(nrow(dat)))
  prepared_pool <- rbind(prepared_pool, data.frame(a = .21, b = .79, x = .35, recipe_id = "new"))
  check_history <- suppressWarnings(svem_score_random(fits[1], goals["y"], candidates = prepared_pool,
                                                 budget = 1L, verbose = FALSE))
  expect_equal(check_history$history_source, "fitted_preparations")
  expect_equal(check_history$worklist$recipe_id, "new")
  interrupted <- fits[1]
  interrupted[[1]]$companions$gaupro$package <- "missingDependencyForTest"
  expect_warning(partial <- svem_score_random(interrupted, goals["y"], candidates = x,
                                               budget = 1L, verbose = FALSE),
                 "gaupro could not score")
  expect_identical(partial$engine_status$status[partial$engine_status$engine == "gaupro"],
                   "prediction_failed")
  expect_true(all(is.na(partial$score_table$y_gaupro_pred)))
  expect_true(all(is.finite(partial$score_table$y_extra_trees_pred)))

  # A missing required response drops just that engine's default optimum
  # slot. A declared budget can still fill the remaining slots with exploration.
  missing <- fits
  missing[[2L]]$companions$extra_trees <- NULL
  score_missing <- function(goals, budget = NULL) {
    svem_score_random(missing, goals, candidates = x, budget = budget, verbose = FALSE)
  }
  expect_warning(automatic <- score_missing(goals), "extra_trees.*z \\(not_requested\\)")
  expect_equal(automatic$nominations$status, c("selected", "selected", "no_finite_score"))
  expect_equal(nrow(automatic$worklist) + automatic$shortfall, 2)
  expect_warning(fixed <- score_missing(goals, budget = 3), "extra_trees could not score")
  expect_equal(nrow(fixed$worklist), 3L)
  expect_true(any(fixed$worklist$role != "optimum"))
  inactive <- goals
  inactive$z$weight <- 0
  expect_no_warning(active_only <- score_missing(inactive))
  expect_true(all(active_only$nominations$status == "selected"))

  # Even with partially finite predictions on both responses, their usable
  # rows must overlap. Infinite means cannot saturate the desirability to 1.
  disjoint <- fits
  for (i in seq_along(disjoint)) {
    disjoint[[i]]$companions$gaupro$model <- list(predict = local({
      parity <- i %% 2L
      function(x) ifelse(seq_len(nrow(x)) %% 2L == parity, Inf, 1)
    }))
  }
  expect_warning(gaps <- svem_score_random(disjoint, goals, candidates = x, verbose = FALSE),
                 "gaupro.*y \\(partial coverage\\), z \\(partial coverage\\)")
  expect_true(all(is.na(gaps$score_table$gaupro_score)))
  expect_identical(gaps$nominations$status[gaps$nominations$engine == "gaupro"], "no_finite_score")
  expect_error(svem_score_random(fits, goals, candidates = transform(x, b = .15),
                                 verbose = FALSE), "mixture bounds or total")
  expect_error(svem_score_random(fits, goals, candidates = x, mixture_groups =
                 list(list(vars = c("b", "a"))), verbose = FALSE), "agree with the geometry")
})

test_that("a supplied pool can cover the shared range of different analysis sets", {
  dat <- expand.grid(x = 0:4, z = c(-1, 0, 1))
  dat$a <- dat$x - dat$z^2 + sin(seq_len(nrow(dat))) / 10
  dat$b <- dat$z + dat$x / 2
  set.seed(19)
  fits <- list(SVEMnet(a ~ x + z, dat, nBoot = 3, glmnet_alpha = 1),
               SVEMnet(b ~ x + z, dat[dat$x < 4, ], nBoot = 3, glmnet_alpha = 1))
  pool <- data.frame(x = c(.5, 1.5), z = c(0, 1))
  tab <- svem_random_table_multi(fits, candidates = pool)
  expect_equal(tab$data, pool)
  expect_true(all(is.finite(as.matrix(tab$pred))))
  expect_error(svem_random_table_multi(fits, candidates = transform(pool, x = 4)),
               "outside its fitted numeric range")
})
