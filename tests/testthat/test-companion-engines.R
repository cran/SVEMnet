test_that("simplex groups are declared once and the original fit stays usable", {
  d <- data.frame(a = seq(.1, .8, length.out = 12))
  d$b <- 1 - d$a
  d$y <- d$a + sin(4 * d$a)
  g <- list(list(vars = c("a", "b")))
  spec <- bigexp_terms(y ~ a + b, d, mixture_groups = g, report = FALSE)
  expect_equal(spec$settings$mixture_vars, c("a", "b"))
  set.seed(9)
  fit <- SVEMnet(spec, d, nBoot = 2)
  expect_identical(fit$sampling_schema$mixture_groups, spec$settings$mixture_groups)
  expect_null(fit$companions)
  expect_equal(predict(fit, d, engine = "all")$svemnet, as.numeric(predict(fit, d)))
  expect_error(SVEMnet(spec, d, mixture_groups = list(list(vars = c("b","a")))),
               "conflicts")
  expect_error(bigexp_terms(y ~ a + b, d, mixture_groups = g, mixture_vars = "a"),
               "same components")
  expect_error(SVEMnet(spec, transform(d, b = b + .01)), "must sum")
  expect_error(SVEMnet(spec, d, companions = "forest"), "companions must")
})

test_that("unsupported companion designs preserve the fitted SVEM branch", {
  d <- data.frame(x = seq_len(18), day = factor(rep(c("A", "B", "C"), 6)))
  d$y <- seq_len(18) + cos(d$x)
  spec <- bigexp_terms(y ~ x, d, blocking = "day", report = FALSE)
  set.seed(7)
  base <- SVEMnet(spec, d, nBoot = 2)
  after <- .Random.seed
  set.seed(7)
  expect_warning(fit <- SVEMnet(spec, d, nBoot = 2, companions = "gaupro"),
                 "requires companion_blocks")
  expect_identical(.Random.seed, after)
  expect_identical(fit$parms, base$parms)
  expect_identical(fit$companions$gaupro$status, "unsupported")
  d$z <- d$y
  set.seed(7)
  expect_warning(other <- SVEMnet(spec, d, response = "z", nBoot = 2,
                                  companions = "gaupro"), "requires companion_blocks")
  expect_false(identical(other$companions$gaupro$seed, fit$companions$gaupro$seed))
  expect_error(predict(fit, d, engine = "gaupro"), "unsupported")
  all <- predict(fit, d, engine = "all")
  expect_equal(all$svemnet, as.numeric(predict(base, d)))
  expect_true(all(is.na(all$gaupro)))
  expect_match(attr(all, "engine_status")[["gaupro"]], "unsupported")
  expect_error(predict(base, d, engine = "extra_trees"), "not requested")
  expect_warning(transformed <- SVEMnet(log(y) ~ x, d, nBoot = 2,
                                         companions = "extra_trees"), "untransformed")
  expect_identical(transformed$companions$extra_trees$status, "unsupported")
})

test_that("worklists honor nominal numeric block labels and preserve legacy midpoint sampling", {
  d <- expand.grid(x = seq(0, 1, length.out = 8), day = 1:2)
  d$y <- d$x^2 + d$day + sin(d$x * 5) / 7
  spec <- bigexp_terms(y ~ x, d, blocking = "day", report = FALSE)
  set.seed(92)
  fit <- suppressWarnings(SVEMnet(spec, d, nBoot = 3, companions = "extra_trees",
                                 companion_blocks = "interacting"))
  goals <- list(y = list(goal = "max", weight = 1))
  new <- suppressWarnings(svem_score_random(fit, goals, n = 20, budget = 1, verbose = FALSE))
  expect_equal(unique(new$score_table$day), 1)
  expect_equal(new$score_table$y_pred, as.numeric(predict(fit, new$score_table)))
  if (fit$companions$extra_trees$status == "ok") expect_true(all(is.finite(new$score_table$y_extra_trees_pred)))
  expect_error(svem_score_random(fit, goals, candidates = data.frame(x = .5, day = 1.5),
                                 verbose = FALSE), "common observed block level")
  fit$companions <- NULL
  legacy <- svem_score_random(fit, goals, n = 20, verbose = FALSE)
  expect_equal(unique(legacy$score_table$day), 1.5)
})

test_that("public companions agree with direct fits on the shared raw rows", {
  skip_if_not_installed("GauPro", minimum_version = "0.2.17")
  skip_if_not_installed("ranger", minimum_version = "0.18.0")
  d <- data.frame(a = rep(seq(.1, .4, length.out = 8), 3),
                  day = factor(rep(c("C","A","B"), each = 8)))
  d$b <- 1 - d$a
  d$y <- 7 * d$a + as.numeric(d$day) + sin(seq_len(nrow(d)))
  d$y[3] <- NA
  g <- list(list(vars = c("a", "b"), lower = c(.1,.6), upper = c(.4,.9)))
  spec <- bigexp_terms(y ~ a + b, d, blocking = "day", mixture_groups = g,
                       report = FALSE)
  set.seed(77)
  base <- SVEMnet(spec, d, nBoot = 2)
  after <- .Random.seed
  set.seed(77)
  fit <- SVEMnet(spec, d, nBoot = 2, companions = TRUE,
                 companion_blocks = "interacting")
  expect_identical(.Random.seed, after)
  expect_identical(fit$parms, base$parms)
  expect_equal(vapply(fit$companions, `[[`, "", "status"),
               c(gaupro = "ok", extra_trees = "ok"))
  expect_identical(attr(fit$companions, "input_map")$omitted, "b")
  prepared <- attr(fit$companions, "prepared_recipes")
  expect_identical(names(prepared), fit$sampling_schema$predictor_vars)
  expect_equal(nrow(prepared), sum(!is.na(d$y)))
  keep <- !is.na(d$y)
  scaled <- function(a) (a - .1) / (.4 - .1)
  train <- data.frame(x1 = scaled(d$a[keep]), x2 = d$day[keep], response = d$y[keep])
  row <- d[c(7, 9, 18, 11), ]
  row$day <- factor(as.character(row$day), levels = rev(levels(d$day)))
  newx <- data.frame(x1 = scaled(row$a), x2 = factor(as.character(row$day), levels = levels(d$day)))
  set.seed(fit$companions$gaupro$seed)
  gp <- GauPro::gpkm(formula = response ~ x1 + x2, data = train, nug.est = TRUE,
                     restarts = 0, normalize = TRUE, parallel = FALSE, verbose = 0)
  expect_identical(fit$companions$gaupro$settings$kernel, "matern52")
  expect_identical(fit$companions$gaupro$settings$fitted_kernel_class, class(gp$kernel)[1L])
  expect_equal(fit$companions$gaupro$deviance, as.numeric(gp$deviance()), tolerance = 1e-7)
  direct <- gp$predict(newx, se.fit = TRUE, mean_dist = TRUE)
  reply <- predict(fit, row, engine = "gaupro", se.fit = TRUE)
  expect_equal(reply$fit, as.numeric(direct$mean), tolerance = 1e-7)
  expect_equal(reply$se.fit, as.numeric(direct$se), tolerance = 1e-7)
  rf <- ranger::ranger(x = train[c("x1","x2")], y = train$response,
                       splitrule = "extratrees", num.trees = 128L, mtry = 2L,
                       num.random.splits = 1L, min.node.size = 1L, min.bucket = 2L,
                       replace = FALSE, sample.fraction = 1,
                       respect.unordered.factors = "partition", oob.error = FALSE,
                       num.threads = 1L, seed = fit$companions$extra_trees$seed)
  all <- predict(fit, row, engine = "all")
  expect_named(all, c("svemnet", "gaupro", "extra_trees"))
  expect_equal(all$extra_trees, predict(rf, newx, num.threads = 1L)$predictions)
  set.seed(77)
  forest_only <- SVEMnet(spec, d, nBoot = 2, companions = "extra_trees",
                         companion_blocks = "interacting")
  expect_identical(forest_only$companions$extra_trees$seed, fit$companions$extra_trees$seed)
  expect_equal(predict(forest_only, row, engine = "extra_trees"), all$extra_trees)
  expect_equal(predict(fit, droplevels(row[1, ]), engine = "gaupro"), all$gaupro[1])
  bad <- row
  bad$day <- as.character(bad$day)
  bad$day[2] <- "new day"
  expect_warning(part <- predict(fit, bad, engine = "gaupro"), "unsupported")
  expect_true(is.na(part[2]))
  expect_equal(part[-2], all$gaupro[-2])
  expect_equal(predict(unserialize(serialize(fit, NULL)), row, engine = "all"), all)
  reopened <- unserialize(serialize(fit, NULL))
  reopened$companions$gaupro$version <- basename(tempfile("previous-"))
  expect_warning(p <- predict(reopened, row, engine = "gaupro"), "fitted with version")
  expect_equal(p, all$gaupro)
  expect_silent(predict(reopened, row, engine = "gaupro"))
  reopened$companions$gaupro$package <- "absentCompanionDependencyForTest"
  expect_warning(partial <- predict(reopened, row, engine = "all"), "install.packages")
  expect_match(attr(partial, "engine_status")[["gaupro"]], "prediction_failed")
  expect_true(all(is.na(partial$gaupro)))
  expect_identical(partial$svemnet, all$svemnet)
  expect_identical(partial$extra_trees, all$extra_trees)
  expect_error(predict(fit, row, engine = "extra_trees", se.fit = TRUE), "only gaupro")
  expect_error(predict(fit, transform(row, b = b - .01), engine = "gaupro"), "must sum")
})

test_that("declared block representations preserve exact IDs and continuous covariates", {
  id <- c(1 + .Machine$double.eps, 1 + 2 * .Machine$double.eps)
  expect_identical(as.character(id[1L]), as.character(id[2L]))
  d <- data.frame(a = rep(seq(.1, .4, length.out = 8), 2), day = rep(id, each = 8),
                  temp = seq(20, 28, length.out = 16), plate = 3.5, fixed = .1)
  d$b <- .9 - d$a
  d$y <- 4 * d$a + sqrt(d$temp) + sin(seq_len(nrow(d)))
  g <- list(list(vars = c("a", "b", "fixed"), lower = c(.1, .5, .1),
                 upper = c(.4, .8, .1)))
  spec <- bigexp_terms(y ~ a + b + fixed, d, blocking = c("day", "temp", "plate"),
                       mixture_groups = g, factorial_order = 1, polynomial_order = 1,
                       include_pc_2way = FALSE, report = FALSE)
  set.seed(103)
  fit <- suppressWarnings(SVEMnet(spec, d, nBoot = 2, companions = "extra_trees",
                                  companion_blocks = "interacting",
                                  companion_continuous_blocks = "temp"))
  map <- attr(fit$companions, "input_map")
  expect_identical(map$levels$day, id)
  expect_named(map$levels, c("day", "plate"), ignore.order = TRUE)
  expect_identical(map$continuous_blocks, "temp")
  expect_identical(map$omitted, "b")
  expect_false(any(c("fixed", "plate", "b") %in% map$active))
  expect_true("plate" %in% names(map$constants))
  expect_identical(map$bounds$temp, c(20, 28))
  both <- d[c(1, 9), ]
  both$temp <- 23.25
  encoded <- SVEMnet:::.svem_companion_inputs(map, both, "error")$x
  expect_equal(as.integer(encoded[[match("day", map$active)]]), 1:2)
  expect_equal(encoded[[match("temp", map$active)]], rep(3.25 / 8, 2))
  expect_error(SVEMnet:::.svem_companion_inputs(map, transform(both, plate = 4), "error"),
               "Unsupported")
  expect_equal(SVEMnet:::.svem_worklist_numeric_blocks(list(fit), both, supplied = TRUE), both)
  expect_error(SVEMnet:::.svem_worklist_numeric_blocks(list(fit), transform(both, day = 1),
                                                        supplied = TRUE), "observed block")
  if (fit$companions$extra_trees$status == "ok") {
    goals <- list(y = list(goal = "max", weight = 1))
    scored <- svem_score_random(fit, goals, candidates = both[1, ], verbose = FALSE)
    expect_identical(scored$score_table$day, id[1L])
    expect_identical(scored$score_table$temp, 23.25)
    expect_true(is.finite(scored$score_table$y_extra_trees_pred))
  }
  expect_warning(SVEMnet:::.svem_companion_map(fit, d, "interacting", TRUE),
                 "unreplicated numeric block.*temp")
  expect_error(SVEMnet:::.svem_companion_map(fit, d, "interacting", TRUE, "a"),
               "continuous numeric blocks")
})
