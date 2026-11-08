.plan_guard_fixture <- function(plain = FALSE) {
  d <- data.frame(x = seq(-1, 1, length.out = 18),
                  g = factor(rep(c("A", "B"), 9)))
  d$y <- 2 + d$x + .05 * sin(seq_len(nrow(d)))
  d$z <- 3 - d$x^2 + .05 * cos(seq_len(nrow(d)))
  d$y[2] <- NA_real_
  d$z[3] <- NA_real_
  spec <- bigexp_terms(y ~ x + g, d, polynomial_order = 2,
                       factorial_order = 1, report = FALSE)
  pool <- expand.grid(x = seq(-1, 1, length.out = 31), g = levels(d$g))
  goals <- list(y = list(goal = "max", weight = .6, lower_acceptable = 0, upper_acceptable = 4),
                z = list(goal = "min", weight = .4, lower_acceptable = 1, upper_acceptable = 4))
  fitted <- withr::with_seed(941, {
    fit_one <- function(response, forward = FALSE) {
      formula <- if (plain) stats::as.formula(paste(response, "~ x + g + I(x^2)")) else spec
      args <- list(formula = formula, data = d, nBoot = 3)
      if (!plain) args$response <- response
      if (forward) do.call(svem_forward, args) else
        do.call(SVEMnet, c(args, list(glmnet_alpha = 1, relaxed = FALSE)))
    }
    list(fits = setNames(lapply(names(goals), fit_one), names(goals)),
          forward = setNames(lapply(names(goals), fit_one, forward = TRUE), names(goals)))
  })
  c(list(data = d, spec = spec, pool = pool, goals = goals), fitted)
}

.plan_guard_call <- function(f, ...) {
  svem_plan_runs(f$spec, f$data, f$goals, f$pool, fits = f$fits,
                 forward = f$forward, score_control = list(verbose = FALSE), ...)
}

test_that("invalid pools and empty reference contexts fail before fitting and RNG", {
  d <- data.frame(x = seq(-1, 1, length.out = 12),
                   day = factor(rep(c("A", "B"), 6)), y = seq_len(12))
  spec <- bigexp_terms(y ~ x, d, blocking = "day", report = FALSE)
  goals <- list(y = list(goal = "max", weight = 1))
  pool <- data.frame(x = seq(-1, 1, length.out = 11), day = "A")
  calls <- 0L
  forbidden <- function(...) { calls <<- calls + 1L; stop("unexpected fit") }
  local_mocked_bindings(svem_fit_responses = forbidden, svem_forward = forbidden,
                         .package = "SVEMnet")
  set.seed(942)
  before <- .Random.seed
  expect_error(svem_plan_runs(spec, d, goals, pool["day"]), "Missing required predictor")
  bad <- pool; bad$day <- "future"
  expect_error(svem_plan_runs(spec, d, goals, bad), "Unseen level.*future.*future day")
  bad <- pool; bad$x[1] <- NA_real_
  expect_error(svem_plan_runs(spec, d, goals, bad), "finite and nonmissing")
  expect_error(svem_plan_runs(spec, d, goals, pool, reference = list(day = "future")),
               "No historical preparations")
  expect_error(svem_plan_runs(spec, d, goals, pool, priority = "unrecognized"), "priority")
  expect_error(svem_plan_runs(spec, d, goals, pool, exploration_fraction = 0), "exploration_fraction")
  expect_error(svem_plan_runs(spec, d, goals, pool, incumbent_summary = "unrecognized"), "arg")
  expect_error(svem_plan_runs(spec, d, goals, pool,
    benchmark = list(missing_id = "control")), "benchmark must")
  d$preparation_id <- paste0("run", seq_len(nrow(d)))
  expect_error(svem_plan_runs(spec, d, goals, pool,
    benchmark = list(preparation_id = "absent")), "exactly one")
  d$preparation_id[2] <- d$preparation_id[1]
  expect_error(svem_plan_runs(spec, d, goals, pool,
    benchmark = list(preparation_id = "run1")), "exactly one.*found 2")
  expect_error(svem_plan_runs(spec, d, goals, pool, budget = 1,
    incumbent = 1L, benchmark = 2L), "budget is too small")
  expect_identical(calls, 0L)
  expect_identical(.Random.seed, before)
})

test_that("stored complete training rows detect stale responses and designs before other fits", {
  f <- .plan_guard_fixture()
  calls <- 0L
  forbidden <- function(...) { calls <<- calls + 1L; stop("unexpected fit") }
  local_mocked_bindings(svem_fit_responses = forbidden, svem_forward = forbidden,
                         .package = "SVEMnet")
  set.seed(943)
  before <- .Random.seed
  changed <- f$data; changed$y[5] <- changed$y[5] + .25
  expect_error(svem_plan_runs(f$spec, changed, f$goals, f$pool, fits = f$fits),
               "Cannot reuse supplied fits.*outcomes or design.*Refit")
  changed <- f$data; changed$x[5] <- changed$x[5] + .01
  expect_error(svem_plan_runs(f$spec, changed, f$goals, f$pool, forward = f$forward),
               "Cannot reuse supplied forward.*outcomes or design.*Refit")
  changed <- f$data; changed$y[5:6] <- rev(changed$y[5:6])
  expect_error(svem_plan_runs(f$spec, changed, f$goals, f$pool, fits = f$fits),
               "Cannot reuse supplied fits")
  changed <- rbind(f$data, f$data[5, , drop = FALSE])
  expect_error(svem_plan_runs(f$spec, changed, f$goals, f$pool, fits = f$fits),
               "complete training rows")
  changed <- f$data
  changed$g <- as.character(changed$g)
  added <- changed[5, , drop = FALSE]; added$g <- "C"
  changed <- rbind(changed, added)
  refreshed <- bigexp_terms(y ~ x + g, changed, polynomial_order = 2,
                            factorial_order = 1, report = FALSE)
  expect_error(svem_plan_runs(refreshed, changed, f$goals, f$pool, fits = f$fits),
               "Cannot reuse supplied fits.*New measured factor levels")
  stale <- f$fits; stale$y$actual_y <- NULL
  expect_error(svem_plan_runs(f$spec, f$data, f$goals, f$pool, fits = stale),
               "Stored training rows.*Older objects")
  expect_error(svem_plan_runs(f$spec, f$data, f$goals, f$pool,
    forward = f$forward, forward_control = list(nBoot = 4)), "cannot modify supplied fits")
  expect_identical(calls, 0L)
  expect_identical(.Random.seed, before)
})

test_that("same-data external fits rescore exactly with missing endpoints and row permutations", {
  for (plain in c(FALSE, TRUE)) {
    f <- .plan_guard_fixture(plain)
    set.seed(944)
    before <- .Random.seed
    a <- .plan_guard_call(f, budget = 2, incumbent = FALSE)
    b <- .plan_guard_call(f, budget = 2, incumbent = FALSE)
    reordered <- f
    reordered$data <- f$data[nrow(f$data):1, , drop = FALSE]
    reordered$data$g <- factor(reordered$data$g, levels = rev(levels(f$data$g)))
    c <- .plan_guard_call(reordered, budget = 2, incumbent = FALSE)
    expect_identical(.Random.seed, before)
    expect_identical(a$worklist, b$worklist)
    expect_identical(a$worklist, c$worklist)
    expect_identical(a$policy$goals, f$goals)
    expect_identical(nrow(f$fits$y$training_X), 17L)
    expect_identical(nrow(f$fits$z$training_X), 17L)
    paths <- c(tempfile(fileext = ".csv"), tempfile(fileext = ".csv"))
    expect_message(out1 <- svem_export_candidates_csv(a, file = paths[1]), "wrote")
    expect_message(out2 <- svem_export_candidates_csv(b$worklist, file = paths[2]), "wrote")
    expect_identical(readBin(paths[1], "raw", file.info(paths[1])$size),
                     readBin(paths[2], "raw", file.info(paths[2])$size))
    expect_identical(out1, out2)
    expect_false(withVisible(svem_export_candidates_csv(a, write_file = FALSE))$visible)
    expect_identical(.Random.seed, before)
  }
})

test_that("plans report absent companions and nominations dropped by budget", {
  f <- .plan_guard_fixture()
  a <- .plan_guard_call(f, budget = 1)
  states <- setNames(a$nominations$status, a$nominations$engine)
  expect_identical(unname(states[c("gaupro", "extra_trees")]), rep("not_requested", 2))
  expect_identical(unname(states[c("svemnet", "forward")]), rep("dropped_by_budget", 2))
  printed <- paste(capture.output(print(a)), collapse = "\n")
  expect_match(printed, "gaupro.*not_requested")
  expect_match(printed, "extra_trees.*not_requested")
  expect_match(printed, "dropped_by_budget")
  expect_match(printed, "Scientific goals and weights")
  expect_match(printed, "High-spread exploration fraction")
  expect_match(printed, "best_ties")
  expect_match(printed, "new_ties")
  malformed <- a; malformed$worklist <- NULL
  expect_error(svem_export_candidates_csv(malformed, write_file = FALSE), "no worklist")
  expect_error(svem_export_candidates_csv(a, a, write_file = FALSE), "All inputs")
})


test_that("training-row sorting treats order argument names as predictors", {
  d <- data.frame(method=seq(-1,1,length.out=16),
    decreasing=sin(seq_len(16)), na.last=cos(seq_len(16)))
  d$y <- 2+d$method+.3*d$decreasing-.2*d$na.last
  fit <- withr::with_seed(982, SVEMnet(y~method+decreasing+na.last,
    data=d,nBoot=3,glmnet_alpha=1,relaxed=FALSE))
  expect_true(all(c("method","decreasing","na.last") %in% colnames(fit$training_X)))
  set.seed(983)
  before <- .Random.seed
  expect_no_error(SVEMnet:::.svem_plan_training_rows(fit,d))
  expect_no_error(SVEMnet:::.svem_plan_training_rows(fit,d[16:1,]))
  changed <- d; changed$method[4] <- changed$method[4]+.01
  expect_error(SVEMnet:::.svem_plan_training_rows(fit,changed),"outcomes or design settings differ")
  expect_identical(.Random.seed,before)
})

test_that("documented original lipid benchmark identifier resolves exactly", {
  example_data <- new.env(parent=emptyenv())
  utils::data("lipid_screen",package="SVEMnet",envir=example_data)
  expect_identical(SVEMnet:::.svem_resolve_benchmark(
    list(RunID="02SEP22-2"),example_data$lipid_screen),2L)
})
