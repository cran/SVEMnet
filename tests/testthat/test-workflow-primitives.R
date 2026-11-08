.workflow_primitive_data <- function() {
  dat <- data.frame(x = seq(-1, 1, length.out = 24), z = rep(c(-1, 1), 12))
  dat$y1 <- 2 + dat$x - .4 * dat$z + .05 * sin(seq_len(nrow(dat)))
  dat$y2 <- 3 - dat$x + .2 * dat$z + .08 * cos(seq_len(nrow(dat)))
  dat
}

test_that("response-list fitting preserves sequential draws and ordinary fit outputs", {
  dat <- .workflow_primitive_data()
  spec <- bigexp_terms(y1 ~ x + z, dat, factorial_order = 1,
                       polynomial_order = 1, report = FALSE)
  set.seed(817)
  manual <- list(
    y2 = SVEMnet(spec, dat, response = "y2", nBoot = 3, objective = "wBIC",
                  glmnet_alpha = 1, relaxed = FALSE, store_member_weights = TRUE),
    y1 = SVEMnet(spec, dat, response = "y1", nBoot = 3, objective = "wBIC",
                  glmnet_alpha = 1, relaxed = FALSE, store_member_weights = TRUE)
  )
  manual_seed <- .Random.seed
  set.seed(817)
  fitted <- svem_fit_responses(spec, dat, c(second = "y2", first = "y1"),
                               nBoot = 3, objective = "wBIC", glmnet_alpha = 1,
                               relaxed = FALSE, store_member_weights = TRUE)
  expect_identical(class(fitted), "list")
  expect_named(fitted, c("y2", "y1"))
  expect_identical(.Random.seed, manual_seed)
  for (response in names(manual)) {
    expect_s3_class(fitted[[response]], "svem_model")
    for (field in c("parms", "coef_matrix", "member_weights", "member_diagnostics",
                    "objective", "y_pred", "training_X", "actual_y")) {
      expect_identical(fitted[[response]][[field]], manual[[response]][[field]])
    }
    expect_identical(predict(fitted[[response]], dat), predict(manual[[response]], dat))
  }
  set.seed(818)
  one <- SVEMnet(spec, dat, nBoot = 2, glmnet_alpha = 1, relaxed = FALSE)
  set.seed(818)
  list_one <- svem_fit_responses(spec, dat, "y1", nBoot = 2,
                                 glmnet_alpha = 1, relaxed = FALSE)
  expect_identical(one$coef_matrix, list_one$y1$coef_matrix)
  expect_named(list_one, "y1")
})

test_that("response-list validation happens before fitting or drawing weights", {
  dat <- .workflow_primitive_data()
  spec <- bigexp_terms(y1 ~ x + z, dat, report = FALSE)
  set.seed(819)
  before <- .Random.seed
  expect_error(svem_fit_responses(spec, dat, c("y1", "y1")), "duplicate responses")
  expect_error(svem_fit_responses(spec, dat, c("y1", "absent")), "Unknown.*absent")
  for (bad in list(character(), NA_character_, "", 1L, list("y1"))) {
    expect_error(svem_fit_responses(spec, dat, bad), "nonempty character")
  }
  expect_error(svem_fit_responses(y1 ~ x, dat, "y1"), "bigexp_spec")
  expect_error(svem_fit_responses(spec, as.matrix(dat), "y1"), "data frame")
  expect_error(svem_fit_responses(spec, dat, "y1", response = "y2"), "separate.*response")
  expect_identical(.Random.seed, before)
})

.workflow_export_fixture <- function() {
  out <- data.frame(recipe_id = c("new-C", "new-A", "new-B"), slot = c(3L, 1L, 2L),
                    role = c("high_spread_medoid", "optimum", "maximum_spread"),
                    nominated_by = c("svemnet", "svemnet + gaupro", "svemnet"),
                    x = c(.3, .1, .2), y_pred = c(4.4, 8.3, 2.7),
                    stringsAsFactors = FALSE)
  attr(out, "budget") <- 5L
  attr(out, "shortfall") <- 2L
  out
}

test_that("budgeted scored results and direct worklists export the same identities", {
  work <- .workflow_export_fixture()
  scored <- list(score_table = data.frame(recipe_id = c("old", work$recipe_id)),
                  original_data_scored = data.frame(x = 0), worklist = work,
                  shortfall = 2L)
  a <- withVisible(svem_export_candidates_csv(scored, write_file = FALSE))
  b <- withVisible(svem_export_candidates_csv(work, write_file = FALSE))
  expect_false(a$visible)
  expect_false(b$visible)
  expect_identical(a$value, b$value)
  expect_named(a$value, names(work))
  expect_identical(a$value$slot, 1:3)
  expect_identical(a$value$recipe_id, c("new-A", "new-B", "new-C"))
  expect_identical(a$value$role, c("optimum", "maximum_spread", "high_spread_medoid"))
  expect_identical(a$value$nominated_by[1], "svemnet + gaupro")
  expect_identical(attr(a$value, "budget"), 5L)
  expect_identical(attr(a$value, "shortfall"), 2L)
  expect_false(any(c("candidate_type", "selection_label") %in% names(a$value)))
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  expect_message(written <- withVisible(svem_export_candidates_csv(scored, file = path)),
                  "wrote candidates")
  expect_false(written$visible)
  readback <- utils::read.csv(path, stringsAsFactors = FALSE)
  expect_named(readback, names(work))
  for (column in names(work)) expect_equal(readback[[column]], a$value[[column]])
  expect_error(svem_export_candidates_csv(work, file = path), "already exists")
  expect_message(svem_export_candidates_csv(work, file = path, overwrite = TRUE), "wrote")
})

test_that("empty shortfalls export headers and arbitrary score pools are rejected", {
  work <- .workflow_export_fixture()[FALSE, ]
  attr(work, "shortfall") <- 5L
  out <- svem_export_candidates_csv(work, write_file = FALSE)
  expect_identical(nrow(out), 0L)
  expect_named(out, names(work))
  expect_identical(attr(out, "shortfall"), 5L)
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  expect_message(svem_export_candidates_csv(work, file = path), "wrote")
  expect_identical(nrow(utils::read.csv(path)), 0L)
  expect_named(utils::read.csv(path), names(work))
  scored <- list(score_table = data.frame(x = 1:2, score = c(.5, .8)),
                  original_data_scored = NULL)
  expect_error(svem_export_candidates_csv(scored, write_file = FALSE),
                "no worklist.*budget.*select rows")
  expect_error(svem_export_candidates_csv(scored$score_table, write_file = FALSE),
                "recipe_id, slot, role, and nominated_by")
  expect_error(svem_export_candidates_csv(work), "file.*CSV path")
  expect_error(svem_export_candidates_csv(work, work, write_file = FALSE), "All inputs")
})

test_that("ambiguous worklist identities fail instead of being relabeled", {
  work <- .workflow_export_fixture()
  bad <- work
  bad$recipe_id[2] <- bad$recipe_id[1]
  expect_error(svem_export_candidates_csv(bad, write_file = FALSE), "recipe_id.*unique")
  for (slots in list(c(1, 1, 3), c(0, 1, 2), c(1, 2.5, 3), c(1, NA, 3))) {
    bad <- work
    bad$slot <- slots
    expect_error(svem_export_candidates_csv(bad, write_file = FALSE), "slot.*whole numbers")
  }
  bad <- work
  bad$role[1] <- NA_character_
  expect_error(svem_export_candidates_csv(bad, write_file = FALSE), "role.*labels")
  bad <- work
  bad$nominated_by[1] <- " "
  expect_error(svem_export_candidates_csv(bad, write_file = FALSE), "nominated_by.*labels")
  bad <- work
  names(bad)[6] <- "slot"
  expect_error(svem_export_candidates_csv(bad, write_file = FALSE), "unique column names")
})

test_that("legacy exports retain rows, columns, labels and invisible return", {
  first <- list(best = data.frame(x = 1, score = .9, Notes = "a"),
                 candidates = data.frame(x = 2, score = .7, Notes = "b"), label = "old label")
  second <- list(best = data.frame(x = 3, score = .8, extra = "c"),
                  candidates = NULL, call = quote(svem_select_from_score_table(label = "literal")))
  expected <- data.frame(candidate_type = c("best", "medoid", "best"),
                         selection_label = c("old label", "old label", "literal"),
                         Notes = c("a", "b", NA), x = c(1, 2, 3),
                         score = c(.9, .7, .8), extra = c(NA, NA, "c"))
  separate <- withVisible(svem_export_candidates_csv(first, second, write_file = FALSE))
  listed <- svem_export_candidates_csv(list(first, second), write_file = FALSE)
  expect_false(separate$visible)
  expect_identical(separate$value, expected)
  expect_identical(listed, expected)
  expect_identical(svem_export_candidates_csv(first, write_file = FALSE),
                    expected[1:2, setdiff(names(expected), "extra")])
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  expect_message(svem_export_candidates_csv(first, second, file = path), "wrote")
  expect_equal(utils::read.csv(path, stringsAsFactors = FALSE), expected)
  expect_error(svem_export_candidates_csv(first, file = path), "already exists")
})

test_that("worklist export respects historical exclusion and merged endorsements", {
  tab <- data.frame(x = 0:3, recipe_id = letters[1:4], y_pred = 4:1,
                     y_gaupro_pred = 4:1, score = c(1, .9, .8, .7),
                     y_lwr = 0, y_upr = c(1, 2, 4, 3),
                     uncertainty_measure = c(.1, .3, .9, .8))
  result <- .svem_worklist(tab, data.frame(x = 0),
                            .svem_recipe_keys(tab, "x", c(x = 1)), c(x = 1),
                            list(predictor_vars = "x"), "gaupro",
                            list(y = list(goal = "max")), c(y = 1), c(y = 1),
                            list(y = c(1, 4)), NULL, 2L)
  scored <- c(list(score_table = tab), result)
  exported <- svem_export_candidates_csv(scored, write_file = FALSE)
  expect_false("a" %in% exported$recipe_id)
  expect_identical(exported$recipe_id, c("b", "c"))
  expect_identical(exported$nominated_by[1], "svemnet + gaupro")
  expect_identical(exported$role, c("optimum", "maximum_spread"))
  expect_equal(nrow(exported), 2L)
})
