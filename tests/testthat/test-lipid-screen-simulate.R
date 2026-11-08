test_that("the lipid example reproduces reference means at new feasible recipes", {
  # Boundary and interior reference evaluations spanning all ingredient labels.
  x <- data.frame(
    PEG = c(.01, .035, .05), Helper = c(.10, .275, .60),
    Ionizable = c(.29, .455, .10), Cholesterol = c(.60, .235, .25),
    Ionizable_Lipid_Type = factor(c("H101", "H102", "H103"),
                                  levels = c("H103", "H101", "H102")),
    N_P_ratio = c(6, 14, 8), flow_rate = c(1, 3, 2.5)
  )
  mu <- lipid_screen_simulate(x)
  expect_equal(mu$Potency_true, c(74.50470489392985, 79.8221770750573, 82.28786349841252))
  expect_equal(mu$Size_true, c(225.06899608733588, 131.902355751172, 48.896546107850455))
  expect_identical(mu$PDI_true, rep(.20, 3))
  x$Operator <- c("B", "A", "C")
  x$Potency <- 123
  expect_identical(lipid_screen_simulate(x), mu)
  expect_equal(nrow(lipid_screen_simulate(x[FALSE, ])), 0L)
  expect_identical(.lipid_screen_round(c(-2.5, -1.5, -.5, .5, 1.5, 2.5)),
                   c(-3, -2, -1, 1, 2, 3))
  bad <- x
  bad$PEG[1L] <- .02
  expect_error(lipid_screen_simulate(bad), "summing to one")
  bad <- x
  bad$Ionizable_Lipid_Type <- rep("unknown", nrow(bad))
  expect_error(lipid_screen_simulate(bad), "lipid label")
})

test_that("new preparations are keyed by ID and do not use the fitting RNG", {
  skip_if_not_installed("digest")
  data("lipid_screen", package = "SVEMnet", envir = environment())
  x <- lipid_screen[1:3, ]
  ids <- c("prep-A", "prep-B", "prep-C")
  set.seed(932)
  original_state <- .Random.seed
  draws <- lipid_screen_simulate(x, ids, 812L)
  expect_identical(.Random.seed, original_state)
  subset_ids <- c(3, 1)
  pr <- lipid_screen_simulate(x[subset_ids, ], ids[subset_ids], 812L)
  expect_equal(unname(as.matrix(pr)), unname(as.matrix(draws[subset_ids, ])))
  expect_false(identical(draws$PDI_draw, lipid_screen_simulate(x, ids, 813L)$PDI_draw))
  repeats <- lipid_screen_simulate(x[c(1, 1), ], c("prep-A", "new-prep-A"), 812L)
  expect_identical(repeats$Potency_true[1L], repeats$Potency_true[2L])
  expect_false(identical(repeats$PDI_draw[1L], repeats$PDI_draw[2L]))
  expect_true(all(draws$PDI >= 0 & draws$PDI <= 1))
  expect_equal(draws$PDI, round(draws$PDI_draw, 2))
  expect_identical(draws$Size, round(draws$Size))
  expect_identical(draws$Potency, round(draws$Potency))
  expect_error(lipid_screen_simulate(x, c("A", "A", "B")), "unique")
  expect_error(lipid_screen_simulate(x, ids, -1), "seed")
})
