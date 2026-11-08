test_that("the planner preserves primary and forward model identities", {
  as_fit <- function(subclass = character()) {
    structure(list(formula = y ~ x, family = "gaussian"),
              class = c(subclass, "svem_model", "SVEMnet"))
  }
  for (kind in c("svem_forward", "svem_forward_ic")) {
    expect_error(.svem_plan_model_list(list(y = as_fit(kind)), "y", "fits"),
                 "Use the forward argument for svem_forward")
  }
  ordinary <- list(y = as_fit())
  expect_identical(.svem_plan_model_list(ordinary, "y", "fits"), ordinary)
  extension <- list(y = as_fit("other_primary_subclass"))
  expect_identical(.svem_plan_model_list(extension, "y", "fits"), extension)
  forward <- list(y = as_fit("svem_forward"))
  expect_identical(.svem_plan_model_list(forward, "y", "forward"), forward)
})
