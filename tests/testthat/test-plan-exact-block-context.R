numeric_day_setup <- function(day_ids) {
  x <- rep(seq(0, 1, length.out=12), length(day_ids))
  data <- data.frame(x=x, day=rep(day_ids,each=12), y=1+x)
  data$y[13:24] <- 100 + 2*(1-data$x[13:24])
  spec <- bigexp_terms(y~x,data,blocking="day",factorial_order=1,
                       polynomial_order=1,report=FALSE)
  list(spec=spec,data=data,goals=list(y=list(goal="max",weight=1)),
       candidates=data.frame(x=seq(0,1,length.out=19),day=day_ids[1]),
       budget=1L)
}

test_that("one numeric block ID can have integer and double storage", {
  args <- numeric_day_setup(c(99999L, 100000L))
  args$candidates$day <- 1e5
  set.seed(9217)
  plan <- do.call(svem_plan_runs, c(args,list(reference=list(day=100000L),
    fit_control=list(nBoot=3,companions=FALSE),forward_control=list(nBoot=3),
    score_control=list(verbose=FALSE))))
  expect_identical(plan$incumbent$source_row,13L)
  expect_equal(plan$incumbent$source_data$day,100000L)
  expect_equal(plan$worklist$x,0)
  expect_equal(plan$worklist$prediction_context,plan$worklist$reference_context)
  # Preserve the accepted text representation of a numeric reference as well.
  again <- do.call(svem_plan_runs, c(args,list(fits=plan$fits,forward=plan$forward,
    reference=list(day="1e5"),score_control=list(verbose=FALSE))))
  expect_identical(again$incumbent$source_row,13L)
})

test_that("nearby distinct days stay separate in selection and provenance", {
  close_day <- 1 + 2*.Machine$double.eps
  args <- numeric_day_setup(c(1, close_day, 3))
  set.seed(1133)
  first <- do.call(svem_plan_runs,c(args,list(
    fit_control=list(nBoot=3,companions=FALSE),forward_control=list(nBoot=3),
    score_control=list(verbose=FALSE))))
  expect_identical(first$incumbent$source_row,12L)
  expect_equal(first$worklist$x,1)
  expect_equal(first$incumbent$recipe_n_preparations,1)
  expect_equal(nrow(first$incumbent$by_context),3L)
  expect_equal(length(unique(first$incumbent$by_context$context_key)),3L)
  expect_equal(length(unique(first$incumbent$by_context$context)),3L)
  expect_true(first$incumbent$higher_raw_in_other_context)
  reuse <- list(fits=first$fits,forward=first$forward,score_control=list(verbose=FALSE))
  other <- do.call(svem_plan_runs,c(args,reuse,list(reference=list(day=close_day))))
  expect_identical(other$incumbent$source_row,13L)
  expect_equal(other$worklist$x,0)
  expect_false(identical(other$worklist$prediction_context, other$worklist$reference_context))
  args$candidates$day <- 1 + 4*.Machine$double.eps
  expect_error(do.call(svem_plan_runs,c(args,reuse)),"no measured preparations")
  expect_identical(.svem_matches_block(factor(c("01","1",NA)),"01"),c(TRUE,FALSE,FALSE))
})

test_that("cross-validated fits cannot be reported as primary SVEM nominees", {
  args <- numeric_day_setup(c(0L,1L))
  cv <- structure(list(formula=y~x,family="gaussian"),class=c("svem_cv","svem_model"))
  expect_error(do.call(svem_plan_runs,c(args,list(fits=list(y=cv)))),"cross-validated fit")
  ordinary <- list(y=structure(list(formula=y~x,family="gaussian"),
                               class=c("custom_primary","svem_model")))
  expect_identical(.svem_plan_model_list(ordinary,"y","fits"),ordinary)
})
