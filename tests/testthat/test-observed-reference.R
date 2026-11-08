test_that("optional observed means can revise a noisy best while literal best persists", {
  d <- data.frame(x = c(0, 0, 1, 1, .5), y = c(10, 2, 7, 7, NA_real_))
  objects <- list(y = list(formula = y ~ x))
  prepared <- list(rows = seq_len(nrow(d)), all_keys = c("a","a","b","b","c"),
    keys = c("a","a","b","b","c"), eligible_history_rows = seq_len(nrow(d)),
    blocks = character(), reference = list(), explicit = FALSE, excluded_rows = integer())
  args <- list(objects = objects, data = d, prepared = prepared,
    parameters = list(y = list()), des = function(y,p) pmin(1,pmax(0,y/10)),
    weights = c(y=1), combine="mean", floor=1e-6,
    goals=list(y=list(goal="max")), pool_keys=c("a","b","c"))
  best <- do.call(SVEMnet:::.svem_observed_choice, args)
  avg <- do.call(SVEMnet:::.svem_observed_choice, c(args,list(summary="mean")))
  expect_identical(best$source_rows, 1L)
  expect_equal(best$observed_responses$y, 10)
  expect_identical(best$summary, "best")
  expect_identical(best$recipe_source_rows, c(1L,2L))
  expect_equal(best$recipe_n_preparations,2)
  expect_equal(best$recipe_n_complete,2)
  expect_equal(best$recipe_mean_responses$y,6)
  expect_equal(best$recipe_mean_score,.6)
  expect_equal(best$n_preparations,1)
  expect_identical(avg$source_rows, c(3L,4L))
  expect_equal(avg$n_preparations, 2)
  expect_equal(avg$observed_responses$y, 7)
  expect_identical(avg$summary, "mean")
  expect_equal(avg$observed_score, .7)
  expect_identical(avg$incomplete_rows, 5L)
  args$prepared$explicit <- TRUE
  args$prepared$rows <- 2L
  forced <- do.call(SVEMnet:::.svem_observed_choice, c(args,list(summary="mean")))
  expect_identical(forced$source_rows, 2L)
  expect_equal(forced$observed_responses$y, 2)
  expect_identical(forced$summary, "explicit_row")
})

test_that("context winners remain separately observed and never pooled", {
  d <- data.frame(x=c(0,1,0,1), day=c("A","A","B","B"), y=c(10,9,1,100))
  objects <- list(y=list(formula=y~x+day))
  prepared <- list(rows=1:2, all_keys=c("a","b","a","b"), keys=c("a","b"),
    eligible_history_rows=1:4, blocks="day", reference=list(day="A"),
    explicit=FALSE, excluded_rows=integer())
  z <- SVEMnet:::.svem_observed_choice(objects,d,prepared,list(y=list()),
    function(y,p) pmin(1,pmax(0,y/100)), c(y=1),"mean",1e-6,
    list(y=list(goal="max")), c("a","b"),"mean")
  expect_identical(z$source_rows, 1L)
  expect_identical(z$source_context, "day=A")
  expect_true(z$higher_raw_in_other_context)
  expect_equal(z$by_context$observed_score, c(.1,1))
  expect_equal(z$by_context$source_row, c(1,4))
  expect_equal(z$by_context$n_preparations,c(1,1))
})

test_that("portfolio export preserves observed and prediction contexts", {
  d <- data.frame(x=seq(-1,1,length.out=12), day=factor(rep(c("A","B"),each=6)),
                  y=seq_len(12))
  spec <- bigexp_terms(y~x,d,blocking="day",report=FALSE)
  pool <- data.frame(x=seq(-1,1,length.out=21),day=factor("A",levels=c("A","B")))
  set.seed(824)
  plan <- svem_plan_runs(spec,d,list(y=list(goal="max",weight=1)),pool,budget=2,
    incumbent=12L, fit_control=list(nBoot=3,companions=FALSE),
    forward_control=list(nBoot=3),score_control=list(verbose=FALSE))
  work <- plan$worklist
  expect_identical(work$prediction_context, rep("day=A",nrow(work)))
  expect_identical(work$reference_context[work$role=="incumbent_repeat"], "day=B")
  expect_identical(work$incumbent_source_rows[work$role=="incumbent_repeat"], "12")
  expect_equal(work$incumbent_n_preparations[work$role=="incumbent_repeat"],1)
  out <- svem_export_candidates_csv(plan,write_file=FALSE)
  expect_identical(out$prediction_context,work$prediction_context)
  expect_identical(out$reference_context,work$reference_context)
  eligibility <- plan$scores$uncertainty_eligibility
  expect_identical(eligibility$recipe_id,plan$scores$score_table$recipe_id)
  expect_identical(eligibility$informative_width,is.finite(eligibility$standardized_width))
  expect_named(plan$scores$uncertainty_anchors,"y")
})

