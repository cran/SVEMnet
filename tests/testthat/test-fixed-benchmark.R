benchmark_fixture <- function() {
  d <- data.frame(RunID=paste0("seed-",1:12), x=seq(-1,1,length.out=12),
                  y=2+seq(-1,1,length.out=12))
  spec <- bigexp_terms(y~x,d,report=FALSE)
  list(d=d,spec=spec,goals=list(y=list(goal="max",weight=1)),
    pool=data.frame(x=seq(-1,1,length.out=31)))
}

test_that("fixed benchmark is mandatory, exportable, and stable by original run ID", {
  f <- benchmark_fixture()
  set.seed(912)
  plan <- svem_plan_runs(f$spec,f$d,f$goals,f$pool,budget=2,
    benchmark=list(RunID="seed-1"),fit_control=list(nBoot=3,companions=FALSE),
    forward_control=list(nBoot=3),score_control=list(verbose=FALSE))
  expect_identical(plan$worklist$role,c("incumbent_repeat","benchmark_repeat"))
  expect_equal(plan$worklist$x,c(1,-1))
  expect_false(any(plan$worklist$new_recipe))
  expect_equal(plan$shortfall,0)
  expect_identical(plan$benchmark$source_row,1L)
  expect_identical(plan$policy$benchmark,list(RunID="seed-1"))
  nom <- plan$nominations[plan$nominations$engine=="benchmark",]
  expect_identical(nom$status,"selected")
  expect_equal(nom$slot,2)
  expect_identical(nom$fit_status,"observed")
  before <- .Random.seed
  again <- svem_plan_runs(f$spec,f$d[12:1,],f$goals,f$pool,budget=2,
    fits=plan$fits,forward=plan$forward,benchmark=list(RunID="seed-1"),
    score_control=list(verbose=FALSE))
  expect_identical(.Random.seed,before)
  expect_equal(again$worklist$x,plan$worklist$x)
  expect_identical(again$benchmark$source_row,12L)
  exported <- svem_export_candidates_csv(plan,write_file=FALSE)
  expect_equal(exported$benchmark_source_row,c(NA,1))
  expect_identical(exported$role,plan$worklist$role)
})

test_that("coincident incumbent and fixed benchmark use one preparation with both roles", {
  f <- benchmark_fixture()
  set.seed(913)
  plan <- svem_plan_runs(f$spec,f$d,f$goals,f$pool,budget=1,
    benchmark=list(RunID="seed-12"),fit_control=list(nBoot=3,companions=FALSE),
    forward_control=list(nBoot=3),score_control=list(verbose=FALSE))
  expect_equal(nrow(plan$worklist),1)
  expect_equal(plan$shortfall,0)
  expect_identical(plan$worklist$role,"incumbent_benchmark_repeat")
  expect_match(plan$worklist$nominated_by,"observed.*benchmark")
  expect_equal(plan$worklist$incumbent_n_preparations,1)
  expect_equal(plan$worklist$benchmark_source_row,12)
  expect_equal(plan$nominations$slot[plan$nominations$engine %in% c("incumbent","benchmark")],c(1,1))
  expect_error(svem_plan_runs(f$spec,f$d,f$goals,f$pool,budget=1,
    fits=plan$fits,forward=plan$forward,incumbent=12L,
    benchmark=list(RunID="seed-1"),score_control=list(verbose=FALSE)),
    "budget|Budget")
  benchmark_only <- svem_plan_runs(f$spec,f$d,f$goals,f$pool,budget=1,
    fits=plan$fits,forward=plan$forward,incumbent=FALSE,
    benchmark=list(RunID="seed-1"),priority="svemnet",score_control=list(verbose=FALSE))
  expect_equal(benchmark_only$worklist$x,-1)
  expect_identical(benchmark_only$worklist$role,"benchmark_repeat")
})

test_that("benchmark need not have a measured outcome and selector is unambiguous", {
  f <- benchmark_fixture()
  f$d$y[1] <- NA_real_
  set.seed(914)
  plan <- svem_plan_runs(f$spec,f$d,f$goals,f$pool,budget=2,
    benchmark=list(RunID="seed-1"),fit_control=list(nBoot=3,companions=FALSE),
    forward_control=list(nBoot=3),score_control=list(verbose=FALSE))
  expect_equal(plan$worklist$x,c(1,-1))
  expect_true(is.na(plan$benchmark$source_data$y))
  expect_error(SVEMnet:::.svem_resolve_benchmark(list(RunID="missing"),f$d),"exactly one|unique|match")
  dup <- rbind(f$d,f$d[1,])
  expect_error(SVEMnet:::.svem_resolve_benchmark(list(RunID="seed-1"),dup),"exactly one|unique|match")
  expect_error(SVEMnet:::.svem_resolve_benchmark(list(no_such_column="x"),f$d),"column|Column")
})


test_that("one mandatory benchmark fits a one-run budget with the incumbent disabled", {
  f <- benchmark_fixture()
  roles <- eval(formals(svem_plan_runs)$priority)
  args <- list(spec=f$spec, data=f$d, candidates=f$pool, incumbent=FALSE,
               reference=NULL, score_control=list(), benchmark=list(RunID="seed-1"),
               budget=1L, priority=roles)
  expect_invisible(do.call(.svem_plan_preflight, args))
  expect_invisible(do.call(.svem_plan_preflight, c(args, list(n_optima=0L, n_explore=0L))))
  expect_error(do.call(.svem_plan_preflight, c(args, list(n_optima=1L))),
               "distinct required references")
})
