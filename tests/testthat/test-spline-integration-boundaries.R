boundary_problem <- function() {
  data <- data.frame(x=seq(-1,1,length.out=16))
  data$y <- 2 + data$x + .12*cos(3*data$x)
  list(spec=bigexp_terms(y~x,data,report=FALSE),data=data,
       goals=list(y=list(goal="max",weight=1)),
       candidates=data.frame(x=seq(-1,1,length.out=41)))
}

test_that("invalid or contradictory spline requests fail before any model fit", {
  args <- boundary_problem()
  calls <- 0L
  forbidden <- function(...) { calls <<- calls+1L; stop("unexpected fit") }
  local_mocked_bindings(svem_fit_responses=forbidden, svem_forward=forbidden,
                        .svem_fit_spline_responses=forbidden, .package="SVEMnet")
  set.seed(4326); initial_rng <- .Random.seed
  for (ctrl in list(list(nBoot=0),list(rank=2),list(lambda=c(1,.1)),
                    list(store_member_weights=NA),list(null_weight=0))) {
    expect_error(do.call(svem_plan_runs,c(args,list(spline=TRUE,spline_control=ctrl))),
                 names(ctrl))
  }
  expect_error(do.call(svem_plan_runs,c(args,list(spline=TRUE,priority="svemnet"))),
               "spline = TRUE requires")
  expect_error(do.call(svem_plan_runs,c(args,list(spline=FALSE,priority=c("svemnet","spline")))),
               "conflicts")
  expect_error(do.call(svem_plan_runs,c(args,list(spline=FALSE,spline_control=list(rank=6)))),
               "spline_control")
  expect_error(do.call(svem_plan_runs,c(args,list(n_explore=1))),"explicit total physical budget")
  expect_error(do.call(svem_plan_runs,c(args,list(n_optima=1))),"explicit total physical budget")
  expect_identical(calls,0L)
  expect_identical(.Random.seed,initial_rng)
})

test_that("an unsupported optional surface is reported before primary fitting", {
  data <- data.frame(u=1:10,v=10:1,y=1:10)
  spec <- bigexp_terms(y~u+v,data,report=FALSE)
  events <- character()
  primary <- function(...) {
    events <<- c(events,"primary")
    stop("test stopped at first primary fit")
  }
  local_mocked_bindings(svem_fit_responses=primary,.package="SVEMnet")
  expect_error(withCallingHandlers(
      svem_plan_runs(spec,data,list(y=list(goal="max")),data[c("u","v")],spline=TRUE),
      warning=function(w){
        if(grepl("Spline candidate unavailable",conditionMessage(w),fixed=TRUE)){
          events <<- c(events,"unavailable")
          invokeRestart("muffleWarning")
        }
      }),"stopped at first primary fit")
  expect_identical(events,c("unavailable","primary"))
})

test_that("counts-only reuse and lower-level scoring consume no new model fits or RNG", {
  args <- boundary_problem()
  set.seed(123)
  fits <- svem_fit_responses(args$spec,args$data,"y",nBoot=3L,companions=FALSE)
  forward <- list(y=svem_forward(args$spec,args$data,nBoot=3L))
  forbidden <- function(...) stop("unexpected refit")
  local_mocked_bindings(svem_fit_responses=forbidden,svem_forward=forbidden,
                        .svem_fit_spline_responses=forbidden,.package="SVEMnet")
  set.seed(541); before <- .Random.seed
  plan <- do.call(svem_plan_runs,c(args,list(budget=7L,fits=fits,forward=forward,
                n_explore=1L,score_control=list(verbose=FALSE))))
  expect_identical(.Random.seed,before)
  expect_identical(plan$fits,fits)
  expect_identical(plan$forward,forward)
  expect_null(plan$spline)
  expect_false("spline" %in% plan$policy$priority)
  expect_gt(plan$allocation_counts$unfilled_model_capacity,0)
  expect_equal(plan$allocation_counts$model_shortfall,0)
  expect_equal(plan$allocation_counts$selected_explore,1)
  expect_equal(plan$shortfall,0)
  scorer <- list(objects=fits,goals=args$goals,candidates=args$candidates,data=args$data,
                 budget=7L,forward=forward,allocation="portfolio",verbose=FALSE)
  scored <- do.call(svem_score_random,c(scorer,list(n_explore=1L)))
  expect_false("spline" %in% scored$allocation_details$priority)
  expect_error(do.call(svem_score_random,c(scorer,list(n_optima=5L))),"prioritized")
  expect_identical(.Random.seed,before)
})

test_that("the scorer preserves legal spline rows on a wider primary pool", {
  skip_if_not_installed("mgcv")
  f <- boundary_problem()
  set.seed(39)
  primary <- svem_fit_responses(f$spec,f$data,"y",nBoot=3,companions=FALSE)
  narrower <- data.frame(x=seq(-.5,.5,length.out=19))
  fit <- svem_spline(spline_terms(y~x,narrower),
    f$data[abs(f$data$x)<=.5,,drop=FALSE],nBoot=3)
  out <- svem_score_random(primary,f$goals,candidates=f$candidates,data=f$data,
      budget=1L,allocation="portfolio",priority="spline",spline=list(y=fit),verbose=FALSE)
  tab <- out$score_table; legal <- abs(tab$x)<=.5+1e-8
  expect_true(any(legal) && any(!legal))
  expect_equal(tab$y_spline_pred[legal],predict(fit,tab[legal,,drop=FALSE]))
  expect_true(all(is.na(tab$y_spline_pred[!legal])))
  expect_true(all(is.finite(tab$y_pred)))
  expect_identical(out$engine_status$status[out$engine_status$engine=="spline"],"partial_coverage")
  expect_true(out$worklist$x %in% tab$x[legal])
})

test_that("the public bridge uses reference-grid geometry and separate physical limits", {
  skip_if_not_installed("mgcv")
  wide <- expand.grid(a=seq(0,1,length.out=6),b=seq(-2,2,length.out=5))
  wide$y <- 2*wide$a-wide$b+.2*wide$a*wide$b
  grid <- expand.grid(a=seq(.1,.9,length.out=7),b=seq(-.3,.3,length.out=7))
  spec <- bigexp_terms(y~a+b,wide,report=FALSE)
  set.seed(59)
  fit <- svem_fit_responses(spec,wide,"y",method="spline",domain_grid=grid,nBoot=3)[[1]]
  direct <- spline_terms(y~a+b,grid)
  expect_identical(fit$spec$lower,direct$lower)
  expect_identical(fit$spec$widths,direct$widths)
  expect_identical(fit$spec$reference_X,direct$reference_X)
  expect_identical(fit$spec$penalty,direct$penalty)
  expect_equal(fit$spec$legal_lower,c(a=0,b=-2))
  expect_true(all(is.finite(predict(fit,wide))))
})

test_that("count-based planning and portfolio scoring match an explicit cap", {
  a <- boundary_problem()
  set.seed(1049)
  fits <- svem_fit_responses(a$spec, a$data, "y", nBoot = 3, companions = FALSE)
  forward <- list(y = svem_forward(a$spec, a$data, nBoot = 3))
  args <- c(a, list(fits = fits, forward = forward, n_optima = 1L, n_explore = 2L,
                    score_control = list(verbose = FALSE)))
  initial_rng <- .Random.seed
  inferred <- do.call(svem_plan_runs, args)
  capped <- do.call(svem_plan_runs, c(args, list(budget = 4L)))
  expect_identical(inferred, capped)
  expect_identical(inferred$policy$budget, 4L)
  expect_output(print(inferred), "Reference repeats: 1 .*New model candidates: 1 .*Exploration runs: 2")
  without <- do.call(svem_plan_runs, c(args, list(incumbent = FALSE)))
  expect_identical(without$policy$budget, 3L)
  expect_identical(without$worklist$x, inferred$worklist$x[-1L])
  score_args <- list(objects = fits, data = a$data, goals = a$goals, candidates = a$candidates,
                     forward = forward, allocation = "portfolio", incumbent = TRUE, n_optima = 1L,
                     n_explore = 2L, verbose = FALSE)
  expect_identical(do.call(svem_score_random, score_args),
                   do.call(svem_score_random, c(score_args, list(budget = 4L))))
  historical <- which.max(a$data$y)
  for (benchmark in c(historical, 1L)) {
    refs <- length(unique(c(historical, benchmark)))
    shared <- do.call(svem_plan_runs, c(args, list(incumbent = historical, benchmark = benchmark)))
    expect_identical(shared$policy$budget, 3L + refs)
    expect_equal(shared$allocation_counts$distinct_repeats, refs)
  }
  args$n_optima <- args$n_explore <- 0L
  repeats <- do.call(svem_plan_runs, args)
  expect_identical(repeats$worklist$role, "incumbent_repeat")
  empty <- do.call(svem_plan_runs, c(args, list(incumbent = FALSE)))
  expect_identical(empty$policy$budget, 0L)
  expect_equal(nrow(empty$worklist), 0L)
  expect_identical(.Random.seed, initial_rng)
})

test_that("an automatically attempted unsupported spline gives one actionable message", {
  d <- data.frame(ingredient=factor(rep(c("A","B"),each=5)),y=1:10)
  spec <- bigexp_terms(y~ingredient,d,report=FALSE)
  local_mocked_bindings(svem_fit_responses=function(...) stop("stopped before fitting"),.package="SVEMnet")
  expect_message(expect_error(svem_plan_runs(spec,d,list(y=list(goal="max")),d["ingredient"],
    n_optima=1,n_explore=0), "stopped before fitting"), "spline = FALSE")
  expect_error(svem_plan_runs(spec,d,list(y=list(goal="max")),d["ingredient"],n_optima=10),
    "maximum 5.*budget = 10")
  missing <- structure(list(message="unsupported",formula=y~1),class="svem_spline_unavailable")
  expect_silent(scores <- .svem_spline_scores(data.frame(x=1:2),data.frame(x=1:2),
    list(y=missing),list(y=NULL),function(x,parameters)x,c(y=1),"mean",1e-6))
  expect_identical(scores$status$status,"fit_failed")
})
