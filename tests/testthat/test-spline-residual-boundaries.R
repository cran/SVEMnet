residual_block_fixture <- function() {
  dat <- expand.grid(x=seq(-1,1,length.out=8), temperature=c(20,25))
  dat$y <- 1+dat$x+.1*(dat$temperature==25)+.02*cos(seq_len(nrow(dat)))
  list(spec=bigexp_terms(y~x,dat,blocking="temperature",report=FALSE),data=dat,
       goals=list(y=list(goal="max",weight=1)),
       candidates=data.frame(x=seq(-1,1,length.out=21),temperature=20),
       budget=3L,fit_control=list(companion_continuous_blocks="temperature"))
}

test_that("a continuous companion block reports spline unavailability before fitting", {
  args <- residual_block_fixture()
  events <- character()
  local_mocked_bindings(svem_fit_responses=function(...) {
    events <<- c(events,"primary")
    stop("stopped at primary fitter")
  }, .package="SVEMnet")
  expect_error(withCallingHandlers(do.call(svem_plan_runs,c(args,list(spline=TRUE))),
    warning=function(w) {
      if (grepl("continuous companion block",conditionMessage(w),fixed=TRUE)) {
        events <<- c(events,"unavailable")
        invokeRestart("muffleWarning")
      }
    }), "stopped at primary fitter")
  expect_identical(events,c("unavailable","primary"))
})

test_that("a supplied block-ID spline conflicts with a declared continuous block", {
  skip_if_not_installed("mgcv")
  args <- residual_block_fixture()
  set.seed(517)
  fit <- svem_spline(spline_terms(y~x,args$candidates,blocking="temperature"),args$data,nBoot=2)
  calls <- 0L
  no_fit <- function(...) { calls <<- calls+1L; stop("unexpected fit") }
  local_mocked_bindings(svem_fit_responses=no_fit,svem_forward=no_fit,.package="SVEMnet")
  before <- .Random.seed
  expect_error(do.call(svem_plan_runs,c(args,list(spline=list(y=fit)))),"continuous companion block")
  expect_identical(calls,0L)
  expect_identical(.Random.seed,before)
})

test_that("different mixture declarations reject a supplied spline before fitting", {
  skip_if_not_installed("mgcv")
  g <- expand.grid(a=seq(.1,.3,length.out=5),b=seq(.2,.4,length.out=5))
  g$c <- 1-g$a-g$b
  d <- transform(g[seq(1,25,by=2),],y=a+b^2)
  declared <- list(list(vars=names(g),lower=c(0,0,0),upper=c(.4,1,1),total=1))
  other <- list(list(vars=names(g),total=1))
  spec <- bigexp_terms(y~a+b+c,d,mixture_groups=declared,report=FALSE)
  set.seed(518)
  fit <- svem_spline(spline_terms(y~a+b+c,g,mixture_groups=other),d,nBoot=2)
  calls <- 0L
  no_fit <- function(...) { calls <<- calls+1L; stop("unexpected fit") }
  local_mocked_bindings(svem_fit_responses=no_fit,svem_forward=no_fit,
                        .svem_fit_spline_responses=no_fit,.package="SVEMnet")
  before <- .Random.seed
  expect_error(svem_plan_runs(spec,d,list(y=list(goal="max",weight=1)),g,
               budget=3L,spline=list(y=fit)), "same mixture constraints")
  expect_identical(calls,0L)
  expect_identical(.Random.seed,before)
})

test_that("a failed forward prediction retains its status and diagnostic", {
  tab <- data.frame(x=c(1,2))
  model <- structure(list(),class="svem_test_no_prediction_method")
  expect_warning(out <- .svem_forward_scores(tab,tab,list(y=model),list(y=list()),
             function(y,p)y,c(y=1),"mean",1e-6),"cannot score the joint goal")
  expect_identical(out$status$status,"prediction_failed")
  expect_true(nzchar(out$status$message))
  expect_true(all(is.na(out$table$y_forward_pred)))
  expect_true(all(is.na(out$table$forward_score)))
})
