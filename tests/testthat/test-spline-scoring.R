spline_score_fixture <- function() {
  pool <- expand.grid(a=0:6/6, b=0:6/6)
  d <- pool[seq(1,49,by=2),]
  d$y <- 2+d$a-d$b+.3*sin(3*d$a)
  d$z <- 1+d$b^2+.1*cos(d$a*4)
  spec <- bigexp_terms(y~a+b, d, factorial_order=1, polynomial_order=1, report=FALSE)
  set.seed(81)
  prim <- lapply(c("y","z"), function(r) SVEMnet(spec,d,response=r,nBoot=4,glmnet_alpha=1,relaxed=FALSE))
  names(prim) <- c("y","z")
  goals <- list(y=list(goal="max",weight=.6,lower_acceptable=0,upper_acceptable=4),
                z=list(goal="min",weight=.4,lower_acceptable=0,upper_acceptable=4))
  smooth <- .svem_fit_spline_responses(spec,d,c("y","z"),pool,
             control=list(nBoot=4,rank=10,null_weight=.02,lambda=c(Inf,1,.01)))
  list(pool=pool, data=d, primary=prim, spline=smooth, goals=goals)
}

test_that("separate spline nominees use the same endpoint goals and weights", {
  skip_if_not_installed("mgcv")
  s <- spline_score_fixture()
  scored <- svem_score_random(s$primary,s$goals,data=s$data,candidates=s$pool,
     spline=s$spline,budget=4,n_optima=1,n_explore=0,combine="geom",allocation="portfolio",verbose=FALSE)
  expect_identical(scored$engine_status$engine,c("svemnet","svemnet","spline","spline"))
  expect_true(all(scored$engine_status$status=="ok"))
  expect_identical(tail(scored$allocation_details$priority,1),"spline")
  expect_true("spline" %in% scored$nominations$engine)
  expect_equal(scored$allocation_counts$requested_optima,1)
  expect_true(is.data.frame(scored$merged_nominations))
  tab <- scored$score_table
  for(r in names(s$primary)) {
    expect_equal(tab[[paste0(r,"_spline_pred")]],predict(s$spline[[r]],tab),tolerance=1e-10)
  }
  p <- scored$desirability_parameters
  clip <- function(x) pmax(0,pmin(1,x))
  dy <- clip((tab$y_spline_pred-p$y$L)/(p$y$U-p$y$L))^p$y$s
  dz <- clip((p$z$U-tab$z_spline_pred)/(p$z$U-p$z$L))^p$z$s
  expected <- exp(.6*log(1e-6+(1-1e-6)*dy)+.4*log(1e-6+(1-1e-6)*dz))
  expect_equal(tab$spline_score,expected,tolerance=1e-10)
  expect_identical(s$spline$y$spec, s$spline$z$spec)
  expect_equal(s$spline$y$spec$rank,10)
  expect_equal(s$spline$z$lambda_grid,c(Inf,1,.01))
})

test_that("supplied models must match the real endpoint and physical schema", {
  skip_if_not_installed("mgcv")
  s <- spline_score_fixture()
  run <- function(fit) svem_score_random(s$primary,s$goals,data=s$data,candidates=s$pool,
                 spline=fit,budget=3,allocation="portfolio",verbose=FALSE)
  expect_error(run(TRUE),"named list")
  expect_error(run(s$spline["y"]),"exactly")
  expect_error(run(setNames(rev(s$spline),names(s$spline))),"response labels")
  conditional <- s$spline
  conditional$y <- svem_spline(spline_terms(y~a,s$pool),s$data,nBoot=3)
  expect_error(run(conditional),"same raw science")
  aliased <- s
  names(aliased$primary) <- names(aliased$spline) <- names(aliased$goals) <- c("first","second")
  expect_warning(result <- svem_score_random(aliased$primary,aliased$goals,candidates=s$pool,
    spline=aliased$spline,budget=2,priority="svemnet",allocation="portfolio",verbose=FALSE), "differ from model response names")
  expect_identical(result$allocation_details$priority,"svemnet")
  expect_true(all(c("first_spline_pred","second_spline_pred") %in% names(result$score_table)))
  expect_equal(result$nominations$status[result$nominations$engine=="spline"],"not_prioritized")
})

test_that("mixture scoring respects physical limits beyond observed component extrema", {
  skip_if_not_installed("mgcv")
  x <- expand.grid(a=seq(.1,.3,length.out=4),b=seq(.2,.4,length.out=4))
  pool <- transform(x,c=1-a-b)
  small <- pool[pool$a>.11 & pool$a<.29,]
  small$y <- 2*small$a+sin(small$b)
  groups <- list(list(vars=c("a","b","c"),total=1))
  s <- bigexp_terms(y~a+b+c,small,mixture_groups=groups,report=FALSE)
  fit <- .svem_fit_spline_responses(s,small,"y",pool,control=list(nBoot=2))
  expect_true(all(is.finite(predict(fit$y,pool))))
  expect_error(.svem_fit_spline_responses(s,small,"y",pool,control=list(unused=3)),"named spline")
})

test_that("quantitative discrete settings keep their eligible physical support", {
  skip_if_not_installed("mgcv")
  pool <- expand.grid(x=c(1,2,4),u=0:5/5)
  data <- transform(pool,y=x+u^2)
  def <- bigexp_terms(y~x+u,data,discrete_numeric=list(x=c(1,2,4)),report=FALSE)
  sm <- .svem_fit_spline_responses(def,data,"y",pool,list(nBoot=2))
  expect_true(all(is.finite(predict(sm$y,pool))))
  expect_error(.svem_fit_spline_responses(def,data,"y",transform(pool,x=3)), "discrete")
})


test_that("disabling the spline excludes it from count-based default priority", {
  standard <- c("incumbent", "svemnet", "gaupro", "extra_trees", "forward")
  disabled <- .svem_resolve_plan_priority(standard, TRUE, FALSE, 5L, 0L)
  expect_identical(disabled, standard)
  expect_error(.svem_validate_allocation_counts(6L, 5L, 0L, disabled), "prioritized")
  expect_identical(.svem_resolve_plan_priority(standard, TRUE, NULL, 5L, 0L),
                   c(standard, "spline"))
  expect_identical(.svem_resolve_plan_priority(standard, FALSE, TRUE, 4L, 0L), standard)
})


test_that("disabled and unavailable spline states are checked before planner fitting", {
  grid <- data.frame(x=0:5)
  data <- transform(grid, y=c(1,4,2,3,5,6))
  spec <- bigexp_terms(y~x,data,report=FALSE)
  goals <- list(y=list(goal="max",weight=1))
  previous <- list(y=structure(list(formula=y~1,message="old fitting failure",
      status="fit_failed"),class="svem_spline_unavailable"))
  # An internally attempted failure remains valid scorer input for that
  # batch. It has no training rows proving reuse against a different batch.
  expect_identical(.svem_spline_model_list(previous,"y"),previous)
  expect_error(svem_plan_runs(spec,data,goals,grid,spline=previous,
      fit_control=list(nBoot=-1)), "previously unavailable spline.*Refit")
  expect_error(svem_plan_runs(spec,data,goals,grid,budget=2,
      spline=FALSE,priority=c("svemnet","spline"),n_optima=2,incumbent=FALSE,
      fit_control=list(nBoot=-1)), "spline = FALSE conflicts")
})

