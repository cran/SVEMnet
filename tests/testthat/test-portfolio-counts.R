count_worklist <- function(budget=6L, n_optima=NULL, n_explore=NULL, incumbent=1L,
                           benchmark=NULL, order=c("incumbent","svemnet","gaupro","extra_trees","forward","spline"),
                           anchor=c(0,1), raw=FALSE) {
  x <- 0:9
  top <- function(k) -abs(x-k)
  tab <- data.frame(recipe_id=letters[1:10],x=x,score=x/10,y_pred=top(1),
      y_gaupro_pred=top(1),y_extra_trees_pred=top(2),y_forward_pred=top(1),
      y_spline_pred=top(3),y_lwr=0,y_upr=x/10,y_ciw_w=x/10,uncertainty_measure=x/10)
  history <- data.frame(x=if(is.null(benchmark)) 0 else unique(c(0,tab$x[benchmark])))
  args <- list(tab=tab, history=history, keys=.svem_recipe_keys(tab,"x"),
      increments=NULL, schema=list(predictor_vars="x",blocking=NULL),
      engines=c("gaupro","extra_trees","forward","spline"), goals=list(y=list(goal="max")),
      weights=c(y=1), uncertainty_weights=c(y=1), anchors=list(y=anchor), gp_endpoint=NULL,
      budget=budget,allocation="portfolio",priority=order,incumbent=incumbent,benchmark=benchmark)
  if (raw) return(args)
  do.call(.svem_worklist,c(args,list(n_optima=n_optima,n_explore=n_explore)))
}

test_that("model quotas count distinct new recipes and merge endorsements", {
  res <- count_worklist(n_optima=3,n_explore=1)
  expect_equal(res$worklist$recipe_id[1:4],letters[1:4])
  expect_equal(res$worklist$role,c("incumbent_repeat",rep("optimum",3),"high_spread_medoid"))
  expect_equal(res$worklist$nominated_by[2],"svemnet + gaupro + forward")
  expect_equal(res$merged_nominations$source,c("gaupro","forward"))
  expect_equal(res$allocation_counts$selected_optima,3)
  expect_equal(res$allocation_counts$selected_explore,1)
  expect_equal(res$allocation_counts$intentionally_unused,1)
  expect_equal(res$allocation_counts$unused_capacity,1)
  expect_equal(res$shortfall,0)
  expect_equal(attr(res$worklist,"shortfall"),0)
  p <- count_worklist(n_optima=2,n_explore=0)
  expect_false("d" %in% p$worklist$recipe_id)
  expect_equal(p$nominations$status[p$nominations$engine=="spline"],"dropped_by_model_count")
  expect_equal(p$worklist$nominated_by[2],"svemnet + gaupro + forward")
})

test_that("references reserve physical space first and can share one run", {
  order <- c("forward","svemnet","incumbent","gaupro","extra_trees","spline")
  both <- count_worklist(budget=3,n_optima=1,n_explore=1,benchmark=1L,order=order)
  expect_equal(both$worklist$role[1],"incumbent_benchmark_repeat")
  expect_equal(both$allocation_counts$distinct_repeats,1)
  expect_equal(both$worklist$nominated_by[1],"observed + benchmark")
  expect_equal(both$worklist$role,c("incumbent_benchmark_repeat","optimum","high_spread_medoid"))
  expect_error(count_worklist(budget=3,n_optima=1,n_explore=1,benchmark=8L),"distinct required references")
  distinct <- count_worklist(budget=3,n_optima=1,n_explore=0,benchmark=8L,order=order)
  expect_equal(distinct$worklist$recipe_id,c("a","h","b"))
  no_refs <- count_worklist(budget=3,n_optima=1,n_explore=1,incumbent=NULL)
  expect_equal(no_refs$allocation_counts$distinct_repeats,0)
  expect_equal(no_refs$allocation_counts$unused_capacity,1)
})

test_that("automatic and exact exploration report role-specific shortages", {
  expect_warning(filled <- count_worklist(budget=8,n_optima=4), "short by 1 model")
  cts <- filled$allocation_counts
  expect_equal(nrow(filled$worklist),8)
  expect_equal(cts$model_shortfall,1)
  expect_equal(cts$unused_capacity,0)
  expect_equal(cts$selected_explore,4)
  expect_equal(cts$requested_shortfall,1)
  expect_warning(fixed <- count_worklist(budget=8,n_optima=4,n_explore=1), "short by 1 model")
  expect_equal(fixed$allocation_counts$selected_explore,1)
  expect_equal(fixed$allocation_counts$unused_capacity,3)
  expect_equal(fixed$allocation_counts$intentionally_unused,2)
  expect_warning(empty <- count_worklist(budget=5,n_optima=0,n_explore=2,anchor=c(0,0)),
                 "0 model and 2 exploration")
  expect_equal(nrow(empty$worklist),1)
  expect_equal(empty$allocation_counts$exploration_shortfall,2)
  expect_equal(empty$allocation_counts$intentionally_unused,2)
  auto <- count_worklist(budget=4,n_explore=1)
  expect_equal(auto$worklist$role,c("incumbent_repeat","optimum","optimum","high_spread_medoid"))
})

test_that("count controls are explicit and omitted controls retain routing", {
  expect_error(count_worklist(budget=NULL,n_optima=2),"explicit total physical budget")
  expect_error(count_worklist(n_optima=6),"prioritized")
  expect_error(count_worklist(n_optima=1.2),"integer")
  expect_error(count_worklist(n_explore=-1),"integer")
  expect_error(count_worklist(budget=2,n_optima=2,n_explore=1),"exceed")
  expect_error(count_worklist(n_optima=2,order=c("incumbent","forward")),"prioritized")
  a <- count_worklist(raw=TRUE)
  old <- do.call(.svem_worklist,a)
  expect_identical(old,do.call(.svem_worklist,c(a,list(n_optima=NULL,n_explore=NULL))))
  expect_null(old$allocation_counts)
  a$allocation <- "legacy"
  expect_error(do.call(.svem_worklist,c(a,list(n_optima=0))),"portfolio")
})

test_that("two exact counts infer the budget from distinct physical references", {
  for (case in list(list(incumbent = NULL), list(incumbent = 1L),
                    list(incumbent = 1L, benchmark = 1L),
                    list(incumbent = 1L, benchmark = 8L))) {
    n_ref <- length(unique(unlist(case, use.names = FALSE)))
    counts <- c(case, list(n_optima = 2L, n_explore = 1L))
    inferred <- do.call(count_worklist, c(counts, list(budget = NULL)))
    capped <- do.call(count_worklist, c(counts, list(budget = n_ref + 3L)))
    expect_identical(inferred, capped)
  }
  expect_error(.svem_count_targets(NULL, 1L, .Machine$integer.max, integer(), "svemnet"),
               "inferred budget")
})

test_that("oversized counts explain the available roster and a total-run alternative", {
  expect_error(count_worklist(budget=NULL, n_optima=10), "maximum 5.*budget = 10")
  expect_error(count_worklist(n_optima=2,order=c("incumbent","forward")), "maximum 1.*forward")
  x <- count_worklist(budget=10,n_explore=2)
  plan <- structure(list(worklist=x$worklist,policy=list(budget=10),shortfall=x$shortfall,
    allocation_counts=x$allocation_counts,nominations=x$nominations,merged_nominations=x$merged_nominations),class="svem_run_plan")
  expect_identical(x$shortfall,0L)
  expect_output(print(plan), "6 of 10 planned preparations; 4 unused places; 0 unmet")
})
