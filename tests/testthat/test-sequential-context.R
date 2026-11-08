test_that("a measured new day is rebuilt and retained across two planned batches", {
  d <- expand.grid(x=seq(-1,1,length.out=12),day=c("A","B"),
                   KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
  d$day <- factor(d$day,levels=c("A","B"))
  d$RunID <- paste0("original-",seq_len(nrow(d)))
  d$y <- 3+d$x+.2*as.integer(d$day)+.05*sin(seq_len(nrow(d)))
  spec <- bigexp_terms(y~x,d,blocking="day",report=FALSE)
  pool <- data.frame(x=seq(-1,1,length.out=41),
                     day=factor("B",levels=levels(d$day)))
  goals <- list(y=list(goal="max",weight=1))
  fit_args <- list(nBoot=3,companions=FALSE,glmnet_alpha=1,relaxed=FALSE)
  set.seed(982)
  first <- svem_plan_runs(spec,d,goals,pool,budget=3,
    benchmark=list(RunID="original-1"),fit_control=fit_args,
    forward_control=list(nBoot=3),score_control=list(verbose=FALSE))
  new <- data.frame(x=first$worklist$x,day="C",
    RunID=paste0("followup-",seq_len(nrow(first$worklist))))
  new$y <- 3+new$x+.6+.03*cos(seq_len(nrow(new)))
  all <- rbind(transform(d,day=as.character(day)),new)
  all$day <- factor(all$day,levels=c("A","B","C"))
  refreshed <- bigexp_terms(y~x,all,blocking="day",report=FALSE)
  next_pool <- pool
  next_pool$day <- factor("C",levels=levels(all$day))
  next_plan <- svem_plan_runs(refreshed,all,goals,next_pool,budget=3,
    benchmark=list(RunID="original-1"),fit_control=fit_args,
    forward_control=list(nBoot=3),score_control=list(verbose=FALSE))
  expect_equal(next_plan$fits$y$nobs,nrow(all))
  expect_equal(next_plan$forward$y$nobs,nrow(all))
  expect_true(all(next_plan$incumbent$source_rows > nrow(d)))
  expect_identical(next_plan$incumbent$source_context,"day=C")
  expect_identical(next_plan$benchmark$source_context,"day=A")
  expect_true(all(next_plan$worklist$prediction_context=="day=C"))
  expect_true(all(as.character(next_plan$worklist$day)=="C"))
  expect_equal(next_plan$benchmark$source_row,1)
  expect_equal(next_plan$worklist$x[next_plan$worklist$recipe_id==next_plan$benchmark$recipe_id],d$x[1])
  unseen <- next_pool
  unseen$day <- "D"
  before <- .Random.seed
  expect_error(svem_plan_runs(refreshed,all,goals,unseen,budget=3,
    fits=next_plan$fits,forward=next_plan$forward), "Unseen level|incompatible")
  expect_identical(.Random.seed,before)
})

