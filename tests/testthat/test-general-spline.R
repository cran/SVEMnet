
general_mixture_fixture <- function(n = 360L) {
  set.seed(12417)
  x <- matrix(rexp(n * 4L), n, 4L)
  x <- .1 + .6 * x / rowSums(x)
  g <- data.frame(A=x[,1],B=x[,2],C=x[,3],D=x[,4],
    ratio=rep(c(6,8,10,12),length.out=n),flow=rep(c(1,2,3,2,1),length.out=n),
    kind=factor(rep(c("K1","K2","K3"),length.out=n)))
  mix <- list(list(vars=c("A","B","C","D"),lower=rep(.1,4),upper=rep(.7,4),total=1))
  d <- g[seq_len(36),]
  d$day <- factor(rep(c("D0","D1"),each=18))
  d$y <- with(d, sin(6*A)+.15*ratio+.2*(kind=="K2")*flow+.3*(day=="D1"))
  list(grid=g,data=d,mix=mix)
}

test_that("general mixture/process/nominal fits solve the weighted Gaussian equations", {
  skip_if_not_installed("mgcv")
  v <- general_mixture_fixture()
  s <- spline_terms(y~A+B+C+D+ratio+flow+kind,v$grid,mixture_groups=v$mix,blocking="day")
  expect_identical(s$version,2L)
  expect_identical(s$dimension,5L)
  expect_identical(s$smoother,"Duchon")
  set.seed(81);f <- svem_spline(s,v$data,nBoot=3,store_member_weights=TRUE,lambda=c(Inf,10,1,.1))
  x <- cbind(1,.svem_spline_design(s,v$data,f$blocks))
  z <- (v$data$y-f$response_center)/f$response_scale
  penalty <- diag(c(0,f$penalty))
  for (i in 1:3) {
    w <- f$member_weights$training[i,];vv <- f$member_weights$validation[i,];wt <- w/sum(w)
    cap <- max(2,min(length(w),sum(vv)^2/(sum(vv^2)+.Machine$double.eps)))
    loss <- rep(Inf,length(f$lambda_grid));co <- matrix(0,ncol(x),length(loss))
    for (j in seq_along(loss)) {
      lam <- f$lambda_grid[j]
      if (is.infinite(lam)) {co[1,j] <- sum(wt*z);df <- 0} else {
        lhs <- crossprod(x,x*wt)+lam*penalty
        co[,j] <- solve(lhs,crossprod(x,wt*z))
        df <- sum(diag(solve(lhs,crossprod(x,x*wt))))-1
      }
      e2 <- as.numeric(z-x%*%co[,j])^2;t <- max(sum(wt*e2),.Machine$double.eps^2)
      if (df<cap) loss[j] <- log(t)+sum(vv*e2)/sum(vv)/t
    }
    expect_equal(f$path_scores[i,],loss,tolerance=1e-7)
    expect_equal(as.numeric(x%*%f$coef_matrix[i,]),as.numeric(x%*%co[,which.min(loss)]),tolerance=1e-7)
    expect_lt(f$member_diagnostics$slope_edf[i],cap)
  }
  ga <- transform(v$grid,day="D0");gb <- transform(v$grid,day="D1")
  expect_lt(diff(range(predict(f,gb)-predict(f,ga))),1e-10)
  expect_equal(predict(f,ga[1,,drop=FALSE]),predict(f,ga)[1],tolerance=1e-9)
  expect_equal(predict(unserialize(serialize(f,NULL)),ga),predict(f,ga),tolerance=0)
  expect_equal(predict(f,ga[FALSE,]),numeric())
  expect_true(all(is.finite(predict(f,ga))))
  expect_true(all(tail(f$penalty,1)==1))
})

test_that("multiple closures and exact fixed components preserve the prediction kernel", {
  skip_if_not_installed("mgcv")
  g0 <- expand.grid(u=0:6/6,v=0:6/6)
  g <- data.frame(A=.03+.04*g0$u,B=.12+.4*g0$v);g$C <- .9-g$A-g$B
  a <- spline_terms(y~A+B+C,g,mixture_groups=list(list(vars=c("A","B","C"),total=.9)))
  z <- transform(g,C=C+.1,Z=0)
  b <- spline_terms(y~A+B+C+Z,z,mixture_groups=list(list(vars=c("A","B","C","Z"),total=1)))
  xa <- .svem_spline_surface(a,g);xb <- .svem_spline_surface(b,z)
  expect_equal(xa%*%(t(xa)/a$penalty),xb%*%(t(xb)/b$penalty),tolerance=1e-7)
  d <- g[seq(1,nrow(g),by=2),];d$y <- sin(5*d$A)+d$B
  dz <- transform(d,C=C+.1,Z=0)
  set.seed(19);fa <- svem_spline(a,d,nBoot=3)
  set.seed(19);fb <- svem_spline(b,dz,nBoot=3)
  expect_equal(predict(fa,g),predict(fb,z),tolerance=1e-7)
  set.seed(56);q <- data.frame(a=runif(180,.1,.4),b=runif(180,.1,.4),c=0,
    e=runif(180,1,2),f=0,p=runif(180))
  q$c <- 1-q$a-q$b;q$f <- 5-q$e
  groups <- list(list(vars=c("a","b","c"),total=1),list(vars=c("e","f"),lower=c(0,0),upper=c(5,5),total=5))
  s <- spline_terms(y~.,q,mixture_groups=groups)
  expect_identical(s$dimension,4L)
  expect_true(all(.svem_spline_legal(s,q)))
  wrong <- q[1,,drop=FALSE];wrong$f <- wrong$f+.1
  expect_false(.svem_spline_legal(s,wrong))
})

test_that("nominal labels have no numeric trend and unmeasured levels are explicit", {
  skip_if_not_installed("mgcv")
  g <- expand.grid(u=0:5/5,v=0:5/5,kind=factor(c("A","B","C")))
  d <- g[seq(1,nrow(g),by=2),];d$y <- with(d,u+.3*v+.5*(kind=="B")*v)
  s <- spline_terms(y~u+v+kind,g)
  set.seed(5);f <- svem_spline(s,d,nBoot=3)
  labels <- c(A="10",B="2",C="1")
  g2 <- g;g2$kind <- factor(labels[as.character(g$kind)],levels=c("2","1","10"))
  d2 <- d;d2$kind <- factor(labels[as.character(d$kind)],levels=c("2","1","10"))
  set.seed(5);f2 <- svem_spline(spline_terms(y~u+v+kind,g2),d2,nBoot=3)
  expect_equal(predict(f,g),predict(f2,g2),tolerance=1e-7)
  set.seed(7);partial <- svem_spline(s,d[d$kind!="C",],nBoot=2)
  expect_error(predict(partial,g[g$kind=="C",]),"unmeasured")
  expect_true(all(!.svem_spline_general_covered(s,partial$nominal_levels,g[g$kind=="C",])))
  expect_true(all(.svem_spline_general_covered(s,partial$nominal_levels,g[g$kind=="A",])))
})

test_that("all-nominal and fixed-coordinate limits are honest Gaussian surfaces", {
  skip_if_not_installed("mgcv")
  g <- expand.grid(a=factor(c("a","b")),b=factor(c("a","b")),c=factor(c("a","b")))
  d <- g[rep(1:8,4),];d$y <- with(d,2+2*(2*(a=="b")-1)*(2*(b=="b")-1)*(2*(c=="b")-1))
  d$y <- d$y+.01*sin(seq_len(nrow(d)))
  s <- spline_terms(y~.,g)
  expect_identical(s$dimension,0L)
  expect_identical(s$smoother,"categorical shrinkage")
  expect_true(any(lengths(s$cells)==3L))
  set.seed(12);f <- svem_spline(s,d,nBoot=3)
  expect_gt(diff(range(predict(f,g))),2)
  expect_true(all(is.finite(predict(f,g))))
  fixed <- data.frame(x=rep(.15,8),y=1:8)
  z <- spline_terms(y~x,data.frame(x=.15))
  set.seed(3);h <- svem_spline(z,fixed,nBoot=2,store_member_weights=TRUE)
  expect_true(all(h$member_diagnostics$intercept_only))
  expect_equal(predict(h,data.frame(x=.15),type="members")[1,],
    drop(h$member_weights$training%*%fixed$y/rowSums(h$member_weights$training)),tolerance=1e-12)
  expect_error(predict(h,data.frame(x=.2)),"domain")
  binary <- spline_terms(y~x,data.frame(x=c(0,1)))
  expect_identical(binary$smoother,"linear limit")
})

test_that("narrow numeric units are not mistaken for fixed coordinates", {
  skip_if_not_installed("mgcv")
  set.seed(84);g <- as.data.frame(matrix(runif(480),120,4));names(g) <- letters[1:4]
  tiny <- g;tiny$a <- tiny$a*1e-12
  s <- spline_terms(y~.,g);ss <- spline_terms(y~.,tiny)
  expect_identical(ss$dimension,4L);expect_length(ss$fixed,0)
  d <- g[1:30,];d$y <- with(d,sin(3*a)+b*c)
  dt <- tiny[1:30,];dt$y <- 10000+100*d$y
  set.seed(16);f <- svem_spline(s,d,nBoot=3)
  set.seed(16);ff <- svem_spline(ss,dt,nBoot=3)
  expect_equal((predict(ff,tiny)-10000)/100,predict(f,g),tolerance=1e-7)
  expect_identical(f$member_diagnostics$grid_index,ff$member_diagnostics$grid_index)
})

test_that("general bases expose resource reductions without a low dimension cutoff", {
  skip_if_not_installed("mgcv")
  set.seed(715);g <- as.data.frame(matrix(runif(1300),100,13));names(g) <- paste0("x",1:13)
  s <- spline_terms(y~.,g)
  expect_identical(s$dimension,13L)
  expect_lte(length(s$penalty),200)
  mixed <- expand.grid(u=0:5/5,v=0:5/5,kind=factor(c("a","b","c")))
  low <- spline_terms(y~.,mixed,max_columns=30)
  expect_lte(length(low$penalty),30)
  expect_gt(length(low$basis_notes),0)
  expect_error(spline_terms(y~.,mixed,max_columns=2),"integer|columns")
  expect_error(spline_terms(y~u+v,mixed,max_columns="200"),"integer")
  expect_error(spline_terms(y~u+v,mixed,interaction_order=4),"integer")
})

test_that("flat general fits do not invent an optimum but constant endpoints can join a goal", {
  skip_if_not_installed("mgcv")
  g <- expand.grid(u=0:5/5,kind=factor(c("a","b")))
  d <- transform(g,y=3,z=u+.3*(kind=="b"))
  s <- spline_terms(y~u+kind,g)
  f <- svem_spline(s,d,nBoot=2)
  h <- svem_spline(s,d,response="z",nBoot=2)
  des <- function(x,p) pmin(1,pmax(.01,x/5))
  one <- .svem_spline_scores(data.frame(id=seq_len(nrow(g))),g,list(y=f),list(y=list()),des,c(y=1),"geom",1e-12)
  expect_true(all(is.na(one$table$spline_score)))
  expect_identical(one$status$status,"no_signal")
  both <- .svem_spline_scores(data.frame(id=seq_len(nrow(g))),g,list(y=f,z=h),list(y=list(),z=list()),des,c(y=.2,z=.8),"geom",1e-12)
  expect_true(all(is.finite(both$table$spline_score)))
  expect_gt(diff(range(both$table$spline_score)),0)
})


test_that("the shared formulation interface fits, scores and reuses general spline nominees", {
  skip_if_not_installed("mgcv")
  v <- general_mixture_fixture(180L)
  pool <- transform(v$grid,day="D1")
  domain <- pool;domain$y <- NA_real_
  dd <- rbind(v$data,domain[names(v$data)])
  spec <- bigexp_terms(y~A+B+C+D+ratio+flow+kind,dd,mixture_groups=v$mix,
    blocking="day",discrete_numeric=c("ratio","flow"),report=FALSE)
  fits <- svem_fit_responses(spec,v$data,"y",method="spline",domain_grid=pool,nBoot=3)
  expect_s3_class(fits$y,"svem_spline")
  expect_identical(fits$y$spec$dimension,5L)
  set.seed(128)
  expect_warning(plan <- svem_plan_runs(spec,v$data,list(y=list(goal="max",weight=1)),pool,n_optima=5,n_explore=2,
    fit_control=list(nBoot=3,companions=FALSE),forward_control=list(nBoot=3),
    spline_control=list(nBoot=3),reference=list(day="D1"),score_control=list(verbose=FALSE)), "Requested allocation is short")
  expect_s3_class(plan$spline$y,"svem_spline")
  expect_identical(subset(plan$scores$engine_status,engine=="spline")$status,"ok")
  expect_equal(nrow(plan$worklist)+plan$shortfall,8)
  expect_false(anyDuplicated(plan$worklist$recipe_id)>0)
  seed <- .Random.seed
  expect_warning(reused <- svem_plan_runs(spec,v$data,list(y=list(goal="max",weight=1)),pool,n_optima=5,n_explore=2,
    fits=plan$fits,forward=plan$forward,spline=plan$spline,reference=list(day="D1"),score_control=list(verbose=FALSE)), "Requested allocation is short")
  expect_identical(.Random.seed,seed)
  expect_identical(reused$worklist,plan$worklist)
  expect_equal(reused$worklist$y_spline_pred,predict(plan$spline$y,reused$worklist),tolerance=1e-9)
  expect_error(.svem_spline_check_training(plan$spline$y,transform(v$data,y=y+.01)),"differ")
  # Known but unmeasured nominal levels retain partial, explicit coverage.
  observed <- v$data[v$data$kind!="K3",]
  partial <- svem_fit_responses(spec,observed,"y",method="spline",domain_grid=pool,nBoot=2)
  primary <- svem_fit_responses(spec,observed,"y",nBoot=2,companions=FALSE)
  scored <- svem_score_random(primary,data=observed,goals=list(y=list(goal="max",weight=1)),candidates=pool,
    spline=partial,allocation="portfolio",budget=3L,verbose=FALSE)
  tab <- scored$score_table
  expect_true(all(is.na(tab$y_spline_pred[tab$kind=="K3"])))
  expect_true(all(is.finite(tab$y_spline_pred[tab$kind!="K3"])))
  expect_identical(subset(scored$engine_status,engine=="spline")$status,"partial_coverage")
})


test_that("rounding noise is fixed without erasing meaningful narrow units", {
  skip_if_not_installed("mgcv")
  g <- expand.grid(u=0:5/5,v=0:5/5)
  g$fixed <- .1+sin(seq_len(nrow(g)))*.Machine$double.eps/8
  sp <- spline_terms(y~u+v+fixed,g)
  expect_identical(sp$dimension,2L)
  expect_identical(sp$fixed,"fixed")
  tiny <- transform(g,small=u*1e-12)
  # The independent small process is intentionally collinear with u here.
  # Detection must retain it so the span check can catch the design deficiency.
  expect_error(spline_terms(y~u+v+small,tiny),"span")
})

test_that("a flat general surface is excluded by the public raw-mean nomination path", {
  skip_if_not_installed("mgcv")
  g <- expand.grid(u=0:5/5,kind=factor(c("a","b")))
  d <- transform(g,y=3)
  sp <- bigexp_terms(y~u+kind,d,report=FALSE)
  primary <- svem_fit_responses(sp,d,"y",nBoot=2,companions=FALSE)
  smooth <- svem_fit_responses(sp,d,"y",method="spline",domain_grid=g,nBoot=2)
  pool <- g;pool$u <- pool$u+.01
  # Keep candidates within the primary legal range, with genuinely new recipes.
  pool <- pool[pool$u<=1,]
  expect_warning(out <- svem_score_random(primary,list(y=list(goal="max",weight=1)),
    data=d,candidates=pool,spline=smooth,budget=1,allocation="portfolio",priority="spline",
    incumbent=FALSE,n_optima=1,n_explore=0,verbose=FALSE),"Requested allocation is short")
  row <- subset(out$nominations,engine=="spline")
  expect_identical(row$status,"no_finite_score")
  expect_true(is.na(row$new_recipe_id))
  expect_equal(nrow(out$worklist),0)
})


test_that("saved marginal curvature is reproducible and EDF allocations sum correctly", {
  skip_if_not_installed("mgcv")
  v <- general_mixture_fixture()
  sp <- spline_terms(y~A+B+C+D+ratio+flow+kind,v$grid,mixture_groups=v$mix,blocking="day")
  expect_gt(length(sp$marginal_specs),0)
  expect_equal(sp$marginal_weight,.25)
  set.seed(17);fit <- svem_spline(sp,v$data,nBoot=4)
  edf <- fit$member_diagnostics[grep("^edf_",names(fit$member_diagnostics))]
  expect_equal(rowSums(edf),fit$member_diagnostics$slope_edf,tolerance=1e-10)
  expect_true("edf_marginal_curvature" %in% names(edf))
  pool <- transform(v$grid,day="D1")
  whole <- predict(fit,pool)
  pieces <- numeric(nrow(pool))
  for (i in split(seq_len(nrow(pool)),rep(1:4,length.out=nrow(pool))))
    pieces[i] <- predict(fit,pool[i,])
  expect_equal(whole,pieces,tolerance=1e-10)
  expect_identical(predict(unserialize(serialize(fit,NULL)),pool),whole)
  lean <- spline_terms(y~A+B+C+D+ratio+flow+kind,v$grid,mixture_groups=v$mix,marginal_modes=0)
  expect_null(lean$marginal_specs)
  expect_lte(length(sp$penalty),sp$max_columns)
})

test_that("numeric slopes can differ across two nominal factors without a numeric label trend", {
  skip_if_not_installed("mgcv")
  g <- expand.grid(u=seq(0,1,length.out=8),a=factor(c("A","B")),b=factor(c("low","high")))
  d <- g;d$y <- with(d,.2*u+2*u*(a=="B")*(b=="high"))+.01*sin(seq_len(nrow(d)))
  sp <- spline_terms(y~u+a+b,g)
  expect_gt(length(sp$pair_modes),0)
  set.seed(12);fit <- svem_spline(sp,d,nBoot=10)
  pred <- predict(fit,g)
  expect_lt(sqrt(mean((pred-d$y)^2)),.15)
  g2 <- g;g2$a <- factor(ifelse(g$a=="A","z","a"));g2$b <- factor(ifelse(g$b=="low","9","1"))
  d2 <- g2;d2$y <- d$y
  set.seed(12);fit2 <- svem_spline(spline_terms(y~u+a+b,g2),d2,nBoot=10)
  expect_equal(pred,predict(fit2,g2),tolerance=1e-8)
})

 test_that("general linear-limit ranks pass through the shared public fitting helper", {
  g <- expand.grid(u=seq(0,1,length.out=8),kind=factor(c("a","b")))
  d <- g;d$y <- d$u+.1*(d$kind=="b")
  sp <- bigexp_terms(y~u+kind,d,report=FALSE)
  set.seed(41)
  f <- svem_fit_responses(sp,d,"y",method="spline",domain_grid=g,rank=2,nBoot=2)
  expect_identical(f$y$spec$smoother,"linear limit")
  expect_true(all(is.finite(predict(f$y,g))))
})
 test_that("categorical and linear limits need no numeric spline dependency", {
  without_dependency <- function(fun) {
    e <- new.env(parent=environment(fun))
    e$requireNamespace <- function(...) stop("unexpected numeric spline dependency")
    environment(fun) <- e
    fun
  }
  construct <- without_dependency(.svem_spline_general_terms)
  design <- without_dependency(.svem_spline_general_design)
  g <- expand.grid(kind=factor(c("a","b")),u=0:3/3)
  categorical <- construct("y","kind",character(),g["kind"],NULL)
  expect_true(all(is.finite(design(categorical,g["kind"]))))
  linear <- construct("y",c("u","kind"),character(),g,NULL,rank=2)
  expect_identical(linear$smoother,"linear limit")
  expect_true(all(is.finite(design(linear,g))))
})