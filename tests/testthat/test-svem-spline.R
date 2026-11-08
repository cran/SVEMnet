test_that("Gaussian smoothing scores and degrees of freedom match weighted normal equations", {
  skip_if_not_installed("mgcv")
  pool <- expand.grid(u = 0:10/10, v = 0:10/10)
  spec <- spline_terms(y ~ u + v, pool)
  d <- pool[c(1, 10, 19, 33, 46, 52, 63, 72, 85, 98, 106, 118), ]
  d$y <- with(d, u - 0.5*v + cos(4*u + 3*v))
  set.seed(48)
  f <- svem_spline(spec, d, nBoot = 4, store_member_weights = TRUE)
  set.seed(48)
  expect_identical(f$member_weights$uniforms[1, ], runif(nrow(d)))
  x <- cbind(1, .svem_spline_design(spec, d, f$blocks))
  z <- (d$y - f$response_center) / f$response_scale
  pen <- diag(c(0, f$penalty))
  for (i in seq_len(f$nBoot)) {
    w <- f$member_weights$training[i, ]
    v <- f$member_weights$validation[i, ]
    cap <- max(2, min(nrow(d), sum(v)^2/(sum(v^2) + .Machine$double.eps)))
    wt <- w/sum(w)
    loss <- rep(Inf, length(f$lambda_grid))
    coef <- matrix(0, ncol(x), length(loss))
    for (j in seq_along(loss)) {
      lambda <- f$lambda_grid[j]
      if (is.finite(lambda)) {
        lhs <- crossprod(x, x*wt) + lambda*pen
        coef[, j] <- solve(lhs, crossprod(x, wt*z))
        df <- sum(diag(solve(lhs, crossprod(x, x*wt)))) - 1
      } else {
        coef[1, j] <- sum(wt*z)
        df <- 0
      }
      e2 <- as.numeric(z - x %*% coef[, j])^2
      t <- max(sum(wt*e2), .Machine$double.eps^2)
      if (df < cap) loss[j] <- log(t) + sum(v*e2)/sum(v)/t
    }
    expect_equal(f$path_scores[i, ], loss, tolerance = 1e-7)
    j <- which.min(loss)
    expect_equal(unname(f$coef_matrix[i, ]), coef[, j], tolerance = 1e-7)
    expect_true(f$member_diagnostics$slope_edf[i] < cap)
  }
  expect_true(any(is.infinite(f$path_scores)))
  expect_true(all(rowMeans(f$member_weights$training) > 1 - 1e-10))
})

test_that("saved spline coordinates and centered responses transport predictions", {
  skip_if_not_installed("mgcv")
  g <- expand.grid(a=0:8/8, b=0:8/8)
  d <- transform(g[c(1,5,9,12,19,25,30,31,44,50,53,62,69,73,80),],
                 y=3*a+cos(2*b))
  spec <- spline_terms(y ~ a + b, g)
  set.seed(71); f <- svem_spline(spec, d, nBoot=4)
  set.seed(71); h <- svem_spline(spec, transform(d, z=10+100*y), response="z", nBoot=4)
  mean <- predict(f,g); members <- predict(f,g,type="members")
  expect_equal(unname(rowMeans(members)), mean, tolerance=1e-12)
  spread <- predict(f,g,se.fit=TRUE)
  expect_equal(unname(spread$fit),mean,tolerance=1e-12)
  expect_equal(unname(spread$se.fit),unname(apply(members,1,sd)),tolerance=1e-12)
  expect_equal((predict(h,g)-10)/100,mean,tolerance=1e-9)
  expect_identical(h$member_diagnostics$grid_index,f$member_diagnostics$grid_index)
  expect_equal(predict(unserialize(serialize(f,NULL)),g[9:2,]),mean[9:2],tolerance=1e-11)
  expect_equal(predict(f,g[FALSE,]),numeric())
  set.seed(1); constant <- svem_spline(spec,transform(d,y=3),nBoot=2)
  expect_true(all(constant$member_diagnostics$intercept_only))
  expect_equal(predict(constant,g),rep(3,nrow(g)))
  expect_equal(unname(predict(constant,g,se.fit=TRUE)$se.fit),rep(0,nrow(g)))
  expect_output(print(f),"effective degrees")
  expect_output(print(spec),"reference recipes")
  expect_invisible(.svem_spline_check_training(f,d[nrow(d):1,]))
  expect_error(.svem_spline_check_training(f,transform(d,y=y+.1)),"differ")
  expect_error(predict(f,g,debias=TRUE),"FALSE")
})

test_that("fixed-sum coordinates and penalized day shifts remain separate", {
  skip_if_not_installed("mgcv")
  d0 <- expand.grid(u=0:6/6,v=0:6/6)
  g <- data.frame(A=.03+.04*d0$u,B=.12+.45*d0$v)
  g$C <- 1 - g$A - g$B
  mix <- list(list(vars=c("A","B","C"),total=1))
  s <- spline_terms(y ~ A + B + C,g,blocking="day",mixture_groups=mix)
  expect_equal(s$dimension,2)
  d <- g[seq(1,nrow(g),by=2),]
  d$day <- factor(rep(c("a","b"),length.out=nrow(d)))
  d$y <- .3*(d$day=="b") + sin(d$B*3)
  set.seed(91); fit <- svem_spline(s,d,nBoot=4)
  p <- predict(fit,transform(g,day="b")) - predict(fit,transform(g,day="a"))
  expect_lt(diff(range(p)),1e-10)
  expect_equal(tail(fit$penalty,1),1)
  expect_error(predict(fit,transform(g,day="c")),"Unknown")
  expect_error(predict(fit,transform(g,A=A+.001,day="a")),"mixture|bounds")
  expect_error(spline_terms(y ~ A + B + C,transform(g,C=C-.01),mixture_groups=mix),"total")
  bad <- transform(d,y=ifelse(seq_len(nrow(d))==1,NA_real_,y))
  expect_equal(nrow(svem_spline(s,bad,nBoot=2)$training_data),nrow(d)-1)
})

test_that("domain, response, optional nuisance and admissible dimension are explicit", {
  skip_if_not_installed("mgcv")
  q <- expand.grid(a=0:3,b=0:3,c=0:3)
  expect_equal(spline_terms(y~a+b+c,q)$rank,26)
  expect_equal(spline_terms(y~a+b,rbind(q[1:16,1:2],q[1:2,1:2]))$rank,16)
  expect_error(spline_terms(y~a+b+c+e,transform(q,e=a+b)),"span")
  expect_identical(spline_terms(y~a+b,transform(q,b=factor(b)))$version,2L)
  expect_error(spline_terms(y~a+b,transform(q,b=as.Date(b,origin="2000-01-01"))),"numeric")
  expect_error(spline_terms(y~a*b,q),"raw predictor")
  expect_error(spline_terms(y~a+b,transform(q,b=a)),"span")
  one <- spline_terms(y~a,q,blocking="day")
  d <- transform(q[1:4,],day=factor(rep("a",4)),y=c(0,2,1,4))
  expect_true(all(is.finite(predict(svem_spline(one,d,nBoot=2),d))))
  expect_error(svem_spline(one,transform(d,y=as.character(y)),nBoot=2),"numeric")
  expect_error(svem_spline(one,d,nBoot=1.5),"integer")
})

test_that("documented smoothing controls and percentile bounds use the chosen model", {
  skip_if_not_installed("mgcv")
  g <- expand.grid(u=0:8/8,v=0:8/8)
  spec <- spline_terms(y~u+v,g,rank=10,null_weight=.02)
  expect_equal(spec$rank,10)
  expect_equal(tail(spec$penalty,2),rep(.02,2))
  expect_error(spline_terms(y~u+v,g,rank=3),"integer")
  expect_error(spline_terms(y~u+v,g,null_weight=0),"null_weight")
  d <- transform(g[seq(2,nrow(g),by=3),],y=sin(3*u)+v)
  fit <- svem_spline(spec,d,nBoot=4,lambda=c(Inf,1,.1))
  expect_identical(fit$lambda_grid,c(Inf,1,.1))
  expect_error(svem_spline(spec,d,lambda=c(Inf,.1,.2)),"decrease")
  m <- predict(fit,g,type="members");int <- predict(fit,g,interval=TRUE,se.fit=TRUE,level=.8)
  expect_equal(unname(int$lwr),unname(apply(m,1,quantile,probs=.1)),tolerance=1e-12)
  expect_equal(unname(int$upr),unname(apply(m,1,quantile,probs=.9)),tolerance=1e-12)
  expect_length(predict(fit,g[1,,drop=FALSE],interval=TRUE)$lwr,1)
})

test_that("the Gaussian training variance floor is finite at interpolation", {
  x <- matrix(c(-1,-1,1,1),ncol=1L)
  y <- as.numeric(x)
  w <- rep(1,4)
  result <- .svem_spline_member(x,1,y,w,w,c(Inf,1e-30))
  expect_equal(result$t[2],0)
  expect_equal(result$score[2],log(.Machine$double.eps^2))
  expect_true(all(is.finite(result$score)))
})

