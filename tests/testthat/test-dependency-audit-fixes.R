# Focused regressions for the dependency-audit repairs.

test_that("SVEM requests the complete path and preserves explicit controls", {
  set.seed(5);X<-matrix(rnorm(240),60,4);y<-3*X[,1]-X[,2]+rnorm(60,sd=.001)
  setup<-.svem_glmnet_setup(X,list(),"gaussian",FALSE)
  defaults<-setup$args$control
  if(is.null(defaults))defaults<-attr(setup$args,"svem_legacy_control")
  expect_equal(defaults$fdev,0);expect_equal(defaults$devmax,1)
  full<-.svem_glmnet_path(setup$design$x,y,rep(1,60),1,setup$args,1,setup$design$sentinel)
  legacy<-.svem_glmnet_setup(X,list(control=list(fdev=1e-5,devmax=.999)),"gaussian",FALSE)
  short<-.svem_glmnet_path(legacy$design$x,y,rep(1,60),1,legacy$args,1,legacy$design$sentinel)
  expect_length(full$fit$lambda,500)
  expect_gt(length(full$fit$lambda),length(short$fit$lambda))
  expect_lt(tail(full$fit$lambda,1),tail(short$fit$lambda,1))
})

test_that("legacy scoped controls restore session state on success and error", {
  old<-glmnet::glmnet.control();on.exit(do.call(glmnet::glmnet.control,old))
  glmnet::glmnet.control(fdev=.02,devmax=.8);before<-glmnet::glmnet.control()
  args<-.glmnet_prepare_call_args(list(family="gaussian",nlambda=20),function(x,y,...)NULL,
    default_control=list(maxit=1e6,fdev=0,devmax=1),scoped_legacy=TRUE)
  expect_null(args$fdev);expect_null(args$devmax)
  expect_equal(attr(args,"svem_legacy_control"),list(fdev=0,devmax=1))
  set.seed(22);X<-matrix(rnorm(80),20,4)
  fit<-suppressWarnings(.svem_glmnet_call(X,X[,1]+rnorm(20),rep(1,20),1,args))
  expect_length(fit$lambda,20);expect_identical(glmnet::glmnet.control(),before)
  bad<-args;bad$family<-"not-a-family"
  expect_error(.svem_glmnet_call(X,X[,1],rep(1,20),1,bad))
  expect_identical(glmnet::glmnet.control(),before)
})

test_that("stored coefficient extraction matches glmnet including relaxed endpoints", {
  set.seed(91);X<-matrix(rnorm(120),30,4);colnames(X)<-paste0("x",1:4)
  for(fam in c("gaussian","binomial"))for(relax in c(FALSE,TRUE))for(exclude in list(NULL,2L)){
    y<-if(fam=="gaussian")X[,1]-X[,2]+rnorm(30)else rep(0:1,15)
    setup<-.svem_glmnet_setup(X,c(list(exclude=exclude),if(relax)list(maxp=2L)),fam,relax)
    b<-.svem_glmnet_path(setup$design$x,y,seq_len(30),.7,setup$args,c(0,.2,.6,1),setup$design$sentinel)
    expect_length(b$paths,4)
    if(relax){
      expect_lt(length(b$fit$relaxed$lambda),length(b$fit$lambda))
      expect_true(any(!b$fit$lambda %in% b$fit$relaxed$lambda))
    }
    for(p in b$paths){
      ref<-if(relax)as.matrix(coef(b$fit,s=b$fit$lambda,gamma=p$gamma))else as.matrix(coef(b$fit,s=b$fit$lambda))
      expect_equal(unname(p$coefficients),unname(ref),tolerance=0)
      expect_identical(dimnames(p$coefficients),dimnames(ref))
    }
  }
  one<-matrix(seq(-1,1,length.out=20),ncol=1,dimnames=list(NULL,"one"))
  setup<-.svem_glmnet_setup(one,list(),"gaussian",TRUE)
  b<-.svem_glmnet_path(setup$design$x,sin(one[,1]),rep(1,20),1,setup$args,c(.2,1),setup$design$sentinel)
  expect_true(all(vapply(b$paths,function(p)nrow(p$coefficients)==2L,logical(1))))
  for(p in b$paths)expect_equal(unname(p$coefficients),
    unname(as.matrix(coef(b$fit,s=b$fit$lambda,gamma=p$gamma))[1:2,,drop=FALSE]),tolerance=0)
  # A singleton relaxed beta must retain its predictor dimension (upstream n=3 edge).
  base<-matrix(1:12,3,4,dimnames=list(c("(Intercept)","a","b"),paste0("s",1:4)))
  fit<-structure(list(lambda=4:1,relaxed=list(lambda=4,a0=2,beta=c(3,4))),class="relaxed")
  expect_equal(.svem_glmnet_blend(fit,base,.2),.2*base+.8*matrix(c(2,3,4),3,4),ignore_attr=TRUE)
})

test_that("GP response normalization preserves means and latent SD across units", {
  skip_if_not_installed("GauPro")
  set.seed(4);d<-data.frame(x=runif(30),z=runif(30),w=runif(30))
  d$y<-with(d,3*sin(3*x)+2*z^2-w*z)+rnorm(30,sd=.1)
  new<-data.frame(x=runif(12),z=runif(12),w=runif(12))
  fits<-lapply(list(c(1,0),c(1e5,0),c(1,1e6)),function(a){
    dd<-d;dd$y<-a[1]*d$y+a[2];set.seed(44)
    m<-SVEMnet(y~x+z+w,dd,nBoot=2,glmnet_alpha=1,relaxed=FALSE,companions="gaupro")
    expect_identical(m$companions$gaupro$status,"ok")
    expect_true(m$companions$gaupro$settings$normalize)
    p<-predict(m,new,engine="gaupro",se.fit=TRUE)
    list(mean=(p$fit-a[2])/a[1],se=p$se.fit/a[1])
  })
  for(p in fits[-1]){
    expect_equal(p$mean,fits[[1]]$mean,tolerance=1e-5)
    # Finite-difference GP optimization is not bitwise affine equivariant.
    expect_equal(p$se,fits[[1]]$se,tolerance=1e-3)
  }
})

test_that("unsupported factor-only GP cannot report a usable fit", {
  skip_if_not_installed("GauPro")
  d<-data.frame(f=factor(rep(letters[1:3],each=5)),y=rep(1:3,each=5)+sin(1:15)/10)
  expect_warning(m<-SVEMnet(y~f,d,nBoot=2,companions="gaupro"),"single factor")
  expect_identical(m$companions$gaupro$status,"unsupported")
  expect_true(all(is.finite(predict(m,d))))
  expect_error(predict(m,d,engine="gaupro"),"unsupported")
})

test_that("logical Gower factors use symmetric matches", {
  d<-data.frame(score=1,x=c(0,.1,.2,.8,1),flag=c(FALSE,FALSE,FALSE,TRUE,TRUE))
  distance<-(abs(outer(d$x,d$x,"-"))+outer(d$flag,d$flag,"!="))/2
  optimum<-which.min(rowSums(distance))
  expect_identical(svem_select_candidates(d,"score",1,1,c("x","flag")),as.integer(optimum))
  d$flag<-factor(d$flag)
  expect_identical(svem_select_candidates(d,"score",1,1,c("x","flag")),as.integer(optimum))
})

test_that("huge PAM shortlists fail before distance allocation", {
  local_mocked_bindings(daisy=function(...)stop("unexpected quadratic allocation"),.package="cluster")
  d<-data.frame(score=rep(1,65537),x=seq_len(65537))
  expect_error(svem_select_candidates(d,"score",1,2,"x"),"65,536")
})


test_that("GP fit messages are retained and constant responses are explicit", {
  skip_if_not_installed("GauPro")
  actual<-GauPro::gpkm
  local_mocked_bindings(gpkm=function(...) {message("fit diagnostic sentinel");actual(...)},.package="GauPro")
  d<-data.frame(x=seq(0,1,length.out=15));d$y<-sin(d$x*4)+cos(seq_len(15))/10
  set.seed(42);m<-SVEMnet(y~x,d,nBoot=2,glmnet_alpha=1,relaxed=FALSE,companions="gaupro")
  expect_identical(m$companions$gaupro$status,"ok")
  expect_true(any(grepl("fit diagnostic sentinel",m$companions$gaupro$warnings)))
  d$y<-rep(2,15)
  expect_warning(m<-SVEMnet(y~x,d,nBoot=2,companions="gaupro"),"must vary")
  expect_identical(m$companions$gaupro$status,"unsupported")
  expect_true(all(is.finite(predict(m,d))))
})


test_that("small-sample relaxed paths retain every gamma despite singleton refits", {
  X<-cbind(a=c(-1,0,1),b=c(0,1,-1));y<-c(0,1,3)
  setup<-.svem_glmnet_setup(X,list(),"gaussian",TRUE)
  set.seed(17)
  b<-.svem_glmnet_path(setup$design$x,y,rep(1,3),1,setup$args,c(.2,.6,1),setup$design$sentinel)
  expect_length(b$errors,0L);expect_length(b$paths,3L)
  expect_true(all(vapply(b$paths,function(p)all(is.finite(p$coefficients)),logical(1))))
})


test_that("coefficient extraction errors keep per-member fallback behavior", {
  set.seed(71);d<-data.frame(x=rnorm(25),z=rnorm(25));d$y<-d$x-d$z+rnorm(25)/10
  actual<-.svem_glmnet_coefficient_matrix;calls<-0L
  local_mocked_bindings(.svem_glmnet_coefficient_matrix=function(fit){
    calls<<-calls+1L;if(calls==1L)stop("extraction sentinel");actual(fit)
  })
  set.seed(81)
  expect_warning(m<-SVEMnet(y~x+z,d,nBoot=3,glmnet_alpha=1,relaxed=FALSE),"1 of 3 bootstrap")
  expect_equal(m$diagnostics$fallback_rate,1/3)
  expect_match(names(m$diagnostics$fallback_reasons),"extraction sentinel")
  expect_true(all(is.finite(predict(m,d))))
})

test_that("solver warnings and effective controls remain inspectable", {
  actual<-.svem_glmnet_call
  local_mocked_bindings(.svem_glmnet_call=function(...) {warning("solver sentinel");actual(...)})
  d<-data.frame(x=seq(-1,1,length.out=20));d$y<-d$x+cos(1:20)/10
  set.seed(7);m<-SVEMnet(y~x,d,nBoot=2,glmnet_alpha=1,relaxed=FALSE,control=list(fdev=.001))
  expect_equal(m$diagnostics$solver_warning_count,2L)
  expect_true(all(grepl("solver sentinel",m$diagnostics$solver_warnings)))
  expect_equal(m$diagnostics$glmnet_control$fdev,.001)
  expect_equal(m$diagnostics$glmnet_control$devmax,1)
  expect_identical(m$diagnostics$glmnet_version,as.character(packageVersion("glmnet")))
})

test_that("GP prediction probe failures are not usable companion fits", {
  skip_if_not_installed("GauPro")
  local_mocked_bindings(gpkm=function(...)list(predict=function(...) {
    message("probe diagnostic");stop("probe sentinel")
  }),.package="GauPro")
  d<-data.frame(x=seq(0,1,length.out=15));d$y<-sin(4*d$x)+cos(1:15)/10
  expect_warning(m<-SVEMnet(y~x,d,nBoot=2,companions="gaupro"),"probe sentinel")
  expect_identical(m$companions$gaupro$status,"failed")
  expect_identical(m$companions$gaupro$warnings,"probe diagnostic")
  expect_output(print(m$companions),"1 fit note")
  expect_true(all(is.finite(predict(m,d))))
})

test_that("only documented relaxed controls pass the SVEMnet whitelist", {
  X<-cbind(a=seq(-1,1,length.out=20),b=cos(1:20))
  expect_warning(setup<-.svem_glmnet_setup(X,list(check.args=FALSE),"gaussian",TRUE),"check.args")
  expect_null(setup$args$check.args)
  setup<-.svem_glmnet_setup(X,list(maxp=1,path=TRUE),"gaussian",TRUE)
  expect_equal(setup$args$maxp,1);expect_true(setup$args$path)
})


test_that("binomial retains its former path defaults and can explicitly opt in", {
  X<-cbind(a=seq(-1,1,length.out=20),b=cos(1:20))
  args<-.svem_glmnet_setup(X,list(),"binomial",FALSE)$args
  control<-if(!is.null(args$control))args$control else attr(args,"svem_legacy_control")
  expect_null(control$fdev);expect_null(control$devmax)
  args<-.svem_glmnet_setup(X,list(control=list(fdev=0,devmax=1)),"binomial",FALSE)$args
  control<-if(!is.null(args$control))args$control else attr(args,"svem_legacy_control")
  expect_equal(control$fdev,0);expect_equal(control$devmax,1)
})


test_that("WMT relaxed controls are filtered on the master by its relaxed flag", {
  allowed<-NULL
  local_mocked_bindings(.svem_unknown_glmnet_args=function(...,extra_allowed=character()) {
    allowed<<-extra_allowed;stop("master preflight sentinel")
  })
  d<-data.frame(x=1:10,y=sin(1:10))
  for(relax in c(FALSE,TRUE)){
    expect_error(svem_significance_test_parallel(y~x,d,relaxed=relax,maxp=1,
      nPerm=20,nPoint=4,nSVEM=1,nBoot=2,verbose=FALSE),"master preflight sentinel")
    expect_identical(allowed,if(relax)c("maxp","path")else character())
  }
})
