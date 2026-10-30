# Conditional Gaussian elastic-net trace. This is not the degrees of freedom
# of the response-adaptive alpha/lambda/gamma/validation search.
.svem_edf_context <- function(X, y, w, args = list(), nvars = ncol(X)) {
  n <- nrow(X); sw <- sum(w)
  ym <- sum(w * y) / sw
  ys <- sqrt(sum(w * (y - ym)^2) / sw)
  mw <- colSums(X * w) / sw
  xc <- sweep(X, 2L, mw, "-")
  sx <- sqrt(colSums(xc^2 * w) / sw)
  sx[!is.finite(sx) | sx <= 0] <- 1
  z <- sweep(xc, 2L, sx, "/") * sqrt(w * n / sw)
  pf <- args$penalty.factor
  if (is.null(pf)) pf <- rep(1, nvars)
  reason <- NULL
  if (length(pf) != nvars || anyNA(pf) || any(pf == -Inf)) {
    reason <- "invalid penalty.factor for conditional EDF"
  } else {
    pf[pf == Inf] <- 1  # glmnet excludes these columns before normalizing
    ex <- args$exclude
    if (is.function(ex)) {
      if (any(pf != 1)) reason <- "conditional EDF with nonuniform penalties requires a fixed exclude vector"
    } else if (length(ex)) {
      if (is.logical(ex)) ex <- which(ex)
      if (any(!is.finite(ex)) || any(ex < 1 | ex > nvars | ex != floor(ex))) {
        reason <- "invalid exclude vector for conditional EDF"
      } else pf[ex] <- 1
    }
    pf <- pmax(pf, 0)
    if (sum(pf) <= 0) reason <- "conditional EDF requires a positive total penalty factor"
    else pf <- pf * nvars / sum(pf)
  }
  bounded <- any(is.finite(args$lower.limits)) || any(is.finite(args$upper.limits))
  if (bounded) reason <- "conditional EDF is not implemented for finite coefficient bounds"
  list(n = n, y_sd = ys, z = z, penalty = pf[seq_len(ncol(X))],
       cache = new.env(parent = emptyenv()), reason = reason)
}

.svem_active <- function(beta, tol = 0) {
  tol <- .svem_numeric_scalar(tol, "tol", lower = 0)
  if (tol == 0) return(.svem_solver_active(beta))
  # Explicit positive tolerances retain the private research harness's
  # historical raw-coefficient convention; production uses exact support.
  slope <- beta[-1L]
  finite <- is.finite(slope)
  m <- if (any(finite)) max(abs(slope[finite])) else 0
  which(finite & abs(slope) > tol * max(1, m))
}

.svem_conditional_edf <- function(context, coefficients, lambda, alpha, gamma = 1,
                                  coef_tol = 0) {
  coef_tol <- .svem_numeric_scalar(coef_tol, "coef_tol", lower = 0)
  if (!is.matrix(coefficients)) coefficients <- matrix(coefficients, ncol = 1L)
  L <- ncol(coefficients)
  if (length(lambda) != L) stop("EDF lambda and coefficient paths are misaligned.")
  if (!is.null(context$reason)) return(rep(NA_real_, L))
  vapply(seq_len(L), function(j) {
    act <- .svem_active(coefficients[, j], coef_tol)
    if (!length(act)) return(1)
    key <- paste(act, collapse = ",")
    cached <- context$cache[[key]]
    if (is.null(cached)) {
      G <- crossprod(context$z[, act, drop = FALSE])
      ev <- pmax(eigen(G, symmetric = TRUE, only.values = TRUE)$values, 0)
      rank <- if (length(ev) && ev[1L] > 0) sum(ev > ev[1L] * 1e-10) else 0
      pf <- context$penalty[act]
      if (all(pf > 0)) {
        H <- G / outer(sqrt(pf), sqrt(pf))
        ridge_ev <- pmax(eigen(H, symmetric = TRUE, only.values = TRUE)$values, 0)
      } else ridge_ev <- NULL
      cached <- list(G = G, rank = rank, pf = pf, ridge_ev = ridge_ev)
      context$cache[[key]] <- cached
    }
    # The conditional lasso projection has the active design's rank, which
    # equals support only when its centered columns are independent.
    if (alpha == 1) return(as.numeric(1L + cached$rank))
    if (!is.finite(context$y_sd) || context$y_sd <= 0) return(NA_real_)
    # glmnet reports Gaussian lambda multiplied by its weighted response SD.
    rho <- context$n * (lambda[j] / context$y_sd) * (1 - alpha)
    if (!is.finite(rho) || rho < 0) return(NA_real_)
    if (rho == 0 || gamma == 0) {
      df <- cached$rank
    } else if (!is.null(cached$ridge_ev)) {
      df <- sum(cached$ridge_ev / (cached$ridge_ev + rho))
    } else {
      H <- cached$G + diag(rho * cached$pf, length(act))
      eg <- eigen(H, symmetric = TRUE)
      ok <- eg$values > max(eg$values) * 1e-12
      if (!any(ok)) df <- 0 else {
        V <- eg$vectors[, ok, drop = FALSE]
        df <- sum(colSums(V * (cached$G %*% V)) / eg$values[ok])
      }
    }
    as.numeric(1 + gamma * df + (1 - gamma) * cached$rank)
  }, numeric(1))
}
