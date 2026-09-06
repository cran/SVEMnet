# SHASHo density and score formulas adapted from gamlss.dist 6.1-1,
# R/SHASHo.R (GPL-2 | GPL-3). Authors: Robert Rigby, Mikis Stasinopoulos,
# and Fiona McElduff; copyright holder: Mikis Stasinopoulos.
# See inst/COPYRIGHTS for provenance. The intercept-only optimizer below
# is new SVEMnet code. No GAMLSS fitting engine is included or required.
# Jones and Pewsey (2009), Sinh-arcsinh distributions, Biometrika 96, 761-780.

.svem_log_cosh <- function(x) {
  ax <- abs(x)
  ax + log1p(exp(-2 * ax)) - log(2)
}

.svem_wmt_shasho_log_density <- function(x, mu, sigma, nu, tau) {
  a <- asinh((x - mu) / sigma)
  b <- tau * a - nu
  log(tau) - log(sigma) - log(2 * pi) / 2 - .svem_log_cosh(a) +
    .svem_log_cosh(b) - sinh(b)^2 / 2
}

# Mean negative log likelihood and its analytic gradient in
# (mu, log(sigma), nu, log(tau)); input data are standardized by the fitter.
.svem_wmt_shasho_objective <- function(par, x, gradient = FALSE) {
  sigma <- exp(par[2L])
  tau <- exp(par[4L])
  a <- asinh((x - par[1L]) / sigma)
  b <- tau * a - par[3L]
  if (any(!is.finite(c(par, sigma, tau, a, b))) || sigma <= 0 || tau <= 0 ||
      any(abs(b) > 350)) {
    return(if (gradient) rep(NA_real_, 4L) else Inf)
  }
  if (!gradient) {
    return(-mean(par[4L] - par[2L] - log(2 * pi) / 2 -
                   .svem_log_cosh(a) + .svem_log_cosh(b) - sinh(b)^2 / 2))
  }
  h <- tanh(b) - sinh(b) * cosh(b)
  za <- tanh(a)
  inv_cosh <- exp(-.svem_log_cosh(a))
  -c(mean((-h * tau + za) * inv_cosh / sigma),
     mean(-1 - h * tau * za + za^2),
     mean(-h),
     mean(1 + h * tau * a))
}

.svem_wmt_fit_shasho <- function(x) {
  if (!is.numeric(x) || length(x) < 5L || any(!is.finite(x)) ||
      length(unique(x)) < 5L) {
    stop("The SHASHo fit requires at least five distinct finite distances.",
         call. = FALSE)
  }
  center <- stats::median(x)
  scale <- stats::IQR(x) / 1.349
  if (!is.finite(scale) || scale <= 0) scale <- stats::sd(x)
  if (!is.finite(scale) || scale <= 0) {
    stop("The SHASHo distance scale must be positive and finite.", call. = FALSE)
  }
  z <- (x - center) / scale

  # Deterministic starts cover either skew direction and light/heavy tails.
  # No random draws or constraints on the statistical parameter space.
  starts <- rbind(c(0, 0, 0, 0), c(0, log(0.2), 0.5, log(0.5)),
                  c(0, log(0.2), -0.5, log(0.5)),
                  c(0, log(2), 0, log(2)))
  fits <- lapply(seq_len(nrow(starts)), function(i) {
    tryCatch(stats::optim(
      starts[i, ], .svem_wmt_shasho_objective,
      gr = function(par, x) .svem_wmt_shasho_objective(par, x, TRUE),
      x = z, method = "BFGS", control = list(maxit = 2000L, reltol = 1e-12)
    ), error = function(e) NULL)
  })
  usable <- vapply(fits, function(fit) {
    if (is.null(fit) || fit$convergence != 0L || !is.finite(fit$value) ||
        any(!is.finite(fit$par))) return(FALSE)
    score <- .svem_wmt_shasho_objective(fit$par, z, TRUE)
    all(is.finite(score)) && max(abs(score)) < 1e-4
  }, logical(1L))
  if (!any(usable)) {
    stop("The SHASHo permutation-distance fit did not converge; no WMT p-value was returned. ",
         "Try increasing nPerm to obtain a better-resolved reference sample.",
         call. = FALSE)
  }
  candidates <- fits[usable]
  best <- candidates[[which.min(vapply(candidates, `[[`, numeric(1L), "value"))]]
  parameters <- c(mu = center + scale * best$par[1L],
                  sigma = scale * exp(best$par[2L]), nu = best$par[3L],
                  tau = exp(best$par[4L]))
  if (any(!is.finite(parameters)) || parameters["sigma"] <= 0 ||
      parameters["tau"] <= 0) {
    stop("The SHASHo fit returned invalid parameters; no WMT p-value was returned.",
         call. = FALSE)
  }
  structure(list(
    parameters = parameters, converged = TRUE,
    loglik = -length(x) * (best$value + log(scale)), n = length(x),
    method = "Maximum likelihood (BFGS, deterministic multiple starts)",
    convergence = best$convergence,
    gradient_max = max(abs(.svem_wmt_shasho_objective(best$par, z, TRUE))),
    successful_starts = sum(usable), attempted_starts = nrow(starts)
  ), class = "svem_shasho_fit")
}

#' @export
#' @noRd
coef.svem_shasho_fit <- function(object, what = c("mu", "sigma", "nu", "tau"), ...) {
  what <- match.arg(what)
  value <- object$parameters[[what]]
  if (what %in% c("sigma", "tau")) value <- log(value)
  stats::setNames(value, "(Intercept)")
}
