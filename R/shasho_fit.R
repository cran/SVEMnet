# Independent implementation derived from Jones and Pewsey (2009),
# "Sinh-arcsinh distributions", doi:10.1093/biomet/asp053, equations (1)-(2).
# Here nu = -epsilon in that paper: T(x) = sinh(tau*asinh((x-mu)/sigma)-nu).
# The paper's scale/tail reparametrization is used with numerical nlminb;
# no third-party implementation or analytic score code is incorporated.

.svem_wmt_log_cosh <- function(x) {
  ax <- abs(x)
  ax + log1p(exp(-2 * ax)) - log(2)
}

.svem_wmt_shasho_standardize <- function(x, mu, sigma) {
  difference <- x - mu
  z <- difference / sigma
  # A finite x can produce an overflowing difference or ratio even though
  # its arcsinh and final probability are representable.
  overflow <- is.infinite(z) & is.finite(x)
  log_ratio <- rep(NA_real_, length(x))
  if (any(overflow)) {
    d <- difference[overflow]
    log_difference <- log(abs(d))
    wide <- is.infinite(d)
    log_difference[wide] <- log(abs(x[overflow][wide] / 2 - mu / 2)) + log(2)
    log_ratio[overflow] <- log_difference - log(sigma)
    z[overflow] <- sign(d) * exp(log_ratio[overflow])
  }
  az <- abs(z)
  log_hypot <- log(pmax(1, az)) + 0.5 * log1p(pmin(az, 1 / az)^2)
  inverse_sinh <- asinh(z)
  remaining <- overflow & is.infinite(z)
  inverse_sinh[remaining] <- sign(z[remaining]) * (log_ratio[remaining] + log(2))
  log_hypot[remaining] <- log_ratio[remaining]
  list(inverse_sinh = inverse_sinh, log_hypot = log_hypot)
}

.svem_wmt_shasho_log_density <- function(x, mu, sigma, nu, tau) {
  pars <- c(mu, sigma, nu, tau)
  if (any(lengths(list(mu, sigma, nu, tau)) != 1L) ||
      !is.numeric(pars) || is.complex(pars) || any(!is.finite(pars)) ||
      sigma <= 0 || tau <= 0) {
    stop("Sinh-arcsinh parameters must be finite, with positive sigma and tau.",
         call. = FALSE)
  }
  standardized <- .svem_wmt_shasho_standardize(x, mu, sigma)
  a <- tau * standardized$inverse_sinh - nu
  value <- -0.5 * log(2 * pi) - 0.5 * sinh(a)^2 + log(tau) -
    log(sigma) + .svem_wmt_log_cosh(a) - standardized$log_hypot
  # These limits also cover overflow in intermediate tail calculations.
  value[is.infinite(x) | is.infinite(a)] <- -Inf
  value
}

.svem_wmt_shasho_upper_tail <- function(q, mu, sigma, nu, tau) {
  pars <- c(mu, sigma, nu, tau)
  if (any(lengths(list(mu, sigma, nu, tau)) != 1L) ||
      !is.numeric(pars) || is.complex(pars) || any(!is.finite(pars)) ||
      sigma <= 0 || tau <= 0) {
    stop("Sinh-arcsinh parameters must be finite, with positive sigma and tau.",
         call. = FALSE)
  }
  standardized <- .svem_wmt_shasho_standardize(q, mu, sigma)
  stats::pnorm(sinh(tau * standardized$inverse_sinh - nu), lower.tail = FALSE)
}

.svem_wmt_shasho_objective <- function(p, x) {
  if (any(!is.finite(p))) return(1e100)
  # p[2] is log(sigma/tau), as motivated by the paper's reparametrization.
  sigma <- exp(p[2L] + p[4L])
  tau <- exp(p[4L])
  if (!is.finite(sigma) || !is.finite(tau) || sigma <= 0 || tau <= 0) {
    return(1e100)
  }
  standardized <- .svem_wmt_shasho_standardize(x, p[1L], sigma)
  a <- tau * standardized$inverse_sinh - p[3L]
  if (any(!is.finite(a)) || any(abs(a) > 350)) return(1e100)
  terms <- 0.5 * sinh(a)^2 - .svem_wmt_log_cosh(a) + standardized$log_hypot
  value <- mean(terms) + 0.5 * log(2 * pi) + p[2L]
  if (is.finite(value)) value else 1e100
}

.svem_wmt_shasho_numerical_gradient <- function(p, objective) {
  vapply(seq_along(p), function(j) {
    h <- .Machine$double.eps^(1 / 3) * max(1, abs(p[j]))
    hi <- lo <- half_hi <- half_lo <- p
    hi[j] <- hi[j] + h
    lo[j] <- lo[j] - h
    half_hi[j] <- half_hi[j] + h / 2
    half_lo[j] <- half_lo[j] - h / 2
    coarse <- (objective(hi) - objective(lo)) / (2 * h)
    fine <- (objective(half_hi) - objective(half_lo)) / h
    # Richardson extrapolation removes the leading truncation error. This
    # matters when the location likelihood is much narrower than data scale.
    (4 * fine - coarse) / 3
  }, numeric(1))
}

.svem_wmt_fit_shasho <- function(x) {
  if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) ||
      any(!is.finite(x)) || length(unique(x)) < 5L) {
    stop("Sinh-arcsinh fitting needs a finite numeric vector with at least five distinct values.",
         call. = FALSE)
  }
  # Normalize before optimizing. The two-stage calculation avoids squaring
  # large observations, and a median keeps isolated extremes from setting the
  # location used by the numerical optimization.
  amplitude <- max(abs(x))
  unit_x <- x / amplitude
  center <- stats::median(unit_x)
  residual <- unit_x - center
  spread <- stats::IQR(residual) / (stats::qnorm(0.75) - stats::qnorm(0.25))
  if (!is.finite(spread) || spread <= 0) spread <- stats::sd(residual)
  if (!is.finite(spread) || spread <= 0) {
    stop("Sinh-arcsinh fitting could not resolve a positive data scale.", call. = FALSE)
  }
  y <- residual / spread
  objective <- function(p) .svem_wmt_shasho_objective(p, y)

  # Each shape start supplies theoretical quantiles; ordinary least squares
  # aligns their location and scale with empirical quantiles. No random draws,
  # finite parameter boxes, or analytic likelihood derivatives are used.
  probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  observed <- as.numeric(stats::quantile(y, probs, names = FALSE))
  normal <- asinh(stats::qnorm(probs))
  shapes <- expand.grid(nu = c(-1, 0, 1), tau = c(0.6, 1, 1.8))
  candidates <- vector("list", nrow(shapes))
  for (i in seq_len(nrow(shapes))) {
    nu0 <- shapes$nu[i]
    tau0 <- shapes$tau[i]
    theoretical <- sinh((normal + nu0) / tau0)
    scale0 <- sum((theoretical - mean(theoretical)) *
                    (observed - mean(observed))) /
      sum((theoretical - mean(theoretical))^2)
    if (!is.finite(scale0) || scale0 <= 0) scale0 <- 1
    mu0 <- mean(observed) - scale0 * mean(theoretical)
    initial <- c(mu0, log(scale0 / tau0), nu0, log(tau0))
    fit <- tryCatch(stats::nlminb(
      initial, objective,
      control = list(eval.max = 1200L, iter.max = 600L,
                     rel.tol = 1e-11, x.tol = 1e-9)
    ), error = function(e) NULL)
    if (is.null(fit) || !is.finite(fit$objective) || fit$objective >= 1e100 ||
        any(!is.finite(fit$par))) next
    gradient <- .svem_wmt_shasho_numerical_gradient(fit$par, objective)
    if (fit$convergence != 0L || any(!is.finite(gradient)) ||
        max(abs(gradient)) > 1e-4) {
      # Retry the same likelihood from this candidate with more accurate
      # numerical derivatives. No distribution or convergence rule changes.
      numerical_gradient <- function(p) .svem_wmt_shasho_numerical_gradient(p, objective)
      fit <- tryCatch(stats::nlminb(
        fit$par, objective, gradient = numerical_gradient,
        control = list(eval.max = 1200L, iter.max = 600L,
                       rel.tol = 1e-11, x.tol = 1e-9)
      ), error = function(e) NULL)
      if (is.null(fit) || !is.finite(fit$objective) || fit$objective >= 1e100 ||
          any(!is.finite(fit$par))) next
      gradient <- numerical_gradient(fit$par)
    }
    if (fit$convergence != 0L || any(!is.finite(gradient)) ||
        max(abs(gradient)) > 1e-4) next
    fit$gradient_max <- max(abs(gradient))
    candidates[[i]] <- fit
  }
  successful <- which(!vapply(candidates, is.null, logical(1)))
  if (!length(successful)) {
    stop("Sinh-arcsinh maximum likelihood optimization did not converge to a finite stationary fit.",
         call. = FALSE)
  }
  best_index <- successful[which.min(vapply(candidates[successful],
                                          function(z) z$objective, numeric(1)))]
  best <- candidates[[best_index]]
  p <- best$par
  log_data_scale <- log(amplitude) + log(spread)
  parameters <- c(mu = (center + spread * p[1L]) * amplitude,
                  sigma = exp(log_data_scale + p[2L] + p[4L]),
                  nu = p[3L], tau = exp(p[4L]))
  if (any(!is.finite(parameters)) || parameters["sigma"] <= 0 ||
      parameters["tau"] <= 0) {
    stop("Sinh-arcsinh optimization produced unrepresentable parameters.", call. = FALSE)
  }
  loglik <- -length(x) * (best$objective + log_data_scale)
  if (!is.finite(loglik)) {
    stop("Sinh-arcsinh optimization produced a nonfinite log likelihood.", call. = FALSE)
  }
  structure(list(parameters = parameters, converged = TRUE,
                 loglik = loglik, n = length(x),
                 method = "Independent numerical maximum likelihood (nlminb)",
                 convergence = best$convergence,
                 gradient_max = best$gradient_max,
                 successful_starts = length(successful),
                 attempted_starts = nrow(shapes)),
            class = "svem_shasho_fit")
}

#' @export
#' @noRd
coef.svem_shasho_fit <- function(object, what = c("mu", "sigma", "nu", "tau"), ...) {
  what <- match.arg(what)
  value <- object$parameters[[what]]
  if (what %in% c("sigma", "tau")) value <- log(value)
  stats::setNames(value, "(Intercept)")
}
