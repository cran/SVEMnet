# The physical budget includes repeated reference recipes. Explicit model and
# exploration counts refer to distinct additional physical preparations.
.svem_validate_allocation_counts <- function(budget, n_optima, n_explore, priority) {
  if (is.null(n_optima) && is.null(n_explore)) return(invisible(NULL))
  if (!is.null(n_optima)) {
    n_optima <- .svem_integer_scalar(n_optima, "n_optima", min = 0L)
    roles <- c("svemnet", "gaupro", "extra_trees", "forward", "spline")
    available <- intersect(priority, roles)
    if (n_optima > length(available)) {
      stop("n_optima exceeds the number of prioritized model roles: requested ", n_optima,
           ", maximum ", length(available), " (", paste(available, collapse = ", "),
           "). Each model nominates at most one recipe, and duplicates can reduce the count. ",
           "For more preparations, use budget = ", n_optima,
           " or request exploration with n_explore.", call. = FALSE)
    }
  }
  if (!is.null(n_explore)) n_explore <- .svem_integer_scalar(n_explore, "n_explore", min = 0L)
  if (is.null(budget) && (is.null(n_optima) || is.null(n_explore))) {
    stop("Supply both n_optima and n_explore, or an explicit total physical budget with one count.", call. = FALSE)
  }
  budget <- .svem_integer_scalar(budget, "budget", min = 0L, allow_null = TRUE)
  if (!is.null(budget) && sum(as.double(n_optima), as.double(n_explore)) > budget) {
    stop("n_optima and n_explore exceed the physical budget before reference repeats.", call. = FALSE)
  }
  invisible(NULL)
}

.svem_count_targets <- function(budget, n_optima, n_explore, required, priority) {
  repeats <- length(required)
  models <- length(setdiff(priority, "incumbent"))
  requested <- sum(as.double(n_optima), as.double(n_explore), repeats)
  if (is.null(budget)) budget <- .svem_integer_scalar(requested, "inferred budget")
  if (requested > budget) {
    stop("The budget cannot accommodate the distinct required references plus the requested model and exploration counts.", call. = FALSE)
  }
  target <- if (!is.null(n_optima)) as.integer(n_optima) else
    as.integer(min(models, budget - repeats - sum(n_explore)))
  list(budget = budget, requested_optima = n_optima, requested_explore = n_explore,
       target_optima = target, distinct_repeats = repeats)
}
