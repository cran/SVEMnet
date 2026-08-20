# Internal scalar validators shared by public entry points.  These reject
# fractional count values instead of silently truncating them with as.integer().

.svem_integer_scalar <- function(x, name, min = 0L,
                                 max = .Machine$integer.max,
                                 allow_null = FALSE) {
  if (allow_null && is.null(x)) return(NULL)
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != floor(x) || x < min || x > max) {
    stop("`", name, "` must be a single integer in [", min, ", ", max, "].",
         call. = FALSE)
  }
  as.integer(x)
}

.svem_numeric_scalar <- function(x, name, lower = -Inf, upper = Inf,
                                 lower_open = FALSE, upper_open = FALSE,
                                 allow_null = FALSE) {
  if (allow_null && is.null(x)) return(NULL)
  valid <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
  if (valid) {
    valid <- if (lower_open) x > lower else x >= lower
  }
  if (valid) {
    valid <- if (upper_open) x < upper else x <= upper
  }
  if (!valid) {
    left  <- if (lower_open) "(" else "["
    right <- if (upper_open) ")" else "]"
    stop("`", name, "` must be a single finite number in ", left,
         lower, ", ", upper, right, ".", call. = FALSE)
  }
  as.numeric(x)
}

.svem_logical_scalar <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
  }
  x
}
