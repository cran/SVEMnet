if (getRversion() >= "2.15.1") utils::globalVariables(
  c("jloop")
)
`%||%` <- function(a, b) if (!is.null(a)) a else b
