if (getRversion() >= "2.15.1") utils::globalVariables(
  c("Order", "Percent of Bootstraps Nonzero", "Variable", "Group", "D", "Source_Type","jloop",".data", "Actual", "Predicted", "Type")
)
`%||%` <- function(a, b) if (!is.null(a)) a else b
