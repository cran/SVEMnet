if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "jloop",
    # plot.svem_model
    "Actual", "Predicted", "Type",
    # plot.svem_binomial
    "p", "y", "p_mean", "y_mean", "n",
    "fpr", "tpr", "recall", "precision",
    # plot.svem_significance_test
    "Group", "D", "Source_Type", "Response",
    # svem_nonzero
    "Order", "PercentNonzero", "Variable"
  ))
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
