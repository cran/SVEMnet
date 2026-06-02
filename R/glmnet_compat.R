.glmnet_direct_control_args <- function() {
  c("thresh", "maxit", "dfmax", "pmax", "trace.it")
}

.glmnet_supports_control <- function(fun) {
  "control" %in% names(formals(fun))
}

.glmnet_check_control_list <- function(control) {
  if (is.null(control)) return(control)
  if (!is.list(control)) {
    stop("'control' must be a named list of glmnet algorithm-control values.")
  }
  nms <- names(control)
  if (is.null(nms) || any(!nzchar(nms))) {
    stop("'control' must be a named list of glmnet algorithm-control values.")
  }
  control
}

.glmnet_prepare_call_args <- function(args, fun, default_control = list(),
                                      direct_control = .glmnet_direct_control_args(),
                                      old_control = .glmnet_direct_control_args()) {
  if (!length(args)) args <- list()

  arg_names <- names(args)
  if (is.null(arg_names)) arg_names <- rep("", length(args))

  user_control <- NULL
  if ("control" %in% arg_names) {
    user_control <- .glmnet_check_control_list(args[["control"]])
    args[["control"]] <- NULL
    arg_names <- names(args)
    if (is.null(arg_names)) arg_names <- rep("", length(args))
  }

  direct_names <- intersect(direct_control, arg_names)
  direct_values <- args[direct_names]
  args <- args[setdiff(seq_along(args), match(direct_names, arg_names))]

  if (.glmnet_supports_control(fun)) {
    control <- default_control
    if (!is.null(user_control)) {
      control <- utils::modifyList(control, user_control, keep.null = TRUE)
    }
    if (length(direct_values)) {
      control <- utils::modifyList(control, direct_values, keep.null = TRUE)
    }
    if (length(control)) args$control <- control
    return(args)
  }

  translated <- list()
  if (!is.null(user_control)) {
    extra <- setdiff(names(user_control), old_control)
    if (length(extra)) {
      stop(
        "The installed glmnet does not support control = list(...). ",
        "Only these control entries can be translated for older glmnet: ",
        paste(old_control, collapse = ", "), ". Unsupported entr",
        if (length(extra) == 1L) "y: " else "ies: ",
        paste(extra, collapse = ", "), "."
      )
    }
    translated <- user_control[intersect(names(user_control), old_control)]
  }

  old_values <- default_control
  if (length(translated)) {
    old_values <- utils::modifyList(old_values, translated, keep.null = TRUE)
  }
  if (length(direct_values)) {
    old_values <- utils::modifyList(old_values, direct_values, keep.null = TRUE)
  }

  if (length(old_values)) {
    for (nm in names(old_values)) {
      if (!nm %in% names(args)) args[[nm]] <- old_values[[nm]]
    }
  }

  args
}
