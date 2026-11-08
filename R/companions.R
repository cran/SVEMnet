# Parameters are canonicalized once, so compatibility checks across response
# fits compare physical domains rather than list names or integer storage.
.svem_mixture_groups <- function(groups) {
  if (is.null(groups) || identical(groups, list())) return(NULL)
  if (!is.list(groups) || !length(groups))
    stop("mixture_groups must be a nonempty list of mixture groups.", call. = FALSE)
  groups <- unname(lapply(groups, function(g) {
    if (!is.list(g)) stop("Each mixture group must be a list.", call. = FALSE)
    vars <- g$vars
    if (!is.character(vars) || !length(vars) || anyNA(vars) || any(!nzchar(vars)))
      stop("Each mixture group must name its components in vars.", call. = FALSE)
    k <- length(vars)
    lower <- if (is.null(g$lower)) rep(0, k) else g$lower
    upper <- if (is.null(g$upper)) rep(1, k) else g$upper
    total <- if (is.null(g$total)) 1 else g$total
    if (!is.numeric(lower) || !is.numeric(upper) || !is.numeric(total) ||
        length(lower) != k || length(upper) != k || length(total) != 1L ||
        any(!is.finite(c(lower, upper, total))) || any(lower > upper) ||
        total < sum(lower) - 1e-12 || total > sum(upper) + 1e-12)
      stop("Invalid or infeasible lower, upper or total in mixture_groups.", call. = FALSE)
    list(vars = unname(vars), lower = as.numeric(lower), upper = as.numeric(upper),
         total = as.numeric(total))
  }))
  if (anyDuplicated(unlist(lapply(groups, `[[`, "vars"), use.names = FALSE)))
    stop("Mixture components cannot repeat within or across groups.", call. = FALSE)
  groups
}

.svem_check_mixture_totals <- function(data, groups) {
  for (g in groups) {
    totals <- rowSums(data[g$vars])
    if (any(abs(totals - g$total) > 1e-8 * max(1, abs(g$total)), na.rm = TRUE))
      stop("Mixture components ", paste(g$vars, collapse = ", "),
           " must sum to ", g$total, ".", call. = FALSE)
  }
}

.svem_companion_names <- function(companions) {
  available <- c("gaupro", "extra_trees")
  if (identical(companions, FALSE)) return(character())
  if (identical(companions, TRUE)) return(available)
  if (!is.character(companions) || !length(companions) ||
      anyNA(companions) || any(!companions %in% available))
    stop("companions must be FALSE, TRUE, or names from 'gaupro' and 'extra_trees'.",
         call. = FALSE)
  available[available %in% companions]
}

# The raw-coordinate interpretation is clear for polynomial and interaction
# formulas. Decline arbitrary transformations rather than silently fitting a
# companion on a different experimental factor or response scale.
.svem_companion_rhs <- function(expr) {
  if (!is.call(expr)) return(TRUE)
  if (!is.symbol(expr[[1L]]) ||
      !as.character(expr[[1L]]) %in% c("~", "+", "-", "*", "/", "^", ":", "I", "("))
    return(FALSE)
  all(vapply(as.list(expr)[-1L], .svem_companion_rhs, logical(1)))
}

.svem_companion_map <- function(object, data, blocks, using_spec, continuous_blocks = character()) {
  if (object$family != "gaussian" || !is.symbol(object$formula[[2L]]))
    stop("Companions require an untransformed Gaussian response.", call. = FALSE)
  spec <- object$sampling_schema
  vars <- spec$predictor_vars
  rhs <- stats::delete.response(object$terms)
  if (!using_spec && (!.svem_companion_rhs(object$formula[[3L]]) ||
                       length(setdiff(all.vars(rhs), names(data)))))
    stop("Supply original data columns for the companion predictors.", call. = FALSE)
  if (any(spec$var_classes == "ordered", na.rm = TRUE))
    stop("Ordered predictors require an explicit numeric or nominal representation.",
         call. = FALSE)
  if (any(!vapply(data[vars], function(x) {
    !is.matrix(x) && (is.factor(x) || is.character(x) || is.logical(x) ||
                        (is.numeric(x) && !is.object(x)))
  }, logical(1))))
    stop("Companion predictors must be numeric or categorical data columns.", call. = FALSE)
  varying <- vapply(data[vars], function(x) length(unique(x)) > 1L, logical(1))
  if (blocks == "respect" && any(spec$blocking %in% vars[varying]))
    stop("A varying declared block requires companion_blocks = 'interacting'.",
         call. = FALSE)

  if (length(continuous_blocks) &&
      (blocks != "interacting" || !all(continuous_blocks %in% spec$blocking) ||
       !all(vapply(data[continuous_blocks], is.numeric, logical(1))) ||
       length(intersect(continuous_blocks, spec$discrete_numeric)))) {
    stop("companion_continuous_blocks must name continuous numeric blocks declared in an interacting fit.",
         call. = FALSE)
  }
  nominal_blocks <- setdiff(spec$blocking, continuous_blocks)
  nominal <- vars[vapply(data[vars], function(x) !is.numeric(x), logical(1)) |
                    vars %in% nominal_blocks]
  unique_ids <- nominal_blocks[vapply(data[nominal_blocks], function(x) {
    is.numeric(x) && length(x) > 1L && !anyDuplicated(x)
  }, logical(1))]
  if (length(unique_ids)) {
    warning("Treating unreplicated numeric block(s) as nominal IDs: ",
            paste(unique_ids, collapse = ", "),
            ". Name measured continuous covariates in companion_continuous_blocks.",
            call. = FALSE)
  }
  levels <- lapply(data[nominal], function(x) {
    if (is.factor(x)) return(base::levels(droplevels(x)))
    if (is.numeric(x)) {
      ids <- unique(x)
      # Preserve the historical ordering of nominal labels and distinguish
      # exact numeric values when their printable labels coincide.
      return(ids[order(as.character(ids), ids)])
    }
    sort(unique(as.character(x)))
  })
  numeric_vars <- setdiff(vars, nominal)
  bounds <- lapply(data[numeric_vars], range)
  for (v in intersect(numeric_vars, colnames(spec$num_ranges))) {
    bounds[[v]] <- as.numeric(spec$num_ranges[, v])
  }
  omitted <- character()
  for (g in spec$mixture_groups) {
    for (i in seq_along(g$vars)) bounds[[g$vars[i]]] <- c(g$lower[i], g$upper[i])
    active <- g$vars[g$vars %in% vars[varying]]
    if (length(active)) omitted <- c(omitted, active[length(active)])
  }
  active <- setdiff(vars[varying], omitted)
  if (!length(active) || length(unique(object$actual_y)) < 2L)
    stop("The response and at least one independent coordinate must vary.", call. = FALSE)
  if (any(vapply(bounds[intersect(active, numeric_vars)], diff, numeric(1)) <= 0))
    stop("A varying numeric predictor has zero range in the declared domain.", call. = FALSE)
  list(vars = vars, active = active, levels = levels, bounds = bounds,
       constants = lapply(data[vars[!varying]], `[[`, 1L),
       omitted = omitted, groups = spec$mixture_groups, continuous_blocks = continuous_blocks)
}

.svem_companion_inputs <- function(map, data, unseen) {
  .svem_check_mixture_totals(data, map$groups)
  x <- data[map$vars]
  invalid <- rep(FALSE, nrow(data))
  for (v in map$vars) {
    value <- x[[v]]
    if (v %in% names(map$levels)) {
      support <- map$levels[[v]]
      missing <- is.na(value)
      if (is.numeric(support)) {
        if (!is.numeric(value)) stop("Expected numeric companion block '", v, "'.", call. = FALSE)
        # Match original doubles before assigning nominal labels. Distinct
        # observed IDs can have the same default as.character() representation.
        value <- factor(match(value, support), levels = seq_along(support))
      } else {
        value <- factor(as.character(value), levels = support)
      }
      invalid <- invalid | (!missing & is.na(value))
    } else {
      if (!is.numeric(value))
        stop("Expected numeric companion predictor '", v, "'.", call. = FALSE)
      if (v %in% names(map$constants)) {
        invalid <- invalid | (!is.na(value) & value != map$constants[[v]])
      } else if (v %in% map$active) {
        bound <- map$bounds[[v]]
        value <- (as.numeric(value) - bound[1L]) / diff(bound)
      }
      value[!is.finite(value)] <- NA_real_
    }
    x[[v]] <- value
  }
  if (any(invalid) && unseen == "error")
    stop("Unsupported companion predictor level or constant in newdata.", call. = FALSE)
  valid <- stats::complete.cases(x) & !invalid
  x <- x[valid, map$active, drop = FALSE]
  # Engine-safe identifiers also handle a user factor named 'response' or
  # containing operators. The reversible order is in map$active.
  names(x) <- paste0("x", seq_along(map$active))
  list(x = x, valid = valid)
}

# The legacy sampler can set a continuous numeric block to its midpoint. In
# an interacting companion the same column is a nominal ID, so its worklist
# requires a level actually observed for each participating response.
.svem_worklist_numeric_blocks <- function(objects, x, supplied) {
  numeric_blocks <- intersect(objects[[1L]]$sampling_schema$blocking,
                              names(x)[vapply(x, is.numeric, logical(1))])
  if (!length(numeric_blocks)) return(x)
  maps <- lapply(objects, function(obj) attr(obj$companions, "input_map", exact = TRUE))
  for (v in numeric_blocks) {
    supports <- Filter(Negate(is.null), lapply(seq_along(maps), function(i) {
      support <- maps[[i]]$levels[[v]]
      if (is.character(support)) {
        # Earlier maps stored rounded labels. Their original numeric IDs
        # are available among the retained preparation settings.
        original <- attr(objects[[i]]$companions, "prepared_recipes", exact = TRUE)[[v]]
        if (is.numeric(original)) return(unique(original))
        return(as.numeric(support))
      }
      support
    }))
    if (!length(supports)) next
    levels <- sort(as.numeric(Reduce(intersect, supports)))
    if (!length(levels)) stop("No common observed block level for '", v, "'.", call. = FALSE)
    if (any(!x[[v]] %in% levels)) {
      if (supplied) stop("Supply a common observed block level for '", v, "'.", call. = FALSE)
      x[[v]][] <- levels[which.min(abs(levels - x[[v]][1L]))]
    }
  }
  x
}

.svem_fit_companions <- function(object, data, engines, blocks, unseen, using_spec,
                                 continuous_blocks = character()) {
  map <- tryCatch(.svem_companion_map(object, data, blocks, using_spec, continuous_blocks),
                  error = identity)
  # Isolate the entire optional branch, including first-time namespace loading.
  # Draw the seeds for both engines in a fixed order so subset requests agree.
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  old_kind <- RNGkind()
  on.exit({
    if (!identical(RNGkind(), old_kind)) do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  })
  # A common caller seed is often reused when fitting several responses.
  # Key the engine streams by the response as well as by their fixed slot.
  # Double arithmetic here is exact for these bounded integer products.
  modulus <- as.double(.Machine$integer.max) - 1
  label_key <- 0
  for (code in utf8ToInt(enc2utf8(paste(deparse(object$formula[[2L]]), collapse = ""))))
    label_key <- (127 * label_key + code) %% modulus
  seeds <- as.integer(1 + ((as.double(sample.int(modulus, 2L)) - 1 + label_key) %% modulus))
  names(seeds) <- c("gaupro", "extra_trees")
  result <- list()
  if (!inherits(map, "error")) input <- .svem_companion_inputs(map, data, "error")$x
  for (engine in engines) {
    pkg <- if (engine == "gaupro") "GauPro" else "ranger"
    entry <- list(status = "ok", message = "", package = pkg, seed = seeds[[engine]])
    if (inherits(map, "error")) {
      entry$status <- "unsupported"
      entry$message <- conditionMessage(map)
    } else if (!requireNamespace(pkg, quietly = TRUE)) {
      entry$status <- "unavailable"
      entry$message <- sprintf("Install the optional dependency with install.packages('%s').", pkg)
    } else if (engine == "gaupro" && ncol(input) == 1L && is.factor(input[[1L]])) {
      entry$status <- "unsupported"
      entry$message <- "GauPro cannot reliably predict with a single factor as its only input; use another companion."
    } else {
      entry$version <- as.character(utils::packageVersion(pkg))
      notes <- character()
      set.seed(entry$seed)
      fit <- tryCatch(withCallingHandlers({
        if (engine == "gaupro") {
          train <- input
          train$response <- object$actual_y
          GauPro::gpkm(formula = stats::reformulate(names(input), "response", env = baseenv()),
                       data = train, kernel = "matern52", nug.est = TRUE, restarts = 0,
                       normalize = TRUE, optimizer = "L-BFGS-B",
                       parallel = FALSE, verbose = 0)
        } else {
          ranger::ranger(x = input, y = object$actual_y, num.trees = 128L,
                         splitrule = "extratrees", mtry = ncol(input),
                         num.random.splits = 1L, min.node.size = 1L, min.bucket = 2L,
                         replace = FALSE, sample.fraction = 1,
                         respect.unordered.factors = "partition", num.threads = 1L,
                         oob.error = FALSE, write.forest = TRUE, seed = entry$seed)
        }
      }, warning = function(w) {
        notes <<- c(notes, trimws(conditionMessage(w)))
      }, message = function(m) {
        notes <<- c(notes, trimws(conditionMessage(m)))
        invokeRestart("muffleMessage")
      }), error = identity)
      if (!inherits(fit, "error") && engine == "gaupro") {
        probe <- tryCatch(withCallingHandlers(as.numeric(fit$predict(input)),
          warning = function(w) notes <<- c(notes, trimws(conditionMessage(w))),
          message = function(m) {
            notes <<- c(notes, trimws(conditionMessage(m)))
            invokeRestart("muffleMessage")
          }), error = identity)
        if (inherits(probe, "error")) fit <- probe
        else if (length(probe) != nrow(input) || any(!is.finite(probe)))
          fit <- simpleError("GauPro returned invalid training predictions.")
      }
      entry$warnings <- unique(notes)
      if (inherits(fit, "error")) {
        entry$status <- "failed"
        entry$message <- conditionMessage(fit)
      } else {
        entry$model <- fit
        if (engine == "gaupro") {
          entry$nugget <- as.numeric(fit$nug)
          entry$settings <- list(kernel = "matern52", fitted_kernel_class = class(fit$kernel)[1L],
                                 restarts = fit$restarts, optimizer = fit$optimizer,
                                 nug.est = fit$nug.est, normalize = fit$normalize)
          entry$deviance <- tryCatch(as.numeric(fit$deviance()), error = function(e) {
            warning("GauPro's fitted deviance was unavailable: ", conditionMessage(e), call. = FALSE)
            NA_real_
          })
        } else {
          entry$settings <- list(splitrule = fit$splitrule, num.trees = fit$num.trees,
                                 mtry = fit$mtry, min.node.size = fit$min.node.size,
                                 min.bucket = 2L, num.random.splits = 1L, replace = FALSE,
                                 sample.fraction = 1, respect.unordered.factors = "partition")
        }
      }
    }
    if (entry$status != "ok") warning(engine, " companion ", entry$status, ": ",
                                      entry$message, call. = FALSE)
    result[[engine]] <- entry
  }
  prepared <- unique(data[object$sampling_schema$predictor_vars])
  row.names(prepared) <- NULL
  structure(result, class = c("svem_companions", "list"),
            input_map = if (!inherits(map, "error")) map,
            prepared_recipes = prepared,
            unseen = unseen, block_mode = blocks)
}

#' @export
print.svem_companions <- function(x, ...) {
  for (engine in names(x)) {
    entry <- x[[engine]]
    cat(engine, ": ", entry$status, if (nzchar(entry$message)) paste0(" - ", entry$message),
        "\n", sep = "")
    if (length(entry$warnings)) cat("  ", length(entry$warnings),
      " fit note(s); inspect $", engine, "$warnings.\n", sep = "")
  }
  invisible(x)
}

.svem_predict_engine <- function(object, newdata, engine, type, debias, se.fit, interval) {
  if (type != "response" || debias || interval || (se.fit && engine != "gaupro"))
    stop("Companions provide ordinary response means; only gaupro supports latent se.fit.",
         call. = FALSE)
  entries <- object$companions
  if (engine == "all") {
    out <- data.frame(svemnet = as.numeric(stats::predict(object, newdata)),
                      row.names = row.names(newdata))
    status <- c(svemnet = "ok")
    for (nm in names(entries)) {
      entry <- entries[[nm]]
      if (entry$status == "ok") {
        pr <- tryCatch(.svem_predict_engine(object, newdata, nm, type, FALSE, FALSE, FALSE),
                       error = identity)
        if (!inherits(pr, "error")) {
          out[[nm]] <- pr
          status[nm] <- "ok"
          next
        }
        entry$status <- "prediction_failed"
        entry$message <- conditionMessage(pr)
        warning(nm, " companion prediction failed: ", entry$message, call. = FALSE)
      }
      out[[nm]] <- rep(NA_real_, nrow(newdata))
      status[nm] <- paste(entry$status, entry$message)
    }
    attr(out, "engine_status") <- status
    return(out)
  }
  entry <- entries[[engine]]
  if (is.null(entry)) stop("Companion '", engine, "' was not requested at fit time.",
                           call. = FALSE)
  if (entry$status != "ok") stop(engine, " companion ", entry$status, ": ", entry$message,
                                 call. = FALSE)
  if (!requireNamespace(entry$package, quietly = TRUE))
    stop("Install the companion dependency with install.packages('", entry$package, "').",
         call. = FALSE)
  current <- as.character(getNamespaceVersion(entry$package))
  if (!is.null(entry$version) && !identical(entry$version, current)) {
    key <- paste(entry$package, entry$version, current, sep = "/")
    if (!exists(key, envir = .svem_companion_version_notices, inherits = FALSE)) {
      assign(key, TRUE, envir = .svem_companion_version_notices)
      warning(entry$package, " companion fitted with version ", entry$version,
              " is being used with version ", current, ". Verify predictions after reloading.",
              call. = FALSE)
    }
  }
  newdata <- .svem_validate_newdata_classes(object, newdata)
  encoded <- .svem_companion_inputs(attr(entries, "input_map"), newdata, attr(entries, "unseen"))
  if (any(!encoded$valid))
    warning(sum(!encoded$valid), " row(s) contain missing or unsupported companion inputs; returning NA.",
            call. = FALSE)
  fit <- sd_mean <- rep(NA_real_, nrow(newdata))
  if (any(encoded$valid)) {
    x <- encoded$x
    if (engine == "extra_trees") {
      fit[encoded$valid] <- stats::predict(entry$model, data = x, num.threads = 1L)$predictions
    } else if (se.fit) {
      pr <- entry$model$predict(x, se.fit = TRUE, mean_dist = TRUE)
      fit[encoded$valid] <- as.numeric(pr$mean)
      sd_mean[encoded$valid] <- as.numeric(pr$se)
    } else {
      fit[encoded$valid] <- as.numeric(entry$model$predict(x))
    }
  }
  if (se.fit) list(fit = fit, se.fit = sd_mean) else fit
}

.svem_companion_version_notices <- new.env(parent = emptyenv())
