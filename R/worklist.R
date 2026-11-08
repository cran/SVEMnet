# Training-time mixture geometry fixes the coordinate map used by companions.
.svem_pool_groups <- function(objects, supplied) {
  fitted <- lapply(objects, function(obj) obj$sampling_schema$mixture_groups)
  declared <- Filter(Negate(is.null), fitted)
  if (!length(declared) && is.null(supplied)) return(NULL)
  has_companions <- vapply(objects, function(obj) length(obj$companions) > 0L, logical(1))
  if (any(has_companions & vapply(fitted, is.null, logical(1)))) {
    stop("Declare mixture_groups when fitting companions, then score using the same groups.")
  }
  if (!length(declared)) return(supplied)  # legacy scoring-only constraints
  groups <- declared[[1L]]
  if (!is.null(supplied)) declared <- c(declared, list(.svem_mixture_groups(supplied)))
  if (!all(vapply(declared, identical, logical(1), groups))) {
    stop("mixture_groups must agree with the geometry in every fitted model.")
  }
  groups
}

.svem_candidate_frame <- function(candidates, schema, discrete, groups, tol) {
  vars <- schema$predictor_vars
  if (!is.data.frame(candidates) || !nrow(candidates) || anyDuplicated(names(candidates)) ||
      !all(vars %in% names(candidates))) {
    stop("candidates must be a nonempty data frame containing each predictor by name.")
  }
  x <- candidates[vars]
  if (anyNA(x)) stop("Missing settings in the supplied candidate recipes.")
  groups <- .svem_mixture_groups(groups)
  mix_vars <- unlist(lapply(groups, "[[", "vars"), use.names = FALSE)
  if (!all(mix_vars %in% colnames(schema$num_ranges)) ||
      length(intersect(mix_vars, c(schema$blocking, names(discrete))))) {
    stop("Mixture components must be continuous numeric predictors outside the blocking set.")
  }
  for (v in vars) {
    if (v %in% colnames(schema$num_ranges)) {
      z <- x[[v]]
      if (!is.numeric(z) || any(!is.finite(z))) {
        stop("Candidate predictor '", v, "' must be finite numeric.")
      }
      range <- schema$num_ranges[, v]
      if (!(v %in% mix_vars) && any(z < range["min"] - tol | z > range["max"] + tol)) {
        stop("Candidate predictor '", v, "' is outside its fitted numeric range.")
      }
      if (length(discrete[[v]]) && any(!z %in% discrete[[v]])) {
        stop("Candidate predictor '", v, "' is outside its discrete support.")
      }
    } else {
      labels <- as.character(x[[v]])
      lev <- schema$factor_levels[[v]]
      if (any(!labels %in% lev)) stop("Unknown candidate level in '", v, "'.")
      cls <- unname(schema$var_classes[v])
      x[[v]] <- switch(cls, logical = as.logical(labels), character = labels,
                       factor(labels, levels = lev, ordered = identical(cls, "ordered")))
    }
  }
  .svem_check_candidate_mixtures(x, groups, tol)
  x
}

.svem_check_candidate_mixtures <- function(data, groups, tol) {
  for (g in groups) {
    if (!all(g$vars %in% names(data)) || !all(vapply(data[g$vars], is.numeric, logical(1)))) {
      stop("Candidates must contain numeric columns for every declared mixture component.")
    }
    part <- as.matrix(data[g$vars])
    if (any(!is.finite(part))) stop("Candidate mixture settings must be finite.")
    if (any(abs(rowSums(part) - g$total) > tol * max(1, abs(g$total))) ||
        any(sweep(part, 2L, g$lower, "-") < -tol) ||
        any(sweep(part, 2L, g$upper, "-") > tol)) {
      stop("Candidate recipe violates the declared mixture bounds or total.")
    }
  }
  invisible(NULL)
}

# Length prefixes distinguish labels containing separators. Recipe labels and
# numeric preparation units are independent of factor codes and row order.
.svem_recipe_keys <- function(x, cols, increments = NULL) {
  parts <- lapply(sort(cols, method = "radix"), function(nm) {
    v <- x[[nm]]
    off_grid <- is.na(v)
    if (nm %in% names(increments)) {
      units <- v / increments[[nm]]
      off_grid <- off_grid | abs(units - round(units)) > 1e-7
      v <- round(units)
    }
    value <- if (is.numeric(v)) sprintf("%.15g", ifelse(v == 0, 0, v)) else enc2utf8(as.character(v))
    value[off_grid] <- NA_character_
    value
  })
  if (!length(parts)) return(rep("", nrow(x)))
  keys <- do.call(paste0, lapply(parts, function(z) paste0(nchar(z, type = "bytes"), ":", z)))
  keys[Reduce("|", lapply(parts, is.na))] <- NA_character_
  keys
}

.svem_pool_index <- function(x, candidates, increments, blocks) {
  if (length(increments) && (is.null(candidates) || !is.numeric(increments) ||
      is.null(names(increments)) || anyDuplicated(names(increments)) ||
      !all(names(increments) %in% names(x)[vapply(x, is.numeric, logical(1))]) ||
      any(!is.finite(increments) | increments <= 0))) {
    stop("increments must be named positive numeric preparation steps for supplied candidates.")
  }
  if (length(blocks) && any(vapply(x[blocks], function(z) length(unique(z)) != 1L, logical(1)))) {
    stop("Supply a candidate pool at one known level of each blocking factor per worklist.")
  }
  keys <- .svem_recipe_keys(x, setdiff(names(x), blocks), increments)
  if (anyNA(keys)) stop("Candidate settings must lie on the declared preparation increments.")
  supplied_ids <- !is.null(candidates) && "recipe_id" %in% names(candidates)
  ids <- if (supplied_ids) candidates$recipe_id else match(keys, sort(unique(keys), method = "radix"))
  if (supplied_ids && (!(is.character(ids) || is.numeric(ids)) || anyNA(ids) ||
                       any(!nzchar(as.character(ids))) || anyDuplicated(ids))) {
    stop("Supplied recipe_id values must be unique, nonmissing numeric or character labels.")
  }
  # For duplicate physical recipes keep the smallest ID regardless of input order.
  ord <- order(keys, ids, method = "radix")
  keep <- sort(ord[!duplicated(keys[ord])])
  list(keep = keep, keys = keys[keep], ids = ids[keep])
}

.svem_worklist_names <- function(responses, engines) {
  c("recipe_id", "slot", "role", "nominated_by", "nearest_distance",
    paste0(responses, "_gp_latent_sd"),
    unlist(lapply(engines, function(e) c(paste0(e, "_score"),
                                          paste0(responses, "_", e, "_pred"))),
           use.names = FALSE))
}

.svem_engine_scores <- function(tab, x, objects, engines, parameters, desirability,
                                weights, combine, floor, gp_endpoint) {
  status <- list(data.frame(engine = "svemnet", response = names(objects),
                            status = "ok", message = ""))
  for (engine in engines) {
    scores <- list()
    incomplete <- character()
    for (r in names(objects)) {
      obj <- objects[[r]]
      use_sd <- identical(engine, "gaupro") && identical(r, gp_endpoint)
      mu <- sd <- rep(NA_real_, nrow(x))
      info <- obj$companions[[engine]]
      state <- if (is.null(info)) "not_requested" else info$status
      msg <- if (is.null(info$message)) "" else info$message
      if (identical(state, "ok")) {
        pr <- tryCatch(predict(obj, x, engine = engine, se.fit = use_sd),
                       error = identity)
        if (inherits(pr, "error")) {
          state <- "prediction_failed"
          msg <- conditionMessage(pr)
          if (weights[[r]] == 0) warning(engine, " could not score ", r, ": ", msg, call. = FALSE)
        } else {
          mu <- if (use_sd) pr$fit else pr
          if (use_sd) sd <- pr$se.fit
          if (!any(is.finite(mu))) {
            state <- "no_finite_predictions"
            msg <- "No finite fitted means on the supplied pool."
          }
        }
      }
      status[[length(status) + 1L]] <- data.frame(engine = engine, response = r,
                                                  status = state, message = msg)
      tab[[paste0(r, "_", engine, "_pred")]] <- mu
      if (use_sd) tab[[paste0(r, "_gp_latent_sd")]] <- sd
      if (weights[[r]] > 0) {
        d <- desirability(mu, parameters[[r]])
        # A saturated desirability can turn an infinite prediction into a
        # plausible finite score. Only finite fitted means may rank.
        d[!is.finite(mu)] <- NA_real_
        if (any(!is.finite(d))) {
          reason <- if (state == "ok") "partial coverage" else state
          incomplete <- c(incomplete, paste0(r, " (", reason, ")"))
        }
        if (combine == "geom") d <- log(floor + (1 - floor) * d)
        scores[[r]] <- weights[[r]] * d
      }
    }
    combined <- Reduce("+", scores)
    if (combine == "geom") combined <- exp(combined)
    tab[[paste0(engine, "_score")]] <- combined
    if (!any(is.finite(combined))) {
      warning(engine, " could not score the joint goal on any candidate. Required responses: ",
              paste(incomplete, collapse = ", "), ". No optimum nomination is available. ",
              "An explicit budget may use exploration instead. See engine_status for details.",
              call. = FALSE)
    }
  }
  list(table = tab, status = do.call(rbind, status))
}

.svem_top_ties <- function(values, eligible) {
  rows <- eligible[is.finite(values[eligible])]
  if (!length(rows)) return(integer())
  top <- max(values[rows])
  rows[values[rows] >= top - 1e-12 * max(1, abs(top))]
}

.svem_stable_top <- function(values, eligible, keys) {
  ties <- .svem_top_ties(values, eligible)
  if (!length(ties)) return(NA_integer_)
  ties[order(keys[ties], method = "radix")[1L]]
}

# Gower-style distance at fixed design/pool ranges, excluding nuisance blocks.
.svem_nearest_recipe <- function(point, reference, spans) {
  if (is.null(reference) || !nrow(reference) || !ncol(point)) return(NA_real_)
  distances <- numeric(nrow(reference))
  for (v in names(point)) {
    z <- reference[[v]]
    if (is.numeric(point[[v]])) {
      span <- spans[[v]]
      if (is.finite(span) && span > 0) distances <- distances + abs(z - point[[v]]) / span
    } else {
      distances <- distances + (as.character(z) != as.character(point[[v]]))
    }
  }
  valid <- is.finite(distances)
  if (any(valid)) min(distances[valid]) / ncol(point) else NA_real_
}

.svem_worklist <- function(tab, history, keys, increments, schema, engines, goals,
                           weights, uncertainty_weights, anchors, gp_endpoint, budget,
                           allocation = c("legacy", "portfolio"),
                           priority = c("incumbent", "svemnet", "gaupro", "extra_trees", "forward"),
                           exploration_fraction = .1, incumbent = NULL, benchmark = NULL,
                           n_optima = NULL, n_explore = NULL) {
  allocation <- match.arg(allocation)
  portfolio <- identical(allocation, "portfolio")
  count_mode <- !is.null(n_optima) || !is.null(n_explore)
  if (count_mode && !portfolio) stop("Explicit allocation counts require portfolio allocation.")
  if (portfolio) {
    allowed <- c("incumbent", "svemnet", "gaupro", "extra_trees", "forward", "spline")
    reserved <- c("new_recipe", "prediction_context", "reference_context",
      "incumbent_source_rows", "incumbent_n_preparations", "incumbent_summary",
      "benchmark_source_row", "benchmark_source_context",
      "incumbent_recipe_n_preparations", "incumbent_recipe_n_complete",
      "incumbent_recipe_mean_score")
    if (any(reserved %in% c(schema$predictor_vars, names(weights)))) {
      stop("Reserved portfolio worklist name: ",
           paste(intersect(reserved, c(schema$predictor_vars, names(weights))), collapse = ", "))
    }
    if (!is.character(priority) || !length(priority) || anyNA(priority) ||
        anyDuplicated(priority) || any(!priority %in% allowed)) {
      stop("priority must contain distinct supported role names.")
    }
    .svem_validate_allocation_counts(budget, n_optima, n_explore, priority)
    exploration_fraction <- .svem_numeric_scalar(exploration_fraction,
      "exploration_fraction", lower = 0, upper = 1, lower_open = TRUE)
    if (!is.null(incumbent)) {
      incumbent <- .svem_integer_scalar(incumbent, "incumbent", min = 1L,
                                        max = nrow(tab))
    }
  }
  if (!is.null(benchmark)) {
    if (!portfolio) stop("A fixed benchmark requires portfolio allocation.")
    benchmark <- .svem_integer_scalar(benchmark, "benchmark", min = 1L, max = nrow(tab))
  }
  predictors <- schema$predictor_vars
  design <- setdiff(predictors, schema$blocking)
  history_keys <- if (!is.null(history)) .svem_recipe_keys(history, design, increments)
  off_grid <- sum(is.na(history_keys))
  if (off_grid) {
    warning(off_grid, " historical row(s) have missing or off-grid settings and could ",
            "not be matched to a lab-ready recipe.", call. = FALSE)
  }
  all_rows <- seq_len(nrow(tab))
  fresh <- all_rows[!keys %in% history_keys]
  if (portfolio && !is.null(incumbent)) fresh <- setdiff(fresh, incumbent)
  if (count_mode && !is.null(benchmark)) fresh <- setdiff(fresh, benchmark)
  responses <- names(weights)[weights > 0]
  raw_endpoint <- length(responses) == 1L && goals[[responses]]$goal %in% c("max", "min")
  rankings <- list()
  for (e in c("svemnet", engines)) {
    score_name <- if (e == "svemnet") "score" else paste0(e, "_score")
    values <- tab[[score_name]]
    if (raw_endpoint) {
      pred_name <- if (e == "svemnet") paste0(responses, "_pred") else
        paste0(responses, "_", e, "_pred")
      values <- tab[[pred_name]] * if (goals[[responses]]$goal == "max") 1 else -1
    }
    # A general spline can retain ordinary mean predictions while explicitly
    # declining an uninformative optimum. Raw endpoint ranking must respect it.
    if (e == "spline" && !is.null(tab[[score_name]]) && !any(is.finite(tab[[score_name]]))) {
      values[] <- NA_real_
    }
    rankings[[e]] <- values
  }
  nominations <- do.call(rbind, lapply(names(rankings), function(e) {
    v <- rankings[[e]]
    best <- .svem_stable_top(v, all_rows, keys)
    new <- .svem_stable_top(v, fresh, keys)
    state <- if (is.na(best)) "no_finite_score" else if (is.na(new)) {
      "no_unprepared_recipe"
    } else "eligible"
    data.frame(engine = e, status = state, best_recipe_id = tab$recipe_id[best],
               best_prepared = if (is.na(best) || is.null(history)) NA else !best %in% fresh,
               best_ties = length(.svem_top_ties(v, all_rows)),
               new_recipe_id = tab$recipe_id[new], new_ties = length(.svem_top_ties(v, fresh)),
               stringsAsFactors = FALSE)
  }))
  required <- unique(c(if ("incumbent" %in% priority) incumbent, benchmark))
  if (count_mode) {
    counts <- .svem_count_targets(budget, n_optima, n_explore, required, priority)
    budget <- counts$budget
  } else if (is.null(budget)) {
    budget <- if (portfolio) {
      sum(nominations$status == "eligible" & nominations$engine %in% priority) +
        as.integer(!is.null(incumbent) && "incumbent" %in% priority) +
        as.integer(!is.null(benchmark))
    } else sum(nominations$status == "eligible")
  }
  if (portfolio && !is.null(benchmark) && budget < length(required)) {
    stop("budget cannot accommodate the distinct required incumbent and fixed benchmark preparations.")
  }
  chosen <- integer()
  roles <- proposers <- character()
  if (count_mode) {
    merges <- data.frame(slot = integer(), recipe_id = tab$recipe_id[integer()],
      role = character(), source = character(), previous_nominated_by = character(),
      stringsAsFactors = FALSE)
  }
  add <- function(i, role, source) {
    if (is.na(i)) return(invisible(NULL))
    same <- match(i, chosen)
    if (!is.na(same)) {
      if (count_mode) {
        if (source %in% strsplit(proposers[same], " + ", fixed = TRUE)[[1L]]) return(invisible(NULL))
        merges <<- rbind(merges, data.frame(slot = same, recipe_id = tab$recipe_id[i],
          role = role, source = source, previous_nominated_by = proposers[same],
          stringsAsFactors = FALSE))
      }
      proposers[same] <<- paste(proposers[same], source, sep = " + ")
    } else if (length(chosen) < budget) {
      chosen <<- c(chosen, i)
      roles <<- c(roles, role)
      proposers <<- c(proposers, source)
    }
  }
  if (portfolio) {
    if (!is.null(benchmark) || count_mode) {
      if (!is.null(incumbent) && "incumbent" %in% priority) {
        add(incumbent, "incumbent_repeat", "observed")
      }
      if (!is.null(benchmark)) add(benchmark, "benchmark_repeat", "benchmark")
      if (!is.null(benchmark) && !is.null(incumbent) && incumbent == benchmark && "incumbent" %in% priority) {
        roles[match(benchmark, chosen)] <- "incumbent_benchmark_repeat"
      }
    }
    for (source in priority) {
      if (source == "incumbent") {
        if (!is.null(incumbent) && !incumbent %in% chosen) add(incumbent, "incumbent_repeat", "observed")
      } else {
        i <- match(source, nominations$engine)
        if (!is.na(i)) {
          row <- match(nominations$new_recipe_id[i], tab$recipe_id)
          if (!count_mode || row %in% chosen || sum(roles == "optimum") < counts$target_optima)
            add(row, "optimum", source)
        }
      }
    }
    nominations$endorses_incumbent <- !is.null(incumbent) &
      !is.na(nominations$best_recipe_id) &
      nominations$best_recipe_id %in% tab$recipe_id[incumbent]
    nominations$best_prepared[nominations$endorses_incumbent] <- TRUE
    if (!is.null(incumbent) && incumbent %in% chosen) {
      for (source in nominations$engine[nominations$endorses_incumbent]) {
        add(incumbent, "incumbent_repeat", source)
      }
    }
  } else {
    for (i in seq_len(nrow(nominations))) {
      add(match(nominations$new_recipe_id[i], tab$recipe_id), "optimum", nominations$engine[i])
    }
  }
  nominations$slot <- match(nominations$new_recipe_id, tab$recipe_id[chosen])
  eligible <- nominations$status == "eligible"
  nominations$status[eligible] <- ifelse(is.na(nominations$slot[eligible]), "dropped_by_budget", "selected")
  if (portfolio) {
    if (count_mode) nominations$status[eligible & is.na(nominations$slot)] <- "dropped_by_model_count"
    nominations$status[eligible & !nominations$engine %in% priority] <- "not_prioritized"
    incumbent_slot <- match(incumbent, chosen)
    incumbent_status <- if (is.null(incumbent)) "not_supplied" else
      if (!"incumbent" %in% priority) "not_prioritized" else
        if (is.na(incumbent_slot)) "dropped_by_budget" else "selected"
    incumbent_row <- nominations[NA_integer_, , drop = FALSE]
    incumbent_row$engine <- "incumbent"
    incumbent_row$status <- incumbent_status
    incumbent_row$best_recipe_id <- tab$recipe_id[if (is.null(incumbent)) NA_integer_ else incumbent]
    incumbent_row$best_prepared <- if (is.null(incumbent)) NA else TRUE
    incumbent_row$best_ties <- if (is.null(incumbent)) 0L else 1L
    incumbent_row$new_ties <- 0L
    incumbent_row$slot <- if (is.null(incumbent)) NA_integer_ else incumbent_slot
    incumbent_row$endorses_incumbent <- FALSE
    nominations <- rbind(incumbent_row, nominations)
    if (!is.null(benchmark)) {
      benchmark_row <- incumbent_row
      benchmark_row$engine <- "benchmark"
      benchmark_row$status <- "selected"
      benchmark_row$best_recipe_id <- tab$recipe_id[benchmark]
      benchmark_row$best_prepared <- TRUE
      benchmark_row$best_ties <- 1L
      benchmark_row$slot <- match(benchmark, chosen)
      benchmark_row$endorses_incumbent <- !is.null(incumbent) && incumbent == benchmark
      nominations <- rbind(nominations[1L, , drop = FALSE], benchmark_row,
        nominations[-1L, , drop = FALSE])
    }
    rownames(nominations) <- NULL
  }

  # Missing/constant intervals may yield a neutral published score. Such a
  # score supplies no information to support a new exploration allocation.
  active <- intersect(names(uncertainty_weights)[uncertainty_weights > 0],
                      names(Filter(Negate(is.null), anchors)))
  if (portfolio) {
    active <- active[vapply(anchors[active], function(a) {
      is.numeric(a) && length(a) == 2L && all(is.finite(a)) && a[2L] > a[1L]
    }, logical(1))]
  }
  width <- rep(NA_real_, nrow(tab))
  if (length(active)) {
    terms <- lapply(active, function(r) {
      w <- tab[[paste0(r, "_upr")]] - tab[[paste0(r, "_lwr")]]
      uncertainty_weights[[r]] * w / diff(anchors[[r]])
    })
    width <- Reduce("+", terms)
  }
  available <- setdiff(fresh, chosen)
  ties <- .svem_top_ties(tab$uncertainty_measure, available[is.finite(width[available])])
  secondary <- .svem_top_ties(width, ties)
  picks <- c(.svem_stable_top(tab$score, ties, keys),
             .svem_stable_top(tab$score, secondary, keys))
  tie_range <- if (length(ties) && any(is.finite(tab$score[ties]))) {
    range(tab$score[ties][is.finite(tab$score[ties])])
  } else c(NA_real_, NA_real_)
  u <- tab$uncertainty_measure[available[is.finite(width[available])]]
  u <- u[is.finite(u)]
  top <- if (length(u)) max(u) else NA_real_
  exact_ties <- sum(u == top)
  display_only <- sum(round(u, 3L) == round(top, 3L)) - length(ties)
  spread_range <- function(x) {
    x <- x[is.finite(x)]
    if (length(x)) range(x) else c(NA_real_, NA_real_)
  }
  tie_endpoint_ranges <- do.call(rbind, lapply(names(weights), function(r) {
    predicted <- spread_range(tab[[paste0(r, "_pred")]][ties])
    raw <- spread_range(tab[[paste0(r, "_upr")]][ties] - tab[[paste0(r, "_lwr")]][ties])
    contribution <- spread_range(tab[[paste0(r, "_ciw_w")]][ties])
    data.frame(response = r, min_prediction = predicted[1L], max_prediction = predicted[2L],
               min_width = raw[1L], max_width = raw[2L],
               min_weighted_spread = contribution[1L], max_weighted_spread = contribution[2L])
  }))
  exploration <- data.frame(rule = c("best_predicted_tie", "widest_spread"),
                            recipe_id = tab$recipe_id[picks], tie_count = length(ties),
                            exact_ties = exact_ties, display_only_ties = max(0L, display_only),
                            min_desirability = tie_range[1L], max_desirability = tie_range[2L],
                            uncertainty_measure = tab$uncertainty_measure[picks],
                            strict_width = width[picks], score = tab$score[picks])
  spread_role <- if (length(ties) == 1L) "maximum_spread" else "best_predicted_tie"
  if (!portfolio) {
    add(picks[1L], spread_role, "svemnet")
    sd_name <- paste0(gp_endpoint, "_gp_latent_sd")
    if (!is.null(gp_endpoint) && sd_name %in% names(tab)) {
      add(.svem_stable_top(tab[[sd_name]], setdiff(fresh, chosen), keys),
          "latent_sd", "gaupro")
    }
  }
  remaining <- setdiff(fresh, chosen)
  remaining <- remaining[is.finite(width[remaining])]
  fraction <- .1
  explore_slots <- if (count_mode && !is.null(n_explore)) as.integer(n_explore) else budget - length(chosen)
  if (portfolio) {
    remaining <- remaining[is.finite(tab$uncertainty_measure[remaining])]
    allocation_details <- list(policy = allocation, priority = priority,
      requested_fraction = exploration_fraction, effective_fraction = 0,
      eligible_count = length(remaining), shortlist_count = 0L, selected_count = 0L,
      cutoff = NA_real_, cutoff_ties = 0L, cutoff_ties_retained = 0L)
    if (explore_slots > 0L && length(remaining) && length(design)) {
      shortlist <- min(length(remaining), max(explore_slots,
        ceiling(exploration_fraction * length(remaining) - 1e-9)))
      fraction <- shortlist / length(remaining)
      ranked <- sort(tab$uncertainty_measure[remaining], decreasing = TRUE)
      cutoff <- ranked[shortlist]
      allocation_details$effective_fraction <- fraction
      allocation_details$shortlist_count <- shortlist
      allocation_details$cutoff <- cutoff
      allocation_details$cutoff_ties <- sum(ranked == cutoff)
      allocation_details$cutoff_ties_retained <- sum(ranked[seq_len(shortlist)] == cutoff)
    }
  }
  if (explore_slots > 0L && length(remaining) && length(design)) {
    # Sorting by recipe key before PAM also makes ties at the top-fraction
    # cutoff independent of input order.
    remaining <- remaining[order(keys[remaining], method = "radix")]
    med <- svem_select_candidates(tab[remaining, , drop = FALSE], by = "uncertainty_measure",
                                  top_frac = fraction, k = explore_slots,
                                  predictor_cols = design)
    for (i in remaining[med]) add(i, "high_spread_medoid", "svemnet")
    if (portfolio) allocation_details$selected_count <- length(med)
  }

  work <- tab[chosen, , drop = FALSE]
  work$slot <- seq_along(chosen)
  work$role <- roles
  work$nominated_by <- proposers
  if (portfolio) work$new_recipe <- chosen %in% fresh
  reference <- if (is.null(history)) tab[integer(), design, drop = FALSE] else history[design]
  spans <- lapply(tab[design], function(x) if (is.numeric(x)) diff(range(x)) else NA_real_)
  work$nearest_distance <- numeric(length(chosen))
  for (j in seq_along(chosen)) {
    row <- work[j, design, drop = FALSE]
    work$nearest_distance[j] <- .svem_nearest_recipe(row, reference, spans)
    reference <- rbind(reference, row)
  }
  selections <- lapply(seq_along(chosen), function(j) {
    svem_select_from_score_table(work[j, , drop = FALSE], target = "slot", k = 0L,
                                 predictor_cols = predictors,
                                 label = paste(work$nominated_by[j], work$role[j]))
  })
  missing <- budget - nrow(work)
  if (count_mode) {
    counts$selected_optima <- sum(roles == "optimum")
    counts$target_explore <- explore_slots
    counts$selected_explore <- sum(roles == "high_spread_medoid")
    counts$unused_capacity <- missing
    counts$unfilled_model_capacity <- counts$target_optima - counts$selected_optima
    counts$model_shortfall <- if (is.null(n_optima)) 0L else counts$unfilled_model_capacity
    counts$exploration_shortfall <- counts$target_explore - counts$selected_explore
    counts$intentionally_unused <- max(0, budget - counts$distinct_repeats -
                                      counts$target_optima - counts$target_explore)
    counts$requested_shortfall <- counts$model_shortfall + counts$exploration_shortfall
    missing <- counts$requested_shortfall
    if (missing > 0L) warning("Requested allocation is short by ", counts$model_shortfall,
      " model and ", counts$exploration_shortfall, " exploration preparation(s).", call. = FALSE)
  } else if (missing > 0) warning("Worklist is short by ", missing, " preparation(s): available ",
                           "nominations and the current eligible uncertainty shortlist ",
                           "do not fill the budget.", call. = FALSE)
  if (portfolio) {
    attr(work, "budget") <- budget
    attr(work, "shortfall") <- missing
  }
  result <- list(worklist = work, nominations = nominations, exploration = exploration,
                 tie_endpoint_ranges = tie_endpoint_ranges, selections = selections,
                 shortfall = missing, historical_unmatched = off_grid)
  if (portfolio) {
    if (count_mode) {
      result$allocation_counts <- counts
      result$merged_nominations <- merges
    }
    result$allocation_details <- allocation_details
    result$uncertainty_eligibility <- data.frame(recipe_id = tab$recipe_id,
      informative_width = is.finite(width), standardized_width = width,
      finite_score = is.finite(tab$uncertainty_measure),
      unprepared = all_rows %in% fresh)
  }
  result
}
