# Shared composition helpers for the default preparation workflow.
.svem_plan_control <- function(x, reserved, label) {
  if (!is.list(x) || (length(x) && (is.null(names(x)) || anyNA(names(x)) ||
      any(!nzchar(names(x))) || anyDuplicated(names(x))))) {
    stop(label, " must be a list with unique argument names.", call. = FALSE)
  }
  bad <- intersect(names(x), reserved)
  if (length(bad)) stop(label, " cannot replace: ", paste(bad, collapse = ", "), call. = FALSE)
  x
}

# Numeric day/plate IDs are compared by value, independent of display format.
.svem_matches_block <- function(observed, reference) {
  if (is.numeric(observed)) {
    key <- if (is.numeric(reference)) reference else
      suppressWarnings(as.numeric(as.character(reference)))
    return(!is.na(observed) & !is.na(key) & observed == key)
  }
  !is.na(observed) & as.character(observed) == as.character(reference)
}

# Distinguish nearby physical block IDs in provenance without changing the
# preparation keys or their rounding conventions for scientific factors.
.svem_block_labels <- function(x) {
  if (!is.numeric(x)) return(as.character(x))
  text <- rep(NA_character_, length(x))
  valid <- !is.na(x)
  text[valid] <- sprintf("%.15g", ifelse(x[valid] == 0, 0, x[valid]))
  longer <- valid & as.numeric(text) != x
  text[longer] <- sprintf("%.17g", x[longer])
  text
}

.svem_observed_pool <- function(objects, data, candidates, groups, increments,
                                 incumbent, reference) {
  obj <- objects[[1L]]
  schema <- obj$sampling_schema
  vars <- schema$predictor_vars
  blocks <- schema$blocking
  design <- setdiff(vars, blocks)
  if (!is.data.frame(data) || !nrow(data) || !all(vars %in% names(data))) {
    stop("An incumbent needs measured `data` containing every predictor.", call. = FALSE)
  }
  if (is.null(candidates)) stop("An incumbent needs an explicit feasible candidate pool.", call. = FALSE)
  discrete <- .svem_discrete_map(schema)
  x <- .svem_candidate_frame(candidates, schema, discrete, groups, 1e-8)
  for (b in blocks) {
    if (length(unique(x[[b]])) != 1L) stop("Set one planned level for blocking factor '", b, "'.", call. = FALSE)
  }
  explicit <- is.numeric(incumbent) && length(incumbent) == 1L
  if (explicit) {
    row <- .svem_integer_scalar(incumbent, "incumbent", min = 1L, max = nrow(data))
    if (!is.null(reference)) stop("Use either a historical row index or reference block values, not both.")
    rows <- row
  } else {
    if (!isTRUE(incumbent)) stop("incumbent must be TRUE, FALSE, or a historical row index.")
    if (is.null(reference)) reference <- as.list(x[1L, blocks, drop = FALSE])
    if (!is.list(reference) || (length(reference) && (is.null(names(reference)) ||
        anyNA(names(reference)) || anyDuplicated(names(reference)) ||
        !setequal(names(reference), blocks)))) {
      stop("reference must name exactly the blocking factors, each with one historical value.")
    }
    if (length(blocks) && !setequal(names(reference), blocks)) {
      stop("reference must name every blocking factor.")
    }
    rows <- seq_len(nrow(data))
    for (b in names(reference)) {
      z <- reference[[b]]
      if (length(z) != 1L || is.na(z)) stop("Each reference block value must be a single nonmissing value.")
      rows <- rows[.svem_matches_block(data[[b]][rows], z)]
    }
  }
  if (!length(rows)) stop("No historical preparations in the reference block context. Supply an incumbent row or set incumbent = FALSE.")
  # Validate all historical physical recipes in the planned execution context.
  # Other contexts are reported, never pooled automatically for selection.
  h_all <- data[, vars, drop = FALSE]
  for (b in blocks) h_all[[b]] <- rep(x[[b]][1L], nrow(h_all))
  valid <- vapply(seq_len(nrow(h_all)), function(i) {
    tryCatch({
      .svem_candidate_frame(h_all[i, , drop = FALSE], schema, discrete, groups, 1e-8)
      !is.na(.svem_recipe_keys(h_all[i, , drop = FALSE], design, increments))
    }, error = function(e) FALSE)
  }, logical(1))
  excluded <- rows[!valid[rows]]
  rows <- rows[valid[rows]]
  h <- h_all[rows, , drop = FALSE]
  if (!length(rows)) stop("No eligible exact historical recipe satisfies the current mixture/support/increment constraints.")
  existing <- .svem_recipe_keys(x, design, increments)
  historical <- .svem_recipe_keys(h, design, increments)
  extra <- which(!historical %in% existing & !duplicated(historical))
  supplied_ids <- "recipe_id" %in% names(candidates)
  if (supplied_ids) {
    ids <- candidates$recipe_id
    if (!(is.numeric(ids) || is.character(ids)) || anyNA(ids) ||
        any(!nzchar(as.character(ids))) || anyDuplicated(ids)) {
      stop("Candidate recipe_id values must be unique and nonmissing numeric or character labels.")
    }
    x$recipe_id <- as.character(ids)
  }
  if (length(extra)) {
    add <- h[extra, , drop = FALSE]
    if (supplied_ids) {
      stable_ids <- paste0("observed_", match(historical[extra],
                                             sort(unique(historical), method = "radix")))
      newids <- make.unique(c(as.character(ids), stable_ids), sep = "_")
      add$recipe_id <- tail(newids, length(extra))
    }
    x <- rbind(x, add)
  }
  list(candidates = x, rows = rows, keys = historical, excluded_rows = excluded,
       reference = reference, explicit = explicit,
       eligible_history_rows = which(valid),
       all_keys = .svem_recipe_keys(h_all, design, increments),
       blocks = blocks)
}

.svem_context_label <- function(data, blocks) {
  if (!length(blocks)) return(rep("unblocked", nrow(data)))
  vapply(seq_len(nrow(data)), function(i) {
    paste(paste0(blocks, "=", vapply(data[i, blocks, drop = FALSE],
      .svem_block_labels, character(1))), collapse = "; ")
  }, character(1))
}

.svem_observed_choice <- function(objects, data, prepared, parameters, des,
                                  weights, combine, floor, goals, pool_keys,
                                  summary = c("best", "mean")) {
  summary <- match.arg(summary)
  active <- names(weights)[weights > 0]
  values <- lapply(active, function(r) {
    y <- eval(objects[[r]]$formula[[2L]], data, environment(objects[[r]]$formula))
    if (!is.numeric(y) || is.matrix(y)) stop("Incumbent scoring requires numeric Gaussian response values.")
    as.numeric(y)
  })
  names(values) <- active
  complete <- Reduce("&", lapply(values, is.finite))
  assess <- function(rows, average) {
    rows <- rows[complete[rows]]
    if (!length(rows)) return(NULL)
    groups <- if (average) split(rows, prepared$all_keys[rows]) else as.list(rows)
    means <- lapply(values, function(v) vapply(groups, function(g) mean(v[g]), numeric(1)))
    contributions <- lapply(active, function(r) {
      d <- des(means[[r]], parameters[[r]])
      if (combine == "geom") d <- log(floor + (1 - floor) * d)
      weights[[r]] * d
    })
    score <- Reduce("+", contributions)
    if (combine == "geom") score <- exp(score)
    rank <- score
    if (length(active) == 1L && goals[[active]]$goal %in% c("max", "min")) {
      rank <- means[[active]] * if (goals[[active]]$goal == "max") 1 else -1
    }
    first <- vapply(groups, "[", integer(1), 1L)
    keys <- prepared$all_keys[first]
    pick <- .svem_stable_top(rank, seq_along(groups), keys)
    if (is.na(pick)) return(NULL)
    list(groups = groups, values = means, score = score, rank = rank,
         keys = keys, pick = pick, first = first)
  }
  use_mean <- summary == "mean" && !prepared$explicit
  selected <- assess(prepared$rows, use_mean)
  if (is.null(selected)) stop("No eligible incumbent has all positively weighted outcomes observed. Supply complete observations or set incumbent = FALSE.")
  pick <- selected$pick
  source_rows <- unname(selected$groups[[pick]])
  index <- match(selected$keys[pick], pool_keys)
  if (is.na(index)) stop("Selected historical recipe was not retained in the feasible candidate pool.")
  contexts <- .svem_context_label(data, prepared$blocks)
  block_data <- data[prepared$blocks]
  block_data[] <- lapply(block_data, .svem_block_labels)
  context_keys <- .svem_recipe_keys(block_data, prepared$blocks)
  # An explicit row can have an unknown source block. Do not pool different
  # unknown contexts into a recipe mean or lose its provenance to NA indexing.
  unknown_context <- which(is.na(context_keys))
  context_keys[unknown_context] <- paste0("unknown-context-row-", unknown_context)
  eligible <- prepared$eligible_history_rows
  by_context <- lapply(split(eligible, context_keys[eligible]), function(rows) {
    z <- assess(rows, summary == "mean")
    if (is.null(z)) return(NULL)
    i <- z$pick
    source <- z$groups[[i]]
    data.frame(context = contexts[source[1L]], context_key = context_keys[source[1L]],
      source_row = source[1L],
      source_rows = paste(source, collapse = ","), n_preparations = length(source),
      observed_score = unname(z$score[i]), ranking_value = unname(z$rank[i]),
      recipe_id = NA_character_, recipe_key = z$keys[i], stringsAsFactors = FALSE)
  })
  by_context <- do.call(rbind, by_context)
  rownames(by_context) <- NULL
  source_context <- contexts[source_rows[1L]]
  higher_elsewhere <- by_context$context_key != context_keys[source_rows[1L]] &
    by_context$ranking_value > selected$rank[pick]
  recipe_rows <- eligible[prepared$all_keys[eligible] == selected$keys[pick] &
    !is.na(context_keys[eligible]) & context_keys[eligible] == context_keys[source_rows[1L]]]
  recipe_mean <- assess(recipe_rows, TRUE)
  list(index = index, source_row = source_rows[1L], source_rows = source_rows,
       n_preparations = length(source_rows),
       recipe_source_rows = recipe_rows,
       recipe_n_preparations = length(recipe_rows),
       recipe_n_complete = sum(complete[recipe_rows]),
       recipe_mean_score = unname(recipe_mean$score[recipe_mean$pick]),
       recipe_mean_responses = lapply(recipe_mean$values,
         function(v) unname(v[recipe_mean$pick])),
       summary = if (prepared$explicit) "explicit_row" else summary,
       observed_score = unname(selected$score[pick]),
       observed_responses = lapply(selected$values, function(v) unname(v[pick])),
       ties = length(.svem_top_ties(selected$rank, seq_along(selected$groups))),
       reference = prepared$reference, source_context = source_context,
       explicit = prepared$explicit, excluded_rows = prepared$excluded_rows,
       incomplete_rows = prepared$rows[!complete[prepared$rows]],
       observations = data.frame(source_row = selected$first,
          source_rows = vapply(selected$groups, paste, character(1), collapse = ","),
          n_preparations = lengths(selected$groups),
          observed_score = unname(selected$score), eligible = is.finite(selected$rank)),
       by_context = by_context, higher_raw_in_other_context = any(higher_elsewhere),
       source_data = data[source_rows, , drop = FALSE],
       rule = if (use_mean) "best observed recipe mean in the reference context" else
         "best observed independent preparation in the reference context")
}

.svem_forward_scores <- function(tab, x, forward, parameters, des, weights, combine, floor) {
  contributions <- list()
  status <- list()
  for (r in names(weights)) {
    pr <- tryCatch(as.numeric(predict(forward[[r]], newdata = x, debias = FALSE)), error = identity)
    state <- if (inherits(pr, "error")) "prediction_failed" else "ok"
    msg <- if (inherits(pr, "error")) conditionMessage(pr) else ""
    if (inherits(pr, "error")) pr <- rep(NA_real_, nrow(x))
    if (state == "ok" && !any(is.finite(pr))) state <- "no_finite_predictions"
    tab[[paste0(r, "_forward_pred")]] <- pr
    status[[r]] <- data.frame(engine = "forward", response = r, status = state, message = msg)
    if (weights[[r]] > 0) {
      d <- des(pr, parameters[[r]])
      d[!is.finite(pr)] <- NA_real_
      if (combine == "geom") d <- log(floor + (1 - floor) * d)
      contributions[[r]] <- weights[[r]] * d
    }
  }
  score <- Reduce("+", contributions)
  if (combine == "geom") score <- exp(score)
  tab$forward_score <- score
  if (!any(is.finite(score))) warning("Forward fits cannot score the joint goal; inspect engine_status.", call. = FALSE)
  list(table = tab, status = do.call(rbind, status))
}
