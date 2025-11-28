#' Select best row and diverse candidates from an SVEM score table
#'
#' @description
#' Given a scored random-search table (e.g. from \code{svem_score_random()}),
#' pick a single "best" row under a chosen objective column and sample a
#' small, diverse set of medoid candidates from the top of that ranking.
#' Any per-response CI columns (e.g. \code{*_lwr} / \code{*_upr}) present in
#' \code{score_table} are carried through unchanged.
#'
#' Optionally, a string \code{label} can be supplied to annotate the returned
#' \code{best} row and \code{candidates} by appending that label to a
#' \code{"Notes_from_SVEMnet"} column. If \code{"Notes_from_SVEMnet"} is missing, it is created. If it
#' exists and is nonempty, the label is appended with \code{"; "} as a
#' separator.
#'
#' @param score_table Data frame with predictors, responses, scores, and
#'   \code{uncertainty_measure}, typically \code{scored$score_table} from
#'   \code{\link{svem_score_random}}. When medoids are requested
#'   (\code{k > 0}), the predictor columns used for clustering are taken
#'   from the \code{"svem_predictor_cols"} attribute by default. If that
#'   attribute is missing, a numeric-column heuristic is used. If you
#'   accidentally pass the full \code{scored} list, a helpful error is
#'   thrown reminding you to use \code{scored$score_table}.
#' @param target Character scalar naming the column in \code{score_table}
#'   to optimize (e.g. \code{"score"}, \code{"wmt_score"},
#'   \code{"uncertainty_measure"}).
#' @param direction Either \code{"max"} or \code{"min"} indicating whether
#'   larger or smaller values of \code{target} are preferred.
#' @param k Integer; desired number of medoid candidates to return. If
#'   \code{k <= 0}, only the \code{best} row is returned and no clustering
#'   is performed.
#' @param top_type Either \code{"frac"} or \code{"n"} specifying whether
#'   \code{top} is a fraction of rows or an integer count.
#' @param top Value for the top set: a fraction in (0,1] if
#'   \code{top_type = "frac"}, or an integer \code{>= 1} if
#'   \code{top_type = "n"}.
#' @param predictor_cols Optional character vector of predictor column names
#'   used to measure diversity in the PAM step when \code{k > 0}. When
#'   \code{NULL} (default), the function first tries
#'   \code{attr(score_table, "svem_predictor_cols")}. If that is unavailable,
#'   it falls back to using numeric, non-meta columns (excluding e.g.
#'   \code{*_pred}, \code{*_des}, \code{*_lwr}, \code{*_upr},
#'   \code{*_ciw_w}, \code{*_p_in_spec_mean}, \code{*_in_spec_point},
#'   \code{score}, \code{wmt_score}, \code{uncertainty_measure},
#'   \code{p_joint_mean}, \code{joint_in_spec_point},
#'   \code{candidate_type}, \code{selection_label}, \code{Notes_from_SVEMnet}). If no
#'   usable predictor columns can be inferred, a warning is issued and only
#'   \code{best} is returned.
#' @param label Optional character scalar. When non-\code{NULL}, this label
#'   is appended into a \code{"Notes_from_SVEMnet"} column for the returned \code{best}
#'   row and \code{candidates}. If \code{"Notes_from_SVEMnet"} is missing, it is created;
#'   if present and nonempty, the label is appended using \code{"; "} as
#'   separator.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{best}{One-row data frame at the optimum of \code{target}
#'     under the specified \code{direction}, including any columns present
#'     in \code{score_table} (e.g. \code{*_lwr} / \code{*_upr}).}
#'   \item{candidates}{Data frame of medoid candidates (possibly
#'     empty or \code{NULL}) drawn from the top \code{top} of the ranking
#'     on \code{target}, with all columns carried through from
#'     \code{score_table}.}
#'   \item{call}{The matched call, including all arguments used to
#'     create this selection object.}
#' }
#'
#' @seealso
#' \code{\link{svem_score_random}},
#' \code{svem_select_candidates()}
#'
#' @examples
#' \donttest{
#' ## ------------------------------------------------------------------------
#' ## Selecting optimal and exploration candidates from a scored SVEM table
#' ## ------------------------------------------------------------------------
#'
#' data(lipid_screen)
#'
#' # Build expansion and fit three SVEM models (Potency, Size, PDI)
#' spec <- bigexp_terms(
#'   Potency ~ PEG + Helper + Ionizable + Cholesterol +
#'     Ionizable_Lipid_Type + N_P_ratio + flow_rate,
#'   data             = lipid_screen,
#'   factorial_order  = 3,
#'   polynomial_order = 3,
#'   include_pc_2way  = TRUE,
#'   include_pc_3way  = FALSE
#' )
#'
#' form_pot <- bigexp_formula(spec, "Potency")
#' form_siz <- bigexp_formula(spec, "Size")
#' form_pdi <- bigexp_formula(spec, "PDI")
#'
#' set.seed(1)
#' fit_pot <- SVEMnet(form_pot, lipid_screen)
#' fit_siz <- SVEMnet(form_siz, lipid_screen)
#' fit_pdi <- SVEMnet(form_pdi, lipid_screen)
#'
#' objs <- list(Potency = fit_pot, Size = fit_siz, PDI = fit_pdi)
#'
#' goals <- list(
#'   Potency = list(goal = "max", weight = 0.6),
#'   Size    = list(goal = "min", weight = 0.3),
#'   PDI     = list(goal = "min", weight = 0.1)
#' )
#'
#' mix <- list(list(
#'   vars  = c("PEG", "Helper", "Ionizable", "Cholesterol"),
#'   lower = c(0.01, 0.10, 0.10, 0.10),
#'   upper = c(0.05, 0.60, 0.60, 0.60),
#'   total = 1.0
#' ))
#'
#' set.seed(3)
#' scored <- svem_score_random(
#'   objects         = objs,
#'   goals           = goals,
#'   n               = 20000,
#'   mixture_groups  = mix,
#'   combine         = "geom",
#'   numeric_sampler = "random",
#'   verbose         = FALSE
#' )
#'
#' # The scored table contains predictors, <resp>_pred, <resp>_des, score,
#' # uncertainty_measure, and per-response CI columns.
#' names(scored$score_table)
#'
#' ## ------------------------------------------------------------------------
#' ## 1) Optimal candidates by multi-response score
#' ## ------------------------------------------------------------------------
#'
#' opt_sel <- svem_select_from_score_table(
#'   score_table = scored$score_table,
#'   target      = "score",    # column to optimize
#'   direction   = "max",      # maximize score
#'   k           = 5,          # 5 medoid candidates
#'   top_type    = "frac",     # sample medoids from top fraction
#'   top         = 0.02,       # top 2% by score
#'   label       = "round1_optimal"
#' )
#'
#' # Single best row (highest score) with predictions and CIs
#' opt_sel$best
#'
#' # Diverse high-score candidates (medoids)
#' head(opt_sel$candidates)
#'
#' ## ------------------------------------------------------------------------
#' ## 2) Exploration candidates: highest model uncertainty
#' ## ------------------------------------------------------------------------
#'
#' explore_sel <- svem_select_from_score_table(
#'   score_table = scored$score_table,
#'   target      = "uncertainty_measure",  # scalar uncertainty
#'   direction   = "max",                  # look for high-uncertainty settings
#'   k           = 5,
#'   top_type    = "frac",
#'   top         = 0.05,                   # top 5% most uncertain
#'   label       = "round1_explore"
#' )
#'
#' explore_sel$best
#' head(explore_sel$candidates)
#'
#' ## ------------------------------------------------------------------------
#' ## 3) Re-ranking by design-space assurance (if p_joint_mean is available)
#' ## ------------------------------------------------------------------------
#'
#' # If svem_score_random() was called with a non-NULL `specs` argument,
#' # the score_table may contain p_joint_mean (joint mean-level in-spec prob).
#' if ("p_joint_mean" %in% names(scored$score_table)) {
#'   inspec_sel <- svem_select_from_score_table(
#'     score_table = scored$score_table,
#'     target      = "p_joint_mean",   # maximize mean-level spec assurance
#'     direction   = "max",
#'     k           = 5,
#'     top_type    = "frac",
#'     top         = 0.10,
#'     label       = "round1_inspec"
#'   )
#'
#'   inspec_sel$best
#'   head(inspec_sel$candidates)
#' }
#'
#' ## ------------------------------------------------------------------------
#' ## 4) Selecting the best existing run from a scored original data table
#' ## ------------------------------------------------------------------------
#'
#' # If svem_score_random() was called with data = lipid_screen,
#' # the original_data_scored component contains scored existing runs.
#' if (!is.null(scored$original_data_scored)) {
#'   best_existing <- svem_select_from_score_table(
#'     score_table = scored$original_data_scored,
#'     target      = "score",
#'     direction   = "max",
#'     k           = 0,          # k <= 0: only return the single best row
#'     top_type    = "frac",
#'     top         = 1.0,
#'     label       = "round1_existing_best"
#'   )
#'
#'   best_existing$best
#' }
#' }
#'
#' @export
svem_select_from_score_table <- function(score_table,
                                         target         = "score",
                                         direction      = c("max", "min"),
                                         k              = 5,
                                         top_type       = c("frac", "n"),
                                         top            = 0.1,
                                         predictor_cols = NULL,
                                         label          = NULL) {

  # Helpful catch: user passed full `scored` object instead of `scored$score_table`
  if (!inherits(score_table, "data.frame") &&
      is.list(score_table) &&
      "score_table" %in% names(score_table) &&
      is.data.frame(score_table$score_table)) {
    stop(
      "It looks like you passed the full `scored` object.\n",
      "Please call:\n",
      "  svem_select_from_score_table(score_table = scored$score_table, ...)",
      call. = FALSE
    )
  }

  if (!is.data.frame(score_table))
    stop("score_table must be a data.frame.")

  direction <- match.arg(direction)
  top_type  <- match.arg(top_type)

  # Capture the call (arguments) for output
  mc <- match.call()

  if (!is.character(target) || length(target) != 1L || !nzchar(target))
    stop("`target` must be a nonempty character scalar naming a column.")
  if (!(target %in% colnames(score_table)))
    stop("Column `", target, "` not found in `score_table`.")

  if (!is.numeric(k) || length(k) != 1L || !is.finite(k))
    stop("`k` must be a single finite numeric value.")
  k <- as.integer(k)

  if (!is.null(label)) {
    if (!is.character(label) || length(label) != 1L) {
      stop("`label` must be a single character string or NULL.")
    }
  }

  # values for ranking
  vals <- score_table[[target]]
  if (!is.numeric(vals))
    stop("Column `", target, "` must be numeric for ranking.")

  decreasing_flag <- (direction == "max")
  ord <- order(vals, decreasing = decreasing_flag, na.last = TRUE)

  if (!length(ord))
    stop("No rows available in score_table for selection.")

  # best index and row
  best_idx <- ord[1L]
  best     <- score_table[best_idx, , drop = FALSE]

  # candidate medoids
  candidates <- NULL
  n_total    <- length(ord)

  if (k > 0L && n_total > 0L) {

    ## --- determine predictor_cols (needed only for PAM) ---

    if (is.null(predictor_cols)) {
      predictor_cols <- attr(score_table, "svem_predictor_cols")
    }

    # ensure predictors actually exist in this table
    if (!is.null(predictor_cols)) {
      predictor_cols <- intersect(predictor_cols, names(score_table))
    }

    # heuristic fallback if attributes / explicit predictors are not usable
    if (is.null(predictor_cols) || !length(predictor_cols)) {
      nm     <- names(score_table)
      is_num <- vapply(score_table, is.numeric, logical(1L))

      is_meta <- grepl("(_pred|_des|_lwr|_upr|_ciw_w|_p_in_spec_mean|_in_spec_point)$", nm) |
        nm %in% c("score",
                  "wmt_score",
                  "uncertainty_measure",
                  "p_joint_mean",
                  "joint_in_spec_point",
                  "candidate_type",
                  "selection_label",
                  "Notes_from_SVEMnet")

      predictor_cols <- nm[is_num & !is_meta]
    }

    if (!length(predictor_cols)) {
      warning(
        "svem_select_from_score_table(): could not infer predictor columns for PAM; ",
        "returning only `best` (no medoid candidates)."
      )
    } else {

      ## --- convert top_type/top into a top_frac in (0, 1] ---

      if (top_type == "frac") {
        if (!(is.numeric(top) && length(top) == 1L &&
              is.finite(top) && top > 0 && top <= 1)) {
          stop("When top_type = 'frac', `top` must be in (0, 1].")
        }
        top_frac <- top
      } else {  # top_type == "n"
        if (!(is.numeric(top) && length(top) == 1L &&
              is.finite(top) && top >= 1)) {
          stop("When top_type = 'n', `top` must be >= 1.")
        }
        n_top   <- min(n_total, as.integer(top))
        top_frac <- n_top / n_total
      }

      ## --- delegate PAM / medoid selection to svem_select_candidates() ---

      med_idx <- svem_select_candidates(
        table          = score_table,
        by             = target,
        top_frac       = top_frac,
        k              = k,
        predictor_cols = predictor_cols,
        direction      = direction,
        metric         = "gower"
      )

      if (length(med_idx)) {
        candidates <- score_table[med_idx, , drop = FALSE]
      } else {
        candidates <- score_table[0L, , drop = FALSE]
      }
    }
  }

  # Optionally append label into a Notes column for best and candidates
  if (!is.null(label)) {
    append_label <- function(df) {
      if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(df)
      if (!("Notes_from_SVEMnet" %in% names(df))) {
        df$Notes_from_SVEMnet <- ""
      }
      old_notes  <- as.character(df$Notes_from_SVEMnet)
      to_replace <- is.na(old_notes) | old_notes == ""
      old_notes[to_replace]  <- label
      old_notes[!to_replace] <- paste0(old_notes[!to_replace], "; ", label)
      df$Notes_from_SVEMnet  <- old_notes
      df
    }

    best       <- append_label(best)
    candidates <- append_label(candidates)
  }

  list(
    best       = best,
    candidates = candidates,
    call       = mc
  )
}
