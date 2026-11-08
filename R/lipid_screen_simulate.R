#' Simulate follow-up preparations for the lipid screening example
#'
#' Returns the underlying responses for recipes from the simulated 23-run
#' lipid screen. Supply preparation IDs to draw new observations as well.
#' Call this function after fixing the follow-up worklist. Only recorded
#' observations belong in a refitted model.
#'
#' @param recipes Data frame containing `PEG`, `Helper`, `Ionizable`,
#'   `Cholesterol` as proportions summing to one, `Ionizable_Lipid_Type`
#'   as an H101, H102 or H103 label, and numeric `N_P_ratio` and
#'   `flow_rate`. Additional columns are ignored. Use the experiment's
#'   declared feasible region when choosing recipes.
#' @param run_id Optional character vector of unique preparation IDs, one per
#'   row. A fresh preparation of the same recipe needs a distinct ID.
#' @param seed Nonnegative integer identifying a set of simulated observations.
#'   For each endpoint the draw is keyed to the seed and preparation ID. Thus
#'   reordering or subsetting a worklist retains its measurements. Generating
#'   observations requires the suggested package `digest` and preserves
#'   the caller's RNG state.
#' @return A data frame in input row order with `Potency_true`, `Size_true`
#'   and `PDI_true`, the latent pre-recording means. With `run_id` supplied,
#'   it also contains `Potency`, `Size` and `PDI` for refitting, plus
#'   `PDI_draw`, the continuous draw before recording.
#' @details Potency and Size use the generating surfaces for the simulated
#'   JoVE example. Fresh observations have independent Normal errors with
#'   standard deviations 2 and 3 respectively and are rounded to whole units,
#'   with half units rounded away from zero. The constant correction inside
#'   the Size generating mean belongs to its archived surface. The archived
#'   observation formula adds Normal noise on the original Size scale.
#'   Operator does not change these
#'   generating surfaces. For PDI the original noise family is unknown.
#'   This demonstration uses an independent Beta(4, 16) draw at every recipe
#'   and operator (mean 0.20, SD 0.0873), recorded to the nearest hundredth.
#'   These choices describe a teaching dataset.
#' @references Karl AT, Essex S, Wisnowski J, Rushing H (2023).
#'   A workflow for lipid nanoparticle (LNP) formulation optimization using
#'   designed mixture-process experiments and Self-Validated Ensemble Models
#'   (SVEM). \doi{10.3791/65200}.
#' @seealso [lipid_screen]
#' @examples
#' data(lipid_screen)
#' lipid_screen_simulate(head(lipid_screen, 2))
#' if (requireNamespace("digest", quietly = TRUE)) {
#'   lipid_screen_simulate(head(lipid_screen, 2), c("prep-24", "prep-25"), seed = 72)
#' }
#' @export
lipid_screen_simulate <- function(recipes, run_id = NULL, seed = 1L) {
  mixture <- c("PEG", "Helper", "Ionizable", "Cholesterol")
  numeric_inputs <- c(mixture, "N_P_ratio", "flow_rate")
  required <- c(numeric_inputs, "Ionizable_Lipid_Type")
  if (!is.data.frame(recipes) || !all(required %in% names(recipes)) ||
      anyDuplicated(names(recipes)[names(recipes) %in% required])) {
    stop("Supply each named lipid-screen recipe column exactly once.", call. = FALSE)
  }
  for (nm in numeric_inputs) {
    x <- recipes[[nm]]
    if (!is.numeric(x) || !is.null(dim(x)) || any(!is.finite(x))) {
      stop("Recipe column '", nm, "' must contain finite numeric values.", call. = FALSE)
    }
  }
  amounts <- as.matrix(recipes[mixture])
  if (any(amounts < 0 | amounts > 1) || any(abs(rowSums(amounts) - 1) > 1e-7)) {
    stop("Mixture inputs must be proportions in [0, 1] summing to one.", call. = FALSE)
  }
  lipid <- as.character(recipes$Ionizable_Lipid_Type)
  types <- c("H101", "H102", "H103")
  lipid_index <- match(lipid, types)
  if (anyNA(lipid_index)) {
    stop("Use an observed lipid label: H101, H102 or H103.", call. = FALSE)
  }

  # Pseudocomponents subtract the lower-bound composition, then scale by
  # the remaining 0.69 mixture fraction.
  ion <- (recipes$Ionizable - 0.10) / 0.69
  helper <- (recipes$Helper - 0.10) / 0.69
  chol <- (recipes$Cholesterol - 0.10) / 0.69
  peg <- (recipes$PEG - 0.01) / 0.69
  np <- recipes$N_P_ratio / 5 - 2
  blend <- 7.6 * ion + 8.15 * helper + 4.43 * chol + 47.73 * peg -
    5.15 * ion * helper + 8.82 * ion * chol + 8.23 * helper * chol -
    93.53 * ion * peg - 100.26 * helper * peg +
    275.55 * ion * helper * peg + 1000 * peg * helper * chol
  lipid_effect <- c(-1, 0.1, -4)[lipid_index] +
    ion * c(0, 0, -3)[lipid_index] + np * c(-1, 0, 1)[lipid_index]
  potency <- (blend + lipid_effect + 0.2 * np + 38) * 10 / 6
  size <- exp(4.5 + (-0.29068009689168 - 0.32) * (recipes$Helper - 0.4) / 0.2 +
                0.0278157917977491^2 / 2)
  result <- data.frame(Potency_true = potency, Size_true = size,
                       PDI_true = rep(0.2, nrow(recipes)))
  if (is.null(run_id)) return(result)
  if (!is.character(run_id) || length(run_id) != nrow(recipes) ||
      anyNA(run_id) || any(!nzchar(run_id)) || anyDuplicated(run_id)) {
    stop("Supply one unique, nonempty character run_id per preparation.", call. = FALSE)
  }
  seed <- .svem_integer_scalar(seed, "seed")
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Install the 'digest' package to simulate new observations.", call. = FALSE)
  }
  # A keyed uniform per preparation and endpoint avoids a dependence on row
  # order, the user's RNG kind, or whether other recipes have already run.
  uniform <- function(endpoint) {
    vapply(run_id, function(id) {
      key <- paste("lipid-screen-v1", seed, endpoint, enc2utf8(id), sep = ":")
      h <- digest::digest(key, algo = "sha256", serialize = FALSE)
      hi <- strtoi(substr(h, 1, 6), 16L)
      lo <- strtoi(substr(h, 7, 12), 16L)
      (hi * 2^24 + lo + 0.5) / 2^48
    }, numeric(1), USE.NAMES = FALSE)
  }
  result$Potency <- .lipid_screen_round(potency + stats::qnorm(uniform("Potency"), sd = 2))
  result$Size <- .lipid_screen_round(size + stats::qnorm(uniform("Size"), sd = 3))
  result$PDI_draw <- stats::qbeta(uniform("PDI"), 4, 16)
  result$PDI <- .lipid_screen_round(100 * result$PDI_draw) / 100
  result
}

.lipid_screen_round <- function(x) sign(x) * floor(abs(x) + 0.5)
