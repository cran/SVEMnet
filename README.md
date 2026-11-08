# SVEMnet: formulation and process optimization

SVEMnet fits self-validated ensemble models and helps turn measured experiments into a feasible follow-up worklist. Use the same factor specification for process-only, mixture, and mixture-process problems, including ingredient categories and blocking factors. Declare scientifically meaningful response goals and use a pool containing only preparations the laboratory can make.

The **default formulation workflow** is `svem_plan_runs()`. Its five priority roles are the best observed reference, elastic-net SVEM, GauPro, ranger Extra Trees, and forward-likelihood SVEM. Shared optima consume one preparation; remaining slots are diverse uncertainty medoids. Inspect the printed tie counts: shared nominees chosen by the same tie key do not demonstrate independent model agreement. The reference is a deliberate repeat, included in the budget. A budget of one with that reference enabled repeats the incumbent. An optional fixed original-design benchmark can be carried alongside the changing incumbent every round. It is mandatory when requested, consumes one slot, and shares a repeat if both refer to the same physical recipe; otherwise the budget must cover both.

```r
library(SVEMnet)
# measured: one row per independent preparation, on the intended response scale
spec <- bigexp_terms(response ~ factor1 + factor2 + ingredient,
                     data = measured, blocking = "day")
goals <- list(response = list(goal = "max", weight = 1,
                              lower_acceptable = 0, upper_acceptable = 100))
# feasible_pool uses allowed settings and one known fitted day as reference.
# Record the actual future execution day separately until it has measurements.
fixed_benchmark <- NULL
# Optional: fixed_benchmark <- list(preparation_id = "original-control")
plan <- svem_plan_runs(spec, measured, goals, feasible_pool, budget = 8L,
                      reference = list(day = "reference_day"),
                      benchmark = fixed_benchmark)
plan$worklist
plan$nominations
plot(plan)
svem_export_candidates_csv(plan, file = "follow-up.csv",
                          planned_date = as.Date("2026-11-02"),
                          preparation_prefix = "studyA-20261102-round1")
```

The optional export labels above distinguish physical `recipe_id` settings from a new `preparation_id` for each execution, including a repeat. Carry supplied stable recipe IDs and historical settings into later pools; automatically generated IDs can change when the pool changes, so otherwise match physical settings. Use a distinct prefix for each batch, even on the same date; the exporter has no global ID registry. `planned_date` is scheduling metadata. Record the actual execution date separately, and declare it as a model block only when scientifically intended. Omitting both labels retains the existing export.

This is a syntax template; replace the factors, measured data, scientific anchors, feasible pool and block levels. The fully runnable example is `vignette("lipid-screen-follow-up", package = "SVEMnet")`. It includes mixture bounds, three outcomes, the core-only five-new-recipe workflow, the reference-inclusive default portfolio, simulated confirmation and a refit.

The planner does not add a round or batch factor automatically. Declare intended blocking factors explicitly from physical conditions before seeing the outcomes. When a new measured day or operator level arrives, rebuild the specification from the augmented data and refit. Primary SVEM uses additive penalized block effects without recipe-by-block interactions; companion models can use interactions. A repeated fixed benchmark helps connect contexts, and exploration broadens recipe coverage, but neither guarantees separation of block and recipe effects. Compare new recipes with repeats prepared in the same batch: an additive batch shift cancels from that contrast. A fixed benchmark chosen in advance avoids selecting the reference for a flattering observed outcome, but repeat measurements are still noisy. Choose any challenge preparation to revisit a poor response explicitly.

The existing `bigexp_terms()` default offers a cubic candidate model, including pure cubics and three-way interactions, for regularized/model-selected fits. Explicitly chosen smaller or larger models remain available. A design may not independently support every offered term. Ensemble spread is an exploration diagnostic, not automatically a calibrated prediction interval.

Install suggested packages `GauPro` and `ranger` for the full portfolio. Keep the scientific inputs and model controls separate: `fit_control`, `forward_control`, `score_control`, `priority` and `exploration_fraction` can be customized, and `fits`, `forward`, `scores`, `incumbent`, `worklist` and `shortfall` remain inspectable. When changing only scientific weights or allocation, reuse `fits` and `forward`. Pass `spline = plan$spline` as well if the earlier plan fitted the spline. The planner verifies their complete training rows. With both primary and forward lists supplied, changing counts alone does not request a new spline fit. Refit on augmented or corrected measurements for the next round. The default `incumbent_summary = "best"` keeps the best observed individual preparation in the reference context, even after a poorer repeat. Use `incumbent_summary = "mean"` explicitly to score recipe means across repeated independent preparations within that same context. WMT and percent-in-spec diagnostics are optional. Existing `SVEMnet()` and selection/export calls remain supported.

For the initial experiment, [A Sphericity Filter for Space-Filling Designs in Nondeterministic Experiments](https://arxiv.org/abs/2610.36012) is a useful starting point. Its low-order design-protection model is separate from the richer analysis basis. This package does not yet implement that initial-design filter.

The [SVEMnet paper](https://arxiv.org/abs/2511.20968) and [formulation workflow paper](https://arxiv.org/abs/2212.11264) provide the methodological background. Their author-maintained arXiv versions may be updated; the journal versions remain static. The portfolio is a custom confirmation workflow, not a full Bayesian optimizer. Its default ordering is practical guidance, not a guarantee that one allocation is best for every application.

## How many preparations?

The counts for new recipes exclude reference repeats. Each model offers at most one nominee. Shared recipes occupy one place, so inspect the attained counts before booking laboratory time.

| Request | Controls | Planned total with distinct nominees |
|---|---|---|
| Two total preparations | `budget = 2` | Best-observed repeat and one new model candidate |
| Two new candidates plus the repeat | `n_optima = 2, n_explore = 0` | Three preparations |
| Three candidates, two explorations, no repeat | `n_optima = 3, n_explore = 2, incumbent = FALSE` | Five preparations |
| Ten preparations with automatic exploration | `budget = 10` | Repeat, available model candidates, then exploration |
| Fixed budget with two exploration runs | `budget = 10, n_explore = 2` | At most eight, or seven without the spline, including the repeat and exactly two exploration runs. Other places remain unused |

A supplied `benchmark` adds a required repeat unless it is the incumbent's recipe. `n_optima = 10` exceeds the roster of four ordinary models and the optional spline. Use a total budget or add exploration for more runs. `allocation_counts$unused_capacity` reports unfilled physical places. With exact counts, `shortfall` reports unmet requested counts and can be zero while budget places remain unused.

## Requesting a spline nominee and reserving exploration

For a supported Gaussian experiment, the optional Gaussian spline can nominate an additional recipe from the same legal pool. It supports general numeric and nominal study factors, mixtures with multiple closure groups, fixed components and additive blocks. The original low-dimensional basis is retained; wider numeric problems use a separate saved general surface with explicit complexity budgets.

Choose the response-free physical reference grid before looking at outcomes. Its coordinates determine the spline basis. Separately, the physical factor bounds determine which recipes are legal. A planner-created spline uses `candidates` as its reference grid. Keep that grid fixed across rounds for comparable fits. To score a smaller candidate list with the original basis, fit using `svem_fit_responses(..., method = "spline", domain_grid = reference_grid)` and pass the fitted list as `spline`. Predicting off its saved grid requires the suggested `mgcv` package, including after reload.

Specify a total physical budget when lab capacity is fixed. Alternatively, specify both new-recipe counts and let the planner add the distinct reference repeats. For example, in a two-factor process study with both optional companions available, the following request totals eight physical preparations, including the incumbent repeat. Request five distinct **new model nominations** and reserve two exploration settings. All eight places are used if the models nominate five distinct new recipes and both exploration settings are available. In the worked example, the five models nominate three distinct new recipes. The repeat and two exploration settings bring the worklist to six, leaving two places unused.

```r
plan <- svem_plan_runs(process_spec, measured, goals, candidates = feasible,
                       incumbent = TRUE, n_optima = 5L, n_explore = 2L,
                       spline = TRUE)
plan$allocation_counts
plan$merged_nominations
plan$nominations
plan$worklist[, c("slot", "role", "nominated_by", "recipe_id")]
```

The ordinary five-role workflow remains the default. On a new fit with explicit counts, the standard priority also tries the optional spline after forward SVEM. Set `spline = FALSE` to omit it. A duplicate endorsement uses one physical preparation and the planner can try the next model for the open model place. Inspect the returned nominations and allocation when fewer distinct settings are available. A separate fixed benchmark consumes a reference place unless it matches the incumbent. Supply `budget` when only one count is given. With both counts, the planner infers the total after merging required reference repeats. An explicit `n_explore` reserves that many exploration settings. With only exploration specified, models can fill the other places and an unfilled model place is unused capacity. With `n_explore = NULL`, exploration can fill the remaining budget. `exploration_fraction` controls the size of the high-spread shortlist used by the medoid calculation.

The executable follow-up vignette includes a two-factor teaching surface and an optional general spline on the full 23-run lipid experiment: four mixture components, two numeric process settings, nominal ingredient type and additive operator blocks. Use `spline = TRUE` or request `n_optima = 5L` to add it after forward SVEM. The saved general surface supports multiple mixture groups and fixed components, uses Duchon bases plus modest one-factor curvature in higher numeric dimensions, and reports any basis-budget reductions. Member diagnostics expose effective degrees of freedom by basis family alongside the validation-size cap. These overlapping basis components are not a unique scientific decomposition. A purely categorical problem uses the penalized categorical component. General factor support is separate from the original low-dimensional simulation evidence and does not establish a universally better optimizer.
