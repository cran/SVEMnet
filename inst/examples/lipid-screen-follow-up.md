---
title: "Choosing and confirming follow-up recipes in a 23-run lipid screen"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Choosing and confirming follow-up recipes in a 23-run lipid screen}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



Start with the 23 measured preparations in `lipid_screen`. The first example
uses only the core package to export **five distinct new preparations**. It
keeps the three response goals, mixture bounds, planned operator and run
budget visible. The advanced example then adds other models, reserves an
incumbent repeat within an eight-preparation budget, simulates measurements
and refits. A separate two-factor process example at the end shows how to
request a spline nominee and reserve distinct exploration preparations.

For another experiment, change the data, factor specification, goals and
feasible candidate settings below. Keep one row per independent preparation;
aggregate technical assay replicates before fitting. A repeated recipe in a
new preparation is a new row with its own preparation ID.

For planning the initial mixture or process experiment, see Karl's
[A Sphericity Filter for Space-Filling Designs in Nondeterministic Experiments](https://arxiv.org/abs/2610.36012).
Its low-order working model protects selected estimation precision while
retaining space-filling coverage. That design choice is separate from the
richer cubic candidate model fitted below; it does not require the final
response surface to be low-order.

## Specify the experiment once


``` r
library(SVEMnet)
data(lipid_screen, package = "SVEMnet")
dat <- lipid_screen
components <- c("PEG", "Helper", "Ionizable", "Cholesterol")
recipe_cols <- c(components, "Ionizable_Lipid_Type", "N_P_ratio", "flow_rate")
predictors <- c(recipe_cols, "Operator")
mix <- list(list(vars = components, lower = c(.01, .10, .10, .10),
                 upper = c(.05, .60, .60, .60), total = 1))
spec <- bigexp_terms(
  Potency ~ PEG + Helper + Ionizable + Cholesterol +
    Ionizable_Lipid_Type + N_P_ratio + flow_rate, data = dat,
  blocking = "Operator", mixture_groups = mix,
  discrete_numeric = c("N_P_ratio", "flow_rate"), report = FALSE
)
goals <- list(
  Potency = list(goal = "max", weight = .6, lower_acceptable = 65, upper_acceptable = 95),
  Size    = list(goal = "min", weight = .3, lower_acceptable = 40, upper_acceptable = 210),
  PDI     = list(goal = "min", weight = .1, lower_acceptable = .04, upper_acceptable = .39)
)
endpoints <- names(goals)

```

The default `bigexp_terms()` expansion is a full cubic candidate model: pure
squares and cubes, two- and three-factor interactions, and partial-cubic
terms such as `x^2:z`. Blocking variables stay additive. The fitting method
reduces this common candidate model; no preliminary choice of a smaller
model is required.

## Five new preparations with the core package

`svem_fit_responses()` returns an ordinary named list, so `core_fits$Potency`
remains a regular fitted model. It makes the same sequential calls and random
draws as fitting the three endpoints individually. With only the response
already named in `spec`, use `SVEMnet(spec, dat)` directly.


``` r
set.seed(2301)
core_fits <- svem_fit_responses(spec, dat, responses = names(goals))

# Change this to the operator who will prepare the next batch.
planned_operator <- "A"
set.seed(2302)
core_pool <- svem_random_table_multi(core_fits, n = 600, mixture_groups = mix)$data
# Include old recipes so exclusion is also exercised in this example.
core_pool <- rbind(dat[predictors], core_pool[predictors])
core_pool$Operator <- factor(planned_operator, levels = levels(factor(dat$Operator)))
core_pool <- unique(core_pool)
core_pool$recipe_id <- sprintf("C%04d", seq_len(nrow(core_pool)))

scored_core <- svem_score_random(
  core_fits, goals = goals, data = dat, candidates = core_pool,
  mixture_groups = mix, budget = 5L, verbose = FALSE
)

# Set this to an explicit file in your project for an actual experiment.
csv_file <- file.path(tempdir(), "lipid-screen-five-new.csv")
svem_export_candidates_csv(scored_core, file = csv_file, overwrite = TRUE,
  planned_date = as.Date("2026-11-02"),
  preparation_prefix = "lipid-screen-20261102-core-example")
new_runs <- read.csv(csv_file, stringsAsFactors = FALSE)
```

The candidate pool respects the declared mixture bounds and total. The
sampler draws proportions continuously; replace this pool with the lab's
allowable settings when preparation increments matter. The advanced example
uses an explicit preparation grid. The operator is set before scoring so all
models predict the context in which the new preparations will be made.

Here `budget = 5L` counts distinct **new** preparations. Historical physical
recipes are excluded even if their recorded operator differs. Inspect
`shortfall` before booking laboratory runs. The core path can leave slots empty
because of its nomination or top-10% shortlist limits, or unavailable
informative spread, even when the pool contains five or more feasible recipes.
For a small pool, use `svem_plan_runs(incumbent = FALSE)` to apply the portfolio
allocator, which widens its shortlist when needed to meet the budget. The exporter writes only the selected
worklist in slot order and retains its role and every endorsing engine. The
explicit optional date and prefix above add scheduling labels. `recipe_id`
identifies physical settings within this pool; `preparation_id` identifies one
new execution, including a repeated recipe. Retain supplied stable recipe IDs
and historical settings in later pools; automatically generated IDs may change
when the pool changes, so otherwise match physical settings. Choose a distinct batch prefix even for two batches
on the same date. Omitting both options leaves the export unchanged.

`planned_date` is not the actual execution date and does not add a model block.
Record the actual date when preparing the runs. The historical data use
date-like `RunID` labels such as `02SEP22-1`; retain those original identifiers.
The exporter does not maintain a global laboratory ID registry.


``` r
stopifnot(scored_core$shortfall == 0L, nrow(new_runs) == 5L,
          identical(new_runs$slot, 1:5), !anyDuplicated(new_runs$recipe_id),
          !anyDuplicated(new_runs$preparation_id),
          all(new_runs$planned_date == "2026-11-02"),
          !anyDuplicated(new_runs[recipe_cols]),
          identical(new_runs$recipe_id, scored_core$worklist$recipe_id),
          identical(new_runs$role, scored_core$worklist$role),
          identical(new_runs$nominated_by, scored_core$worklist$nominated_by),
          all(new_runs$Operator == planned_operator),
          all(abs(rowSums(new_runs[components]) - 1) < 1e-10),
          nrow(merge(new_runs[recipe_cols], dat[recipe_cols], by = recipe_cols)) == 0L)
for (i in seq_along(components)) {
  stopifnot(all(new_runs[[components[i]]] >= mix[[1]]$lower[i]),
            all(new_runs[[components[i]]] <= mix[[1]]$upper[i]))
}
knitr::kable(new_runs[c("slot", "recipe_id", "preparation_id", "planned_date")],
             row.names = FALSE)
```



| slot|recipe_id |preparation_id                         |planned_date |
|----:|:---------|:--------------------------------------|:------------|
|    1|C0294     |lipid-screen-20261102-core-example-001 |2026-11-02   |
|    2|C0523     |lipid-screen-20261102-core-example-002 |2026-11-02   |
|    3|C0575     |lipid-screen-20261102-core-example-003 |2026-11-02   |
|    4|C0553     |lipid-screen-20261102-core-example-004 |2026-11-02   |
|    5|C0272     |lipid-screen-20261102-core-example-005 |2026-11-02   |



``` r
knitr::kable(new_runs[c("slot", "role", "nominated_by", "recipe_id", predictors)],
             digits = 4, row.names = FALSE)
```



| slot|role               |nominated_by |recipe_id |    PEG| Helper| Ionizable| Cholesterol|Ionizable_Lipid_Type | N_P_ratio| flow_rate|Operator |
|----:|:------------------|:------------|:---------|------:|------:|---------:|-----------:|:--------------------|---------:|---------:|:--------|
|    1|optimum            |svemnet      |C0294     | 0.0469| 0.5603|    0.1026|      0.2902|H103                 |        14|       1.5|A        |
|    2|maximum_spread     |svemnet      |C0523     | 0.0105| 0.5670|    0.1840|      0.2386|H103                 |         6|       3.0|A        |
|    3|high_spread_medoid |svemnet      |C0575     | 0.0305| 0.5636|    0.1211|      0.2848|H103                 |        12|       2.0|A        |
|    4|high_spread_medoid |svemnet      |C0553     | 0.0345| 0.1546|    0.5209|      0.2900|H102                 |        12|       2.0|A        |
|    5|high_spread_medoid |svemnet      |C0272     | 0.0495| 0.2326|    0.5369|      0.1810|H102                 |        10|       2.5|A        |



For an existing workflow, keep using `svem_select_from_score_table()`. Its
`k` is the requested number of medoids **in addition to** its best row;
`k = 4` can therefore produce five rows. The budgeted worklist has a different
meaning: its budget is a total number of distinct new preparations.

## Advanced: several models and a carried-forward reference

The SVEMnet paper starts with the same data, candidate terms and scientific
weights. Its workflow fits each endpoint, scores feasible settings, chooses
diverse high-score and high-uncertainty recipes, and refits after measurement.
Here we also ask forward-selection SVEM, GauPro and Extra Trees to nominate
recipes. The fixed publication did not benchmark these additions. The journal
articles remain the versions of record; the author occasionally updates the
[SVEMnet arXiv manuscript](https://arxiv.org/abs/2511.20968) and the earlier
[workflow arXiv manuscript](https://arxiv.org/abs/2212.11264).

## What the Gaussian process contributes

This package implements a custom confirmation workflow. GauPro proposes a
recipe using its fitted mean (or the desirability of several fitted means).
A separate exploration slot can use the largest latent standard deviation
for a chosen response. The scientist allocates preparations, measures them
and decides when to refit. A full Bayesian optimizer specifies an acquisition
policy for those sequential decisions. For example, expected improvement
uses both the predicted improvement and uncertainty to choose a point.
[Jones, Schonlau and Welch (1998)](https://doi.org/10.1023/A:1008306431147)
and [Frazier (2018)](https://arxiv.org/abs/1807.02811).

SVEMnet offers the full candidate expansion and shrinks the additive operator
effect. The two optional companions use the original factors and free mixture
coordinates. Here we deliberately let those models learn operator-by-recipe
patterns. Both operators handled every lipid type, but some combinations
have little support. We specify operator A for every new preparation.


``` r
limits <- list(Potency = list(lower = 78), Size = list(upper = 100),
               PDI = list(upper = .25))
with(dat, table(Operator, Ionizable_Lipid_Type))
#>         Ionizable_Lipid_Type
#> Operator H101 H102 H103
#>        A    2    5    5
#>        B    7    2    2
# Installing GauPro and ranger enables both companions. Either works alone.
packages <- c(gaupro = "GauPro", extra_trees = "ranger")
enabled <- names(packages)[vapply(packages, requireNamespace,
                                  quietly = TRUE, FUN.VALUE = logical(1))]
used <- c("SVEMnet", unname(packages[enabled]))
knitr::kable(data.frame(package = used, version = vapply(used, function(p)
  as.character(packageVersion(p)), character(1))), row.names = FALSE)
```



|package |version |
|:-------|:-------|
|SVEMnet |4.0.0   |
|GauPro  |0.2.17  |
|ranger  |0.18.0  |



With both packages installed, `companions = TRUE` is the shortcut for this
request. `predict(fits$Potency, newdata)` returns the ordinary SVEM response.
Use `engine = "all"` to inspect the models side by side. We use the fitting
defaults throughout and ordinary SVEM predictions with `debias = FALSE`.

## Prepare the candidate table before optimizing

The proportions below are built in integer basis points. The fourth component
is determined by the total. The seed recipes are included in the candidate
table so that an endorsement of a previous preparation can be recognized.
The paper draws 25,000 feasible settings and derives default desirability
anchors from predicted response quantiles. Here a smaller discrete pool
and fixed anchors make the models and both rounds comparable.
The resulting recipe IDs and scores are specific to this example. In practice, supply
the full set of recipes the lab would prepare or draw a larger feasible pool.


``` r
composition <- expand.grid(PEG = seq(100, 500, 100),
                           Helper = seq(1000, 6000, 500),
                           Ionizable = seq(1000, 6000, 500))
composition$Cholesterol <- 10000 - rowSums(composition)
composition <- subset(composition, Cholesterol >= 1000 & Cholesterol <= 6000)
settings <- expand.grid(Ionizable_Lipid_Type = levels(factor(dat$Ionizable_Lipid_Type)),
                        N_P_ratio = sort(unique(dat$N_P_ratio)),
                        flow_rate = sort(unique(dat$flow_rate)), stringsAsFactors = FALSE)
grid <- merge(composition / 10000, settings, by = NULL)
set.seed(2302)
pool <- unique(rbind(dat[recipe_cols], grid[sample.int(nrow(grid), 600), recipe_cols]))
pool$Operator <- factor("A", levels = levels(factor(dat$Operator)))
pool$recipe_id <- sprintf("R%04d", seq_len(nrow(pool)))
increments <- setNames(rep(.0001, length(components)), components)
stopifnot(all(abs(rowSums(pool[components]) - 1) < 1e-10))

# Map each preparation to a physical recipe ID, regardless of operator.
# Quantize only the mixture columns and reject off-grid inputs.
physical_settings <- function(rows) {
  key <- rows[recipe_cols]
  for (nm in names(increments)) {
    units <- key[[nm]] / increments[[nm]]
    stopifnot(all(abs(units - round(units)) <= 1e-7))
    key[[nm]] <- round(units)
  }
  key
}
lookup <- cbind(recipe_id = pool$recipe_id, physical_settings(pool))
stopifnot(!anyDuplicated(lookup[recipe_cols]))
screen_rows <- cbind(screen_row = seq_len(nrow(dat)), physical_settings(dat))
matched <- merge(screen_rows, lookup, by = recipe_cols, all.x = TRUE)
dat$recipe_id <- matched$recipe_id[match(screen_rows$screen_row, matched$screen_row)]
stopifnot(nrow(matched) == nrow(dat), !anyNA(dat$recipe_id))
```

## Allocate eight preparations, including an incumbent repeat

The convenience function earns its place when several models and a reference
are needed: it fits the endpoints, collects model nominees and delegates the
allocation to the scorer. All intermediates remain inspectable. `plan$fits`
and `plan$forward` are ordinary model lists; `plan$scores` contains the full
candidate table, and `plan$worklist` and `plan$shortfall` describe the actual
preparation plan. It does not write a file.

Unlike the earlier new-recipe-only scorer call, `svem_plan_runs()` counts the
incumbent repeat **within** its total budget. Its default five slots can hold
one reference repeat plus up to four distinct model optima. This example
uses eight total preparations, leaving room for exploration. Optionally retain
a fixed original-design benchmark as well: set `fixed_benchmark` to its unique
preparation ID below. That benchmark is mandatory and counts within the total,
reducing the slots available for model optima or exploration. It shares one
repeat if its physical recipe is also the incumbent; otherwise both repeats
must fit in the budget. Keep the same ID selector on subsequent refits so
adding or reordering rows does not change which benchmark is carried forward.


``` r
fit_controls <- list(companions = if (length(enabled)) enabled else FALSE,
                     companion_blocks = "interacting")
score_controls <- list(increments = increments, specs = limits,
                       gp_endpoint = "Potency", verbose = FALSE)
# Choose the intended incumbent by stable preparation ID, not row position.
incumbent_row <- match("02SEP22-20", dat$RunID)
stopifnot(!is.na(incumbent_row))
fixed_benchmark <- NULL
# Optional fixed control from the original design, retained in every round:
# fixed_benchmark <- list(RunID = "02SEP22-2")
set.seed(2301)
plan <- svem_plan_runs(
  spec, data = dat, goals = goals, candidates = pool, budget = 8L,
  fit_control = fit_controls, score_control = score_controls,
  incumbent = incumbent_row, benchmark = fixed_benchmark
)
fits <- plan$fits
forward <- plan$forward
scored <- plan$scores
batch <- plan$worklist
stopifnot(plan$shortfall == 0L, nrow(batch) == 8L,
          !anyDuplicated(batch$recipe_id), sum(!batch$new_recipe) %in% 1:2)
knitr::kable(batch[c("slot", "role", "nominated_by", "recipe_id", "new_recipe", recipe_cols)],
             digits = 4, row.names = FALSE)
```



| slot|role               |nominated_by                    |recipe_id |new_recipe |    PEG| Helper| Ionizable| Cholesterol|Ionizable_Lipid_Type | N_P_ratio| flow_rate|
|----:|:------------------|:-------------------------------|:---------|:----------|------:|------:|---------:|-----------:|:--------------------|---------:|---------:|
|    1|incumbent_repeat   |observed                        |R0019     |FALSE      | 0.0474| 0.4297|    0.1396|      0.3833|H101                 |         6|       2.5|
|    2|optimum            |svemnet + extra_trees + forward |R0432     |TRUE       | 0.0500| 0.5000|    0.1000|      0.3500|H102                 |        14|       2.0|
|    3|optimum            |gaupro                          |R0259     |TRUE       | 0.0400| 0.4500|    0.1000|      0.4100|H102                 |        14|       3.0|
|    4|high_spread_medoid |svemnet                         |R0171     |TRUE       | 0.0300| 0.1500|    0.5500|      0.2700|H102                 |        12|       2.0|
|    5|high_spread_medoid |svemnet                         |R0600     |TRUE       | 0.0100| 0.5000|    0.1000|      0.3900|H103                 |         8|       2.0|
|    6|high_spread_medoid |svemnet                         |R0444     |TRUE       | 0.0400| 0.6000|    0.1500|      0.2100|H103                 |        12|       1.0|
|    7|high_spread_medoid |svemnet                         |R0384     |TRUE       | 0.0500| 0.4000|    0.1000|      0.4500|H103                 |        12|       1.5|
|    8|high_spread_medoid |svemnet                         |R0466     |TRUE       | 0.0500| 0.2500|    0.5500|      0.1500|H102                 |         6|       2.0|



The default order is the incumbent, elastic-net SVEM, GauPro, Extra Trees,
then forward likelihood. A smaller budget truncates that order; with one
slot and an enabled incumbent, the plan is confirmation only. Set
`incumbent = FALSE` when every slot must be a new recipe. The priority is a
practical allocation rule, not a claim that its ordering is best for every
experiment. For this multi-response decision each model maximizes the same
weighted desirability of its fitted endpoint means. Desirability stops improving
at its declared limits. When several multi-response recipes reach the same
combined desirability, they remain exact utility ties and the stable recipe
key breaks the tie. Extra predicted improvement beyond a saturated goal is
not assigned additional utility. Change the scientific goals explicitly if
those differences should matter.

Models endorsing the same new recipe spend one slot and retain all their
attributions. Inspect `best_ties` and `new_ties` in the nominations. When exact
ties are resolved by the same recipe key, a shared nomination is not evidence
of independent agreement among the fitted models. Remaining slots are
distinct high-spread medoids. An
unavailable engine is visible in the diagnostics; it is not silently replaced
by another fitted model. The code above uses whichever suggested companion
packages are installed.

For this first batch we explicitly carry forward preparation `02SEP22-20`.
Its original operator was B; the planned repeat uses A. This is a declared
incumbent choice across contexts, not an automatic comparison of raw A and B
outcomes. Matching `RunID` keeps the choice stable if the table is reordered.
Original measurements and source context remain visible for review. Selecting
a noisy winner does not establish that it is the best true recipe; the repeat
supplies confirmation evidence. Without an explicit row, the planner chooses
the best observed preparation within the requested reference context. A separate fixed benchmark keeps its
original selector in `plan$policy$benchmark` and its source row and context in
`plan$benchmark`, regardless of which preparation becomes the incumbent.


``` r
reference_row <- plan$incumbent$source_data
reference_row$observed_score <- plan$incumbent$observed_score
knitr::kable(reference_row[c("recipe_id", recipe_cols, "Operator", endpoints,
                             "observed_score")], digits = 4, row.names = FALSE)
```



|recipe_id |    PEG| Helper| Ionizable| Cholesterol|Ionizable_Lipid_Type | N_P_ratio| flow_rate|Operator | Potency| Size| PDI| observed_score|
|:---------|------:|------:|---------:|-----------:|:--------------------|---------:|---------:|:--------|-------:|----:|---:|--------------:|
|R0019     | 0.0474| 0.4297|    0.1396|      0.3833|H101                 |         6|       2.5|B        |      95|   85| 0.3|         0.7961|




``` r
knitr::kable(plan$nominations, digits = 3, row.names = FALSE)
```



|engine      |status   |best_recipe_id |best_prepared | best_ties|new_recipe_id | new_ties|endorses_incumbent | slot|fit_status |
|:-----------|:--------|:--------------|:-------------|---------:|:-------------|--------:|:------------------|----:|:----------|
|incumbent   |selected |R0019          |TRUE          |         1|NA            |        0|FALSE              |    1|observed   |
|svemnet     |selected |R0432          |FALSE         |         1|R0432         |        1|FALSE              |    2|ok         |
|gaupro      |selected |R0259          |FALSE         |         1|R0259         |        1|FALSE              |    3|ok         |
|extra_trees |selected |R0432          |FALSE         |         1|R0432         |        1|FALSE              |    2|ok         |
|forward     |selected |R0432          |FALSE         |         1|R0432         |        1|FALSE              |    2|ok         |



**Exploration choices.** The portfolio uses diverse medoids from the top
tenth of the eligible SVEM spread scores for its remaining slots. SVEM
summarizes widths of member-percentile intervals across the three responses.
Widths are scaled between their 2nd and 98th percentiles and clipped, which
can produce exact ties. Ensemble spread is useful for exploration; it is not
automatically calibrated predictive uncertainty.

The diagnostics below also show alternatives: choose the largest combined
spread and break ties by fitted desirability, or choose the widest underlying
spread. They do not allocate extra preparations in this portfolio. The
separate scorer can also nominate the largest latent GP standard deviation
for a chosen response; no such slot is reserved by the default portfolio.
Extra Trees supplies a mean nominee, with trees varying through randomized
splits of the full training rows.


``` r
knitr::kable(scored$exploration, digits = 4, row.names = FALSE)
```



|rule               |recipe_id | tie_count| exact_ties| display_only_ties| min_desirability| max_desirability| uncertainty_measure| strict_width|  score|
|:------------------|:---------|---------:|----------:|-----------------:|----------------:|----------------:|-------------------:|------------:|------:|
|best_predicted_tie |R0250     |         1|          1|                 0|           0.0092|           0.0092|              0.8381|       1.1334| 0.0092|
|widest_spread      |R0250     |         1|          1|                 0|           0.0092|           0.0092|              0.8381|       1.1334| 0.0092|



``` r
knitr::kable(scored$tie_endpoint_ranges, digits = 4, row.names = FALSE)
```



|response | min_prediction| max_prediction| min_width| max_width| min_weighted_spread| max_weighted_spread|
|:--------|--------------:|--------------:|---------:|---------:|-------------------:|-------------------:|
|Potency  |        78.8114|        78.8114|   13.0967|   13.0967|              0.5135|              0.5135|
|Size     |       213.2066|       213.2066|   21.6754|   21.6754|              0.2245|              0.2245|
|PDI      |         0.2388|         0.2388|    0.3717|    0.3717|              0.1000|              0.1000|



``` r
tab <- scored$score_table
unprepared <- tab[!tab$recipe_id %in% dat$recipe_id, , drop = FALSE]
summary(data.frame(PDI_member_width = tab$PDI_upr - tab$PDI_lwr,
                   PDI_weighted_spread = tab$PDI_ciw_w))
#>  PDI_member_width  PDI_weighted_spread
#>  Min.   :0.06012   Min.   :0.00000    
#>  1st Qu.:0.09726   1st Qu.:0.01331    
#>  Median :0.12562   Median :0.02648    
#>  Mean   :0.13728   Mean   :0.03143    
#>  3rd Qu.:0.16015   3rd Qu.:0.04252    
#>  Max.   :0.40903   Max.   :0.10000
```

An exact tie is possible even when fitted responses and underlying widths
differ substantially. For example, consider a single response with widths
1 through 98, then 105 and 120. The 98th percentile falls below 105, so both
of the last two recipes receive a scaled width of 1. If their desirabilities
are 0.85 and 0.20, respectively, the practical tie-break chooses the former.
The raw-width rule chooses the latter. The tie counts printed above refer to
the actual multi-response pool.

For diverse representatives of the high-spread region, the established
selector can return feasible medoids. A scientist following the paper's
SVEM-only allocation can also select distinct representatives of the
high-score region. Both selections use physical recipe coordinates.


``` r
diverse <- svem_select_from_score_table(unprepared, target = "uncertainty_measure",
                                        top = .10, k = 3, predictor_cols = recipe_cols)
high_score <- svem_select_from_score_table(unprepared, target = "score",
                                           top = .05, k = 3, predictor_cols = recipe_cols)
```

The paper also nominates a recipe for which the modeled means appear likely
to meet specifications. The product of the three marginal probabilities is a
ranking proxy. It addresses uncertainty about the means; the measurements
from an individual future preparation also include assay and preparation
variation. Rank only unprepared recipes. Choosing this alternative would
use one of the available slots.


``` r
finite <- which(is.finite(unprepared$p_joint_mean))
stopifnot(length(finite) > 0L)
by_spec <- finite[which.max(unprepared$p_joint_mean[finite])]
knitr::kable(unprepared[by_spec, c("recipe_id", "score", "p_joint_mean")],
             digits = 4, row.names = FALSE)
```



|recipe_id |  score| p_joint_mean|
|:---------|------:|------------:|
|R0027     | 0.8858|        0.995|



Whole-model permutation tests provide another optional diagnostic. With WMT
enabled, compare `score` and `wmt_score` on the same recipes. A weak global
signal alone is insufficient to discard an endpoint. Here the scientist's
nominal weights govern every model and both rounds.


``` r
forms <- setNames(lapply(endpoints, function(y) bigexp_formula(spec, y)), endpoints)
wmt <- svem_wmt_multi(forms, data = dat, mixture_groups = mix,
                      wmt_control = list(seed = 123, nCore = 1L), plot = FALSE)
svem_specification_scores <- svem_score_random(
  fits, goals = goals, data = dat, candidates = pool,
  wmt = wmt, specs = limits)
head(svem_specification_scores$score_table[c("recipe_id", "score", "wmt_score",
                                             "p_joint_mean")])
```

## Change the budget without refitting

To request specific kinds of runs, set `n_optima = 3L, n_explore = 2L` in the planner call below and omit `budget`. This requests three new model candidates, two exploration recipes and the best-observed repeat, for six planned preparations. In this worked run, these fitted models nominate only two distinct new recipes, so the call returns five preparations and reports one unmet model request. Setting `incumbent = FALSE` removes the repeat from the planned total. Inspect selected counts and unused places before preparing the runs.

Forward likelihood and the relaxed elastic net can emphasize different terms
in the same candidate model. Forward's validation score evaluates each
training-fitted Gaussian mean and variance on the validation-weighted copy.
The wrapper fits this method for each endpoint as well as the primary SVEM
models. The same physical pool, scientific weights and planned operator
apply to every nominee.

For a smaller batch, reuse the fitted models and rescore. These five total
preparations include a reference repeat; they are separate from the first
example's five **new** preparations. We still simulate and refit only the
eight-preparation batch stored in `batch`.


``` r
five_total <- svem_plan_runs(
  spec, data = dat, goals = goals, candidates = pool, budget = 5L,
  fits = fits, forward = forward, score_control = score_controls,
  incumbent = incumbent_row, benchmark = fixed_benchmark
)
five_total
#> SVEMnet default formulation workflow
#> 5 of 5 preparations; 4 new recipes; 0 unfilled slots
#> Reference repeats: 1 | New model candidates: 2 | Exploration runs: 2
#>  slot recipe_id               role                    nominated_by
#>     1     R0019   incumbent_repeat                        observed
#>     2     R0432            optimum svemnet + extra_trees + forward
#>     3     R0259            optimum                          gaupro
#>     4     R0037 high_spread_medoid                         svemnet
#>     5     R0551 high_spread_medoid                         svemnet
#> Nomination status (including unavailable and budget-dropped roles):
#>       engine   status fit_status best_ties new_ties new_recipe_id slot
#>    incumbent selected   observed         1        0          <NA>    1
#>      svemnet selected         ok         1        1         R0432    2
#>       gaupro selected         ok         1        1         R0259    3
#>  extra_trees selected         ok         1        1         R0432    2
#>      forward selected         ok         1        1         R0432    2
#> Priority: incumbent > svemnet > gaupro > extra_trees > forward 
#> High-spread exploration fraction: 0.1 
#> Observed incumbent rule: best 
#> Scientific goals and weights (full settings retained in policy):
#> $Potency
#> $Potency$goal
#> [1] "max"
#> 
#> $Potency$weight
#> [1] 0.6
#> 
#> $Potency$lower_acceptable
#> [1] 65
#> 
#> $Potency$upper_acceptable
#> [1] 95
#> 
#> 
#> $Size
#> $Size$goal
#> [1] "min"
#> 
#> $Size$weight
#> [1] 0.3
#> 
#> $Size$lower_acceptable
#> [1] 40
#> 
#> $Size$upper_acceptable
#> [1] 210
#> 
#> 
#> $PDI
#> $PDI$goal
#> [1] "min"
#> 
#> $PDI$weight
#> [1] 0.1
#> 
#> $PDI$lower_acceptable
#> [1] 0.04
#> 
#> $PDI$upper_acceptable
#> [1] 0.39
```

If scientific priorities change, pass a revised `goals` list and the same
fits and forward fits to inspect their consequences. Keeping the weights
visible matters: agreement among fitted models does not settle the tradeoff
between potency, size and PDI. `exploration_fraction = .25` is an optional
wider shortlist for medoids, not a demonstrated improvement over .10.


``` r
plot(plan)
```

<img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAfgAAAEgCAMAAABvm5EPAAABcVBMVEUAAAAAADUAADoAAGYAOjoAOmYAOpAAVJUAZmYAZrYAujg1AAA6AAA6ADo6AGY6OgA6Ojo6OmY6OpA6ZmY6ZrY6kLY6kNtNTU1NTW5NTY5NbqtNjshUldFeqOthnP9lKABmAABmADpmAGZmOgBmOjpmOpBmZjpmZmZmZpBmkJBmtttmtv9tpsFuTU1uTY5ubqtujshuq6tuq+R2MDCEyeuOTU2OTY6ObquOjk2OjsiOq+SOyP+QOgCQOjqQOmaQZgCQZjqQkJCQ29uQ2/+RkZGVlZWampqbm5uhoaGioqKoXgCpqamrbk2rq6ur5P+0tLS2ZgC2Zjq2kGa2tma2tra225C2/7a2/9u2///BwcHExMTIjk3Ijm7Ijo7I///R0dHV1dXbkDrbkGbbtmbb25Db2//b/7bb/9vb///kq27kq47k/8jk///l5eXryYTr6+v4dm3/tmb/yI7/25D/29v/5Kv//7b//8j//9v//+T///8wALKEAAAACXBIWXMAAAsSAAALEgHS3X78AAAgAElEQVR4nO1dj58ct1VfYjh6xjaxG/ZIC+aSkNJy0EKgDRAnBa41F5btXOEIt5PGuiUhPUOOzfri87F/PfoxGv160kjzSzO7+tqf29mZN5rZ+Y6enp7ekyabhJ3EJPYNJMRBIn5HkYjfUSTidxR1iF8lyVYloyARH18yChLx8SWjIBEfXzIKEvHxJaMgER9fMgoS8fEloyARH18yChLx8SWjoB7xCPlK+pe5u5JRUIt4dHXlx/xYHn0i3guJ+HYloyCp+viSUZCMu/iSUZCIjy8ZBYn4+JJRkIiPLxkF20x87mmCxr7PKNhi4tGFZ6cz9i+KgkR8R1cfI/H/9y8H3y0+v/MVcHwkxCdV7wJE/Iv3Np9+wD6ffwAcHwvxY5GMAoj4//4F4VwmfpXQFnpl1wGI+F8VxG+eH7wHnZNcthuPJzAYimG4avzz98h/E9s7SLOcTPa9JC1PQH4dRkg8b+Oxnn/xXeD4gIk/nUzunNcu8/rR5e2TIx9J+Akoe0dIPLXqX73/1asfH7z5DDg+XFV/++FxkzKXU++rg09g7MRXYRhqGQKusk3KPP3BPU9Vb8HIVX0VhkS8WvPWvykzF0783uXm1EvV+5c5VIyceL2t/eZyszwCJX3KJKRXqPtE/DAkASNLMBdc5np/k2q8AwMiXlf1mLkGNZ505yrsu0T8MCVPZeYGfJ/xsTXEA/2rQd7nULAtxEMelfauLt6qRPzAJDslXipcH+rFVoHFLkjE9yPZpaoXxEPBHdcPzvVdYVePgq0hvh3Jr7/+GjpevlUQ8ZYOYCK+V0mZuap6bPrVv/7awnwJM6pHdRMLJOKbSypq3CmpMFfhV5esgpV0NjnfPghllnlq6fkn4htKIqQabq4yS+Y8ru4g3jHsbJR58zbYwifim0piErJQ4pHX1U1VX4d4m6ZPxDeUJCT4qnpUEk/PWPmGDZQtN9cXAap+bRvHTcQ3lNRJcBF/xXmnVTanfz3YF7Z6lW2XHDiDkhS9LVTUWEG8T5hY5NSLKNgG4mVui3dAqHqv+MC4qRdRMG7iGc2c2+Vk8tuLhW6y+TT0aZDGC4MhnjPOuCUxsn/0Byjr6+qtSUbBqOfAUdU4If6v/syo8Z1dvYDNHtwK4rU0oPziIo+Sf6QjV27jf3/nN/467/vGaD8iQL5jPr2xRdOdrfe0dIg+rq76CuuVGQWC+JvDyZ8q+QhWrDZZVi3FJH3RgiSJslTcKYl4B0ricWU5PbL6HxUMt8arIyaJeAdK4m/ePj89so44KFghudPklPS+D0tcDXAZe5mnWspjL6+ddSR3LMTTGk8qTTV6Ix7ULK4y1RclWfUOKG28lGnqQm+qPpR4TT714x2QVb3vOf3148NUfSI+APWIH6pkDFXvlLx9UjnVQiQIVV+RNCZhAA90KJIO1Ucll/ubqqkWIkHU+EMSIe7VyPf26MHHurJLaHZW2/dJLqZKuoydUtJjuoUIGPIgDfxYV1YJvWfV8n3Si4UTH9CG9olaVr136b0Sb/hSuie+UtXjx/qWlze0d8j9eOb8qsZAVX3XxAOqvrrM6/vD5F2z6j09d96lN5YMSYzqnPgakpb0qgFg0DU+MBWy4za+juSpNacyNobcxgfnwHZs1XclGQVDtuo3o57uYDTEYzXvaYmM5YEyyRRsCUK08SQIw3M83rv0JpI2v0dgmV4DSv45N9tH/NCs+ut7PRKf+w43biHxA7Pqbz/8XjvEk2zbSlI7JP7x48fep/SJoVr1y2kdVQ+OzlXVeXxSd6r+8eOBMj9Qqx4bGzWIh8fjK4gnhzv7RY8fD5V5yarfx3XeawSxB+Ltk0n5Ec9y5zbSF+dJu0w8tuqX+8yqH8YqVE1UPXsDPK+OjKEXB7aPeBJlO2VWvbRChWUVqj5Cr8KIV+8oiHhrmW1IDpV3qcZ/9Mlbx6zG8zVpPvsbXuO1NKDBpFCV0O+o90wqKzDvyvdIPBuQ2vjJ3pds7JivQoWr/Qt4MaK+4uq9JYMDsYckGQWuVag+ezYa4kPjcQclGQVVq1CBbfwwU6hGKxkFrlWo7Fb9QPLjt0QyCgbqwHHMP+VbJuqud96yZBQUxLPg6gGFV9tnlfXtnV9d5bWv3q9kFKgJFdbZ+hT08Jjs2ceJ+HZQa1i2+9Hr5sRLqh4h9bN+mdtIfMCwbA+j1y0QbwzShMxQW1Vmq5JRoA7LevHeS9iCvY0XkxFWKZ5EvAP1Jj/qIVCJzEkLXUZMP0rmtXbeSFL1Dgjil4Oy6jfWpbol4hcLp+pJ3TkHhHH31vHpka9V34sDx+IYllR95m5zEvEOyFb9cupr1ffisuXEa2+ZLOl+ARPxDsiBGOv9YRHPrqJfrCOSPJXY9hGPu3K/fuKX5tWXrx6Bxng3xPu+y1tIvD96VowOVR9cpuWN7YR4z9pxfW8yea0yg+nm3UtHtotfIoyCWlOaRs07aSJpLHAjJFtX9b7eDkJatePMTW0T4gOmNI2bd9KQeLAn0MXVQ4i/fvTFw3funK/JFFnXD98hE2Uti+0752Tr9slk74tHl3R9++tvH7KZtPjRKSljHTi9Vp0pTcdJPK3S+iJ2XV7dV9XTGv/F/WO6hSvf/WNcCfE2JuOa7SVbj7Cq/+LBOe51X+O/lCd6lO169AU91/vm6k5pOkJVXxI+tCAt1sZT+mkyAdlaTslutrfcun70+T4ZRCULMtA2mZ5T7Po8dAaGoaZQWVA/Bt5pvcXsIrLmlZPItpbT9R1Wk/FevuUmPnQaxWFPaaqjSdaL6567WYXKr6cgiC+UPP2LCV0XeoBvlapeJr5U9fQs75urSXy0YMuO0p3QmWfUMGgwWk5dVXmUCwjiN5Jxh/Xvtw6psc22bg4L4+5ooxC/4btqG3chU5rGi7LtJt0py6W1qwLLtD6LvGrw0IIaXbM6qDWlKfL9RX1ZA3Ct81bgTWq8i/g6vPdOvBNqVtDgMqjgGwq4Ta+MK4sQ3f3Pk8kf+sl3zKc3hhqBEySppUfzve2uGOtKxMQGmGEjDYZiGLWmQhlcP15Nj+a7dFXfLAKnIgN3tMQPK8q2pqRU9U37v2nMnTLZgo6lMafEWIinK3f5ReAMT9Wb8ZcdEF8hqfeKxkJ8QC5N7pss2+OYm7nT7Ph1G2ypT+QwFuIDMArie7s6AXOc1i0zCupMfhSdeL3qWsMrOrk6hFNzhGQsxEuTH1WhW+KrJiSMFyYFL1Rgddl6Xz0KoMmPqtCpcecRLd8D8RCb0NIkDpMxkHhrWnhHgCY/qkKn3TmfNInOI2Jhc7FL4o2UsfWU/yk/nFjCrXSxm5agiECTH1Wh27azWtV3efXiHuAknu5UvZkkKpPdnHhTZARRtkGSMg0tq/qKMvVTmhL/7Xt3ztfTm8PXvndEt+lebH/jvfTQ63c+OZzsXdL4HLz79aPytH/DMjfkGN1NwvUuSUFchGGwSZP1JJWq2ufVHz9+rCmJpsTvXS6P1tPlEYnKINt0L/64/a81IRL3IPGX9fSbS/yJt4reJBHdZ3vIyXj3ckoCeqZChKFWP35onjuBWMTTCSwbEA+28YT4Xx7TD66mcb197Rh/vX4D/yErHZH428kRleKnUWl6jO7Gf4i0EGGoQ/x8NiAHjqZ9WlL1gZLQlLUr/Xbc0K36gnhaUQXx9AAJxpryfZjQU7nGF8TTY3413juhAn18MqtPfMvOFq/VXQPLDJe0Ee8dqmSiIJ42zYL4JfGpr18nf6bUyz5dTyZ/IjXgnHh6zKuN90+oyOaL+sS37V7dWuIL6BWVwMfIr0SNhAp0Nqsdoda+X716dVcPhdtOGw9I+qt6C6h1TjZI+H2RYwcTLwn4oE5CxdyX9x5UvY9kdb2reXXpt6i2XdgqtFFQI6Gi5ZCm+pKekx91Rjy6WogqoAR9gT6+YaGGVT8U4iv5hFS9iNRAkKT/1UkRCwR3b0ZEfNCUprl3vHjPxFdn0osMOuXkmveZ2UYVRqXq/ZHLKs6J5sSb6tyu6j3mzmiXeGdjMxriA8KrZ1mThQrAp2V5TIA6hyVJoT6TprSi6r1s9bEQL4VX81WoNpuXP4Lmq29GPNw0NySeybU0aUpVPfbrnY+FeCm8mq9Qgd8AcKGC+aKJqg8i3qnqqwqt20mrSqgG37IGV6e4e/du2AkNAU1bztek2Tz/y/eBVajyi7O5R8ZRaCpSM7RYqCv1il2G/K2bR2bh4e5djXmre04/8I27ebYdhqx6vgrVyx/9z/tAjW+UVAzAmahQs8xGkva6jK4uUCHh35n0wd27BvM2aMTfvOsk3nrYtQrV84ODA2AVqpaJH8LigL55tQXx9JZbVfUm8Wsy/HLn/PYJG5BZHvEQCzpEQz25VGBZpMXzwAwapXFThmgsbVnzrlWoNmRFIvNwj8mIdcusllRo8/9FKKd/fay75sTvs6H49T4jvgixoIOyNM6CCvAqzQMzaJSGCNHwqPGnparnq1BZiG974pDeVb2hqINfZZ/+XJhxB7TxYly9GKItP6/foHEWbPyVE18EZojDNESjmvigOXA6kOwvnIvOcqjWWfoqt9w7b2jVF8QTHh99PmXRGGxPUeOL90Einu0sD9MQDY8a/0vvEb0uSCrDuVBV62krMzO6mC7fgHqNlnrnNZ1CIHiMzSmPpSyJL9r4CRO4fVK08TwwQxwmIRr8sIE6CxV0EXPHy7TNPVlZZoaMeWzCfAPNia/rBm4VvlEadVR9JylUqPShZ+6hAAfxdT131GRzTuxTlqyVqTy1SMSrERjBxAeo+nzWIOauShI/Yuugl7PMLAvy1cv7SCfNWeHFQbXMm0NFRbap6rtGnVmv5qZWtaBeTxrgQHqk9jJDRufUnTWJv/3oE7ueHAvxAeigxqsdKq72y10yK44ykWqcr/RD1tyofAPqAsTddOKgtviZo4EcDfH+89W33cZjexzqScszGVURL1gj7UTxTZJ0TmFsHeplgzHyu6Isfoa2gvhI89WTJp20HIBTSH7iFapemvRogRa8WwARH3Kfu0F8wHz1LfbjyeOjxEMBDsCcRpYy5dmuFqWRB6j6wDAQJP4W2D5VHzJfvXfpXsRT14vZk5YjJCsz4uR3hMEiCc+FWHs+7i0gPmS++hY9d7woTrwyX9kG2DbL1PVE6ZF1Eq+WaX3ndHTmsu0btZIm3VNWSAiOZVNnKCyPGqpeUQeqv4e2HgvriB/37DqI7zsti8LIxOkYNWa9QmeWcHITwY9JnZPUKqmoAx7dXs5riop3IVTVAxOjivOKPV0saUAA5GCV+OaylWw5FTVmvUIXsxPP0sMdOEYXHixToaZw9Ek1vzgcEPGnBumaV+fHOkkm2cBZlxwVMTY1UWfWq3ZUvSsUolLZqiczThZZ6Vdi+rzCDJS+0Q6q3mYocoL49ufq30DEs+xmEk6zJIPsZfyNHIwTcgUddWa9amViBGcoRAXxsD44yZDafnuYlsV1mGfCNUZQqnp17NARJ9eQeDafAc1sf/eSBtYs9/nkGPxLyBV01Jj1qtnECGUpQCiEdFS7gNSjts5zhxazusSzforPwKxqBroiJBu28cUMJvgPI17E34hgnKAraKhh1Wfz2dO2VX3FY5J9aPYJDkWNZWVDhhjsI1qBByEoPX5nbGxDq57VeBJOoxN/2jLx3oEYATXe3wZ2S6rE2wjSXgho4MfyzoT5G0q0SbwO1sa/TgNt9wXxSjBOk/KFcec3uSGBt3Fns4EBq7kiZFtR9RVS4NULy81CvPLawUP2rOeAeiOeooNuHEe9dec8K7KFeODpgwuAhYyk8RPkc5Qx/qK3D5apvSLgkH2hcnJld3ttPIw+iHeuO1c7GwnMMuLZR9LR/OzMECViIZlR+fziYm5NbJqfORKj5LMsS1rlYBmYd/8bpOiMyUDUWneu4rjwt+cuAaVqwcOyUNylNYSSzH6snKBKutRULrv+rIKkudhpX33FT5LGQS5cbULlmBu85redeOJbqRrAhbEy7sdyhS4WOI2DWhE41c+n2DiDXWFA1XKOiCtb1geqTD/kMAOBGu1HPDlz+4j3j8Cp9FeXTzY78zXuwEAM+By7aanH6MmRu3D8nnL1lhOjxkJ8wASH/imGvlY96ya5qpw4Zg38MmL0VuIbb/tx6xFiN5jYPuIDInAau2X0OKiim+SOcOaHHMSLQRZax03imQHo3dAA2D7iQyJwvEun9dhCppSzyvvHTl2LKlW9Msgi97lLVS+te1IdzgViC4n3R1Dola0aKzmrCLKWNT8cl7U/UD10K4fKBKMtnGMKiobYaeKD4updxEueWPMx6Z73auLV/h/4MinCyn1axbjLz3fEz6PMIaBWzF1QCpVD1Rcb5JFWEi8Nk1ijsH3rsXFLiXgvtJxJYyHe9sasbKrAOAmMqzfOKCTtncSk6jnanhgBVvX2q6sZNva4GanMiihan0UMt5N4Ml9Ka776jiWlRHZkZtPDyVblsKouwyW9/LVh9+ktGQW8xi8nfgsKE7T548EAB5eg0o83x3HU9ErN42tl16rqDYff1hG/odNe+WbSVEeaeipG2FYHSFA8d6LVpf9t7bduDWjES5Pm2O5Tbka2U9UXWPvNgWOboF/AUj+MGMoMIl7Vy8q+DX3tFrKdbbfuVuhKC4aWi5QnzbEN55g/VG4+3M9gNMSTZt5z2vK6xGtmGalQGZJqkvDgIrMfz7+sOGUQ8cpwzGoDOuYLWImX1bsRaS8kq6yCkRC/DFjAqLaq14ZRaIWSXhHJkUf+WDqNq/KIyJmSLiBRrbiPzPqpqnrVrWCKCElADr7PQaOWVV+7O6faSrRCQcQzCSvx9gmqDOLVfDirLymX3xB9kE+5+S1U9f5oZUYMmq0mmnJzPB6uedQRq/bKtQZBVvVqaXbv8YXq64VuQR/jr0QiPlTS1OCiOWZDL0qrbr8ZteXOHMSfZZbwW6O/sG3EkwiM0yPPGOtQVV8proQ0QcHXGvFyiXA0tGlhSEN3xv1YIoWgQneb+GoRJHXS4IqG1MkvdOJVcjRVD11K2QX0KYTlaFr6vjFFWxhz1zLx/CnbiZdCIvToVWT25dxXN2w9sDNZSMkXLuEXU7SNUbZdEg802trzN2uSVD91ToyrA24WdpZGPPPmGE4d/3q8lcTLa9Lw5cde/fjgzWfAOQGqHklaW+8ciSgMoEzeHDt9ZxubDGSIFa8aqB48SdpGVX8oOe340iTPP9g8B9akcQcqKY81V+KrAPWdnSi6wSzORryiPQytwHrtuSpJG3cgwrejejwO4snk16Xzplx+DL8DbG0a78wwLfGMpr+Vu8xEuHk+/9i9ilc+d6bfSXtyeZP8Z6ukyZJ0r3G9TpZDs6Jvgm2Q+vHlpFe/KomH16RxuWXUCqWoelOY9K1PTFUvd9bMcwuDUfOim/42tMiJstAdgODd7HKNpziltb6s8a9+AjXxAare+eMJQayfhk/KlegK2Ge6kWJ1Tmb6iKs85xUtjPXObZ5f5W53nPhlUeN5G//yhyDvQT/e4bwp3eeEKWnuPIkqbYBWROfJU1zR70gaOuc9dbbEkMPN6BO5q8hvn3FHFT3PwufLj30KLzhYNz6OQ86jYeRbidcCZ6zE4/p+Il2J9+PLsQAkvz7iLDAawI7t7M4FzJrWjHhz5INY5bCq5zu4GqeqnhKohljSycqlt4VxnHMfIH2VgGJZk+Bbj1GW8uP9JU1Vb4x16uHV3MUmpDI54S2H4nNOFkZ9L6ctJBOjZycibkNuR4qhG59fhOXwVXIPSYaREH9zSJY46yJ3zoCw9Iu/KvG8XVdU92JR7lTGBssSMkWLsG9EN7DOO13wXroaPxm3F3RZA3nk3gY21p/7/fDNaIg/PdqQaRQ7ma/e4mARdTCrIL5sjWntNIfmtM4dD8pblS+BbPnJL9Qiw28EU/UVATU8sHPbVD1x0i+n7Y3OyZJWl6ocNGeOx6s1UOxcSTsKba0nVRTDuIJ4i2uA2ITcb0jsAHulD7X/gySjQBB/++Fx68SzDhUQm4lkOnxaWT2upuzoyYF7pXBZ453hlko/Hs1mDtGtJf72ydH1A1zpLQuRqvD+SeiCzSqYWXR9+c1pVyvDtFzZih30tdLqqhSIYcuwUs5YUS1x4nhHunT1RAFv46/vTY5un/iFV/sTf6ZkmspH9D0rKXMdlfYft87KMVVhXkljL0ijkcdtkD6FbS5y9QYqVb0s6YmREB+C6qV7+PH5wjTUVAm+nZfj80WkbBEdmRX7i6G1xUKb58Iw1TciUmuFwJnykHlLPla9kPTDNhJfFVIuWkTTUAP2EGk78VdZyS2t9ytxWiZ388TlS+KJmjdC42U3cPLVB8GXeKjzA3TAWD8eVPWZFP3K95nEI72NL8ieL+gadNpkR9AY0EoWcP92D5lQySjoUNUjyN1RjpkixSGvTH5aUnSyABrokvisUPXSoqLqpc6ovW8oBfGeKMRXvc3q1duUjIIOXbaUeCPertjQfe1nhleGuFdOsoXBFu/H86EceVFR5VKsTyH7fvQrKKo+EV8J35+EVb3so9EPSl8ylIuqXZ4zy372VIzLaSvGlrH2pNXgC1GpDOcKu6bDT/tFSdVXocboXDn2DgMtPpZ0ehlo+XRWDssZxIvB27LPpvgDJd+Afl1ktAxtT+4SJhkFfY3OIWBcVAKUy6JHYAtVrw+vFyWzeo9O+EuwUo7SzQyZ9V2M+Hn9Il/sOPESnLEwUC4LrZlGXYTmq0Gip4CbfLIYFXkJDOILL582pKOP+Pn/opYko6BH4p2NqC2XxfS8ceLNgBq2m424kdOEF6FsPGhfQjMstRE/NxLxPUmqsTWFJAIU9mbDR35IUEbhFJJDtugLkIkmQZzWTSRdIr4KtLZZJbOFmuuG5lofvCxCibFhbmCFeGFoOiY4qUIi3leyUokyPuRRdnVLHlglo7A577+Lj8LCL4y7knyhwJV+HjBel4j3QkiLWO0WKYjX+uCKv0VysWYl8ZIc69MXxp2o3XrvHNkqPNciXr+ofcko8CNeyzmyrvFlwkOWpjAxOSndim+Ro7nYOZ/P5YwolqQ1p4uXyUWwpcvAlClgz1luLn7WFTrm0xudT4XiqR2M2FnEwq5J/ZT0gBzuyI37chxGDushu88WV9B4vHwNpmc0LVLxi1qXjIJ6qt5vDfWg0ETzBaGdr+xpJjzxlORcH2ZDmvHOT6fEu9YGFwM4SdV7oXJYlgH5LDhSlgmcT3pyi5//vIyMoUPyaG6sKiO9CXKQF8ot3T5xg4vZQnYD17vPxpJR0D3x1XHLtEzoEDbaZouTE6lARAZbZ7ZKjDLa61fG3FzxVDwOOxHvh+rxeIbqRcV4LbY9JpLyIIXC0diMM5ompw/zsm1M/Ez0+1fVb5068OOBHSc+SNL5lvBhNyBWBxVvhTx6jq07Oicd1K4Tyy/LslJ/K8RbbsKhbxy/qGXJKIjuuStINCSLGDvF44Z78bPFBYKJR2VY5oKqe0Q98EqsDQnXUcZn+NuQiPdC2z9eD71iO8UIfjEKi/ecIHSSs4BbUc8F8SfcQ4Mxm13NhQsXFTHaJ/I6SuKlScR7of0fX6w7rcTKqMoZ4WqM6zE2xnJ1ojrBXplTS/KZZ1m2mEtBV4WBn4jnGBLxSCJKFzjBZh7jL1cnqjN1PbXkib5f6UH3SdULDIN4purtxKOrpz+bFcG0uqtHcsRlsrceCr3SdUmKqw9CZ4/JDIssQLvzNjNQAp0SyR7VjxazTNURcly9/322LBkFsYg3wh3dBZ1kdPyN2OsOSTYBAis748lW0tuE2/2ZPHt5Ij4MLfx4OMDZAdo24xq7EL56QIhMeVJcYIEWOd0n5VOg2clTJfsiqfogxCCegkTjLP6pSGsFGgTqgRWj7tRvqCwNm81gd28i3gt9qnrFI0tTJz7OSGqV6s+j7wHrzknxV6viHFEQMqaurrrPPiSjYFDGnQnZI1sk0J5hQ49ETUr+PLp/MaN1WzoDcBhbBpQT8V7om/gyKovW7BwhJFS6RHyRYYmQ3nK7kqdaus9GklEwTOLVQBzFY18MtiqyTNXzhl+31cF0yXbusx3JKBgk8erEwyIZzlqmEpTjIt4ySpeI90K7Px6ZI+Ly6p8bnSytTKWOF7Z7pkpWRw8k4r3Q7rAsI81R451lqnUcVPXt3GeHklEwSOKddRQgXo2+S8T7IDrxkKoPKROFWQNeZfYsGQUQ8XwVKv6pY7APVLPZB3ufQwBEPF+hgn/qGO4DdZqBnV+9rmQUQMTzNWnEalR95RftAHpl1wGIeL4KlViNSsVYatJYJKPAr8arGMsDHYtkFGxXGz9OySiwWvWv3v9qdFb9SCWjIH4/PklGQSI+vmQUJOLjS0ZBIj6+ZBTUIT5hC5CI31Ek4ncUifgdRSJ+R5GI31Ek4ncUAcQrgTk2P34geDGvfnzw5jP85fd/0V6ZtLiW7/P5wcHBey3dZ1wEEK8M2tlG7gLBi3n+web5e6/+tnmBUpm0uJbvE+OzZy3dZ1wEEK8M09vG6gMhFfPig5d/fvCdr9orkxbX+n3ij5buMy4CiFcCc2zROYEQxbx6/yvyTP+uvTJpca3f52fPNi3dZ1wMpca/+skz8tG2Fmn9Pl/+RVFw8zLjYiBt/MsfPqMN/Yt27YYXbd8nZbyl+4yLUKueB+a0aS3jMj8trOU2y+ziPgnpbZUZF6kfv6NIxO8oEvE7ikT8jiIRv6NIxO8oBkn8ejKZHN0+mZLt5d7n+Ntrx9f3jvCXO+dk3+2TyQ8effnu5c27l+U5t0+OtA2G60eX2vb1oy/wvy+lc3cRQyT++sH55vr+8XrvktK43qc7X39jc/vhQ0o8JRD/kUm1Ea8UXMiTT/ncXcQgiWec3H54TLcL4h/9w+XN9z8ixN8cTu78+8N3Jr/1zmTvEmsHfHyJv1C+2QbdifXCZErq98N37pwvyZ5rfNI+q/F4a++nWKWwwk6u2KAAAAFESURBVHcQQyR+czqhKn2JSTmdUsW/h2v3vx6vf+8jXuMJe5dFpT89wjqCNgUbtsF2Yk5v3j7HgvePyRlk+x5uQI6YqienPzh3qIctxyCJJ5UaU4+pwv83vMZ/fvSf/2gQT96KyXSJ6+4poZBtsJ3X94guoIL44x6xE3AbspyWxNOXYFcV/kCJ5xSSSs+J//L7H/2HSfy+EC43uP7GbB8x4td3zm/e0onHJS+nkX5edAyReMwSM9GuH/7xsSD+13//xq2p6u8fk4r74PzmsFD1eIPtZC1FQfz+Zk17BrKq39y8xYzFXcQQiccWGjHLiIFODHuquIlSXk5l4r883MP/S+Putd8tjDuywXaeMtuAUnw4+RZ+IYRxR06/pFbEjmKQxPcF2m/YUewy8evJ7lb4nSZ+p5GI31Ek4ncUifgdRSJ+R/H/HIGlWFkrwi8AAAAASUVORK5CYII=" alt="Follow-up recipe predictions" style="max-width:100%;" />

``` r
plot(plan, type = "responses")
```

<img src="data:image/png;base64,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" alt="Follow-up recipe predictions" style="max-width:100%;" />

The plots show fitted values available before measurement. Exporting the
plan preserves every role and nomination, including the incumbent repeat.
The physical recipe ID can recur in a later round. The optional prefix gives
each planned preparation its own label; replace it with the lab's actual batch
identifier and reconcile labels when executing the rows. Use a different prefix
for the next batch even if its planned date is unchanged. These labels and the
planned date are export metadata; neither changes the fitted model or adds a
day effect.


``` r
eight_csv <- file.path(tempdir(), "lipid-screen-eight-total.csv")
svem_export_candidates_csv(plan, file = eight_csv, overwrite = TRUE,
  planned_date = "2026-11-02",
  preparation_prefix = "lipid-screen-20261102-round1")
eight_readback <- read.csv(eight_csv, stringsAsFactors = FALSE)
stopifnot(nrow(eight_readback) == 8L,
          identical(eight_readback$recipe_id, batch$recipe_id),
          identical(eight_readback$role, batch$role),
          !anyDuplicated(eight_readback$preparation_id),
          all(eight_readback$planned_date == "2026-11-02"))
```

## Freeze predictions, then reveal simulated outcomes


``` r
engines <- c("svemnet", enabled, "forward")
predictions <- function(models, rows, forward_models) {
  setNames(lapply(engines, function(e) {
    if (e == "forward") return(as.data.frame(lapply(forward_models, predict, newdata = rows)))
    as.data.frame(lapply(models, function(f) {
      if (e != "svemnet" && !identical(f$companions[[e]]$status, "ok"))
        return(rep(NA_real_, nrow(rows)))
      as.numeric(predict(f, rows, engine = e))
    }))
  }), engines)
}
before <- predictions(fits, batch, forward)
```

Until this point, only the 23 observed outcomes entered the analysis. We now
use the example simulator. Potency and Size use the original generating
surfaces and independent Normal errors with SDs 2 and 3, recorded to whole
units. This teaching simulation has no additional calendar-day shift.
Its response means are also the same for both operators and contain no
operator-by-recipe interaction.
The fitted companions were allowed to learn such interactions from sparse
operator-by-lipid cells. This difference matters when interpreting the
comparison. One batch cannot establish what caused a particular prediction
error or which model will give better nominees in another study.
PDI is a negative control. Its example distribution is Beta(4,16) for
every recipe and operator, with mean 0.20 and SD 0.0873. Those moments match
the initial PDI measurements reasonably well. The historical PDI noise law is
unknown. The simulator records PDI to hundredths. Each newly prepared run
receives one set of measurements, however many models nominated it. The
explicit preparation IDs key the teaching simulator's reproducible random
measurements; changing those IDs changes this illustrative realization even
with the same seed. They do not change the generating mean surfaces.


``` r
preparation_id <- eight_readback$preparation_id
observed <- lipid_screen_simulate(batch, run_id = preparation_id, seed = 731)
latent <- setNames(observed[paste0(endpoints, "_true")], endpoints)
comparison <- do.call(rbind, lapply(endpoints, function(y) {
  out <- data.frame(recipe_id = batch$recipe_id, endpoint = y)
  for (e in engines) out[[e]] <- before[[e]][[y]]
  out$latent_mean <- latent[[y]]
  out$recorded <- observed[[y]]
  out
}))
knitr::kable(comparison, digits = 3, row.names = FALSE)
```



|recipe_id |endpoint | svemnet|  gaupro| extra_trees| forward| latent_mean| recorded|
|:---------|:--------|-------:|-------:|-----------:|-------:|-----------:|--------:|
|R0019     |Potency  |  91.795|  94.649|      86.894|  92.309|      93.689|    97.00|
|R0432     |Potency  |  94.411|  90.989|      85.367|  95.191|      96.487|    98.00|
|R0259     |Potency  |  91.768|  93.555|      84.109|  92.033|      93.387|    91.00|
|R0171     |Potency  |  78.844|  77.881|      73.092|  77.610|      77.398|    77.00|
|R0600     |Potency  |  75.322|  76.234|      76.766|  73.596|      70.186|    70.00|
|R0444     |Potency  |  81.193|  76.932|      79.210|  81.399|      77.469|    76.00|
|R0384     |Potency  |  92.468|  82.761|      81.376|  94.045|      91.851|    94.00|
|R0466     |Potency  |  79.605|  73.568|      73.680|  78.154|      76.585|    77.00|
|R0019     |Size     |  88.103|  95.505|      86.161|  83.691|      82.245|    84.00|
|R0432     |Size     |  68.766|  69.175|      70.365|  64.985|      66.357|    64.00|
|R0259     |Size     |  79.969|  79.429|      81.053|  76.256|      77.302|    76.00|
|R0171     |Size     | 190.319| 185.555|     189.901| 193.133|     193.202|   196.00|
|R0600     |Size     |  71.044|  72.294|      73.972|  66.987|      66.357|    69.00|
|R0444     |Size     |  46.126|  53.431|      56.758|  48.680|      48.897|    48.00|
|R0384     |Size     |  94.634|  95.341|      96.303|  91.601|      90.052|    91.00|
|R0466     |Size     | 143.548| 143.021|     148.217| 146.549|     142.365|   144.00|
|R0019     |PDI      |   0.198|   0.298|       0.223|   0.210|       0.200|     0.12|
|R0432     |PDI      |   0.152|   0.208|       0.086|   0.165|       0.200|     0.23|
|R0259     |PDI      |   0.166|   0.180|       0.088|   0.172|       0.200|     0.26|
|R0171     |PDI      |   0.230|   0.211|       0.165|   0.252|       0.200|     0.16|
|R0600     |PDI      |   0.195|   0.131|       0.196|   0.198|       0.200|     0.10|
|R0444     |PDI      |   0.164|   0.239|       0.164|   0.193|       0.200|     0.18|
|R0384     |PDI      |   0.171|   0.208|       0.173|   0.200|       0.200|     0.20|
|R0466     |PDI      |   0.212|   0.192|       0.207|   0.243|       0.200|     0.15|



To compare models with the generating surface we evaluate the original
desirability at the three **latent means**. The **realized desirability** uses
the three recorded outcomes from one preparation. A noisy winner is a
candidate for independent confirmation. It supplies little evidence for
ranking the engines.


``` r
utility <- function(y) {
  parts <- lapply(endpoints, function(nm) {
    g <- goals[[nm]]
    z <- (y[[nm]] - g$lower_acceptable) / (g$upper_acceptable - g$lower_acceptable)
    if (g$goal == "min") z <- 1 - z
    g$weight * log(1e-6 + (1 - 1e-6) * pmin(1, pmax(0, z)))
  })
  exp(Reduce(`+`, parts) / sum(vapply(goals, function(g) g$weight, numeric(1))))
}
stopifnot(isTRUE(all.equal(utility(before$svemnet), batch$score, tolerance = 1e-7,
                           check.attributes = FALSE)))
pool_truth <- lipid_screen_simulate(pool)
best_feasible <- max(utility(setNames(pool_truth[paste0(endpoints, "_true")], endpoints)))
result <- data.frame(recipe_id = batch$recipe_id, latent = utility(latent),
                     recorded = utility(observed))
result$gap_to_best_in_pool <- best_feasible - result$latent
for (e in engines) {
  result[[paste0(e, "_pred")]] <- utility(before[[e]])
  fixed_pdi <- before[[e]]
  fixed_pdi$PDI <- rep(.20, nrow(fixed_pdi))
  result[[paste0(e, "_PDI_fixed")]] <- utility(fixed_pdi)
}
knitr::kable(result, digits = 3, row.names = FALSE)
```



|recipe_id | latent| recorded| gap_to_best_in_pool| svemnet_pred| svemnet_PDI_fixed| gaupro_pred| gaupro_PDI_fixed| extra_trees_pred| extra_trees_PDI_fixed| forward_pred| forward_PDI_fixed|
|:---------|------:|--------:|-------------------:|------------:|-----------------:|-----------:|----------------:|----------------:|---------------------:|------------:|-----------------:|
|R0019     |  0.841|    0.891|               0.054|        0.796|             0.796|       0.771|            0.830|            0.699|                 0.708|        0.809|             0.813|
|R0432     |  0.894|    0.883|               0.000|        0.899|             0.879|       0.812|            0.816|            0.737|                 0.703|        0.912|             0.897|
|R0259     |  0.845|    0.774|               0.049|        0.824|             0.811|       0.852|            0.844|            0.692|                 0.661|        0.834|             0.822|
|R0171     |  0.276|    0.262|               0.618|        0.304|             0.310|       0.315|            0.317|            0.230|                 0.226|        0.271|             0.280|
|R0600     |  0.312|    0.317|               0.582|        0.468|             0.467|       0.505|            0.490|            0.503|                 0.502|        0.422|             0.422|
|R0444     |  0.547|    0.513|               0.348|        0.654|             0.643|       0.516|            0.528|            0.593|                 0.582|        0.647|             0.645|
|R0384     |  0.793|    0.828|               0.102|        0.806|             0.794|       0.608|            0.610|            0.588|                 0.580|        0.828|             0.828|
|R0466     |  0.403|    0.418|               0.491|        0.458|             0.461|       0.337|            0.335|            0.329|                 0.330|        0.416|             0.427|



The paired prediction columns show how much each model's score changes when
only its fitted PDI is replaced with 0.20. This is an after-selection check.
The known PDI mean contributed no information to the original nominations.
For context, simulate the range of outcomes one preparation might produce.
These extra draws are kept out of the follow-up data.


``` r
ref_interval <- t(vapply(seq_len(nrow(batch)), function(i) {
  repeats <- batch[rep(i, 500L), , drop = FALSE]
  ids <- paste0("interval-", batch$recipe_id[i], "-", seq_len(nrow(repeats)))
  sim <- lipid_screen_simulate(repeats, run_id = ids, seed = 903)
  unname(quantile(utility(sim), c(.05, .95)))
}, numeric(2)))
knitr::kable(data.frame(recipe_id = batch$recipe_id, lower_90 = ref_interval[, 1],
                       recorded = result$recorded, upper_90 = ref_interval[, 2]),
             digits = 3, row.names = FALSE)
```



|recipe_id | lower_90| recorded| upper_90|
|:---------|--------:|--------:|--------:|
|R0019     |    0.704|    0.891|    0.900|
|R0432     |    0.746|    0.883|    0.937|
|R0259     |    0.680|    0.774|    0.907|
|R0171     |    0.203|    0.262|    0.327|
|R0600     |    0.118|    0.317|    0.423|
|R0444     |    0.380|    0.513|    0.638|
|R0384     |    0.676|    0.828|    0.865|
|R0466     |    0.304|    0.418|    0.472|



## Refit using the eight new measurements


``` r
analysis_cols <- c("RunID", "recipe_id", predictors, endpoints)
followup <- cbind(RunID = preparation_id, batch[c("recipe_id", predictors)],
                  observed[endpoints])
augmented <- rbind(dat[analysis_cols], followup[analysis_cols])
stopifnot(!anyDuplicated(augmented$RunID))
set.seed(2301)
next_plan <- svem_plan_runs(
  spec, data = augmented, goals = goals, candidates = pool, budget = 8L,
  fit_control = fit_controls, score_control = score_controls,
  reference = list(Operator = "A"), benchmark = fixed_benchmark
)
updated <- next_plan$fits
updated_forward <- next_plan$forward
after <- predictions(updated, batch, updated_forward)
refit_scores <- data.frame(recipe_id = batch$recipe_id)
for (e in engines) {
  refit_scores[[paste0(e, "_before")]] <- utility(before[[e]])
  refit_scores[[paste0(e, "_after")]] <- utility(after[[e]])
}
knitr::kable(refit_scores, digits = 3, row.names = FALSE)
```



|recipe_id | svemnet_before| svemnet_after| gaupro_before| gaupro_after| extra_trees_before| extra_trees_after| forward_before| forward_after|
|:---------|--------------:|-------------:|-------------:|------------:|------------------:|-----------------:|--------------:|-------------:|
|R0019     |          0.796|         0.840|         0.771|        0.889|              0.699|             0.887|          0.809|         0.851|
|R0432     |          0.899|         0.899|         0.812|        0.883|              0.737|             0.881|          0.912|         0.897|
|R0259     |          0.824|         0.832|         0.852|        0.814|              0.692|             0.746|          0.834|         0.828|
|R0171     |          0.304|         0.264|         0.315|        0.277|              0.230|             0.277|          0.271|         0.257|
|R0600     |          0.468|         0.351|         0.505|        0.342|              0.503|             0.393|          0.422|         0.354|
|R0444     |          0.654|         0.569|         0.516|        0.514|              0.593|             0.573|          0.647|         0.580|
|R0384     |          0.806|         0.806|         0.608|        0.818|              0.588|             0.814|          0.828|         0.820|
|R0466     |          0.458|         0.415|         0.337|        0.404|              0.329|             0.406|          0.416|         0.407|



``` r
next_batch <- next_plan$worklist
stopifnot(next_plan$shortfall == 0L,
          !any(next_batch$recipe_id[next_batch$new_recipe] %in% augmented$recipe_id))
knitr::kable(next_batch[c("slot", "role", "recipe_id", "new_recipe", recipe_cols)],
             digits = 4, row.names = FALSE)
```



| slot|role               |recipe_id |new_recipe |    PEG| Helper| Ionizable| Cholesterol|Ionizable_Lipid_Type | N_P_ratio| flow_rate|
|----:|:------------------|:---------|:----------|------:|------:|---------:|-----------:|:--------------------|---------:|---------:|
|    1|incumbent_repeat   |R0019     |FALSE      | 0.0474| 0.4297|    0.1396|      0.3833|H101                 |         6|       2.5|
|    2|optimum            |R0071     |TRUE       | 0.0500| 0.5000|    0.1000|      0.3500|H101                 |        12|       3.0|
|    3|optimum            |R0509     |TRUE       | 0.0500| 0.5500|    0.1000|      0.3000|H101                 |         8|       1.5|
|    4|high_spread_medoid |R0503     |TRUE       | 0.0500| 0.1500|    0.5500|      0.2500|H101                 |        12|       1.0|
|    5|high_spread_medoid |R0374     |TRUE       | 0.0400| 0.1000|    0.3000|      0.5600|H101                 |         8|       2.5|
|    6|high_spread_medoid |R0506     |TRUE       | 0.0400| 0.1000|    0.5000|      0.3600|H103                 |         6|       1.0|
|    7|high_spread_medoid |R0355     |TRUE       | 0.0200| 0.1000|    0.6000|      0.2800|H102                 |        10|       2.5|
|    8|high_spread_medoid |R0117     |TRUE       | 0.0500| 0.5500|    0.2500|      0.1500|H103                 |        12|       3.0|



The next round returns to automatic selection within operator A. Inspect
`next_plan$incumbent$by_context` to see the best raw observed scores in each
context; a larger score elsewhere does not make raw block outcomes directly
comparable. An explicit incumbent override remains available when the
scientist wants a deliberate bridge across contexts.

Within the reference context, the default `incumbent_summary = "best"` retains
the best observed **individual preparation**. An older high observation can
therefore remain the incumbent after a poorer repeat; the rule does not erase
that older measurement or average confirmation results automatically. Inspect
`incumbent$recipe_n_preparations` for all attempts at that physical recipe in
its source context, `recipe_n_complete` for complete preparations, and
`recipe_mean_score` for desirability evaluated at their mean responses. These
diagnostics do not change the selected statistic: `n_preparations` is still
one for the default best individual preparation. To compare recipe means over
repeated independent preparations in the same context instead, rescore the
updated fits explicitly:


``` r
mean_reference <- svem_plan_runs(
  spec, augmented, goals, pool, budget = 8L,
  fits = updated, forward = updated_forward, score_control = score_controls,
  reference = list(Operator = "A"), incumbent_summary = "mean",
  benchmark = fixed_benchmark
)
```

This optional rule averages measured response values for each recipe before
applying the goal maps. Keep technical readings summarized within preparation
and keep distinct preparation IDs. Reusing fitted lists is appropriate for
this choice or changed scientific weights on unchanged data. Adding or
correcting measurements requires refitting; the planner checks the stored
training rows before accepting supplied models.

Review optional-engine diagnostics as well as the nominations. A diagnostic
such as a fitted nugget reaching its optimization bound warrants inspecting
the affected fit; it is not evidence that any nominated recipe has been
confirmed. The prospective comparisons still
use the predictions frozen before these new measurements were revealed.

The after columns are fitted on data that include the follow-up preparations.
Only the before columns are prospective predictions for this batch. With
real measurements, replace `observed[endpoints]` with the laboratory results,
record the calendar date and operator with each unique preparation ID, and
refit with the same declared goals and candidate domain. Operator is a model
factor here; the teaching simulation has no additional calendar-day shift. The
planner does not add a round or batch factor automatically. Declare intended
blocks from physical conditions, such as a changed lot, operator or instrument,
before seeing the results. Do not add or remove a block because a round
looked good or bad. Once a new level has measurements,
rebuild the specification from the augmented data and refit. Primary SVEM then
uses additive penalized block effects without recipe-by-block interactions.

The portfolio already reserves the selected incumbent repeat; the optional fixed
benchmark can connect successive contexts even when the incumbent changes.
Exploration also broadens recipe coverage. Neither guarantees separation of
block effects from recipe effects, so plan bridging preparations deliberately.
Compare new recipes with the repeats prepared in the same batch, rather than
only with their earlier measurements. Under an additive batch shift, that
within-batch contrast removes the shift without relying on a fitted block
term. A fixed benchmark chosen in advance avoids choosing the reference for
an unusually high observed response; it still has preparation and measurement
noise. A deliberate challenge run to revisit a poor response remains a manual
choice.

For the **first follow-up round**, suppose the operator had been unknown and
A and B equally likely. Using the original 23 preparations, score each
physical recipe twice and average its two scores before ranking. This answers
a different operational question from preparing every recipe with A. The
same recipe IDs align the two tables.


``` r
pool_b <- pool
pool_b$Operator <- factor("B", levels = levels(pool$Operator))
by_b <- svem_score_random(fits, goals = goals, data = dat, candidates = pool_b,
                          increments = increments, budget = 3L, verbose = FALSE)
idx <- match(tab$recipe_id, by_b$score_table$recipe_id)
stopifnot(!anyNA(idx))
cols <- intersect(c("score", "gaupro_score", "extra_trees_score"), names(tab))
eligible <- !tab$recipe_id %in% dat$recipe_id
avg <- (tab[cols] + by_b$score_table[idx, cols, drop = FALSE]) / 2
vapply(avg, function(s) {
  use <- which(eligible & is.finite(s))
  if (length(use)) as.character(tab$recipe_id[use[which.max(s[use])]]) else NA_character_
}, character(1))
#>             score      gaupro_score extra_trees_score 
#>           "R0432"           "R0259"           "R0112"
```

## A spline nominee on a two-factor process

The full lipid experiment above has four mixture components, two numeric
process settings and a nominal ingredient type. The optional general spline
can now represent that factor specification, including its five independent
numeric coordinates and additive operator shifts. A final optional comparison
below adds its nomination to the same initial data. This separate teaching
experiment keeps two continuous coordinates so its fitted ridge is easy to
display and compare with the known simulated response.

Both coded process settings run from zero to one. A response ridge follows a
line through the interior of the feasible grid. We know its noiseless response
for teaching and can check the nominations against it after fitting. In the
laboratory the measured follow-up supplies that check.


``` r
process_grid <- expand.grid(rate = seq(0, 1, length.out = 9),
                            time = seq(0, 1, length.out = 9))
process_grid$recipe_id <- sprintf("P%03d", seq_len(nrow(process_grid)))
process_mean <- function(rate, time) {
  along <- exp(-0.5 * ((rate - 0.65) / 0.35)^2)
  across <- exp(-0.5 * ((time - 0.2 - 0.6 * rate) / 0.13)^2)
  0.15 * rate + along * across
}
set.seed(7312)
process_data <- process_grid[sample.int(nrow(process_grid), 24L), ]
process_data$yield <- with(process_data, process_mean(rate, time)) +
  rnorm(nrow(process_data), sd = 0.08)
process_spec <- bigexp_terms(yield ~ rate + time, data = process_data,
                             report = FALSE)
process_goals <- list(yield = list(goal = "max", weight = 1))
```

Specify both counts when you know how many new candidates to test and how
many runs to spend on exploration. The planner adds the best-observed repeat
and any distinct fixed benchmark. With both optional companions available,
five requested model candidates and two exploration runs plus the incumbent
give eight planned preparations. Set `incumbent = FALSE` to request the same
numbers of new recipes without that repeat. To impose a lab limit, supply `budget`.
The worklist can be shorter if model nominees coincide.
The standard order tries elastic-net SVEM, GauPro, Extra Trees, forward SVEM
and then the optional spline. When two models nominate the same recipe, their
endorsements share one physical place. The planner considers later model
nominees before declaring a model shortfall. An explicit exploration request
reserves those places for exploration.


``` r
if (requireNamespace("mgcv", quietly = TRUE)) {
  n_models <- 3L + length(enabled)  # Elastic net, forward, spline, and enabled companions
  set.seed(7313)
  process_plan <- svem_plan_runs(
    process_spec, process_data, process_goals, candidates = process_grid,
    incumbent = TRUE, n_optima = n_models, n_explore = 2L,
    spline = TRUE,
    fit_control = list(companions = if (length(enabled)) enabled else FALSE),
    score_control = list(verbose = FALSE)
  )
  counts <- process_plan$allocation_counts
  print(data.frame(type = c("Reference repeats", "New model candidates", "Exploration"),
    requested = c(1L, n_models, 2L),
    selected = c(counts$distinct_repeats, counts$selected_optima, counts$selected_explore)))
  print(process_plan$nominations[, c("engine", "status", "new_recipe_id", "slot")])
  if (nrow(process_plan$merged_nominations)) print(process_plan$merged_nominations)
  process_runs <- process_plan$worklist[, c("slot", "role", "nominated_by",
    "recipe_id", "rate", "time", "yield_spline_pred")]
  process_runs$example_mean <- with(process_runs, process_mean(rate, time))
  print(process_runs)
} else {
  message("Install the suggested mgcv package to run the spline example.")
}
#> Warning: Requested allocation is short by 2 model and 0 exploration
#> preparation(s).
#>                   type requested selected
#> 1    Reference repeats         1        1
#> 2 New model candidates         5        3
#> 3          Exploration         2        2
#>        engine   status new_recipe_id slot
#> 1   incumbent selected          <NA>    1
#> 2     svemnet selected          P054    2
#> 3      gaupro selected          P042    3
#> 4 extra_trees selected          P054    2
#> 5     forward selected          P054    2
#> 6      spline selected          P061    4
#>   slot recipe_id             role      source previous_nominated_by
#> 1    2      P054          optimum extra_trees               svemnet
#> 2    2      P054          optimum     forward svemnet + extra_trees
#> 3    1      P051 incumbent_repeat      gaupro              observed
#> 4    1      P051 incumbent_repeat extra_trees     observed + gaupro
#>    slot               role                    nominated_by recipe_id  rate
#> 51    1   incumbent_repeat observed + gaupro + extra_trees      P051 0.625
#> 54    2            optimum svemnet + extra_trees + forward      P054 1.000
#> 42    3            optimum                          gaupro      P042 0.625
#> 61    4            optimum                          spline      P061 0.750
#> 74    5 high_spread_medoid                         svemnet      P074 0.125
#> 81    6 high_spread_medoid                         svemnet      P081 1.000
#>     time yield_spline_pred example_mean
#> 51 0.625        0.68821923   1.02008854
#> 54 0.625        0.56499915   0.39510586
#> 42 0.500        0.63176544   0.93828328
#> 61 0.750        0.63999195   0.82664139
#> 74 1.000        0.09754001   0.01875006
#> 81 1.000        0.44148405   0.33573545
```

This run requested 5 distinct model candidates and obtained 3. Overlapping nominations share a preparation. Counting reference repeats (1) and exploration settings (2), the worklist contains 6 distinct preparations. The number of unused places is 2.

The `yield_spline_pred` column is the fitted ensemble mean for the recipe in
that row. `example_mean` is the noiseless response used to generate this
teaching experiment. Neither column is a new assay. Check the actual number of
distinct preparations and obtain independent measurements before adopting a
model nominee. With `n_explore = NULL`, the existing high-spread medoid policy
fills spare places. With an explicit count, inspect the returned allocation for
any unused capacity or unavailable nominations.

To fit the same spline independently, the shared response-fitting helper
accepts the original factor specification and the response-free feasible grid.
Set this reference grid before looking at outcomes and keep it fixed across
refits. Its physical coordinates determine the basis. The declared factor
bounds separately determine which recipes are legal. If the planner fits the
spline itself, its `candidates` table supplies the reference grid. To score a
different pool on the existing basis, supply the saved spline fits. Prediction
at recipes off the saved grid requires `mgcv`, including after reload.


``` r
spline_fits <- svem_fit_responses(process_spec, process_data, "yield",
                                  method = "spline", domain_grid = process_grid)
predict(spline_fits[["yield"]], newdata = process_grid[1:3, ])
```

The spline regularizes curvature and chooses a smoothing strength for each
paired training and validation weighting. A member scores its training-fitted
Gaussian mean and residual variance on the validation-weighted copy. The
capacity limit uses the smoother's effective degrees of freedom. The basis follows the
thin-plate regression spline construction of [Wood (2003)](https://doi.org/10.1111/1467-9868.00374).


## An optional general spline on the full lipid experiment

Return to the original 23 preparations for this comparison. The specification,
goals, mixture bounds, ingredient types, process settings and operator remain
visible in the earlier code. Reuse the four fitted learners and add the spline.
Five requested model optima, two exploration medoids and one incumbent repeat
give a nominal budget of eight preparations. Distinct model nominees can be
fewer than five when learners agree; the returned shortfall and endorsements
show what can actually be prepared.


``` r
if (requireNamespace("mgcv", quietly = TRUE)) {
  set.seed(2401)
  spline_plan <- svem_plan_runs(
    spec, data = dat, goals = goals, candidates = pool,
    fits = plan$fits, forward = plan$forward, spline = TRUE,
    incumbent = incumbent_row, benchmark = fixed_benchmark,
    n_optima = 5L, n_explore = 2L,
    spline_control = list(nBoot = 200L), score_control = score_controls
  )
  capacity <- do.call(rbind, lapply(names(spline_plan$spline), function(r) {
    f <- spline_plan$spline[[r]]
    data.frame(response = r, numeric_coordinates = f$spec$dimension,
               model_columns = length(f$penalty),
               median_effective_df = median(f$member_diagnostics$effective_df),
               intercept_member_fraction = mean(f$member_diagnostics$intercept_only))
  }))
  print(knitr::kable(capacity, digits = 3, row.names = FALSE))
  print(knitr::kable(spline_plan$worklist[c("slot", "role", "nominated_by", "recipe_id")],
               row.names = FALSE))
  notes <- spline_plan$scores$engine_status
  notes <- notes[notes$engine == "spline" & nzchar(notes$message), , drop = FALSE]
  if (nrow(notes)) print(knitr::kable(notes, row.names = FALSE))
  cc <- spline_plan$allocation_counts
  print(knitr::kable(data.frame(reference_repeats = cc$distinct_repeats,
    model_candidates = cc$selected_optima, exploration = cc$selected_explore,
    unused_capacity = cc$unused_capacity, shortfall = spline_plan$shortfall),
    row.names = FALSE))
}
#> Warning: Requested allocation is short by 2 model and 0 exploration
#> preparation(s).
```



|response | numeric_coordinates| model_columns| median_effective_df| intercept_member_fraction|
|:--------|-------------------:|-------------:|-------------------:|-------------------------:|
|Potency  |                   5|            75|                6.87|                      0.21|
|Size     |                   5|            75|               10.51|                      0.11|
|PDI      |                   5|            75|                1.00|                      0.70|


| slot|role               |nominated_by                    |recipe_id |
|----:|:------------------|:-------------------------------|:---------|
|    1|incumbent_repeat   |observed                        |R0019     |
|    2|optimum            |svemnet + extra_trees + forward |R0432     |
|    3|optimum            |gaupro                          |R0259     |
|    4|optimum            |spline                          |R0112     |
|    5|high_spread_medoid |svemnet                         |R0037     |
|    6|high_spread_medoid |svemnet                         |R0551     |


|engine |response |status |message                                                                                                      |
|:------|:--------|:------|:------------------------------------------------------------------------------------------------------------|
|spline |PDI      |ok     |70% of members selected the intercept-only fit; the directional nomination comes from the remaining members. |


| reference_repeats| model_candidates| exploration| unused_capacity| shortfall|
|-----------------:|----------------:|-----------:|---------------:|---------:|
|                 1|                3|           2|               2|         2|

The stored general surface uses a shared Duchon basis and modest one-factor
curvature here. Ingredient levels
add penalized shifts and deviations; the operator remains an additive block.
The smoothing strength is selected separately for each self-validated member
using its training-fitted Gaussian mean and variance on the validation copy.
The self-validation score and effective-df cap can favor a simple fit in a small experiment. A flat
general surface supplies predictions but no directional optimum. This extends
factor coverage; it does not transfer the earlier two-coordinate simulation
results to every formulation problem. Confirm the proposed recipes experimentally.

For a new experiment, the same call uses its own shared `bigexp_terms()`
specification and response-free candidate grid. Multiple disjoint mixture
groups and fixed components are supported. With no varying numeric coordinate,
the general smoother reduces to penalized categorical effects or an intercept.
The default four-learner workflow remains available without a spline request.
