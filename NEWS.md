# woven 0.99.6

## New features

* New `woven_tune()`: chooses `gamma_y`, `ridge_w`, `screen_top` (and
  optionally `lambdas`) by the held-out balanced error rate of
  class-stratified cross-validation over the anchors, refitting
  standardization, graphs and projections in every fold, then refits all
  subjects with the best setting. Within a fold the eigendecomposition is
  solved once per (`gamma_y`, `lambda`) and reused for every (`ridge_w`,
  `screen_top`), giving the same predictions as refitting `woven()` for
  each setting at a fraction of the cost. Returns a `wovenTune` object
  (`tune$best`, `tune$results`, `tune$fit`; `plot(tune)` shows the error
  over the grid).
* `woven()` accepts one table per modality covering different subjects:
  when every table has row names (subject IDs), woven matches them and
  treats a subject absent from a table as missing that modality. Tables may
  be data frames, and `Y` may be a vector named by subject ID.
* `woven()` reports progress: a summary of the input (features and subjects
  per modality, anchors), one line per step with timings, and a closing
  summary of coverage against a complete-case analysis and the top features
  per modality. `plot = TRUE` (the default in interactive sessions) draws
  the latent space when the fit finishes.
* The fit stores the class labels, so `plot(fit)`, `summary(fit)` and
  `woven_metrics(fit)` no longer need them passed again.
* `plot(fit, type = )` is one entry point for every plot: `"scores"`
  (default), `"features"`, `"loadings"`, `"coverage"` and `"variance"`.
* New `woven_plot_coverage()`: subjects by their pattern of observed
  modalities, split by class, against what a complete-case analysis keeps.
* New `woven_top_features()`: the top features per modality as a table,
  with each feature's VIP and the class in which it is highest.
* `woven_plot_vip()` now shows every modality by default (one panel each)
  and colours each feature by the class in which it is highest
  (`color_by = "modality"` for the previous single colour).
* `woven_scores()` and `woven_predict()` accept new data in the same forms:
  data frames, tables matched by row name, and named lists that leave out
  modalities (scored without them).

## Changes

* VIP scores normalize each dimension's loading vector to unit length, as
  in PLS, so the mean squared VIP is 1 within a modality and VIP = 1 marks
  the average feature. Rankings can differ slightly from 0.99.5.
* Features without column names are named `<modality>_<j>`; unnamed
  modalities are named `Modality 1`, `Modality 2`, ...
* Plots adapt to the data: point size and transparency scale with the
  number of subjects; more than 8 classes get a non-repeating palette;
  ellipses are drawn only for classes with at least 5 anchors; a fit with
  `K = 1` is drawn as the distribution of dimension 1 per class; numeric
  labels with many values (e.g. age) are drawn as a gradient; long feature
  names are shortened; features removed by `screen_top` are left out; the
  coverage plot combines rare patterns of modalities (`max_patterns`).
* In the latent-space plot, subjects scored from a subset of modalities are
  drawn as open points.
* `woven_plot_loadings()` gives each modality its own scale and keeps the
  modality order. `woven_plot_variance()` is relabelled as the share of the
  squared canonical correlations (it is not variance explained) and drops
  the second axis.
* Fits are numerically unchanged.

## Minor

* The `woven()` MultiAssayExperiment example no longer wraps itself in a
  `requireNamespace()` check, since MultiAssayExperiment is in Imports.
* Source and roxygen lines wrapped to 80 characters. No change in results
  or error messages.
* README, DESCRIPTION and the vignette describe `woven_tune()`, and
  `?woven_top_features` notes how `screen_top` affects the VIP scale.

# woven 0.99.5

* `Depends: R (>= 4.6.0)`, matching the R version the package is built and
  checked with on Bioconductor.
* `woven_precompute()` accepts a `MultiAssayExperiment`, converted the same
  way as in `woven()`, so MAE users can reuse precomputed graphs for tuning.
  The vignette's Tuning section shows this on the ACC data.
* `MultiAssayExperiment` and `SummarizedExperiment` moved from Suggests to
  Imports, since a `MultiAssayExperiment` is now the main input; the
  `requireNamespace()` guard is gone.
* `woven()` split into input checking (`.woven_check_input()`), the fit, and
  result assembly (`.woven_assemble()`) helpers. Results are unchanged.
* Documentation and vignette now say that precomputed graphs only apply to
  fits on the same subjects; cross-validation folds are fit without
  `precomp`. The vignette describes the fold-wise tuning procedure.
* Code restyled with `styler` (4-space indent).

# woven 0.99.4

## Changes in behaviour and API

* Accessors renamed to snake_case to match the rest of the API, and `Z()`
  renamed so it is not a one-letter name: `latent_scores()` (was `Z()`),
  `projection_matrices()` (was `WList()`), `anchor_indices()` (was
  `anchorIdx()`), `singular_values()` (was `singularValues()`). They are
  documented together at `?"woven-accessors"`.
* A user-supplied `anchor_idx` is now checked in `woven()` and
  `woven_mcca_dual()`: indices must be distinct integers within `1..n`, and
  every anchor must be observed in every modality. Offending subjects are
  named in the error (previously this failed later with "non-conformable
  arrays").
* `woven_precompute()` records the `scale`, `k_nn` and number of subjects it
  used. `woven()` stops if its own `scale` differs from the precomputed
  graphs, stops if the number of subjects differs, warns if `k_nn` differs
  (the precomputed graphs' `k_nn` is used and recorded in the fit), and warns
  if a hand-built `precomp` carries no scale record.
* `woven()` and `woven_mcca_dual()` warn when `K` exceeds a modality's
  number of features (each projection matrix has rank at most
  `min(n_a, p_v)`).
* `woven_scores()` and `woven_predict()` accept a `MultiAssayExperiment`;
  experiments are matched to the training modalities by name and features by
  row name.
* `woven_predict()` builds its class centroids (or k-NN reference set) from
  the anchors' latent scores computed by the same projection as the new
  subjects, so reference and new subjects lie in the same space. Previously
  it used the anchors' eigenvector scores, which differ in scale from
  projected scores when a modality has fewer features than anchors.
  `method` is now matched with `match.arg()`.
* `woven_nystrom_error()` now standardizes `X_list` with the fit's stored
  center and scale, refits with all of the fit's hyperparameters (including
  `k_nn`, `ridge_w`, `screen_top`), aligns the sign of each refitted
  dimension with the full fit, compares the held-out subject's projected
  score with its score in `latent_scores(fit)`, and checks that `n_loo` is
  between 1 and the number of anchors. The unused `sigma_proj` argument was
  removed.
* `.mae_to_xlist()` (used by `woven()` for MultiAssayExperiment input) stops
  with an error naming the subject and experiment when an experiment has
  replicate samples for one subject, instead of silently keeping the first.
* `woven_metrics()` no longer takes `...` (extra arguments were computed and
  then discarded).
* Laplacians in `woven_precompute()` are built in a serial loop inside one
  internal helper. BiocParallel and parallel are no longer dependencies: one
  graph per modality takes about 0.5 s for a 300 x 5,000 modality, so
  parallel workers saved little.

## Documentation

* The vignette now uses BiocStyle, has a floating table of contents, and its
  main workflow is a real block-missing dataset in a MultiAssayExperiment
  (TCGA adrenocortical carcinoma, `miniACC`), including mapping scores back
  to `colData`, a held-out prediction, a glossary, comparison with DIABLO,
  MOFA2 and missRows, and a Limitations section.
* The package help page (`?woven-package`) lists the main entry points.
* DESCRIPTION defines SUMCOR, block-missing ("ragged") data, ESS retention
  and the Nystrom extension; biocViews no longer include Sequencing or
  MultipleComparison.
* `?woven_nmi`, `?woven_nystrom_error` and `?woven_metrics` explain that
  results depend on the random number generator and show `set.seed()` in
  the examples.
* The `Y_labels` slot is documented as it behaves: class names
  `levels(as.factor(Y))`, which for numeric `Y` are the numbers as strings.
* `?woven_example` gives the data source and license; the generating script
  moved from `data-raw/` to `inst/scripts/make_woven_example.R` and
  reproduces `data/woven_example.rda` exactly.
* README corrected: MOFA2 accepts samples with whole views missing.

## Bug fixes and minor improvements

* `woven_all_metrics()` (internal) passed two matrices to
  `woven_nystrom_error()` positionally, so the second bound to an unused
  argument and the Nystrom error was computed from one modality only. It now
  takes `X_list`.
* `woven_mcca_dual()` is split into named internal steps (anchor
  Laplacians, metric square roots, label kernel, block dual matrix,
  eigendecomposition, anchor scores, feature screening, W recovery). Results
  are bit-for-bit identical to 0.99.3.
* Shared internal helpers for detecting block-missing rows
  (`rowSums(is.na(X)) == ncol(X)`), for projecting subjects (used by
  `woven()` and `woven_scores()`), for checking new data, and for checking
  that ggplot2 is installed (now called by every plotting function).
* The plotting functions give informative errors for an out-of-range `dim`
  or an unknown `modality`.
* Removed unused internal code (`compute_B()`, `compute_B_kernel()`, the
  truncated RSpectra path in `mat_sqrt()`); RSpectra is no longer in
  Suggests. Internal helpers no longer have Rd pages.
* Every `importFrom` is used; other packages are called with `::`.
* Removed the stale `data-raw/test_solver.R` and comments that referred to
  files outside the package.
* Removed `Author:` and `Maintainer:` from DESCRIPTION (`Authors@R` is the
  single source).
* New tests: solver results against hand-computed eigenvalues and scores
  (V = 2 and V = 3), input validation, precompute/scale checks,
  reproducibility, `woven_scores()`, `woven_predict()`, the plotting
  functions, `woven_effect_bias()`, `woven_all_metrics()` and internal
  helpers.

# woven 0.99.3

## Changes in behaviour and API

* DESCRIPTION no longer claims sparse projection matrices; `woven()` recovers
  dense projection matrices, optionally ridge-regularized and
  feature-screened via `ridge_w` and `screen_top`.

## Bug fixes and minor improvements

* `.scale_fit()` (internal) rewritten without `<<-`.
* `summary()` prints through an internal `wovenSummary` class with its own
  `show` method; `summary()` still returns the metrics vector invisibly and
  the printed output is unchanged.
* Escaped literal braces in a roxygen comment that were parsed as Rd markup.
* `.Rbuildignore` excludes check-output directories.
* Added tests for `woven_precompute()` reuse, a malformed `precomp`,
  `summary()`, `plot()` and `woven_nystrom_error()`.

# woven 0.99.2

## Changes in behaviour and API

* The object returned by `woven()` is an S4 class (`woven`) instead of an S3
  list. `show`, `summary` and `plot` are S4 methods, and accessor generics
  were added for the most used slots.
* `woven()` accepts a `MultiAssayExperiment` as `X_list`. Each experiment
  becomes one modality, subjects are aligned through `sampleMap`, a subject
  absent from an experiment gets an all-NA row, and `Y` may be the name of a
  `colData` column.
* `plot()` on a fit returns the ggplot object visibly instead of printing it
  and returning it invisibly.
* `woven_nystrom_error()` referred to a non-existent field and always
  returned `NaN`; it now computes the error. It also passed a too-short `Y`
  to the internal refit; fixed.
* `woven_precompute()` uses BiocParallel instead of `parallel::mclapply`,
  which does not fork on Windows.

## Bug fixes and minor improvements

* `?woven` no longer describes a removed two-solver design.
* Fixed the "Nystrm" typo.
* Vignette introduction rewritten: defines "assay block", explains the
  Nystrom method, and compares with DIABLO concretely.
* (First available in 0.99.1) `woven()` exposes `ridge_w` and `screen_top`
  for ridge-regularized and feature-screened projection matrices.

# woven 0.99.1

* `woven()` standardizes each modality to zero mean and unit variance by
  default (`scale = TRUE`). Stored centers and scales are reapplied to new
  subjects in `woven_scores()` and `woven_predict()`. `scale = FALSE` keeps
  the previous raw-feature behaviour.

# woven 0.99.0

* Initial Bioconductor submission.
* Closed-form dual SUMCOR MCCA solver (`woven_mcca_dual()`) for any number
  of modalities V >= 2.
* Nystrom (linear out-of-sample) projection of block-missing subjects.
* Graph Laplacian regularization built on observed rows only, without
  imputation.
* Label-augmented cross-covariance for supervised integration.
* The method was benchmarked against DIABLO, MOFA2, IntegrAO and
  imputation followed by DIABLO, and applied to ADNI data, in the companion
  manuscript; the benchmark code is distributed with the manuscript, not
  with this package.
