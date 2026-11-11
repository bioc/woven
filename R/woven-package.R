#' woven: supervised multi-omics integration for block-missing data
#'
#' WOVEN (Weighted Omics View Embedding via Nystrom) learns a shared,
#' supervised latent space across several omics modalities when some
#' subjects are missing whole modalities. The model is fitted on the anchor
#' subjects (those observed in every modality), and every other subject is
#' then scored from the modalities it does have, without imputing the
#' missing ones.
#'
#' @section Main functions:
#' \itemize{
#'   \item [woven()]: fit a model from a list of matrices or a
#'     MultiAssayExperiment.
#'   \item [woven_precompute()]: build the per-modality graphs once for
#'     repeated fits (cross-validation, tuning).
#'   \item [woven_scores()]: latent scores for new subjects.
#'   \item [woven_predict()]: class predictions for new subjects.
#'   \item [woven_metrics()]: silhouette, Davies-Bouldin, NMI and ESS
#'     retention of a fit.
#'   \item [woven-accessors]: [latent_scores()], [projection_matrices()],
#'     [anchor_indices()], [singular_values()].
#' }
#'
#' @section Plots:
#' \code{plot()} on a fit ([plot,woven,missing-method]) shows the latent
#' space; [woven_plot_vip()], [woven_plot_loadings()] and
#' [woven_plot_variance()] show feature importance, loadings and variance
#' per dimension.
#'
#' @section Getting started:
#' See \code{vignette("woven-intro", package = "woven")}.
#'
#' @importFrom methods setClass setGeneric setMethod setValidity show
"_PACKAGE"
