#' Simulated three-modality example dataset
#'
#' A small simulated dataset (90 subjects, three modalities, three classes)
#' used in the help-page examples. All values are synthetic: three latent
#' factors with class-specific means are mapped to each modality by random
#' loadings, Gaussian noise is added, and each feature is standardized. In
#' the block-missing version, 30 subjects each lose one whole modality,
#' chosen at random (missing completely at random). For a real block-missing
#' example see the package vignette, which uses the TCGA adrenocortical
#' carcinoma data \code{miniACC} from MultiAssayExperiment.
#'
#' @format A list with three components:
#' \describe{
#'   \item{X_complete}{Named list of three numeric matrices with subjects in
#'     rows: \code{RNA} (90 x 25), \code{Methylation} (90 x 20) and
#'     \code{Proteomics} (90 x 15). No missing values.}
#'   \item{X_missing}{The same three matrices with 30 subjects (33\%) each
#'     missing one entire modality block (all-NA row). 60 subjects are
#'     observed in all three modalities.}
#'   \item{Y}{Character vector of 90 class labels, \code{"CN"},
#'     \code{"MCI"} and \code{"AD"} (30 each). The names are illustrative
#'     only; no clinical data were used.}
#' }
#' @source Simulated with \code{inst/scripts/make_woven_example.R}
#'   (\code{system.file("scripts", "make_woven_example.R", package =
#'   "woven")}), which reproduces \code{data/woven_example.rda} exactly.
#'   The data are released under the package license (MIT).
#' @examples
#' data(woven_example)
#' fit <- woven(woven_example$X_missing,
#'     Y = woven_example$Y, K = 3L,
#'     verbose = FALSE
#' )
#' fit # all 90 subjects are scored, including the 30 block-missing ones
"woven_example"
