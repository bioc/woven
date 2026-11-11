# accessors.R -- S4 accessor generics for the "woven" class
#
# Named accessors are the supported way to read a fit in user code. Names are
# snake_case to match the rest of the API (woven(), woven_scores(), ...). The
# full slot list is documented at ?`woven-class`.

#' Accessors for a fitted WOVEN model
#'
#' Extract the main results from a \code{\linkS4class{woven}} fit.
#'
#' \itemize{
#'   \item \code{latent_scores()}: n x K matrix of consensus latent scores
#'     for every subject (anchors and block-missing subjects). Rows of
#'     subjects with no observed modality are \code{NA}. Row names come from
#'     the row names of the input matrices (or \code{colData} row names for
#'     a MultiAssayExperiment).
#'   \item \code{projection_matrices()}: named list of V projection
#'     matrices, each p_v x K (features by latent dimensions).
#'   \item \code{anchor_indices()}: integer row indices of the anchor
#'     subjects (observed in every modality) used to fit the model.
#'   \item \code{singular_values()}: length-K vector of supervised canonical
#'     correlations (eigenvalues of the dual problem), largest first.
#' }
#'
#' @param x a \code{\linkS4class{woven}} object.
#' @return See Details.
#' @name woven-accessors
#' @examples
#' data(woven_example)
#' fit <- woven(woven_example$X_complete,
#'     Y = woven_example$Y, K = 3L,
#'     verbose = FALSE
#' )
#' dim(latent_scores(fit))
#' names(projection_matrices(fit))
#' length(anchor_indices(fit))
#' singular_values(fit)
NULL

#' @rdname woven-accessors
#' @export
setGeneric("latent_scores", function(x) standardGeneric("latent_scores"))

#' @rdname woven-accessors
#' @export
setMethod("latent_scores", "woven", function(x) x@Z)

#' @rdname woven-accessors
#' @export
setGeneric(
    "projection_matrices",
    function(x) standardGeneric("projection_matrices")
)

#' @rdname woven-accessors
#' @export
setMethod("projection_matrices", "woven", function(x) x@W_list)

#' @rdname woven-accessors
#' @export
setGeneric("anchor_indices", function(x) standardGeneric("anchor_indices"))

#' @rdname woven-accessors
#' @export
setMethod("anchor_indices", "woven", function(x) x@anchor_idx)

#' @rdname woven-accessors
#' @export
setGeneric("singular_values", function(x) standardGeneric("singular_values"))

#' @rdname woven-accessors
#' @export
setMethod("singular_values", "woven", function(x) x@singular_values)
