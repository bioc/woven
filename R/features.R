# features.R -- which features drive the WOVEN latent space

# Variable importance in projection (VIP) for every feature of modality v:
# the feature's squared loadings (each dimension's loading vector scaled to
# unit length, as in PLS), weighted by each dimension's share of the squared
# singular values and scaled by the number of features, so the mean squared
# VIP is 1 within the modality. Comparable within a modality, not across
# modalities of different sizes.
.vip_scores <- function(fit, v) {
    W <- fit@W_list[[v]]
    nrm <- sqrt(colSums(W^2))
    nrm[nrm == 0] <- 1
    W <- sweep(W, 2L, nrm, "/")
    svals <- fit@singular_values
    wts <- svals^2 / sum(svals^2)
    vip <- sqrt(nrow(W) * rowSums(sweep(W^2, 2L, wts, "*")))
    names(vip) <- if (!is.null(rownames(W))) {
        rownames(W)
    } else {
        paste0("Feature_", seq_len(nrow(W)))
    }
    vip
}

# For each feature of modality v, the class with the highest mean, or NA
# when the fit carries no class means (fits made before woven 0.99.6).
.highest_class <- function(fit, v) {
    cm <- if (length(fit@class_means) >= v) fit@class_means[[v]] else NULL
    if (is.null(cm) || ncol(cm) == 0L) {
        return(rep(NA_character_, nrow(fit@W_list[[v]])))
    }
    out <- colnames(cm)[max.col(replace(cm, is.na(cm), -Inf),
        ties.method = "first"
    )]
    out
}

#' Top features of a WOVEN fit
#'
#' Ranks the features of each modality by their variable importance in
#' projection (VIP) and returns them as a table, with the class in which
#' each feature is highest. This is the tabular counterpart of
#' \code{plot(fit, type = "features")}.
#'
#' VIP combines a feature's loadings over all latent dimensions, weighting
#' each dimension by its share of the squared singular values, and is scaled
#' so its mean square is 1 within a modality (so values above 1 are above
#' average). It ranks features within a modality; it is not comparable
#' across modalities with different numbers of features. Features removed
#' by \code{screen_top} have zero loadings and a VIP of 0; because the mean
#' square is taken over all of a modality's features, the features kept
#' score higher when screening removes many of them.
#'
#' @param fit a woven object from [woven()]
#' @param n integer: number of features per modality (default 10). Use
#'   \code{Inf} for all.
#' @param modality \code{NULL} (default, all modalities), or one or more
#'   modality names or indices.
#' @return a data.frame with one row per feature and columns
#'   \code{modality}, \code{rank} (within modality), \code{feature},
#'   \code{vip}, and \code{highest_in} (the class with the highest mean
#'   value of the feature, on the standardized scale the model was fit on).
#' @examples
#' data(woven_example)
#' fit <- woven(woven_example$X_missing,
#'     Y = woven_example$Y, K = 3L,
#'     verbose = FALSE
#' )
#' woven_top_features(fit, n = 5)
#' woven_top_features(fit, n = 3, modality = "RNA")
#' @seealso [woven_plot_vip()], [woven_plot_loadings()]
#' @export
woven_top_features <- function(fit, n = 10L, modality = NULL) {
    stopifnot(methods::is(fit, "woven"))
    if (length(n) != 1L || is.na(n) || n < 1) {
        stop("n must be a single number >= 1.", call. = FALSE)
    }
    mods <- if (is.null(modality)) {
        seq_len(fit@V)
    } else {
        vapply(modality, function(m) .resolve_modality(m, fit), integer(1L))
    }
    rows <- lapply(mods, function(v) {
        vip <- .vip_scores(fit, v)
        hi <- .highest_class(fit, v)
        k <- min(n, length(vip))
        o <- order(vip, decreasing = TRUE)[seq_len(k)]
        data.frame(
            modality = fit@mod_names[v],
            rank = seq_len(k),
            feature = names(vip)[o],
            vip = unname(vip[o]),
            highest_in = hi[o],
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
}
