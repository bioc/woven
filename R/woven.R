# woven.R  -- Top-level user API for WOVEN
#
# Typical workflow:
#   precomp <- woven_precompute(X_list)       # build Laplacians once (optional)
#   fit  <- woven(X_list, Y, K = 5, precomp = precomp)
#   plot(fit, labels = Y)                     # plot latent space
#   scores  <- woven_scores(fit, X_list_new)  # latent scores for new subjects
#   pred    <- woven_predict(fit, X_list_new) # class predictions

#' Pre-compute Laplacian graphs for reuse across multiple woven() calls
#'
#' Builds k-NN RBF Laplacians for each modality from the observed data.
#' Pass the result to \code{woven(..., precomp = precomp)} to avoid
#' rebuilding the graph on every call when fitting the same subjects
#' several times, for example to compare hyperparameter settings or run a
#' sensitivity analysis. The graphs cover exactly the subjects in
#' \code{X_list}, so they cannot be reused for fits on a subset of subjects
#' (such as cross-validation training folds); fit those without
#' \code{precomp}.
#'
#' @param X_list list of V matrices (n x p_v), or a
#'   \code{MultiAssayExperiment} (converted exactly as in [woven()], so the
#'   graphs line up with a later \code{woven(mae, ...)} call). Block-missing
#'   rows (all NA) are automatically excluded from the k-NN graph.
#' @param k_nn integer, number of nearest neighbours (default 10)
#' @param scale logical, standardize features to unit variance before building
#'   Laplacians (default TRUE); use the same value passed to [woven()].
#' @return list of V sparse Laplacian matrices (n x n), one per modality.
#'   Pass directly to the \code{precomp} argument of [woven()]. The list
#'   records the \code{scale}, \code{k_nn} and \code{n} it was built with
#'   (as attributes \code{"woven_scale"}, \code{"woven_k_nn"},
#'   \code{"woven_n"}), and [woven()] stops if its own \code{scale} or the
#'   number of subjects does not match.
#' @examples
#' set.seed(1)
#' n <- 40
#' K <- 2L
#' X1 <- matrix(rnorm(n * 8), n, 8)
#' X2 <- matrix(rnorm(n * 6), n, 6)
#' Y <- rep(1:2, each = n / 2)
#' miss <- matrix(runif(n * 2) < 0.3, n, 2)
#' for (i in which(rowSums(miss) == 2)) miss[i, sample(2, 1)] <- FALSE
#' X1[miss[, 1], ] <- NA
#' X2[miss[, 2], ] <- NA
#' precomp <- woven_precompute(list(X1, X2), k_nn = 10L)
#' fit <- woven(list(X1, X2), Y = Y, K = K, precomp = precomp)
#' @seealso [woven()]
#' @export
woven_precompute <- function(X_list, k_nn = 10L, scale = TRUE) {
    if (methods::is(X_list, "MultiAssayExperiment")) {
        X_list <- .mae_to_xlist(X_list)
    }
    if (!is.list(X_list)) {
        stop("X_list must be a list of matrices or a MultiAssayExperiment.",
            call. = FALSE
        )
    }
    # Build the graph on the same features woven() will fit on: with
    # scale = TRUE, woven() standardizes each modality before building its
    # Laplacians, and .scale_fit() is deterministic, so the graphs match.
    X_fit <- if (isTRUE(scale)) .scale_fit(X_list)$X else X_list
    out <- .build_laplacians(X_fit, k_nn)
    names(out) <- names(X_list)
    attr(out, "woven_scale") <- isTRUE(scale)
    attr(out, "woven_k_nn") <- as.integer(k_nn)
    attr(out, "woven_n") <- nrow(X_list[[1L]])
    out
}

# One Laplacian per modality. A serial loop: each graph is one k-NN search on
# n_v x p_v data, measured at ~0.5 s for a 300 x 5,000 modality, so parallel
# workers would save little and add startup and memory overhead. Kept in one
# helper so a parallel backend could be swapped in without API changes.
.build_laplacians <- function(X_list, k_nn) {
    lapply(X_list, function(X) build_laplacian(X, k = k_nn))
}

# Check a woven_precompute() result against the data and settings of a
# woven() call.
.check_precomp <- function(precomp, V, n, scale, k_nn) {
    if (!is.list(precomp) || length(precomp) != V) {
        stop(sprintf(
            "precomp must be the output of woven_precompute() -- a list of %d Laplacians.",
            V
        ), call. = FALSE)
    }
    dims_ok <- vapply(precomp, function(L) all(dim(L) == c(n, n)), logical(1L))
    if (!all(dims_ok)) {
        stop(
            sprintf(
                "precomp Laplacians must be %d x %d (one row per subject); ",
                n, n
            ), "rebuild them with woven_precompute() on the same X_list.",
            call. = FALSE
        )
    }
    pc_scale <- attr(precomp, "woven_scale")
    if (is.null(pc_scale)) {
        warning("precomp does not record how it was scaled; make sure it was ",
            "built with woven_precompute(scale = ", isTRUE(scale), ").",
            call. = FALSE
        )
    } else if (!identical(pc_scale, isTRUE(scale))) {
        fmt <- paste0(
            "precomp was built with scale = %s but woven() was called ",
            "with scale = %s. Rebuild it with ",
            "woven_precompute(X_list, scale = %s)."
        )
        stop(sprintf(
            fmt,
            pc_scale, isTRUE(scale), isTRUE(scale)
        ), call. = FALSE)
    }
    pc_k <- attr(precomp, "woven_k_nn")
    if (!is.null(pc_k) && !identical(as.integer(k_nn), pc_k)) {
        warning(sprintf(
            "k_nn = %d is ignored: precomp was built with k_nn = %d.",
            as.integer(k_nn), pc_k
        ), call. = FALSE)
        return(pc_k)
    }
    as.integer(k_nn)
}

# Convert new-subject input for woven_scores(): a MultiAssayExperiment is
# turned into a list ordered like the training modalities, with features
# matched by name when both sides have names.
.as_new_xlist <- function(fit, X_new) {
    if (!methods::is(X_new, "MultiAssayExperiment")) {
        return(.as_new_list(fit, X_new))
    }
    xl <- .mae_to_xlist(X_new)
    miss <- setdiff(fit@mod_names, names(xl))
    if (length(miss) > 0L) {
        stop(sprintf(
            "The MultiAssayExperiment lacks experiment(s) used in the fit: %s.",
            paste(miss, collapse = ", ")
        ), call. = FALSE)
    }
    xl <- xl[fit@mod_names]
    for (v in seq_along(xl)) {
        f_fit <- rownames(fit@W_list[[v]])
        if (!is.null(f_fit) && !is.null(colnames(xl[[v]])) &&
            all(f_fit %in% colnames(xl[[v]]))) {
            xl[[v]] <- xl[[v]][, f_fit, drop = FALSE]
        }
    }
    xl
}

# New data given as a list: accept data frames, match blocks by row name
# (subject ID) when their subjects differ, and when the list is named with
# the fit's modality names, put the blocks in fit order and treat any
# modality left out as missing for every new subject. Unnamed lists are
# taken positionally and checked by .check_new_data().
.as_new_list <- function(fit, X_new) {
    if (!is.list(X_new) || is.data.frame(X_new)) return(X_new)
    nm <- names(X_new)
    xl <- lapply(seq_along(X_new), function(v) {
        lab <- if (!is.null(nm) && nzchar(nm[v])) {
            sprintf("'%s'", nm[v])
        } else {
            sprintf("[[%d]]", v)
        }
        .as_block_matrix(X_new[[v]], lab)
    })
    names(xl) <- nm
    xl <- .align_blocks(xl)$X_list
    by_name <- !is.null(nm) && all(nzchar(nm)) && !anyDuplicated(nm)
    if (!by_name) return(xl)
    unknown <- setdiff(nm, fit@mod_names)
    if (length(unknown) > 0L) {
        stop(sprintf(
            "Unknown modality name(s) in new data: %s. The fit has: %s.",
            paste(unknown, collapse = ", "),
            paste(fit@mod_names, collapse = ", ")
        ), call. = FALSE)
    }
    n_new <- nrow(xl[[1L]])
    ids <- rownames(xl[[1L]])
    full <- lapply(seq_len(fit@V), function(v) {
        m <- fit@mod_names[v]
        if (m %in% names(xl)) return(xl[[m]])
        matrix(NA_real_, n_new, nrow(fit@W_list[[v]]),
            dimnames = list(ids, rownames(fit@W_list[[v]]))
        )
    })
    names(full) <- fit@mod_names
    for (v in seq_along(full)) {
        f_fit <- rownames(fit@W_list[[v]])
        if (!is.null(f_fit) && !is.null(colnames(full[[v]])) &&
            all(f_fit %in% colnames(full[[v]]))) {
            full[[v]] <- full[[v]][, f_fit, drop = FALSE]
        }
    }
    full
}

#' Fit a supervised WOVEN model
#'
#' Learns a shared supervised latent space across V omics modalities, handling
#' block-missing data via anchor-restricted alignment and Nystrom projection.
#' Labels Y are required  -- WOVEN is a supervised method (cf. DIABLO).
#'
#' All V use the same closed-form dual SUMCOR MCCA solver: a single
#' eigendecomposition of a block matrix built from the per-modality
#' Laplacians and the label kernel, with no iterative optimization. For V=2
#' this is equivalent to dual supervised CCA; for V>=3 it is the global
#' optimum of the relaxed SUMCOR trace objective (see \code{woven_mcca_dual}).
#'
#' @param X_list the data, one block per modality (subjects in rows,
#'   features in columns), as either
#'   \itemize{
#'     \item a named list of numeric matrices or data frames. When every
#'       block has row names (subject IDs), the blocks may cover different
#'       subjects in any order: woven matches them by ID and treats a
#'       subject absent from a block as missing that modality. Without row
#'       names, all blocks must have the same rows in the same order, with a
#'       subject missing a modality given an all-NA row; or
#'     \item a \code{MultiAssayExperiment}.
#'   } When a MultiAssayExperiment is
#'   supplied, each experiment becomes one modality, subjects are aligned to
#'   \code{colData(X_list)} via \code{sampleMap(X_list)}, and a subject not
#'   present in a given experiment gets an all-NA row for it automatically.
#' @param Y class labels for all subjects (factor, character or integer).
#'   Plots and tables list the classes in the order of
#'   \code{levels(as.factor(Y))}; pass a factor to choose the order.
#'   Either in subject order (the row order of the blocks, or of
#'   \code{colData(X_list)} for a MultiAssayExperiment) or a named vector
#'   whose names are the subject IDs, in which case order does not matter.
#'   A single string is instead treated as the name of a \code{colData}
#'   column to use as labels, but only when \code{X_list} is a
#'   MultiAssayExperiment. Only anchor subjects' labels enter the supervised
#'   objective.
#' @param anchor_idx integer vector  -- indices of fully-observed subjects
#'   (observed in all V modalities). Must have length >= K. If NULL (default),
#'   anchors are detected automatically as subjects with no block-missing
#'   modalities. A supplied \code{anchor_idx} is checked: indices must be
#'   distinct, within \code{1..n}, and observed in every modality; offending
#'   subjects are named in the error.
#' @param K integer  -- number of latent dimensions (default 5). Must be
#'   smaller than n and no larger than the number of anchors. Each
#'   modality's projection matrix has rank at most \code{min(n_a, p_v)}, so
#'   a warning is given if K exceeds a modality's number of features.
#' @param lambdas numeric scalar or length-V vector  -- Laplacian regularization
#'   strength per modality (default 0.1 for all)
#' @param gamma_y numeric >= 0  -- supervision strength. 0 = unsupervised CCA.
#'   Default 1.0 (equal weight to cross-modal alignment and label alignment).
#'   Tune via cross-validation on anchor set if labels are noisy.
#' @param k_nn integer  -- k-nearest-neighbors for Laplacian graph (default 10).
#' @param ridge_w numeric >= 0  -- ridge shrinkage applied to the per-modality
#'   Gram matrix during projection-matrix (W) recovery (default 1e-8, a purely
#'   numerical-stability nudge with no regularizing effect). Increase (e.g.
#'   0.1-5.0) to control capacity when p_v >> n_a and the true signal is
#'   diffuse across many features -- select it with [woven_tune()] (held-out
#'   classification error); silhouette is insensitive to this parameter even
#'   where it materially changes classification accuracy.
#' @param screen_top integer or Inf  -- restricts W recovery to the top
#'   \code{screen_top} features per modality by ANOVA F-statistic against Y,
#'   computed on anchor subjects only (default Inf = no screening, uses all
#'   features). A genuine re-fit on the reduced feature set, not a post-hoc
#'   truncation. Complementary to \code{ridge_w}; both target overfitting when
#'   p_v is large relative to n_a. Select jointly with \code{ridge_w} via
#'   cross-validation, not independently (the two interact).
#' @param scale logical  -- standardize each modality's features to zero mean
#'   and unit variance before fitting (default TRUE), so no modality
#'   dominates by measurement scale. Apply biological transforms
#'   (CLR/log/VST) upstream.
#'   The center and scale are stored and reapplied to new subjects in
#'   [woven_scores()]. When \code{precomp} is supplied it must have been
#'   built with the same \code{scale} value; a mismatch is an error.
#' @param precomp optional output of [woven_precompute()] -- pre-built
#'   Laplacians for the same \code{X_list}. Pass this when calling
#'   \code{woven()} multiple times on the same subjects (e.g. comparing
#'   hyperparameter settings) to avoid rebuilding the graph each time. Not
#'   usable for fits on a subset of subjects, such as cross-validation
#'   training folds.
#'   When supplied, \code{k_nn} is taken from \code{precomp}.
#' @param verbose logical  -- report progress (default TRUE): a summary of
#'   the input, one line per step with timings, and a closing summary of
#'   coverage and top features. All output is via \code{message()}.
#' @param plot logical  -- draw the latent space when the fit finishes
#'   (default: in interactive sessions only). Needs ggplot2.
#'
#' @return an object of S4 class \code{\linkS4class{woven}}. Key slots:
#'   \describe{
#'     \item{Z}{n x K matrix of consensus latent scores for ALL n subjects
#'       (anchors and block-missing). The primary output for downstream
#'       analysis.}
#'     \item{W_list}{list of V projection matrices, each p_v x K}
#'     \item{Z_anchors}{list of V anchor latent score matrices, each n_a x K}
#'     \item{singular_values}{K supervised canonical correlations}
#'     \item{anchor_idx}{indices of anchor (fully-observed) subjects}
#'     \item{Y_levels}{integer codes of the classes used during fitting}
#'     \item{Y_labels}{class names, \code{levels(as.factor(Y))}; for
#'       numeric \code{Y} these are the numbers as strings}
#'     \item{K, V, n}{dimensions}
#'     \item{lambdas, gamma_y}{hyperparameters}
#'   }
#'   See \code{\linkS4class{woven}} for the complete slot list.
#'
#' @examples
#' set.seed(1)
#' n <- 60
#' p1 <- 20
#' p2 <- 15
#' K <- 2
#' Y <- rep(1:2, each = n / 2)
#' X1 <- matrix(rnorm(n * p1), n, p1)
#' X2 <- matrix(rnorm(n * p2), n, p2)
#' # 30% block missingness; enforce >= 1 view per subject
#' miss <- matrix(runif(n * 2) < 0.3, n, 2)
#' for (i in which(rowSums(miss) == 2)) miss[i, sample(2, 1)] <- FALSE
#' X1[miss[, 1], ] <- NA
#' X2[miss[, 2], ] <- NA
#' # anchor_idx auto-detected from NA pattern -- no need to specify
#' fit <- woven(list(X1, X2), Y = Y, K = K)
#' dim(latent_scores(fit)) # 60 x 2 -- all subjects scored
#'
#' # Most often each assay is its own table covering its own subjects.
#' # Give the tables row names (subject IDs) and woven matches them:
#' ids <- paste0("S", seq_len(n))
#' rna <- matrix(rnorm(50 * p1), 50, p1, dimnames = list(ids[1:50], NULL))
#' prot <- matrix(rnorm(45 * p2), 45, p2, dimnames = list(ids[16:60], NULL))
#' labels <- setNames(Y, ids)
#' fit_tab <- woven(list(RNA = rna, Protein = prot), Y = labels, K = K)
#'
#' # Equivalently, from a MultiAssayExperiment: each experiment is one
#' # modality, and subjects missing an experiment need no special handling.
#' library(MultiAssayExperiment)
#' ids <- paste0("S", seq_len(n))
#' rownames(X1) <- ids
#' rownames(X2) <- ids
#' present1 <- ids[!is.na(X1[, 1])]
#' present2 <- ids[!is.na(X2[, 1])]
#' exps <- ExperimentList(list(
#'     Mod1 = t(X1[present1, , drop = FALSE]),
#'     Mod2 = t(X2[present2, , drop = FALSE])
#' ))
#' cd <- DataFrame(row.names = ids, group = Y)
#' mae <- MultiAssayExperiment(exps, colData = cd)
#' fit_mae <- woven(mae, Y = "group", K = K)
#' dim(latent_scores(fit_mae)) # same 60 x 2 -- all subjects scored
#'
#' @seealso [woven_scores()], [woven_predict()], [woven_metrics()],
#'   [woven_precompute()], [woven-accessors]
#' @export
woven <- function(
        X_list, Y, anchor_idx = NULL, K = 5L, lambdas = 0.1, gamma_y = 1.0,
        k_nn = 10L, ridge_w = 1e-8, screen_top = Inf, scale = TRUE,
        precomp = NULL, verbose = TRUE, plot = interactive()) {
    t_start <- proc.time()[["elapsed"]]
    inp <- .woven_check_input(X_list, Y, anchor_idx, K, verbose)
    X_list <- inp$X_list
    anchor_idx <- inp$anchor_idx
    V <- inp$V
    n <- inp$n

    Y_fct <- as.factor(inp$Y)
    Y_labels <- levels(Y_fct) # original label names e.g. "CN","MCI","Dementia"
    Y <- as.integer(Y_fct)

    if (length(lambdas) == 1L) lambdas <- rep(lambdas, V)
    stopifnot(length(lambdas) == V)

    if (!is.null(precomp)) {
        k_nn <- .check_precomp(precomp, V, n, scale, k_nn)
    }
    if (verbose) .woven_header(X_list, Y_fct, anchor_idx, inp$aligned)

    #  Standardize features to unit variance so no modality dominates by
    #  measurement scale. Stored for out-of-sample scaling in woven_scores().
    #  The transform is deterministic, so a precomp built with the same
    #  scale value matches the graph woven() would build itself.
    scale_center <- NULL
    scale_scale <- NULL
    if (isTRUE(scale)) {
        .sf <- .scale_fit(X_list)
        X_list <- .sf$X
        scale_center <- .sf$center
        scale_scale <- .sf$scale
    }

    #  Step 1: per-modality similarity graphs (built on all subjects observed
    #  in that modality; the solver uses the anchor rows). Building them here
    #  rather than inside the solver gives the same graphs and lets each
    #  modality report progress.
    if (is.null(precomp)) {
        precomp <- .build_laplacians_verbose(X_list, k_nn, verbose)
    } else if (verbose) {
        message("[1/3] Similarity graphs: using precomputed graphs")
    }
    La_list_precomp <- lapply(precomp, function(L) {
        as.matrix(L[anchor_idx, anchor_idx, drop = FALSE])
    })

    #  Step 2: closed-form dual SUMCOR MCCA for all V
    t0 <- proc.time()[["elapsed"]]
    raw <- woven_mcca_dual(
        X_list          = X_list,
        anchor_idx      = anchor_idx,
        Y               = Y,
        K               = K,
        lambdas         = lambdas,
        gamma_y         = gamma_y,
        k_nn            = k_nn,
        ridge_w         = ridge_w,
        screen_top      = screen_top,
        La_list_precomp = La_list_precomp,
        verbose         = FALSE
    )
    if (verbose) {
        message(sprintf(
            "[2/3] Closed-form alignment on %d anchors, K = %d (%s)",
            length(anchor_idx), raw$K, .elapsed(t0)
        ))
    }

    #  Step 3: score every subject from its observed modalities
    t0 <- proc.time()[["elapsed"]]
    fit <- .woven_assemble(
        raw, X_list, Y, Y_labels, anchor_idx, lambdas, gamma_y,
        k_nn, ridge_w, screen_top, scale, scale_center, scale_scale,
        Y_fct, verbose = FALSE
    )
    if (verbose) {
        message(sprintf(
            "[3/3] Scored subjects from their available modalities (%s)",
            .elapsed(t0)
        ))
        .woven_footer(fit, t_start)
    }
    if (isTRUE(plot) && requireNamespace("ggplot2", quietly = TRUE)) {
        methods::show(.plot_scores(fit))
    }
    fit
}

# Input handling for woven(): convert a MultiAssayExperiment or a list of
# matrices/data frames (aligned by row names when their subjects differ),
# validate the modalities and labels, and settle the anchor set. Returns the
# aligned list of matrices with n, V, the (checked) anchor indices, the
# labels and whether row-name alignment was needed.
.woven_check_input <- function(X_list, Y, anchor_idx, K, verbose) {
    #  MultiAssayExperiment input: convert to the internal X_list
    #  representation, then proceed exactly as with a plain list.
    if (methods::is(X_list, "MultiAssayExperiment")) {
        mae <- X_list
        X_list <- .mae_to_xlist(mae)
        if (is.character(Y) && length(Y) == 1L) {
            Y <- MultiAssayExperiment::colData(mae)[[Y]]
            names(Y) <- rownames(MultiAssayExperiment::colData(mae))
        }
    }

    #  Input validation
    if (!is.list(X_list) || is.data.frame(X_list)) {
        stop(
            "X_list must be a list of matrices, ",
            "e.g. list(RNA = X_rna, Methyl = X_meth), ",
            "or a MultiAssayExperiment."
        )
    }
    V <- length(X_list)
    if (V < 2L) {
        stop(sprintf("X_list must contain at least 2 modalities; got %d.", V))
    }
    nm_in <- names(X_list)
    X_list <- lapply(seq_len(V), function(v) {
        lab <- if (!is.null(nm_in) && nzchar(nm_in[v])) {
            sprintf("'%s'", nm_in[v])
        } else {
            sprintf("[[%d]]", v)
        }
        .as_block_matrix(X_list[[v]], lab)
    })
    names(X_list) <- nm_in
    al <- .align_blocks(X_list)
    X_list <- al$X_list
    n_rows <- vapply(X_list, nrow, integer(1L))
    if (length(unique(n_rows)) > 1L) {
        nm <- if (!is.null(names(X_list))) {
            names(X_list)
        } else {
            paste0("[[", seq_len(V), "]]")
        }
        fmt <- paste0(
            "All matrices must have the same number of rows (subjects), ",
            "or row names (subject IDs) so woven can match them.\n",
            "  Row counts: %s"
        )
        stop(sprintf(fmt, paste(sprintf("%s=%d", nm, n_rows), collapse = ", ")))
    }
    n <- n_rows[1L]
    X_list <- .name_blocks(X_list)
    Y <- .align_labels(Y, rownames(X_list[[1L]]), n)
    if (anyNA(Y)) {
        stop("Y contains NA values. Every subject must have a class label.")
    }
    if (K < 1L) stop("K must be >= 1.")
    if (K >= n) stop(sprintf("K=%d must be less than n=%d.", K, n))

    #  Auto-detect anchors (subjects with all modalities observed)
    if (is.null(anchor_idx)) {
        block_missing <- vapply(X_list, .block_missing_rows, logical(n))
        anchor_idx <- which(rowSums(block_missing) == 0L)
        if (length(anchor_idx) == 0L) {
            stop(
                "No fully-observed subjects found. Every subject is ",
                "block-missing in at least one modality."
            )
        }
    }

    anchor_idx <- .check_anchor_idx(X_list, anchor_idx)
    if (length(anchor_idx) < K) {
        stop(sprintf(
            "Need at least K=%d anchor subjects; got %d.", K, length(anchor_idx)
        ))
    }

    list(
        X_list = X_list, Y = Y, anchor_idx = anchor_idx, V = V, n = n,
        aligned = al$aligned
    )
}

# Build the S4 fit from the solver output: name the projection matrices,
# score every subject from its observed modalities and fill the slots.
.woven_assemble <- function(
        raw, X_list, Y, Y_labels, anchor_idx, lambdas, gamma_y, k_nn, ridge_w,
        screen_top, scale, scale_center, scale_scale, Y_fct = NULL,
        verbose = FALSE) {
    V <- length(X_list)
    n <- nrow(X_list[[1L]])
    W_list <- raw$W_list
    Z_anchors <- raw$Za_list
    K <- raw$K
    svals <- raw$singular_values
    fit_mcca <- raw

    #  Propagate feature names and modality names onto W_list
    mod_names <- if (!is.null(names(X_list))) {
        names(X_list)
    } else {
        paste0("Modality ", seq_len(V))
    }
    for (v in seq_len(V)) {
        rownames(W_list[[v]]) <- colnames(X_list[[v]])
    }
    names(W_list) <- mod_names

    #  Latent scores for all n subjects from their observed modalities
    if (verbose) message("  Projecting all subjects...")
    Z_all <- .project_views(X_list, W_list, K)
    if (verbose) {
        n_scored <- sum(!is.na(Z_all[, 1L]))
        message(sprintf(
            "  Done. %d / %d subjects scored (%.0f%% retained).",
            n_scored, n, 100 * n_scored / n
        ))
    }

    #  Bookkeeping for plots and tables: which modalities each subject has,
    #  each subject's label, and each feature's mean per class (on the scale
    #  the model was fit on) for colouring features by their highest class.
    observed <- vapply(X_list, function(X) !.block_missing_rows(X), logical(n))
    observed <- matrix(observed, nrow = n,
        dimnames = list(rownames(X_list[[1L]]), mod_names)
    )
    if (is.null(Y_fct)) Y_fct <- factor(Y_labels[Y], levels = Y_labels)
    names(Y_fct) <- rownames(X_list[[1L]])
    class_means <- lapply(seq_len(V), function(v) {
        .class_means(X_list[[v]], Y_fct, observed[, v])
    })
    names(class_means) <- mod_names

    #  Return unified S4 object
    methods::new("woven",
        Z               = Z_all,
        W_list          = W_list,
        Z_anchors       = Z_anchors,
        singular_values = svals,
        anchor_idx      = as.integer(anchor_idx),
        Y_anchor        = Y[anchor_idx],
        Y_levels        = sort(unique(Y)),
        Y_labels        = Y_labels,
        K               = as.integer(K),
        V               = as.integer(V),
        n               = as.integer(n),
        mod_names       = mod_names,
        lambdas         = lambdas,
        gamma_y         = gamma_y,
        k_nn            = as.integer(k_nn),
        ridge_w         = ridge_w,
        screen_top      = screen_top,
        scaled          = isTRUE(scale),
        scale_center    = scale_center,
        scale_scale     = scale_scale,
        fit_mcca        = fit_mcca,
        Y               = Y_fct,
        observed        = observed,
        class_means     = class_means
    )
}

#' Show method for WOVEN fit
#'
#' Called automatically when a \code{\linkS4class{woven}} object is printed
#' or auto-printed at the console.
#'
#' @param object a woven object from [woven()]
#' @return Invisibly returns \code{NULL}; called for its side effect.
#' @examples
#' set.seed(1)
#' n <- 20
#' K <- 2L
#' X1 <- matrix(rnorm(n * 5), n, 5)
#' X2 <- matrix(rnorm(n * 4), n, 4)
#' Y <- rep(1:2, each = n / 2)
#' miss <- matrix(FALSE, n, 2)
#' miss[c(15, 17, 19), 1] <- TRUE
#' miss[c(16, 18, 20), 2] <- TRUE
#' X1[miss[, 1], ] <- NA
#' X2[miss[, 2], ] <- NA
#' anchor_idx <- which(rowSums(miss) == 0)
#' fit <- woven(list(X1, X2), Y = Y, anchor_idx = anchor_idx, K = K)
#' fit
#' @export
setMethod("show", "woven", function(object) {
    x <- object
    n_anchor <- length(x@anchor_idx)
    n_scored <- sum(!is.na(x@Z[, 1L]))
    cat("WOVEN fit\n")
    cat(sprintf(
        "  Modalities : %d    Subjects: %d    Dimensions: %d\n",
        x@V, x@n, x@K
    ))
    cat(sprintf(
        "  Anchors    : %d (%.0f%%)    Scored: %d (%.0f%%)\n",
        n_anchor, 100 * n_anchor / x@n,
        n_scored, 100 * n_scored / x@n
    ))
    cat(sprintf(
        "  Solver     : mcca_dual (closed-form, globally optimal)\n"
    ))
    cat(sprintf(
        "  gamma_y    : %.2f    lambda: %s    k_nn: %d\n",
        x@gamma_y,
        paste(round(x@lambdas, 3), collapse = "/"),
        x@k_nn
    ))
    cat(sprintf(
        "  Singular values: %s\n",
        paste(round(x@singular_values[seq_len(min(5L, x@K))], 3),
            collapse = ", "
        )
    ))
    if (!is.null(x@Y_labels)) {
        cat(sprintf("  Classes    : %s\n", paste(x@Y_labels, collapse = ", ")))
    }
    if (!is.null(x@mod_names)) {
        cat(sprintf("  Modalities : %s\n", paste(x@mod_names, collapse = ", ")))
    }
    cat("\n")
    cat("  -- Next steps --\n")
    cat("  plot(fit)                       # latent space by class\n")
    cat("  plot(fit, type = \"features\")    # top features per modality\n")
    cat("  plot(fit, type = \"coverage\")    # who has which modality\n")
    cat("  woven_top_features(fit)         # top features as a table\n")
    cat("  woven_tune(X_list, Y)           # choose settings by CV\n")
    cat("  woven_metrics(fit)              # silhouette, NMI, ESS\n")
    cat("  woven_predict(fit, new_data)    # classify new subjects\n")
    invisible(NULL)
})

#' Summarise a WOVEN fit
#'
#' Prints a compact metrics table: silhouette, Davies-Bouldin, NMI, and ESS
#' for the scored subjects, using the class labels stored in the fit unless
#' \code{labels} is given.
#'
#' @param object a woven object from [woven()]
#' @param labels character, factor, or integer vector of length n.
#'   Defaults to the labels the model was fit with.
#' @param ... unused
#' @return Invisibly returns a named numeric vector of metrics.
#' @export
#' @examples
#' data(woven_example)
#' fit <- woven(woven_example$X_complete, Y = woven_example$Y, K = 3L)
#' summary(fit)
setMethod("summary", "woven", function(object, labels = NULL, ...) {
    n_scored <- sum(!is.na(object@Z[, 1L]))
    if (is.null(labels) && length(object@Y) == object@n) labels <- object@Y
    m <- if (!is.null(labels)) woven_metrics(object, labels) else NULL
    show(methods::new("wovenSummary",
        V = length(object@W_list),
        K = object@K,
        n = object@n,
        mod_names = object@mod_names,
        y_labels = object@Y_labels,
        n_scored = n_scored,
        singular_values = object@singular_values,
        metrics = m
    ))
    invisible(m)
})

#' Print a summary(woven) report
#'
#' Prints the report built by \code{\link{summary,woven-method}}. All
#' \code{cat()} calls for the summary display live here, inside a
#' \code{show} method, rather than inside \code{summary()} itself.
#'
#' @param object a \code{\linkS4class{wovenSummary}} object, as returned
#'   internally by \code{summary()} on a \code{\linkS4class{woven}} fit.
#' @return Invisibly returns \code{NULL}; called for its printed side effect.
#' @keywords internal
setMethod("show", "wovenSummary", function(object) {
    cat(sprintf(
        "WOVEN fit  [V=%d  K=%d  n=%d]\n",
        object@V, object@K, object@n
    ))
    cat(sprintf(
        "  Modalities : %s\n",
        paste(object@mod_names %||% seq_len(object@V), collapse = ", ")
    ))
    cat(sprintf(
        "  Classes    : %s\n",
        paste(object@y_labels %||% "unknown", collapse = ", ")
    ))
    cat(sprintf(
        "  Scored     : %d / %d  (ESS = %.2f)\n",
        object@n_scored, object@n, object@n_scored / object@n
    ))
    cat(sprintf(
        "  Singular values: %s\n",
        paste(round(object@singular_values[seq_len(min(5L, object@K))], 3),
            collapse = ", "
        )
    ))
    if (!is.null(object@metrics)) {
        m <- object@metrics
        cat(sprintf("\n  Silhouette : %6.3f\n", m["Silhouette"]))
        cat(sprintf("  Davies-Bouldin : %6.3f\n", m["Davies-Bouldin"]))
        cat(sprintf("  NMI        : %6.3f\n", m["NMI"]))
        cat(sprintf("  ESS        : %6.3f\n", m["ESS"]))
    } else {
        cat("  (Pass labels = Y to compute silhouette / NMI / ESS)\n")
    }
    invisible(NULL)
})

#' Plot a WOVEN fit
#'
#' One entry point for every WOVEN plot. The default draws the latent space:
#' two latent dimensions with subjects coloured by class (the labels stored
#' in the fit, unless \code{labels} is given), anchors (subjects observed in
#' every modality, used to learn the projections) as filled points, and
#' subjects scored from a subset of modalities as open points. The other
#' types call the matching \code{woven_plot_*()} function.
#'
#' @param x a woven object from [woven()]
#' @param y unused; present for S4 `plot` signature compatibility
#' @param type which plot: \code{"scores"} (default, the latent space),
#'   \code{"features"} (top features per modality, [woven_plot_vip()]),
#'   \code{"loadings"} ([woven_plot_loadings()]), \code{"coverage"} (which
#'   subjects have which modalities, [woven_plot_coverage()]), or
#'   \code{"variance"} ([woven_plot_variance()]).
#' @param labels for \code{type = "scores"}: a vector of length n to colour
#'   points by. Defaults to the class labels used to fit the model. A
#'   numeric vector with more than 10 distinct values (for example age) is
#'   drawn as a colour gradient.
#' @param dims for \code{type = "scores"}: integer vector of length 2, the
#'   latent dimensions to plot (default \code{c(1, 2)}). A fit with
#'   \code{K = 1} is drawn as the distribution of dimension 1 per class.
#' @param highlight_anchors for \code{type = "scores"}: logical, draw
#'   anchors filled and other subjects open (default TRUE)
#' @param ... passed to the \code{woven_plot_*()} function for the other
#'   types, e.g. \code{n_top} or \code{modality}.
#' @return a ggplot object, printed automatically when called at the top
#'   level (standard R plot-method behaviour).
#' @examples
#' data(woven_example)
#' fit <- woven(woven_example$X_missing,
#'     Y = woven_example$Y, K = 3L,
#'     verbose = FALSE
#' )
#' plot(fit)
#' plot(fit, type = "features", n_top = 8)
#' plot(fit, type = "coverage")
#' @export
setMethod("plot", signature(x = "woven", y = "missing"), function(
        x, y,
        type = c("scores", "features", "loadings", "coverage", "variance"),
        labels = NULL, dims = c(1L, 2L), highlight_anchors = TRUE, ...) {
    type <- match.arg(type)
    switch(type,
        scores = .plot_scores(x, labels, dims, highlight_anchors),
        features = woven_plot_vip(x, ...),
        loadings = woven_plot_loadings(x, ...),
        coverage = woven_plot_coverage(x, ...),
        variance = woven_plot_variance(x, ...)
    )
})

# Latent-space scatter for plot(fit, type = "scores").
.plot_scores <- function(
        x, labels = NULL, dims = c(1L, 2L), highlight_anchors = TRUE) {
    .require_ggplot2()
    Z <- x@Z
    if (x@K == 1L) dims <- c(1L, 1L)
    if (length(dims) != 2L || any(dims < 1L) || any(dims > x@K)) {
        stop(sprintf(
            "dims must be two integers in 1..%d (the fit has K=%d dimensions).",
            x@K, x@K
        ), call. = FALSE)
    }
    if (is.null(labels) && length(x@Y) == x@n) labels <- x@Y
    if (!is.null(labels) && length(labels) != x@n) {
        stop(sprintf(
            "labels must have one entry per subject (%d); got %d.",
            x@n, length(labels)
        ), call. = FALSE)
    }
    # A numeric label with many distinct values (age, a score) is drawn as
    # a colour gradient; anything else as classes.
    continuous <- is.numeric(labels) && length(unique(labels)) > 10L

    is_anchor <- seq_len(x@n) %in% x@anchor_idx
    df <- data.frame(
        z1 = Z[, dims[1L]],
        z2 = Z[, dims[2L]],
        anchor = is_anchor
    )
    df$group <- if (is.null(labels)) {
        factor("all")
    } else if (continuous) {
        as.numeric(labels)
    } else {
        as.factor(labels)
    }
    df <- df[!is.na(df$z1), , drop = FALSE]

    # Point size and transparency scale with the number of subjects.
    n_pts <- nrow(df)
    size <- if (n_pts <= 200L) 2.6 else if (n_pts <= 1000L) 1.8 else 1.1
    alpha <- if (n_pts <= 200L) 0.9 else if (n_pts <= 1000L) 0.7 else 0.45
    hollow <- highlight_anchors && any(!df$anchor)
    df$shape <- if (hollow) {
        factor(ifelse(df$anchor, "Anchor", "Partial data"),
            levels = c("Anchor", "Partial data")
        )
    } else {
        factor("Anchor")
    }

    if (x@K == 1L) {
        # One dimension: show its distribution per class.
        p_out <- ggplot2::ggplot(df, ggplot2::aes(
            x = z1, y = if (is.factor(group)) group else "",
            color = group, shape = shape
        )) +
            ggplot2::geom_boxplot(
                data = df[df$anchor, , drop = FALSE], color = "gray55",
                fill = NA, outlier.shape = NA, width = 0.5,
                show.legend = FALSE, inherit.aes = FALSE,
                mapping = ggplot2::aes(
                    x = z1, y = if (is.factor(group)) group else ""
                )
            ) +
            ggplot2::geom_jitter(
                size = size, alpha = alpha, stroke = 0.7,
                height = 0.15, width = 0
            ) +
            ggplot2::labs(x = "WOVEN Dimension 1", y = NULL)
    } else {
        p_out <- ggplot2::ggplot(df, ggplot2::aes(
            x = z1, y = z2, color = group, shape = shape
        )) +
            ggplot2::geom_point(size = size, alpha = alpha, stroke = 0.7)
        # Ellipses for classes with enough anchors to estimate one.
        if (is.factor(df$group) && !is.null(labels)) {
            cnt <- table(df$group[df$anchor])
            ok <- names(cnt)[cnt >= 5L]
            df_ell <- df[df$anchor & df$group %in% ok, , drop = FALSE]
            if (nrow(df_ell) > 0L && length(ok) <= 8L) {
                p_out <- p_out + ggplot2::stat_ellipse(
                    data = df_ell, ggplot2::aes(group = group),
                    type = "norm", level = 0.68, linewidth = 0.7,
                    linetype = "dashed", show.legend = FALSE,
                    inherit.aes = TRUE
                )
            }
        }
        p_out <- p_out + ggplot2::labs(
            x = sprintf("WOVEN Dimension %d", dims[1L]),
            y = sprintf("WOVEN Dimension %d", dims[2L])
        )
    }

    p_out <- p_out + ggplot2::scale_shape_manual(
        values = c("Anchor" = 16L, "Partial data" = 1L),
        name = NULL, guide = if (hollow) "legend" else "none"
    )
    if (is.null(labels)) {
        p_out <- p_out +
            ggplot2::scale_color_manual(values = .pal_n(1L), guide = "none")
    } else if (continuous) {
        p_out <- p_out + ggplot2::scale_color_viridis_c(name = NULL)
    } else {
        G <- nlevels(df$group)
        p_out <- p_out + ggplot2::scale_color_manual(
            values = stats::setNames(.pal_n(G), levels(df$group)),
            name = "Class"
        )
    }

    n_scored <- nrow(df)
    n_anchor <- sum(df$anchor)
    p_out +
        ggplot2::labs(
            title = "WOVEN Latent Space",
            subtitle = sprintf(
                "%d of %d subjects scored (%.0f%%); %d anchors, %d %s",
                n_scored, x@n, 100 * n_scored / x@n, n_anchor,
                n_scored - n_anchor, "scored from a subset of modalities"
            )
        ) +
        .theme_woven() +
        ggplot2::theme(legend.position = "right")
}

#' Extract latent scores for new subjects
#'
#' Projects new subjects into the trained WOVEN latent space and returns an
#' n_new x K score matrix: each observed modality's features (standardized
#' with the training center and scale when the fit used \code{scale = TRUE})
#' are multiplied by that modality's projection matrix, and the results are
#' averaged over the subject's observed modalities. This linear map is the
#' Nystrom out-of-sample extension for WOVEN's linear embedding; no kernel
#' evaluation against anchor subjects is needed.
#'
#' For class predictions on new subjects, use [woven_predict()] instead.
#'
#' @param fit woven object from [woven()]
#' @param X_list_new the new subjects, as either
#'   \itemize{
#'     \item a list named with the training modality names, of matrices or
#'       data frames (subjects in rows). Blocks are matched by name, so
#'       their order does not matter and a modality may be left out (every
#'       new subject is then scored without it); blocks with row names are
#'       matched by subject ID; features are matched by column name when
#'       all training features are present;
#'     \item an unnamed list of V matrices in training modality and feature
#'       order, with all-NA rows for subjects missing a modality; or
#'     \item a \code{MultiAssayExperiment} whose experiments have the
#'       training modality names.
#'   }
#' @return numeric matrix n_new x K of consensus latent scores, with rownames
#'   from \code{X_list_new[[1]]} and colnames \code{Dim1}, \code{Dim2}, etc.
#'   Subjects with no observed data in any view receive a row of NA.
#' @examples
#' set.seed(1)
#' n <- 30
#' K <- 2L
#' X1 <- matrix(rnorm(n * 8), n, 8)
#' X2 <- matrix(rnorm(n * 6), n, 6)
#' Y <- rep(1:2, each = n / 2)
#' miss <- matrix(runif(n * 2) < 0.3, n, 2)
#' for (i in which(rowSums(miss) == 2)) miss[i, sample(2, 1)] <- FALSE
#' X1[miss[, 1], ] <- NA
#' X2[miss[, 2], ] <- NA
#' fit <- woven(list(X1, X2), Y = Y, K = K)
#' Z_new <- woven_scores(fit, list(X1[1:5, ], X2[1:5, ]))
#' dim(Z_new)
#' @seealso [woven_predict()] for class predictions, [woven()] for model fitting
#' @export
woven_scores <- function(fit, X_list_new) {
    stopifnot(methods::is(fit, "woven"))
    X_list_new <- .as_new_xlist(fit, X_list_new)
    .check_new_data(fit, X_list_new)
    if (isTRUE(fit@scaled) && !is.null(fit@scale_center)) {
        X_list_new <- .scale_apply(
            X_list_new, fit@scale_center, fit@scale_scale
        )
    }
    .project_views(X_list_new, fit@W_list, fit@K)
}

#' Predict class probabilities for new subjects
#'
#' Projects new subjects into the WOVEN latent space with [woven_scores()]
#' and returns soft class assignments from a nearest-centroid (or k-NN)
#' classifier built on the anchor subjects' latent scores (computed by the
#' same projection). Works the same
#' way for complete and block-missing subjects.
#'
#' @param fit woven object from [woven()]
#' @param X_list_new new subjects, in any form accepted by
#'   [woven_scores()] (for example \code{list(RNA = new_rna)} to classify
#'   subjects that have RNA only).
#' @param method \code{"centroid"} (default): inverse-squared-distance
#'   weights to the class centroids of the anchor scores.
#'   \code{"knn"}: k-NN vote among anchor subjects.
#' @param k_pred integer  -- number of neighbors for knn method (default 5)
#'
#' @return data.frame with n_new rows and columns \code{predicted_class}
#'   (class name, as in \code{levels(as.factor(Y))}), \code{confidence}
#'   (probability of the predicted class), and one \code{p_<class>} column
#'   per class with soft probabilities. Subjects with no observed modality
#'   get uniform probabilities.
#' @examples
#' set.seed(1)
#' n <- 40
#' K <- 2L
#' X1 <- matrix(rnorm(n * 5), n, 5)
#' X2 <- matrix(rnorm(n * 4), n, 4)
#' Y <- rep(1:2, each = n / 2)
#' miss <- matrix(FALSE, n, 2)
#' miss[c(31, 33, 35), 1] <- TRUE
#' miss[c(32, 34, 36), 2] <- TRUE
#' X1[miss[, 1], ] <- NA
#' X2[miss[, 2], ] <- NA
#' anchor_idx <- which(rowSums(miss) == 0)
#' fit <- woven(list(X1, X2), Y = Y, anchor_idx = anchor_idx, K = K)
#' pred <- woven_predict(fit, list(X1[1:5, ], X2[1:5, ]))
#' pred$predicted_class
#' @export
woven_predict <- function(
        fit, X_list_new, method = c("centroid", "knn"), k_pred = 5L) {
    stopifnot(methods::is(fit, "woven"))
    method <- match.arg(method)
    Z_new <- woven_scores(fit, X_list_new)
    pred_rownames <- rownames(Z_new)

    #  Classify in latent space
    Y_a <- fit@Y_anchor
    levels_Y <- fit@Y_levels
    C <- length(levels_Y)
    n_new <- nrow(Z_new)

    uniform_p <- rep(1 / C, C) # fallback for unscored subjects
    centroids <- NULL
    Z_ref_knn <- NULL

    classify_row <- function(z) {
        if (any(is.na(z))) {
            return(uniform_p)
        }
        if (method == "centroid") {
            dists <- sqrt(rowSums(sweep(centroids, 2, z)^2))
            dists <- pmax(dists, 1e-10)
            w <- 1 / dists^2
            w / sum(w)
        } else {
            dists <- sqrt(rowSums(sweep(Z_ref_knn, 2, z)^2))
            nn_idx <- order(dists)[seq_len(min(k_pred, length(dists)))]
            nn_y <- Y_a[nn_idx]
            tabulate(match(nn_y, levels_Y), nbins = C) / length(nn_idx)
        }
    }

    # Reference set: the anchors' latent scores computed the same way as the
    # new subjects' (projection of observed features), so both live in the
    # same space.
    Z_ref <- fit@Z[fit@anchor_idx, , drop = FALSE]

    if (method == "centroid") {
        centroids <- do.call(rbind, lapply(levels_Y, function(g) {
            idx <- which(Y_a == g)
            if (length(idx) == 0L) {
                return(rep(NA_real_, fit@K))
            }
            colMeans(Z_ref[idx, , drop = FALSE])
        }))
        Z_ref_knn <- NULL # not used
    } else {
        centroids <- NULL
        Z_ref_knn <- Z_ref
    }

    probs <- do.call(rbind, lapply(seq_len(n_new), function(i) {
        classify_row(Z_new[i, ])
    }))

    display_labels <- fit@Y_labels
    colnames(probs) <- paste0("p_", display_labels)
    pred_idx <- apply(probs, 1, which.max)

    data.frame(
        predicted_class = display_labels[pred_idx],
        confidence = probs[cbind(seq_len(n_new), pred_idx)],
        as.data.frame(probs),
        row.names = pred_rownames,
        stringsAsFactors = FALSE
    )
}
