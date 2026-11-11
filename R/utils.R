# utils.R -- internal helpers shared across the package

# Logical vector: TRUE for rows in which every entry is NA, i.e. subjects with
# no data for that modality (a block-missing row).
.block_missing_rows <- function(X) {
    rowSums(is.na(X)) == ncol(X)
}

# Logical vector: TRUE for columns in which every entry is NA.
.all_na_cols <- function(X) {
    colSums(is.na(X)) == nrow(X)
}

# Label used to name a subject in error messages: its row name if present,
# otherwise its integer row index.
.subject_ids <- function(X_list, idx) {
    rn <- rownames(X_list[[1L]])
    if (is.null(rn)) as.character(idx) else rn[idx]
}

# Standardize each modality's features to zero mean and unit variance using
# observed-value statistics. Block-missing (all-NA) rows stay NA. Returns the
# standardized X_list plus per-feature center/scale so out-of-sample subjects
# can be scaled identically in woven_scores()/woven_predict().
.scale_fit <- function(X_list) {
    fitted <- lapply(X_list, function(X) {
        ctr <- colMeans(X, na.rm = TRUE)
        scl <- apply(X, 2L, stats::sd, na.rm = TRUE)
        ctr[!is.finite(ctr)] <- 0
        scl[!is.finite(scl) | scl == 0] <- 1
        list(
            X = sweep(sweep(X, 2L, ctr, "-"), 2L, scl, "/"),
            center = ctr, scale = scl
        )
    })
    Xs <- lapply(fitted, `[[`, "X")
    names(Xs) <- names(X_list)
    list(
        X = Xs,
        center = lapply(fitted, `[[`, "center"),
        scale = lapply(fitted, `[[`, "scale")
    )
}

# Apply stored center/scale (from .scale_fit) to new modality matrices.
.scale_apply <- function(X_list, centers, scales) {
    out <- lapply(seq_along(X_list), function(v) {
        sweep(sweep(X_list[[v]], 2L, centers[[v]], "-"), 2L, scales[[v]], "/")
    })
    names(out) <- names(X_list)
    out
}

# Drop block-missing rows and all-NA columns, then impute any remaining
# feature-level NAs (scattered within otherwise observed rows) with the
# column median. Never used to fill a whole missing modality block.
na_impute_median <- function(X) {
    Xc <- X[!.block_missing_rows(X), , drop = FALSE]
    Xc <- Xc[, !.all_na_cols(Xc), drop = FALSE]
    col_med <- MatrixGenerics::colMedians(Xc, na.rm = TRUE)
    for (j in which(colSums(is.na(Xc)) > 0L)) {
        Xc[is.na(Xc[, j]), j] <- col_med[j]
    }
    Xc[!is.finite(Xc)] <- 0
    Xc
}

# Symmetric matrix square root and inverse square root via eigendecomposition.
# Eigenvalues are floored at tol * max|eigenvalue| before taking roots.
mat_sqrt <- function(B, tol = 1e-10) {
    eig <- eigen(B, symmetric = TRUE)
    vals <- pmax(eig$values, tol * max(abs(eig$values)))
    vecs <- eig$vectors
    list(
        sqrt = vecs %*% diag(sqrt(vals), length(vals)) %*% t(vecs),
        inv_sqrt = vecs %*% diag(1 / sqrt(vals), length(vals)) %*% t(vecs)
    )
}

# Validate anchor_idx against X_list. Anchors must be distinct integer row
# indices in 1..n, and every anchor must be observed (not block-missing) in
# every modality. Offending subjects are named in the error message.
.check_anchor_idx <- function(X_list, anchor_idx) {
    n <- nrow(X_list[[1L]])
    if (!is.numeric(anchor_idx) || anyNA(anchor_idx) ||
        any(anchor_idx != round(anchor_idx))) {
        stop("anchor_idx must be a vector of integer row indices without NA.",
            call. = FALSE
        )
    }
    anchor_idx <- as.integer(anchor_idx)
    bad <- anchor_idx[anchor_idx < 1L | anchor_idx > n]
    if (length(bad) > 0L) {
        stop(sprintf(
            "anchor_idx contains indices outside 1..%d: %s.",
            n, paste(utils::head(bad, 10L), collapse = ", ")
        ), call. = FALSE)
    }
    dup <- unique(anchor_idx[duplicated(anchor_idx)])
    if (length(dup) > 0L) {
        stop(sprintf(
            "anchor_idx contains duplicated subjects: %s.",
            paste(.subject_ids(X_list, utils::head(dup, 10L)), collapse = ", ")
        ), call. = FALSE)
    }
    mod_nm <- names(X_list)
    if (is.null(mod_nm)) mod_nm <- paste0("modality ", seq_along(X_list))
    msgs <- character(0)
    for (v in seq_along(X_list)) {
        miss <- anchor_idx[.block_missing_rows(
            X_list[[v]][anchor_idx, , drop = FALSE]
        )]
        if (length(miss) > 0L) {
            msgs <- c(msgs, sprintf(
                "%s: %s", mod_nm[v],
                paste(.subject_ids(X_list, utils::head(miss, 10L)),
                    collapse = ", "
                )
            ))
        }
    }
    if (length(msgs) > 0L) {
        stop(
            "Anchor subjects must be observed in every modality, but these ",
            "anchors are missing a whole modality block:\n  ",
            paste(msgs, collapse = "\n  "),
            "\nDrop them from anchor_idx, or leave anchor_idx = NULL to ",
            "detect anchors automatically.",
            call. = FALSE
        )
    }
    invisible(anchor_idx)
}

# Shared input checks for woven() and woven_mcca_dual().
check_woven_inputs <- function(X_list, anchor_idx, K) {
    V <- length(X_list)
    if (V < 2L) stop("WOVEN requires at least 2 modalities.", call. = FALSE)
    ns <- vapply(X_list, nrow, integer(1L))
    if (length(unique(ns)) > 1L) {
        stop("All modality matrices must have the same number of rows.",
            call. = FALSE
        )
    }
    .check_anchor_idx(X_list, anchor_idx)

    n_a <- length(anchor_idx)
    if (n_a < K) {
        stop(sprintf(
            "Anchor set (n_a=%d) must be >= K=%d. Use more anchors or reduce K.",
            n_a, K
        ), call. = FALSE)
    }
    if (n_a < 5L * K) {
        fmt <- paste0(
            "Anchor set (n_a=%d) is small relative to K=%d ",
            "(recommend n_a >= 5K). Results may be unstable."
        )
        warning(sprintf(
            fmt,
            n_a, K
        ), call. = FALSE)
    }
    # Each modality's projection matrix W_v (p_v x K) has rank at most
    # min(n_a, p_v), so a modality with fewer features than K cannot span
    # K independent latent directions on its own.
    p_v <- vapply(X_list, ncol, integer(1L))
    small <- which(p_v < K)
    if (length(small) > 0L) {
        nm <- names(X_list)
        if (is.null(nm)) nm <- paste0("modality ", seq_len(V))
        fmt <- paste0(
            "K=%d exceeds the number of features in %s. Scores from ",
            "such a modality span at most p_v of the K dimensions ",
            "(each W_v has rank <= min(n_a, p_v))."
        )
        warning(sprintf(
            fmt,
            K, paste(sprintf("%s (p=%d)", nm[small], p_v[small]),
                collapse = ", "
            )
        ), call. = FALSE)
    }
    invisible(TRUE)
}

# Project subjects into the latent space: for each subject, apply W_v to every
# modality it has observed and average over those modalities. Feature-level
# NAs within an observed row contribute zero. Subjects with no observed
# modality get a row of NA. Used by woven() and woven_scores().
.project_views <- function(X_list, W_list, K) {
    n <- nrow(X_list[[1L]])
    Z_acc <- matrix(0, n, K)
    obs_cnt <- integer(n)
    for (v in seq_along(X_list)) {
        obs <- which(!.block_missing_rows(X_list[[v]]))
        if (length(obs) == 0L) next
        Xv_obs <- X_list[[v]][obs, , drop = FALSE]
        Xv_obs[is.na(Xv_obs)] <- 0
        Z_acc[obs, ] <- Z_acc[obs, ] + Xv_obs %*% W_list[[v]]
        obs_cnt[obs] <- obs_cnt[obs] + 1L
    }
    none <- obs_cnt == 0L
    obs_cnt[none] <- 1L
    Z <- Z_acc / obs_cnt
    Z[none, ] <- NA_real_
    dimnames(Z) <- list(rownames(X_list[[1L]]), paste0("Dim", seq_len(K)))
    Z
}

# Check that X_list_new matches the modalities and feature counts of a fit.
.check_new_data <- function(fit, X_list_new) {
    if (!is.list(X_list_new) || length(X_list_new) != fit@V) {
        stop(sprintf(
            "X_list_new must be a list of %d matrices (one per modality).",
            fit@V
        ), call. = FALSE)
    }
    for (v in seq_len(fit@V)) {
        p_new <- ncol(X_list_new[[v]])
        p_fit <- nrow(fit@W_list[[v]])
        if (!identical(as.integer(p_new), as.integer(p_fit))) {
            stop(sprintf(
                "X_list_new[[%d]] (%s) has %d features but model was trained on %d.",
                v, fit@mod_names[v], p_new, p_fit
            ), call. = FALSE)
        }
    }
    invisible(TRUE)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Suppress R CMD check notes for variables used inside closures
utils::globalVariables(c("centroids", "Z_ref_knn"))
