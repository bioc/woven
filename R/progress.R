# progress.R -- what woven() reports while it runs
#
# A fit prints a short header describing the data as woven understood it
# (modalities, features, how many subjects each modality covers, anchors),
# one line per step with timings, and a closing summary: coverage against a
# complete-case analysis and the top features per modality. Everything goes
# through message(), so verbose = FALSE or suppressMessages() silences it.

# Seconds since t0, formatted for a progress line.
.elapsed <- function(t0) {
    s <- proc.time()[["elapsed"]] - t0
    if (s < 60) sprintf("%.1fs", s) else sprintf("%.1f min", s / 60)
}

# Describe the input before fitting.
.woven_header <- function(X_list, Y_fct, anchor_idx, aligned) {
    n <- nrow(X_list[[1L]])
    V <- length(X_list)
    tab <- table(Y_fct)
    message(sprintf(
        "WOVEN | %d subjects | %d modalities | %d classes (%s)",
        n, V, length(tab),
        paste(sprintf("%s %d", names(tab), as.integer(tab)), collapse = ", ")
    ))
    if (aligned) {
        message(sprintf(
            "  Matched subjects across modalities by row name (%d unique IDs).",
            n
        ))
    }
    w <- max(nchar(names(X_list)))
    for (v in seq_len(V)) {
        n_obs <- sum(!.block_missing_rows(X_list[[v]]))
        message(sprintf(
            "  %-*s %6d features  %5d subjects (%3.0f%%)",
            w, names(X_list)[v], ncol(X_list[[v]]), n_obs, 100 * n_obs / n
        ))
    }
    n_a <- length(anchor_idx)
    message(sprintf(
        "  Anchors (observed in every modality): %d of %d (%.0f%%)",
        n_a, n, 100 * n_a / n
    ))
}

# Build one Laplacian per modality, reporting each as it finishes.
.build_laplacians_verbose <- function(X_list, k_nn, verbose) {
    out <- vector("list", length(X_list))
    names(out) <- names(X_list)
    for (v in seq_along(X_list)) {
        t0 <- proc.time()[["elapsed"]]
        out[[v]] <- build_laplacian(X_list[[v]], k = k_nn)
        if (verbose) {
            message(sprintf(
                "[1/3] Similarity graph %d/%d: %s (%s)",
                v, length(X_list), names(X_list)[v], .elapsed(t0)
            ))
        }
    }
    out
}

# Close the run: coverage against complete-case, top features, next steps.
.woven_footer <- function(fit, t_start) {
    n_scored <- sum(!is.na(fit@Z[, 1L]))
    n_a <- length(fit@anchor_idx)
    message(sprintf(
        "Done in %s. Scored %d of %d subjects (%.0f%%); a complete-case %s",
        .elapsed(t_start), n_scored, fit@n, 100 * n_scored / fit@n,
        sprintf("method would keep %d (%.0f%%).", n_a, 100 * n_a / fit@n)
    ))
    top <- woven_top_features(fit, n = 3L)
    w <- max(nchar(fit@mod_names))
    message("Top features:")
    for (m in fit@mod_names) {
        f <- top$feature[top$modality == m]
        message(sprintf("  %-*s %s", w, m, paste(f, collapse = ", ")))
    }
    message(
        "Next: plot(fit), plot(fit, type = \"features\"), ",
        "woven_top_features(fit), woven_predict(fit, new_data)"
    )
    if (fit@ridge_w <= 1e-8 && is.infinite(fit@screen_top)) {
        message(
            "Tip: woven_tune() picks ridge_w and screen_top by ",
            "cross-validation; the defaults leave both off."
        )
    }
}

# Mean of each feature in each class, over subjects observed in the
# modality. Returns a p x G matrix (rows = features, columns = classes).
.class_means <- function(X, Y_fct, obs) {
    lv <- levels(Y_fct)
    M <- vapply(lv, function(g) {
        rows <- obs & Y_fct == g
        if (!any(rows)) return(rep(NA_real_, ncol(X)))
        colMeans(X[rows, , drop = FALSE], na.rm = TRUE)
    }, numeric(ncol(X)))
    M <- matrix(M, nrow = ncol(X), dimnames = list(colnames(X), lv))
    M
}
