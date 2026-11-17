# tune.R -- choose woven()'s settings by cross-validated classification error
#
# woven_tune() splits the anchor subjects into class-stratified folds. For
# each fold it refits everything on the remaining subjects (scaling, graphs,
# eigendecomposition, projection matrices), scores the held-out anchors from
# their features alone, and classifies them by linear discriminant analysis
# (LDA) on the training anchors' latent scores. The balanced error rate
# (BER) of those held-out predictions picks the settings; the full data are
# then refit with them.
#
# Within a fold the eigendecomposition depends only on gamma_y and lambda,
# not on ridge_w or screen_top, so it is solved once per (gamma_y, lambda)
# and the projection matrices are recovered for every (ridge_w, screen_top)
# pair from that one solution. This gives exactly the fits woven() would
# give for each setting, at a fraction of the cost.

#' S4 class for a woven_tune() result
#'
#' @slot best named list: the selected \code{gamma_y}, \code{ridge_w},
#'   \code{screen_top}, \code{lambda} and their held-out \code{ber}.
#' @slot results data.frame, one row per setting tried, sorted by
#'   \code{ber}: the settings, the held-out balanced error rate \code{ber}
#'   and its standard deviation across folds \code{ber_sd}.
#' @slot fit the \code{\linkS4class{woven}} fit on all data with the best
#'   settings, or \code{NULL} when \code{refit = FALSE}.
#' @slot folds list of integer vectors (one per fold and repeat): the
#'   held-out subjects (row indices) of each fold.
#' @slot chance numeric, the balanced error rate of guessing
#'   (\code{1 - 1/G} for G classes).
#' @seealso [woven_tune()]
#' @exportClass wovenTune
setClass("wovenTune", representation(
    best = "list",
    results = "data.frame",
    fit = "ANY",
    folds = "list",
    chance = "numeric"
))

#' Tune woven() by cross-validation
#'
#' Picks \code{gamma_y}, \code{ridge_w}, \code{screen_top} (and optionally
#' \code{lambdas}) by the balanced error rate of held-out anchor subjects,
#' then fits the model on all data with the best settings. This is the
#' protocol of the WOVEN paper and is what lets the model control its
#' capacity when a modality has many more features than there are anchors.
#'
#' The anchor subjects are split into \code{n_folds} class-stratified folds.
#' In each fold everything is refit on the other subjects, including
#' standardization and the similarity graphs, so no held-out subject
#' influences the model that predicts it. Held-out anchors are scored from
#' their features with the fold's projection matrices and classified by
#' linear discriminant analysis on the training anchors' scores (nearest
#' class centroid if the within-class covariance is singular). The balanced
#' error rate (BER) is the mean, over classes, of the fraction of that
#' class misclassified; guessing gives \code{1 - 1/G}.
#'
#' Subjects missing a modality are kept in every training fold (they shape
#' the similarity graphs and are scored like any other subject) but are not
#' used for evaluation, because a complete-case comparison needs the same
#' subjects for every setting. Fold assignment is random; call
#' \code{set.seed()} first for reproducible results. The BER of the selected
#' setting is optimistic, since it was the minimum over the grid; report a
#' separate held-out estimate if you need an unbiased one.
#'
#' @param X_list,Y the data and labels, in any form accepted by [woven()].
#' @param K integer, number of latent dimensions (default 5).
#' @param gamma_y numeric vector of supervision strengths to try.
#' @param ridge_w numeric vector of ridge strengths to try.
#' @param screen_top vector of feature-screening depths to try (\code{Inf}
#'   means no screening). Values at or above every modality's feature count
#'   are the same as \code{Inf} and are skipped.
#' @param lambdas numeric vector of Laplacian strengths to try (each used
#'   for all modalities; default 0.1 only, since results are usually
#'   insensitive to it).
#' @param k_nn,scale as in [woven()].
#' @param n_folds integer, number of cross-validation folds (default 3). It
#'   is reduced, with a message, if a class has fewer anchors than folds.
#' @param n_repeats integer, number of times to repeat the cross-validation
#'   with new folds (default 1); the BER is averaged over repeats.
#' @param refit logical, fit the full data with the best settings (default
#'   TRUE).
#' @param verbose logical, report progress (default TRUE).
#' @param plot logical, draw the tuning results when done (default: in
#'   interactive sessions only). Needs ggplot2.
#' @return a \code{wovenTune} object. \code{tune$best} holds the selected
#'   settings, \code{tune$results} every setting's held-out BER,
#'   \code{tune$fit} the final model; \code{plot(tune)} shows the BER over
#'   the grid.
#' @examples
#' data(woven_example)
#' set.seed(1)
#' tune <- woven_tune(woven_example$X_missing, woven_example$Y,
#'     K = 3L,
#'     gamma_y = c(1, 5), ridge_w = c(1e-8, 1), screen_top = c(Inf, 10),
#'     verbose = FALSE
#' )
#' tune
#' tune$best
#' fit <- tune$fit
#' @seealso [woven()]
#' @export
woven_tune <- function(
        X_list, Y, K = 5L,
        gamma_y = c(0.5, 1, 5, 10),
        ridge_w = c(1e-8, 0.1, 0.5, 1, 5),
        screen_top = c(Inf, 500, 200, 50),
        lambdas = 0.1, k_nn = 10L, scale = TRUE,
        n_folds = 3L, n_repeats = 1L, refit = TRUE,
        verbose = TRUE, plot = interactive()) {
    t_start <- proc.time()[["elapsed"]]
    inp <- .woven_check_input(X_list, Y, NULL, K, FALSE)
    X_all <- inp$X_list
    anchors <- inp$anchor_idx
    Y_fct <- as.factor(inp$Y)
    G <- nlevels(Y_fct)
    if (G < 2L) stop("Y must have at least two classes.", call. = FALSE)

    p_max <- max(vapply(X_all, ncol, integer(1L)))
    screen_top <- unique(ifelse(screen_top >= p_max, Inf, screen_top))
    grid <- expand.grid(
        gamma_y = unique(gamma_y), ridge_w = unique(ridge_w),
        screen_top = screen_top, lambda = unique(lambdas),
        KEEP.OUT.ATTRS = FALSE
    )
    n_folds <- .tune_n_folds(Y_fct[anchors], n_folds, verbose)
    folds <- unlist(lapply(seq_len(n_repeats), function(r) {
        .stratified_folds(anchors, Y_fct[anchors], n_folds)
    }), recursive = FALSE)

    if (verbose) {
        fmt <- paste0(
            "Tuning WOVEN | %d settings x %d folds%s on %d anchors ",
            "(%d classes; guessing gives BER %.2f)"
        )
        message(sprintf(
            fmt,
            nrow(grid), n_folds,
            if (n_repeats > 1L) sprintf(" x %d repeats", n_repeats) else "",
            length(anchors), G, 1 - 1 / G
        ))
    }

    pred <- matrix(NA_integer_, nrow(grid), length(Y_fct))
    pred_by_rep <- vector("list", n_repeats)
    fold_ber <- matrix(NA_real_, nrow(grid), length(folds))
    seen <- new.env(parent = emptyenv())
    seen$warns <- character(0)
    for (f in seq_along(folds)) {
        t0 <- proc.time()[["elapsed"]]
        r <- (f - 1L) %/% n_folds + 1L
        if ((f - 1L) %% n_folds == 0L) {
            pred <- matrix(NA_integer_, nrow(grid), length(Y_fct))
        }
        out <- withCallingHandlers(
            .tune_fold(
                X_all, Y_fct, folds[[f]], grid, K, k_nn, scale
            ),
            warning = function(w) {
                seen$warns <- c(seen$warns, conditionMessage(w))
                invokeRestart("muffleWarning")
            }
        )
        pred[, folds[[f]]] <- out
        fold_ber[, f] <- apply(out, 1L, .ber,
            truth = as.integer(Y_fct[folds[[f]]]), G = G
        )
        if (f %% n_folds == 0L) pred_by_rep[[r]] <- pred
        if (verbose) {
            message(sprintf(
                "  Fold %d/%d%s done (%s)", (f - 1L) %% n_folds + 1L, n_folds,
                if (n_repeats > 1L) sprintf(", repeat %d", r) else "",
                .elapsed(t0)
            ))
        }
    }
    for (w in unique(seen$warns)) warning(w, call. = FALSE)

    truth <- as.integer(Y_fct[anchors])
    ber_rep <- vapply(pred_by_rep, function(P) {
        apply(P[, anchors, drop = FALSE], 1L, .ber, truth = truth, G = G)
    }, numeric(nrow(grid)))
    ber <- rowMeans(matrix(ber_rep, nrow = nrow(grid)))
    res <- grid
    res$ber <- ber
    res$ber_sd <- if (ncol(fold_ber) > 1L) {
        apply(fold_ber, 1L, stats::sd)
    } else {
        NA_real_
    }
    res <- res[order(res$ber, -res$ridge_w), , drop = FALSE]
    rownames(res) <- NULL
    best <- as.list(res[1L, c("gamma_y", "ridge_w", "screen_top", "lambda",
        "ber")])

    if (verbose) {
        fmt <- paste0(
            "Best: gamma_y = %s, ridge_w = %s, screen_top = %s, ",
            "lambda = %s -> held-out BER %.3f (guessing %.2f). %s"
        )
        message(sprintf(
            fmt,
            format(best$gamma_y), format(best$ridge_w),
            format(best$screen_top), format(best$lambda), best$ber,
            1 - 1 / G, .elapsed(t_start)
        ))
    }

    fit <- NULL
    if (isTRUE(refit)) {
        if (verbose) message("Refitting on all subjects with these settings:")
        fit <- woven(X_all, Y_fct,
            K = K, lambdas = best$lambda, gamma_y = best$gamma_y,
            k_nn = k_nn, ridge_w = best$ridge_w,
            screen_top = best$screen_top, scale = scale,
            verbose = verbose, plot = FALSE
        )
    }
    out <- methods::new("wovenTune",
        best = best, results = res, fit = fit, folds = folds,
        chance = 1 - 1 / G
    )
    if (isTRUE(plot) && requireNamespace("ggplot2", quietly = TRUE)) {
        methods::show(.plot_tune(out))
    }
    out
}

# Number of folds that every class can fill with at least one anchor.
.tune_n_folds <- function(y_anchor, n_folds, verbose) {
    smallest <- min(table(y_anchor))
    if (smallest < 2L) {
        stop(
            "Every class needs at least 2 anchor subjects for ",
            "cross-validation; class '",
            names(which.min(table(y_anchor))), "' has ", smallest, ".",
            call. = FALSE
        )
    }
    if (smallest < n_folds) {
        if (verbose) {
            message(sprintf(
                "Using %d folds: the smallest class has %d anchors.",
                smallest, smallest
            ))
        }
        n_folds <- smallest
    }
    as.integer(n_folds)
}

# Class-stratified folds over the anchor indices.
.stratified_folds <- function(anchors, y_anchor, n_folds) {
    fold_of <- integer(length(anchors))
    for (g in levels(droplevels(y_anchor))) {
        idx <- which(y_anchor == g)
        idx <- idx[sample.int(length(idx))]
        fold_of[idx] <- rep_len(sample.int(n_folds), length(idx))
    }
    lapply(seq_len(n_folds), function(f) anchors[fold_of == f])
}

# Balanced error rate of integer predictions against integer truth.
.ber <- function(pred, truth, G) {
    errs <- vapply(seq_len(G), function(g) {
        in_g <- truth == g
        if (!any(in_g)) return(NA_real_)
        mean(pred[in_g] != g | is.na(pred[in_g]))
    }, numeric(1L))
    mean(errs, na.rm = TRUE)
}

# Fit every setting on one fold; return predictions (settings x held-out).
.tune_fold <- function(X_all, Y_fct, test, grid, K, k_nn, scale) {
    n <- nrow(X_all[[1L]])
    train <- setdiff(seq_len(n), test)
    X_tr <- lapply(X_all, function(X) X[train, , drop = FALSE])
    X_te <- lapply(X_all, function(X) X[test, , drop = FALSE])
    y_tr <- as.integer(Y_fct[train])
    if (isTRUE(scale)) {
        sf <- .scale_fit(X_tr)
        X_tr <- sf$X
        X_te <- .scale_apply(X_te, sf$center, sf$scale)
    }
    bm <- vapply(X_tr, .block_missing_rows, logical(length(train)))
    anc <- which(rowSums(matrix(bm, nrow = length(train))) == 0L)
    L_full <- .build_laplacians(X_tr, k_nn)
    La <- lapply(L_full, function(L) as.matrix(L[anc, anc, drop = FALSE]))
    X_tr_a <- lapply(X_tr, function(X) X[anc, , drop = FALSE])
    y_a <- y_tr[anc]

    out <- matrix(NA_integer_, nrow(grid), length(test))
    solves <- unique(grid[, c("gamma_y", "lambda")])
    for (s in seq_len(nrow(solves))) {
        raw <- tryCatch(
            woven_mcca_dual(
                X_tr, anc, y_tr,
                K = K, lambdas = solves$lambda[s],
                gamma_y = solves$gamma_y[s], k_nn = k_nn,
                La_list_precomp = La, verbose = FALSE
            ),
            error = function(e) NULL
        )
        if (is.null(raw)) next
        rows <- which(grid$gamma_y == solves$gamma_y[s] &
            grid$lambda == solves$lambda[s])
        for (i in rows) {
            W <- lapply(seq_along(X_tr), function(v) {
                .recover_W(
                    raw$Xa_list[[v]], raw$Za_list[[v]], y_a,
                    grid$ridge_w[i], grid$screen_top[i]
                )
            })
            Z_a <- .project_views(X_tr_a, W, raw$K)
            Z_te <- .project_views(X_te, W, raw$K)
            out[i, ] <- .lda_predict(.lda_fit(Z_a, y_a), Z_te)
        }
    }
    out
}

# Linear discriminant analysis with a pooled covariance; falls back to the
# nearest class centroid when that covariance is singular.
.lda_fit <- function(Z, y) {
    lv <- sort(unique(y))
    mu <- t(vapply(lv, function(g) {
        colMeans(Z[y == g, , drop = FALSE])
    }, numeric(ncol(Z))))
    mu <- matrix(mu, nrow = length(lv))
    R <- Z - mu[match(y, lv), , drop = FALSE]
    S <- crossprod(R) / max(1L, nrow(Z) - length(lv))
    S <- S + diag(1e-8 * max(mean(diag(S)), 1e-12), ncol(Z))
    S_inv <- tryCatch(solve(S), error = function(e) NULL)
    prior <- as.numeric(table(factor(y, levels = lv))) / length(y)
    list(classes = lv, mu = mu, S_inv = S_inv, log_prior = log(prior))
}

.lda_predict <- function(m, Z) {
    if (is.null(m$S_inv)) {
        d <- vapply(seq_along(m$classes), function(g) {
            rowSums(sweep(Z, 2L, m$mu[g, ], "-")^2)
        }, numeric(nrow(Z)))
        d <- matrix(d, nrow = nrow(Z))
        return(m$classes[max.col(-d, ties.method = "first")])
    }
    sc <- vapply(seq_along(m$classes), function(g) {
        mg <- m$mu[g, ]
        drop(Z %*% m$S_inv %*% mg) - 0.5 * drop(mg %*% m$S_inv %*% mg) +
            m$log_prior[g]
    }, numeric(nrow(Z)))
    sc <- matrix(sc, nrow = nrow(Z))
    m$classes[max.col(sc, ties.method = "first")]
}

#' Show a woven_tune() result
#'
#' @param object a \code{wovenTune} object from [woven_tune()]
#' @return Invisibly returns \code{NULL}; called for its printed output.
#' @examples
#' data(woven_example)
#' set.seed(1)
#' tune <- woven_tune(woven_example$X_missing, woven_example$Y,
#'     K = 3L, gamma_y = c(1, 5), ridge_w = c(1e-8, 1), screen_top = Inf,
#'     verbose = FALSE
#' )
#' tune
#' @export
setMethod("show", "wovenTune", function(object) {
    b <- object@best
    cat(sprintf(
        "WOVEN tuning: %d settings, held-out BER (guessing %.2f)\n",
        nrow(object@results), object@chance
    ))
    cat(sprintf(
        "  Best: gamma_y = %s, ridge_w = %s, screen_top = %s, lambda = %s\n",
        format(b$gamma_y), format(b$ridge_w), format(b$screen_top),
        format(b$lambda)
    ))
    cat(sprintf("  Held-out BER: %.3f\n\n", b$ber))
    top <- utils::head(object@results, 5L)
    top$ber <- round(top$ber, 3L)
    top$ber_sd <- round(top$ber_sd, 3L)
    cat("  Top settings:\n")
    print(top, row.names = FALSE)
    cat(if (is.null(object@fit)) {
        "\n  tune$fit: not refit (refit = FALSE)\n"
    } else {
        "\n  tune$fit: model refit with the best settings\n"
    })
    cat("  plot(tune): BER over the grid\n")
    invisible(NULL)
})

#' Access parts of a woven_tune() result
#'
#' \code{tune$best}, \code{tune$results}, \code{tune$fit}, \code{tune$folds}
#' and \code{tune$chance} return the slots of the same name.
#'
#' @param x a \code{wovenTune} object
#' @param name one of \code{best}, \code{results}, \code{fit},
#'   \code{folds}, \code{chance}
#' @return the requested slot
#' @examples
#' data(woven_example)
#' set.seed(1)
#' tune <- woven_tune(woven_example$X_missing, woven_example$Y,
#'     K = 3L, gamma_y = 1, ridge_w = c(1e-8, 1), screen_top = Inf,
#'     verbose = FALSE
#' )
#' tune$best
#' @aliases $,wovenTune-method
#' @export
setMethod("$", "wovenTune", function(x, name) {
    if (!name %in% methods::slotNames(x)) {
        stop(sprintf(
            "'%s' is not part of a woven_tune() result. Use one of: %s.",
            name, paste(methods::slotNames(x), collapse = ", ")
        ), call. = FALSE)
    }
    methods::slot(x, name)
})

#' Plot woven_tune() results
#'
#' Held-out balanced error rate for every setting tried: one panel per
#' \code{screen_top} value, \code{ridge_w} across, \code{gamma_y} down,
#' darker for lower error. The selected setting is outlined. When several
#' \code{lambda} values were tried, the best over \code{lambda} is shown.
#'
#' @param x a \code{wovenTune} object from [woven_tune()]
#' @param y unused
#' @param ... unused
#' @return a ggplot object
#' @examples
#' data(woven_example)
#' set.seed(1)
#' tune <- woven_tune(woven_example$X_missing, woven_example$Y,
#'     K = 3L, gamma_y = c(1, 5), ridge_w = c(1e-8, 1), screen_top = Inf,
#'     verbose = FALSE
#' )
#' plot(tune)
#' @export
setMethod("plot", signature(x = "wovenTune", y = "missing"),
    function(x, y, ...) .plot_tune(x)
)

.plot_tune <- function(x) {
    .require_ggplot2()
    r <- x@results
    r <- stats::aggregate(ber ~ gamma_y + ridge_w + screen_top,
        data = r, FUN = min
    )
    fmt <- function(v) ifelse(is.infinite(v), "all features",
        paste("top", vapply(v, format, character(1L)))
    )
    fmt_ridge <- function(v) ifelse(v <= 1e-6, "off",
        vapply(v, format, character(1L))
    )
    r$screen <- factor(fmt(r$screen_top),
        levels = fmt(sort(unique(r$screen_top), decreasing = TRUE))
    )
    r$ridge <- factor(fmt_ridge(r$ridge_w),
        levels = fmt_ridge(sort(unique(r$ridge_w)))
    )
    r$gamma <- factor(vapply(r$gamma_y, format, character(1L)),
        levels = vapply(sort(unique(r$gamma_y)), format, character(1L))
    )
    # Dark text on the light (low-error) tiles, white on the dark ones.
    r$ink <- ifelse(r$ber <= stats::median(r$ber), "gray10", "white")
    b <- x@best
    r$best <- r$gamma_y == b$gamma_y & r$ridge_w == b$ridge_w &
        r$screen_top == b$screen_top
    ggplot2::ggplot(r, ggplot2::aes(x = ridge, y = gamma, fill = ber)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.6) +
        ggplot2::geom_tile(
            data = r[r$best, , drop = FALSE], fill = NA,
            color = "black", linewidth = 1.1
        ) +
        ggplot2::geom_text(ggplot2::aes(
            label = sprintf("%.2f", ber), color = ink
        ), size = 3.2) +
        ggplot2::scale_color_identity() +
        ggplot2::facet_wrap(~screen, nrow = 1L) +
        ggplot2::scale_fill_viridis_c(
            name = "Held-out\nBER", direction = -1, option = "D",
            end = 0.9
        ) +
        ggplot2::labs(
            x = "ridge_w", y = "gamma_y",
            title = "WOVEN tuning",
            subtitle = sprintf(
                paste0(
                    "Held-out balanced error rate; guessing gives %.2f. ",
                    "Outlined: selected."
                ),
                x@chance
            )
        ) +
        .theme_woven() +
        ggplot2::theme(panel.grid = ggplot2::element_blank())
}
