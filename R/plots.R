# plots.R -- ggplot2-based publication-ready plots for WOVEN fits
# ggplot2 is in Suggests: every plotting function calls .require_ggplot2().
#
# woven_plot_vip()      -- top features per modality by VIP, coloured by class
# woven_plot_coverage() -- which subjects have which modalities, by class
# woven_plot_loadings() -- Feature loadings for one latent dimension
# woven_plot_variance() -- share of squared canonical correlation per dimension

# Suppress R CMD check notes for ggplot2 aesthetic variable names
utils::globalVariables(c(
    "feature", "vip_score",
    "loading", "sign_str", "modality_label",
    "dim_int", "pct_var", "cumulative",
    "z1", "z2", "group", "point_type", "anchor", "shape",
    "abs_loading", "key", "highest_in", "pattern", "Freq", "n", "ber", "ridge",
    "gamma", "ink"
))

.require_ggplot2 <- function() {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop(
            "Package 'ggplot2' is needed for WOVEN plots.\n",
            "Install it with: install.packages(\"ggplot2\")",
            call. = FALSE
        )
    }
    invisible(TRUE)
}

# Resolve a modality argument (index or name) to a single valid index.
.resolve_modality <- function(modality, fit) {
    if (length(modality) != 1L || is.na(modality)) {
        stop("modality must be a single modality index or name.", call. = FALSE)
    }
    if (is.character(modality)) {
        idx <- match(modality, fit@mod_names)
        if (is.na(idx)) {
            stop(sprintf(
                "Modality '%s' not found. Available: %s",
                modality, paste(fit@mod_names, collapse = ", ")
            ), call. = FALSE)
        }
        return(idx)
    }
    if (!is.numeric(modality) || modality != round(modality) ||
        modality < 1L || modality > fit@V) {
        stop(sprintf(
            "modality must be an integer in 1..%d or one of: %s",
            fit@V, paste(fit@mod_names, collapse = ", ")
        ), call. = FALSE)
    }
    as.integer(modality)
}

# Check that dim is a single latent dimension in 1..K.
.check_dim <- function(dim, fit) {
    if (length(dim) != 1L || !is.numeric(dim) || is.na(dim) ||
        dim != round(dim) || dim < 1L || dim > fit@K) {
        stop(sprintf(
            "dim must be a single integer in 1..%d (the fit has K=%d dimensions).",
            fit@K, fit@K
        ), call. = FALSE)
    }
    as.integer(dim)
}

# Colorblind-safe palette (Wong 2011, Nature Methods)
.pal_woven <- c(
    "#0072B2", "#E69F00", "#009E73", "#D55E00",
    "#CC79A7", "#56B4E9", "#F0E442", "#999999"
)

# n colours: the colour-blind-safe palette for up to 8 classes, an evenly
# spaced HCL palette beyond that (so colours never repeat).
.pal_n <- function(n) {
    if (n <= length(.pal_woven)) .pal_woven[seq_len(n)] else
        grDevices::hcl.colors(n, "Dark 3")
}

# Axis labels: drop the "<panel>___" key prefix used to keep identically
# named features in different panels apart, and shorten very long names.
.axis_labels <- function(width = 32L) {
    function(x) {
        x <- sub("^[0-9]+___", "", x)
        long <- nchar(x) > width
        x[long] <- paste0(substr(x[long], 1L, width - 3L), "...")
        x
    }
}

# Shared minimal theme for all WOVEN plots
.theme_woven <- function(base_size = 12) {
    ggplot2::theme_bw(base_size = base_size) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                face = "bold", hjust = 0.5, size = base_size + 1
            ),
            plot.subtitle = ggplot2::element_text(
                hjust = 0.5, color = "gray40", size = base_size - 1
            ),
            panel.grid.minor = ggplot2::element_blank(),
            strip.background = ggplot2::element_rect(
                fill = "gray95", color = "gray80"
            ),
            strip.text       = ggplot2::element_text(face = "bold")
        )
}

#' Plot the top features of a WOVEN fit by VIP score
#'
#' Bar chart of the features with the highest Variable Importance in
#' Projection (VIP) in each modality, one panel per modality. VIP combines
#' a feature's loadings over all K latent dimensions, weighted by each
#' dimension's share of the squared singular values (see
#' [woven_top_features()] for the same ranking as a table). By default each
#' bar is coloured by the class in which that feature has its highest mean,
#' as in DIABLO's \code{plotLoadings(contrib = "max")}, so the plot shows
#' both which features matter and which class they mark. A dashed line at
#' VIP = 1 marks above-average importance.
#'
#' @param fit woven object from [woven()]
#' @param modality \code{NULL} (default) for every modality, or one or more
#'   modality names or indices (1..V).
#' @param n_top integer: number of top features per modality (default 15)
#' @param feature_names optional character vector of feature labels, one per
#'   feature; only with a single \code{modality}. If NULL, uses the feature
#'   (column) names of the data.
#' @param color_by \code{"class"} (default): colour each bar by the class
#'   with the highest mean for that feature. \code{"modality"}: one colour
#'   per modality.
#' @param main character: plot title. If NULL, a default is generated.
#' @return a ggplot object (printed automatically; add layers with \code{+})
#' @examples
#' data(woven_example)
#' fit <- woven(woven_example$X_missing,
#'     Y = woven_example$Y, K = 3L,
#'     verbose = FALSE
#' )
#' woven_plot_vip(fit)
#' woven_plot_vip(fit, modality = "RNA", n_top = 10)
#' @seealso [woven_top_features()], [woven_plot_loadings()],
#'   [woven_plot_variance()]
#' @export
woven_plot_vip <- function(
        fit, modality = NULL, n_top = 15L, feature_names = NULL,
        color_by = c("class", "modality"), main = NULL) {
    .require_ggplot2()
    stopifnot(methods::is(fit, "woven"))
    color_by <- match.arg(color_by)
    mods <- if (is.null(modality)) {
        seq_len(fit@V)
    } else {
        vapply(modality, .resolve_modality, integer(1L), fit = fit)
    }
    if (!is.null(feature_names) && length(mods) != 1L) {
        stop("feature_names can only be given with a single modality.",
            call. = FALSE
        )
    }

    rows <- lapply(mods, function(v) {
        vip <- .vip_scores(fit, v)
        p <- length(vip)
        fn <- if (is.null(feature_names)) names(vip) else feature_names
        if (length(fn) != p) {
            stop(sprintf(
                "feature_names length %d != %d features in modality %d.",
                length(fn), p, v
            ), call. = FALSE)
        }
        # Features removed by screen_top have zero loadings; leave them out.
        k <- min(as.integer(n_top), sum(vip > 0))
        o <- order(vip, decreasing = TRUE)[seq_len(k)]
        data.frame(
            modality_label = fit@mod_names[v],
            feature = fn[o],
            key = paste0(v, "___", fn[o]),
            vip_score = unname(vip[o]),
            highest_in = .highest_class(fit, v)[o],
            stringsAsFactors = FALSE
        )
    })
    df <- do.call(rbind, rows)
    df$key <- factor(df$key, levels = rev(unique(df$key)))
    df$modality_label <- factor(df$modality_label,
        levels = fit@mod_names[mods]
    )

    by_class <- color_by == "class" && !all(is.na(df$highest_in))
    if (by_class) {
        lv <- fit@Y_labels
        df$highest_in <- factor(df$highest_in, levels = lv)
        p_out <- ggplot2::ggplot(
            df, ggplot2::aes(x = vip_score, y = key, fill = highest_in)
        ) +
            ggplot2::geom_col(width = 0.7) +
            ggplot2::scale_fill_manual(
                values = stats::setNames(.pal_n(length(lv)), lv),
                name = "Highest in"
            )
    } else {
        p_out <- ggplot2::ggplot(
            df, ggplot2::aes(x = vip_score, y = key, fill = modality_label)
        ) +
            ggplot2::geom_col(width = 0.7) +
            ggplot2::scale_fill_manual(
                values = stats::setNames(
                    .pal_n(fit@V)[mods], fit@mod_names[mods]
                ),
                guide = "none"
            )
    }

    if (is.null(main)) {
        main <- if (length(mods) == 1L) {
            sprintf("Top features - %s", fit@mod_names[mods])
        } else {
            "Top features per modality"
        }
    }
    subtitle <- if (by_class) {
        "VIP score; colour = class in which the feature is highest"
    } else {
        "VIP score"
    }

    p_out <- p_out +
        ggplot2::scale_y_discrete(labels = .axis_labels()) +
        ggplot2::scale_x_continuous(
            expand = ggplot2::expansion(mult = c(0, 0.05))
        ) +
        ggplot2::labs(x = "VIP score", y = NULL, title = main,
            subtitle = subtitle
        ) +
        .theme_woven() +
        ggplot2::theme(
            panel.grid.major.y = ggplot2::element_blank(),
            legend.position = "bottom"
        )
    if (length(mods) > 1L) {
        p_out <- p_out + ggplot2::facet_wrap(~modality_label,
            scales = "free", ncol = min(length(mods), 3L)
        )
    }
    if (max(df$vip_score) > 0.9) {
        p_out <- p_out + ggplot2::geom_vline(
            xintercept = 1, linetype = "dashed",
            color = "gray50", linewidth = 0.5
        )
    }
    p_out
}

#' Plot which subjects have which modalities
#'
#' One bar per pattern of observed modalities (for example "RNA +
#' Proteomics"), showing how many subjects have exactly that pattern,
#' split by class. The subtitle compares the subjects a complete-case
#' method can use (those observed in every modality, the anchors) with the
#' subjects WOVEN scores. The class split shows whether the subjects a
#' complete-case analysis would drop differ in their class mix.
#'
#' @param fit woven object from [woven()]
#' @param max_patterns integer: show at most this many bars (default 12).
#'   With many modalities there can be dozens of patterns; the least common
#'   are combined into one "Other patterns" bar.
#' @param main character: plot title (default "Modality coverage")
#' @return a ggplot object (printed automatically; add layers with \code{+})
#' @examples
#' data(woven_example)
#' fit <- woven(woven_example$X_missing,
#'     Y = woven_example$Y, K = 3L,
#'     verbose = FALSE
#' )
#' woven_plot_coverage(fit)
#' @seealso [woven()], [woven_ess_retention()]
#' @export
woven_plot_coverage <- function(fit, max_patterns = 12L,
                                main = "Modality coverage") {
    .require_ggplot2()
    stopifnot(methods::is(fit, "woven"))
    obs <- fit@observed
    if (length(obs) == 0L) {
        stop("This fit does not record which modalities each subject has; ",
            "refit it with woven >= 0.99.6.",
            call. = FALSE
        )
    }
    pattern <- apply(obs, 1L, function(r) {
        if (!any(r)) "(no data)" else paste(colnames(obs)[r], collapse = " + ")
    })
    n_mod <- rowSums(obs)
    is_anchor <- n_mod == ncol(obs)
    pattern[is_anchor] <- "All modalities (anchors)"
    # Anchors first, then the most common patterns. With many modalities
    # there can be dozens of patterns, so beyond max_patterns the smallest
    # are lumped together.
    tab <- table(pattern)
    lv <- names(tab)[order(
        -(names(tab) %in% pattern[is_anchor]), -as.integer(tab), names(tab)
    )]
    if (length(lv) > max_patterns) {
        other <- sprintf("Other patterns (%d)", length(lv) - max_patterns + 1L)
        pattern[pattern %in% lv[max_patterns:length(lv)]] <- other
        tab <- table(pattern)
        lv <- c(lv[seq_len(max_patterns - 1L)], other)
    }
    group <- factor(as.character(fit@Y), levels = fit@Y_labels)
    df <- as.data.frame(table(pattern = pattern, group = group),
        responseName = "n", stringsAsFactors = FALSE
    )
    df <- df[df$n > 0L, , drop = FALSE]
    df$pattern <- factor(df$pattern, levels = rev(lv))
    df$group <- factor(df$group, levels = fit@Y_labels)
    totals <- data.frame(
        pattern = factor(names(tab), levels = rev(lv)),
        Freq = as.integer(tab)
    )
    n_a <- sum(is_anchor)
    n_scored <- sum(!is.na(fit@Z[, 1L]))
    G <- length(fit@Y_labels)

    p_out <- ggplot2::ggplot(
        df, ggplot2::aes(x = n, y = pattern, fill = group)
    ) +
        ggplot2::geom_col(
            width = 0.7, position = ggplot2::position_stack(reverse = TRUE)
        ) +
        ggplot2::geom_text(
            data = totals, inherit.aes = FALSE,
            ggplot2::aes(x = Freq, y = pattern, label = Freq),
            hjust = -0.3, size = 3.5, color = "gray25"
        ) +
        ggplot2::scale_fill_manual(
            values = stats::setNames(.pal_n(G), fit@Y_labels),
            name = "Class"
        ) +
        ggplot2::scale_y_discrete(labels = .axis_labels(60L)) +
        ggplot2::scale_x_continuous(
            expand = ggplot2::expansion(mult = c(0, 0.1))
        ) +
        ggplot2::labs(
            x = "Subjects", y = NULL, title = main,
            subtitle = sprintf(
                paste0(
                    "Complete-case keeps %d of %d (%.0f%%); ",
                    "WOVEN scores %d (%.0f%%)"
                ),
                n_a, fit@n, 100 * n_a / fit@n, n_scored,
                100 * n_scored / fit@n
            )
        ) +
        .theme_woven() +
        ggplot2::theme(
            panel.grid.major.y = ggplot2::element_blank(),
            legend.position = "bottom"
        )
    p_out
}

#' Plot feature loadings for one WOVEN latent dimension
#'
#' For a given latent dimension, shows the top features by absolute loading
#' for each modality (or a selected subset), colored by loading sign.
#' Positive loadings are blue; negative loadings are red-orange.
#' Equivalent to DIABLO's \code{plotLoadings()}.
#'
#' @param fit woven object from [woven()]
#' @param dim integer: which latent dimension to plot (1..K, default 1)
#' @param n_top integer: number of top features per modality (default 15)
#' @param feature_names optional list of V character vectors (one per modality).
#'   If a plain character vector is passed for a single-modality call, it is
#'   used for that modality. If NULL, uses rownames of W or "Feature_j".
#' @param modality NULL, or one or more modality indices or names: plot only
#'   those modalities.
#'   If NULL (default), all V modalities are shown in faceted panels.
#' @param main character: plot title. If NULL, a default is used.
#' @return a ggplot object (printed automatically; add layers with \code{+})
#' @examples
#' set.seed(1)
#' n <- 60
#' p1 <- 30
#' p2 <- 20
#' K <- 3
#' Y <- rep(1:3, each = n / 3)
#' X1 <- matrix(rnorm(n * p1), n, p1)
#' X2 <- matrix(rnorm(n * p2), n, p2)
#' miss <- matrix(runif(n * 2) < 0.3, n, 2)
#' for (i in which(rowSums(miss) == 2)) miss[i, sample(2, 1)] <- FALSE
#' X1[miss[, 1], ] <- NA
#' X2[miss[, 2], ] <- NA
#' fit <- woven(list(X1, X2), Y = Y, K = K)
#' woven_plot_loadings(fit, dim = 1L)
#' @seealso [woven_plot_vip()], [woven_plot_variance()]
#' @export
woven_plot_loadings <- function(fit, dim = 1L, n_top = 15L,
                                feature_names = NULL, modality = NULL,
                                main = NULL) {
    .require_ggplot2()
    stopifnot(methods::is(fit, "woven"))
    dim <- .check_dim(dim, fit)

    mod_seq <- if (is.null(modality)) {
        seq_len(fit@V)
    } else {
        vapply(modality, .resolve_modality, integer(1L), fit = fit)
    }

    # Build feature name list (length V)
    if (is.null(feature_names)) {
        feature_names <- lapply(seq_len(fit@V), function(v) {
            W <- fit@W_list[[v]]
            if (!is.null(rownames(W))) {
                rownames(W)
            } else {
                paste0("Feature_", seq_len(nrow(W)))
            }
        })
    } else if (is.character(feature_names)) {
        fn <- feature_names
        feature_names <- lapply(seq_len(fit@V), function(v) fn)
    }

    col_pos <- "#0072B2"
    col_neg <- "#D55E00"

    mod_names <- if (!is.null(fit@mod_names)) {
        fit@mod_names
    } else {
        paste0("Modality ", seq_len(fit@V))
    }

    rows <- lapply(mod_seq, function(v) {
        W <- fit@W_list[[v]]
        w <- W[, dim]
        p <- length(w)
        nm <- feature_names[[v]]
        if (length(nm) != p) {
            stop(sprintf(
                "feature_names[[%d]] length %d != %d.", v, length(nm), p
            ))
        }

        n_show <- min(as.integer(n_top), sum(w != 0))
        top_idx <- order(abs(w), decreasing = TRUE)[seq_len(n_show)]

        data.frame(
            modality_label = mod_names[[v]],
            feature = nm[top_idx],
            loading = w[top_idx],
            sign_str = ifelse(w[top_idx] >= 0, "positive", "negative"),
            abs_loading = abs(w[top_idx]),
            stringsAsFactors = FALSE
        )
    })
    df <- do.call(rbind, rows)

    # Order features by loading within each panel; the key keeps features
    # with the same name in different modalities apart.
    df$key <- paste0(match(df$modality_label, mod_names), "___", df$feature)
    df$key <- factor(df$key, levels = df$key[order(
        match(df$modality_label, mod_names), df$loading
    )])
    df$modality_label <- factor(df$modality_label,
        levels = mod_names[mod_seq]
    )

    if (is.null(main)) {
        main <- sprintf("Feature Loadings - Dimension %d", dim)
    }

    n_mod <- length(mod_seq)
    ncol_facets <- min(n_mod, 3L)

    p_out <- ggplot2::ggplot(
        df, ggplot2::aes(x = loading, y = key, fill = sign_str)
    ) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::geom_vline(xintercept = 0, color = "gray30", linewidth = 0.4) +
        ggplot2::facet_wrap(
            ~modality_label, scales = "free", ncol = ncol_facets
        ) +
        ggplot2::scale_y_discrete(labels = .axis_labels()) +
        ggplot2::scale_fill_manual(
            values = c("positive" = col_pos, "negative" = col_neg),
            labels = c("positive" = "Positive", "negative" = "Negative")
        ) +
        ggplot2::labs(
            x = sprintf("Loading (Dimension %d)", dim),
            y = NULL,
            fill = NULL,
            title = main,
            subtitle = "Each panel has its own scale"
        ) +
        .theme_woven() +
        ggplot2::theme(
            panel.grid.major.y = ggplot2::element_blank(),
            legend.position    = "bottom"
        )
    p_out
}

#' Plot the strength of each WOVEN latent dimension
#'
#' Bar chart of each latent dimension's share of the squared supervised
#' canonical correlations (\code{singular_values(fit)^2}), with the
#' cumulative share as a line. A dimension with a small share adds little
#' shared, class-related signal, which helps when choosing \code{K}. The
#' shares describe the cross-modal alignment, not the variance of the data.
#'
#' @param fit woven object from [woven()]
#' @param main character: plot title
#' @return a ggplot object (printed automatically; add layers with \code{+})
#' @examples
#' data(woven_example)
#' fit <- woven(woven_example$X_missing,
#'     Y = woven_example$Y, K = 4L,
#'     verbose = FALSE
#' )
#' woven_plot_variance(fit)
#' @seealso [woven_plot_vip()], [woven_plot_loadings()]
#' @export
woven_plot_variance <- function(fit,
                                main = "Strength of each latent dimension") {
    .require_ggplot2()
    stopifnot(methods::is(fit, "woven"))
    svals <- fit@singular_values
    K <- fit@K
    share <- 100 * svals^2 / sum(svals^2)
    df <- data.frame(
        dim_int = seq_len(K),
        pct_var = share,
        cumulative = cumsum(share)
    )
    col_bar <- "#0072B2"
    col_cum <- "#D55E00"

    p_out <- ggplot2::ggplot(df, ggplot2::aes(x = dim_int)) +
        ggplot2::geom_col(ggplot2::aes(y = pct_var),
            fill = col_bar, width = 0.6
        ) +
        ggplot2::geom_text(
            data = if (K <= 12L) df else df[0, , drop = FALSE],
            ggplot2::aes(y = pct_var, label = sprintf("%.0f%%", pct_var)),
            vjust = -0.5, size = 3.5, color = col_bar, fontface = "bold"
        ) +
        ggplot2::geom_line(
            data = if (K > 1L) df else df[0, , drop = FALSE],
            ggplot2::aes(y = cumulative), color = col_cum, linewidth = 0.9
        ) +
        ggplot2::geom_point(ggplot2::aes(y = cumulative),
            color = col_cum, size = 2.5
        ) +
        ggplot2::scale_x_continuous(
            breaks = if (K <= 15L) seq_len(K) else pretty(seq_len(K)),
            labels = function(b) paste0("Dim ", b)
        ) +
        ggplot2::scale_y_continuous(
            limits = c(0, 108), breaks = seq(0, 100, 25),
            labels = function(x) paste0(x, "%")
        ) +
        ggplot2::labs(
            x = NULL, y = "Share of squared canonical correlation",
            title = main,
            subtitle = "Bars: each dimension; line: cumulative"
        ) +
        .theme_woven() +
        ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
    p_out
}
