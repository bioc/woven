# metrics.R  -- evaluation metrics for WOVEN and other latent embeddings
#
# All functions take a latent matrix Z (n x K) and return a scalar or named vector.
# Designed to work on: WOVEN output, DIABLO variate scores, MOFA+ factor scores.
#
# Metrics:
#   woven_silhouette       -- average silhouette width (cluster geometry)
#   woven_davies_bouldin   -- Davies-Bouldin index (cluster compactness/separation)
#   woven_nmi              -- normalized mutual information (label recovery)
#   woven_rv               -- RV coefficient vs ground-truth factor matrix
#   woven_ess_retention    -- effective sample size retention (N_used / N_total)
#   woven_effect_bias      -- subgroup treatment-effect estimate bias
#   woven_nystrom_error    -- leave-anchor-out Nystrom projection error
#   woven_all_metrics      -- compute full battery, returns named list

#  Silhouette

#' Average silhouette width
#'
#' Works on any n x K latent score matrix, so the same metric can be
#' computed for WOVEN and for embeddings from other methods (for example
#' DIABLO variates) when comparing them. [woven_metrics()] computes it
#' directly from a fit.
#'
#' @param Z numeric matrix n x K (latent scores)
#' @param labels integer or factor of length n (subgroup labels)
#' @return scalar in [-1, 1], higher is better
#' @examples
#' set.seed(1)
#' Z <- matrix(rnorm(20 * 2), 20, 2)
#' labels <- rep(1:2, each = 10)
#' woven_silhouette(Z, labels)
#' @export
woven_silhouette <- function(Z, labels) {
    labels <- as.integer(as.factor(labels))
    if (length(unique(labels)) < 2L) {
        return(NA_real_)
    }
    d <- stats::dist(Z)
    sw <- cluster::silhouette(labels, d)
    mean(sw[, "sil_width"])
}

#  Davies-Bouldin

#' Davies-Bouldin index
#'
#' Works on any n x K latent score matrix, so the same metric can be
#' computed for WOVEN and for embeddings from other methods (for example
#' DIABLO variates) when comparing them. [woven_metrics()] computes it
#' directly from a fit.
#'
#' DB = (1/K) sum_i max_\{j != i\} (s_i + s_j) / d(c_i, c_j)
#' where s_i = mean intra-cluster distance, d(c_i, c_j) = centroid distance.
#'
#' @param Z numeric matrix n x K
#' @param labels integer or factor of length n
#' @return scalar >= 0, lower is better
#' @examples
#' set.seed(1)
#' Z <- matrix(rnorm(20 * 2), 20, 2)
#' labels <- rep(1:2, each = 10)
#' woven_davies_bouldin(Z, labels)
#' @export
woven_davies_bouldin <- function(Z, labels) {
    labels <- as.integer(as.factor(labels))
    K_cl <- length(unique(labels))
    if (K_cl < 2L) {
        return(NA_real_)
    }

    centroids <- do.call(rbind, lapply(sort(unique(labels)), function(g) {
        colMeans(Z[labels == g, , drop = FALSE])
    }))

    s <- vapply(sort(unique(labels)), function(g) {
        Zg <- Z[labels == g, , drop = FALSE]
        if (nrow(Zg) < 2L) {
            return(0)
        }
        mean(sqrt(rowSums(sweep(Zg, 2, centroids[g, ])^2)))
    }, numeric(1L))

    db_vals <- vapply(seq_len(K_cl), function(i) {
        others <- setdiff(seq_len(K_cl), i)
        max(vapply(others, function(j) {
            d_ij <- sqrt(sum((centroids[i, ] - centroids[j, ])^2))
            if (d_ij < 1e-12) {
                return(0)
            }
            (s[i] + s[j]) / d_ij
        }, numeric(1L)))
    }, numeric(1L))

    mean(db_vals)
}

#  NMI

#' Normalized mutual information between cluster assignments and true labels
#'
#' Works on any n x K latent score matrix, so the same metric can be
#' computed for WOVEN and for embeddings from other methods (for example
#' DIABLO variates) when comparing them. [woven_metrics()] computes it
#' directly from a fit.
#'
#' Clusters Z with k-means and computes the normalized mutual information
#' between the cluster assignments and \code{labels}. k-means uses
#' \code{n_start} random starts, so the result depends on R's random number
#' generator: call \code{set.seed()} first for a reproducible value (the
#' package itself never sets the seed). [woven_metrics()] and
#' \code{summary()} call this function and inherit the same behaviour.
#'
#' @param Z numeric matrix n x K
#' @param labels integer or factor of length n (true labels)
#' @param n_cl integer, number of clusters (default = number of unique labels)
#' @param n_start integer, k-means random starts
#' @return scalar in [0, 1], higher is better
#' @examples
#' Z <- rbind(matrix(rnorm(40), 20, 2), matrix(rnorm(40, mean = 3), 20, 2))
#' labels <- rep(1:2, each = 20)
#' set.seed(1) # k-means is random; set the seed for a reproducible result
#' woven_nmi(Z, labels)
#' @export
woven_nmi <- function(Z, labels, n_cl = NULL, n_start = 10L) {
    labels <- as.integer(as.factor(labels))
    if (is.null(n_cl)) n_cl <- length(unique(labels))
    if (n_cl < 2L) {
        return(NA_real_)
    }

    km <- stats::kmeans(Z, centers = n_cl, nstart = n_start, iter.max = 100L)
    pred <- km$cluster

    # NMI via entropy decomposition (no external package needed)
    .entropy <- function(x) {
        px <- tabulate(x) / length(x)
        px <- px[px > 0]
        -sum(px * log(px))
    }
    .joint_entropy <- function(x, y) {
        n <- length(x)
        tbl <- table(x, y)
        pxy <- tbl / n
        pxy <- pxy[pxy > 0]
        -sum(pxy * log(pxy))
    }

    H_true <- .entropy(labels)
    H_pred <- .entropy(pred)
    H_joint <- .joint_entropy(labels, pred)
    MI <- H_true + H_pred - H_joint

    denom <- (H_true + H_pred) / 2
    if (denom < 1e-12) {
        return(NA_real_)
    }
    MI / denom
}

#  RV coefficient

#' RV coefficient between latent scores and ground-truth factor matrix
#'
#' Works on any n x K latent score matrix, e.g. to compare an embedding
#' with known factors in a simulation.
#'
#' RV(X, Y) = trace(X X' Y Y') / sqrt(trace(X X' X X') * trace(Y Y' Y Y'))
#' Measures similarity of two cross-product matrices; 1 = identical subspace.
#'
#' @param Z numeric matrix n x K (inferred latent scores)
#' @param Z_true numeric matrix n x K_true (e.g. ground-truth factor scores
#'   in a simulation)
#' @return scalar in [0, 1], higher is better
#' @examples
#' set.seed(1)
#' Z <- matrix(rnorm(20 * 2), 20, 2)
#' Z_true <- matrix(rnorm(20 * 3), 20, 3)
#' woven_rv(Z, Z_true)
#' @export
woven_rv <- function(Z, Z_true) {
    stopifnot(nrow(Z) == nrow(Z_true))
    # Center columns
    Z <- scale(Z, center = TRUE, scale = FALSE)
    Z_true <- scale(Z_true, center = TRUE, scale = FALSE)

    S <- tcrossprod(Z) # n x n
    T_ <- tcrossprod(Z_true) # n x n

    num <- sum(S * T_)
    denom <- sqrt(sum(S * S) * sum(T_ * T_))
    if (denom < 1e-12) {
        return(NA_real_)
    }
    num / denom
}

#  Effective sample size retention

#' Effective sample size (ESS) retention
#'
#' The fraction of the cohort that a method can score: \code{n_used /
#' n_total}. For WOVEN, every subject with at least one observed modality is
#' scored. A complete-case method such as DIABLO scores only subjects
#' observed in every modality, so its retention equals the anchor fraction.
#'
#' @param n_used integer, number of subjects with a latent score
#' @param n_total integer, total subjects in dataset
#' @return scalar in [0, 1], higher is better
#' @examples
#' woven_ess_retention(n_used = 80, n_total = 100)
#' @export
woven_ess_retention <- function(n_used, n_total) {
    stopifnot(n_used >= 0, n_total > 0, n_used <= n_total)
    n_used / n_total
}

#  Subgroup effect estimate bias

#' Subgroup treatment-effect estimate bias
#'
#' For simulation studies where the true treatment effect in each subgroup is
#' known. Within each subgroup defined by \code{labels}, fits a linear model
#' of a continuous outcome on a binary treatment indicator, adjusting for the
#' latent scores \code{Z}, and compares the estimated treatment effect with
#' the true effect. Subgroups with fewer than 5 subjects are skipped.
#'
#' bias_g = |estimated_g - true_g| / |true_g|   (relative)
#' Returns mean bias across subgroups.
#'
#' @param Z numeric matrix n x K (latent scores; used as covariates)
#' @param outcome numeric vector of length n (simulated continuous outcome)
#' @param treatment integer/logical vector of length n (0/1 treatment indicator)
#' @param labels integer or factor of length n (subgroup labels)
#' @param true_effects numeric vector, true treatment effect per subgroup, in
#'   the order of \code{levels(as.factor(labels))}
#' @return scalar >= 0, lower is better
#' @examples
#' set.seed(1)
#' n <- 200
#' Z <- matrix(rnorm(n * 3), n, 3)
#' treatment <- rep(0:1, n / 2)
#' labels <- rep(1:2, each = n / 2)
#' true_eff <- c(0.5, 1.0)
#' outcome <- true_eff[labels] * treatment + Z[, 1] + rnorm(n, sd = 0.5)
#' woven_effect_bias(Z, outcome, treatment, labels, true_eff)
#' @export
woven_effect_bias <- function(Z, outcome, treatment, labels, true_effects) {
    labels <- as.integer(as.factor(labels))
    treatment <- as.numeric(treatment)
    groups <- sort(unique(labels))

    biases <- vapply(seq_along(groups), function(gi) {
        g <- groups[gi]
        idx <- which(labels == g)
        if (length(idx) < 5L) {
            return(NA_real_)
        }

        df <- data.frame(y = outcome[idx], trt = treatment[idx], Z[idx, , drop = FALSE])
        fit <- tryCatch(stats::lm(y ~ ., data = df), error = function(e) NULL)
        if (is.null(fit)) {
            return(NA_real_)
        }

        est <- stats::coef(fit)[["trt"]]
        true <- true_effects[gi]
        if (is.na(true) || abs(true) < 1e-12) {
            return(NA_real_)
        }
        abs(est - true) / abs(true)
    }, numeric(1L))

    mean(biases, na.rm = TRUE)
}

#  Nystrom leave-anchor-out error

#' Leave-anchor-out Nystrom projection error
#'
#' A diagnostic for how stable the projection of an unseen subject is. For
#' each of \code{n_loo} randomly chosen anchor subjects, WOVEN is refitted
#' without that subject (same hyperparameters as \code{fit}), the held-out
#' subject is projected with the refitted projection matrices (the linear
#' Nystrom out-of-sample extension used for block-missing subjects), and the
#' Euclidean distance to its latent score in the full fit
#' (\code{latent_scores(fit)}) is recorded. The function
#' returns the mean distance. Small values mean that a subject's latent score
#' does not depend much on whether it was in the training set.
#'
#' Because latent dimensions are identified only up to sign, each refitted
#' dimension is sign-aligned to the full fit before the distance is taken.
#' The anchors to hold out are drawn with \code{sample()}, so the result
#' depends on R's random number generator: call \code{set.seed()} first for
#' a reproducible value (the package itself never sets the seed).
#'
#' @param fit a woven object from [woven()]
#' @param X_list the same list of modality matrices passed to [woven()],
#'   on the original (unscaled) feature scale
#' @param n_loo integer, number of anchors to hold out, between 1 and the
#'   number of anchors (default \code{min(20, n_a)})
#' @return scalar >= 0, lower is better (mean distance per held-out anchor)
#' @examples
#' data(woven_example)
#' fit <- woven(woven_example$X_complete,
#'     Y = woven_example$Y, K = 3L,
#'     verbose = FALSE
#' )
#' set.seed(1) # held-out anchors are sampled at random
#' woven_nystrom_error(fit, woven_example$X_complete, n_loo = 5L)
#' @export
woven_nystrom_error <- function(fit, X_list, n_loo = NULL) {
    if (!methods::is(fit, "woven")) stop("fit must be a woven object.")
    if (!is.list(X_list) || length(X_list) != fit@V ||
        nrow(X_list[[1L]]) != fit@n) {
        stop(sprintf(
            "X_list must be the list of %d matrices (n = %d) used to fit `fit`.",
            fit@V, fit@n
        ), call. = FALSE)
    }
    anchor_idx <- fit@anchor_idx
    n_a <- length(anchor_idx)
    if (is.null(n_loo)) n_loo <- min(20L, n_a)
    if (length(n_loo) != 1L || !is.numeric(n_loo) || is.na(n_loo) ||
        n_loo < 1L || n_loo > n_a || n_loo != round(n_loo)) {
        stop(sprintf(
            "n_loo must be a single integer between 1 and the number of anchors (%d).",
            n_a
        ), call. = FALSE)
    }
    if (isTRUE(fit@scaled) && !is.null(fit@scale_center)) {
        X_list <- .scale_apply(X_list, fit@scale_center, fit@scale_scale)
    }
    loo_set <- sample(seq_len(n_a), n_loo)

    errors <- vapply(loo_set, function(j) {
        held_out <- anchor_idx[j]
        remain_idx <- anchor_idx[-j]
        X_loo <- lapply(X_list, function(X) {
            X[held_out, ] <- NA_real_
            X
        })
        # Only anchor positions of Y are read by the solver.
        Y_full <- rep(NA_integer_, fit@n)
        Y_full[remain_idx] <- fit@Y_anchor[-j]
        mini <- tryCatch(
            # Input warnings (e.g. a small anchor set) were already given
            # when `fit` was made; muffle the repeats from each refit.
            withCallingHandlers(
                woven_mcca_dual(X_loo,
                    anchor_idx = remain_idx, Y = Y_full,
                    K = fit@K, lambdas = fit@lambdas, gamma_y = fit@gamma_y,
                    k_nn = fit@k_nn, ridge_w = fit@ridge_w,
                    screen_top = fit@screen_top, verbose = FALSE
                ),
                warning = function(w) invokeRestart("muffleWarning")
            ),
            error = function(e) NULL
        )
        if (is.null(mini) || mini$K != fit@K) {
            return(NA_real_)
        }
        # Sign-align each refitted dimension to the full fit, using the
        # remaining anchors' projected scores under both fits.
        X_rem <- lapply(X_list, function(X) X[remain_idx, , drop = FALSE])
        Z_rem <- .project_views(X_rem, mini$W_list, fit@K)
        sgn <- sign(colSums(fit@Z[remain_idx, , drop = FALSE] * Z_rem))
        sgn[sgn == 0] <- 1
        Xh <- lapply(X_list, function(X) X[held_out, , drop = FALSE])
        Z_proj <- sweep(.project_views(Xh, mini$W_list, fit@K), 2L, sgn, "*")
        sqrt(sum((Z_proj - fit@Z[held_out, , drop = FALSE])^2))
    }, numeric(1L))

    mean(errors, na.rm = TRUE)
}

#  Full metric battery

#' Compute the full metric battery on a latent score matrix
#'
#' @param Z numeric matrix n x K (consensus latent scores)
#' @param labels integer or factor of length n
#' @param n_total integer, total subjects before any filtering (for ESS)
#' @param Z_true optional numeric matrix n x K_true (ground-truth factors for RV)
#' @param outcome optional numeric vector (for effect bias)
#' @param treatment optional 0/1 vector (for effect bias)
#' @param true_effects optional named numeric (for effect bias)
#' @param fit optional WOVEN fit object (for the Nystrom LOO error)
#' @param X_list optional list of modality matrices used to fit \code{fit}
#'   (for the Nystrom LOO error)
#' @param n_loo integer, anchors to hold out for the Nystrom LOO error
#'   (default \code{min(20, n_a)})
#' @details \code{woven_nmi()} and \code{woven_nystrom_error()} are random;
#'   set the seed before calling for reproducible values.
#' @return named list of metric values
#' @noRd
woven_all_metrics <- function(
        Z, labels, n_total, Z_true = NULL, outcome = NULL, treatment = NULL,
        true_effects = NULL, fit = NULL, X_list = NULL, n_loo = NULL) {
    stopifnot(nrow(Z) == length(labels))

    out <- list(
        silhouette      = woven_silhouette(Z, labels),
        davies_bouldin  = woven_davies_bouldin(Z, labels),
        nmi             = woven_nmi(Z, labels),
        ess_retention   = woven_ess_retention(nrow(Z), n_total)
    )
    # Classification error is deliberately not computed here: a supervised
    # embedding encodes the labels into Z, so error measured on the same Z
    # is optimistic by construction. It must be estimated by refitting the
    # model within each cross-validation fold.

    if (!is.null(Z_true)) {
        out$rv_coefficient <- woven_rv(Z, Z_true)
    }

    if (!is.null(outcome) && !is.null(treatment) && !is.null(true_effects)) {
        out$effect_bias <- woven_effect_bias(Z, outcome, treatment, labels, true_effects)
    }

    if (!is.null(fit) && !is.null(X_list)) {
        out$nystrom_error <- woven_nystrom_error(fit, X_list, n_loo = n_loo)
    }

    out
}

#' Convenience wrapper: compute core metrics directly from a woven fit
#'
#' Computes silhouette width, Davies-Bouldin index, NMI and ESS retention on
#' the subjects the fit scored, using \code{latent_scores(fit)} and the fit's
#' subject count. NMI uses k-means, so call \code{set.seed()} first for
#' reproducible values.
#'
#' @param fit woven object from [woven()]
#' @param labels integer, factor, or character vector of length n with subgroup
#'   labels for all subjects (same Y passed to [woven()])
#' @return named numeric vector with elements \code{Silhouette},
#'   \code{Davies-Bouldin}, \code{NMI} and \code{ESS}
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
#' fit <- woven(list(X1, X2), Y = Y, K = K)
#' set.seed(1) # NMI uses k-means
#' woven_metrics(fit, Y)
#' @seealso [woven_silhouette()], [woven_davies_bouldin()], [woven_nmi()],
#'   [woven_ess_retention()]
#' @export
woven_metrics <- function(fit, labels) {
    stopifnot(methods::is(fit, "woven"))
    stopifnot(length(labels) == fit@n)

    scored <- !is.na(fit@Z[, 1L])
    res <- woven_all_metrics(
        Z       = fit@Z[scored, , drop = FALSE],
        labels  = labels[scored],
        n_total = fit@n
    )

    out <- unlist(res[c("silhouette", "davies_bouldin", "nmi", "ess_retention")])
    names(out) <- c("Silhouette", "Davies-Bouldin", "NMI", "ESS")
    out
}
