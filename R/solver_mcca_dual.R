# solver_mcca_dual.R -- closed-form dual SUMCOR MCCA for WOVEN (V >= 2)
#
# SUMCOR ("sum of correlations") multi-set CCA maximizes the sum of pairwise
# correlations between the latent scores of all modality pairs. In the dual
# (subject) space this has a closed-form solution from one eigendecomposition.
#
# Derivation (anchor subjects only, n_a of them):
#   M_v  = I + lambda_v L_a^(v) / n_a       graph-regularized metric, n_a x n_a
#   IKY  = I + gamma_y K_Y                     label-augmented kernel, n_a x n_a
#   objective: max sum_{v<u} tr(Z_v^T IKY Z_u)  s.t.  Z_v^T M_v Z_v = I_K
#
#   Block dual matrix P, (V n_a) x (V n_a), symmetric:
#     P_vu = M_v^{-1/2} IKY M_u^{-1/2}   (v != u),   P_vv = 0
#   The top-K eigenvectors phi of P give Z_v = M_v^{-1/2} phi_v.
#   Projection matrices: W_v = X_a^T (X_a X_a^T + ridge)^{-1} Z_v (min-norm).
#
# Cost: building P is O(V^2 n_a^2); the eigendecomposition is O((V n_a)^3);
# W recovery is O(V n_a^2 p_v).
#
# The solver is split into one internal helper per step:
#   .anchor_laplacians()  -> n_a x n_a anchor Laplacians
#   .metric_roots()       -> M_v, M_v^{1/2}, M_v^{-1/2}, chol(M_v)
#   .label_kernel()       -> IKY
#   .block_dual_matrix()  -> P
#   .top_eigen()          -> top-K positive eigenpairs of P
#   .anchor_scores()      -> M-orthonormal anchor scores Z_v
#   .anova_screen()       -> optional supervised feature screening
#   .recover_W()          -> projection matrix W_v

#' Closed-form dual SUMCOR MCCA solver (low-level)
#'
#' The solver that \code{\link{woven}} calls internally. Most users should
#' call \code{woven()}, which adds input standardization, automatic anchor
#' detection, MultiAssayExperiment input, latent scores for every subject,
#' and an S4 result object. \code{woven_mcca_dual()} is exported for
#' method developers who want the raw solver output (anchor scores, metric
#' matrices, anchor Laplacians) without that wrapping.
#'
#' The fit is a single eigendecomposition of a (V n_a) x (V n_a) block
#' matrix, with no iterations, random restarts or local optima. SUMCOR
#' ("sum of correlations") is the multi-set CCA criterion that maximizes the
#' sum of pairwise correlations between the modalities' latent scores.
#'
#' @param X_list list of V matrices, each n x p_v. Rows that are entirely
#'   \code{NA} mark subjects missing that modality (block-missing).
#' @param anchor_idx integer vector of subjects observed in every modality.
#'   Validated: indices must be distinct, within 1..n, and not block-missing
#'   in any modality.
#' @param Y vector of length n of class labels. Only \code{Y[anchor_idx]} is
#'   used.
#' @param K integer, number of latent dimensions. Each projection matrix has
#'   rank at most \code{min(n_a, p_v)}; a warning is given when \code{K}
#'   exceeds a modality's feature count.
#' @param lambdas numeric scalar or length-V vector, Laplacian regularization.
#' @param gamma_y numeric >= 0, label supervision strength.
#' @param k_nn integer, k-NN graph size (ignored if \code{La_list_precomp}
#'   is supplied).
#' @param La_list_precomp optional list of V pre-extracted n_a x n_a anchor
#'   Laplacians.
#' @param ridge_w numeric >= 0, relative ridge strength in W recovery
#'   (\code{ridge_w * max(diag(Ka_v))} is added to the diagonal of the anchor
#'   Gram matrix \code{Ka_v}). The default 1e-8 only stabilizes the solve.
#' @param screen_top integer or \code{Inf}; if finite, W recovery uses only
#'   the \code{screen_top} features per modality with the largest ANOVA
#'   F-statistic against \code{Y} on anchor subjects. Screened-out features
#'   get a zero row in W.
#' @param verbose logical, print progress messages.
#'
#' @return a list with components \code{W_list} (V matrices, p_v x K),
#'   \code{Za_list} (V anchor score matrices, n_a x K), \code{Xa_list}
#'   (anchor data used), \code{singular_values} (length K),
#'   \code{M_list}, \code{L_a_list}, \code{anchor_idx}, \code{K}, \code{V},
#'   \code{lambdas}, \code{gamma_y}, \code{k_nn}, \code{n}, \code{p_v}.
#' @seealso [woven()]
#' @examples
#' set.seed(1)
#' n <- 20
#' X1 <- matrix(rnorm(n * 5), n, 5)
#' X2 <- matrix(rnorm(n * 4), n, 4)
#' X3 <- matrix(rnorm(n * 3), n, 3)
#' Y <- rep(1:2, each = n / 2)
#' fit <- woven_mcca_dual(list(X1, X2, X3),
#'     anchor_idx = seq_len(14L), Y = Y,
#'     K = 2L, verbose = FALSE
#' )
#' fit$singular_values
#' @export
woven_mcca_dual <- function(X_list, anchor_idx, Y,
                            K = 5L,
                            lambdas = 0.1,
                            gamma_y = 1.0,
                            k_nn = 10L,
                            La_list_precomp = NULL,
                            ridge_w = 1e-8,
                            screen_top = Inf,
                            verbose = TRUE) {
    V <- length(X_list)
    n <- nrow(X_list[[1L]])
    check_woven_inputs(X_list, anchor_idx, K)
    anchor_idx <- as.integer(anchor_idx)
    n_a <- length(anchor_idx)
    if (length(Y) != n) {
        stop(sprintf("length(Y) = %d but X_list has %d rows.", length(Y), n),
            call. = FALSE
        )
    }
    if (length(lambdas) == 1L) lambdas <- rep(lambdas, V)
    if (length(lambdas) != V) {
        stop(sprintf("lambdas must have length 1 or V = %d.", V), call. = FALSE)
    }
    if (verbose) {
        message(sprintf(
            "WOVEN MCCA-dual | V=%d, n=%d, n_a=%d, K=%d, gamma_y=%.2f",
            V, n, n_a, K, gamma_y
        ))
    }

    Xa_raw <- lapply(X_list, function(X) X[anchor_idx, , drop = FALSE])
    col_ok <- lapply(Xa_raw, function(X) which(!.all_na_cols(X)))
    Xa_list <- lapply(Xa_raw, na_impute_median)
    Y_a <- Y[anchor_idx]

    L_a_list <- .anchor_laplacians(X_list, anchor_idx, k_nn, La_list_precomp)
    roots <- .metric_roots(L_a_list, lambdas, n_a)
    IKY <- .label_kernel(Y_a, gamma_y)
    P <- .block_dual_matrix(roots$inv_sqrt, IKY)
    eig <- .top_eigen(P, K)
    K_use <- length(eig$values)

    Za_list <- vector("list", V)
    W_list <- vector("list", V)
    for (v in seq_len(V)) {
        rows <- ((v - 1L) * n_a + 1L):(v * n_a)
        Za_list[[v]] <- .anchor_scores(
            eig$vectors[rows, , drop = FALSE],
            roots$inv_sqrt[[v]], roots$chol[[v]]
        )
        W_list[[v]] <- .recover_W(
            Xa_list[[v]], Za_list[[v]], Y_a, ridge_w, screen_top
        )
    }

    if (verbose) {
        message(sprintf(
            "  Done. Top eigenvalues: %s",
            paste(round(eig$values[seq_len(min(5L, K_use))], 4),
                collapse = ", "
            )
        ))
    }

    list(
        W_list = W_list,
        Z_list = Za_list,
        Za_list = Za_list,
        Xa_list = Xa_list,
        col_ok_list = col_ok,
        singular_values = eig$values,
        M_list = roots$M,
        L_a_list = L_a_list,
        anchor_idx = anchor_idx,
        K = K_use, V = V,
        lambdas = lambdas, gamma_y = gamma_y, k_nn = k_nn,
        n = n, p_v = vapply(X_list, ncol, integer(1L)),
        dual = TRUE
    )
}

# Step 1: n_a x n_a anchor Laplacians, either extracted from precomputed
# full-cohort Laplacians or built from each modality's observed rows.
.anchor_laplacians <- function(X_list, anchor_idx, k_nn, La_list_precomp) {
    if (!is.null(La_list_precomp)) {
        return(lapply(La_list_precomp, as.matrix))
    }
    lapply(X_list, function(X) {
        as.matrix(build_laplacian(X, k = k_nn)[anchor_idx, anchor_idx])
    })
}

# Step 2: per-modality metric M_v = I + lambda_v L_a / n_a (with a tiny
# diagonal nudge for numerical stability), its square roots and Cholesky
# factor.
.metric_roots <- function(L_a_list, lambdas, n_a) {
    M <- lapply(seq_along(L_a_list), function(v) {
        Mv <- diag(n_a) + lambdas[v] * L_a_list[[v]] / n_a
        Mv + diag(1e-8 * max(diag(Mv)), n_a)
    })
    sq <- lapply(M, mat_sqrt)
    list(
        M = M,
        sqrt = lapply(sq, `[[`, "sqrt"),
        inv_sqrt = lapply(sq, `[[`, "inv_sqrt"),
        chol = lapply(M, chol)
    )
}

# Step 3: label-augmented kernel IKY = I + gamma_y * Ytilde Ytilde^T / n_a,
# where Ytilde is the column-centered one-hot label matrix.
.label_kernel <- function(Y_a, gamma_y) {
    n_a <- length(Y_a)
    Y_onehot <- stats::model.matrix(~ 0 + factor(Y_a))
    Y_tilde <- scale(Y_onehot, center = TRUE, scale = FALSE)
    diag(n_a) + gamma_y * tcrossprod(Y_tilde) / n_a
}

# Step 4: symmetric block dual matrix with off-diagonal blocks
# M_v^{-1/2} IKY M_u^{-1/2} and zero diagonal blocks.
.block_dual_matrix <- function(inv_sqrt_list, IKY) {
    V <- length(inv_sqrt_list)
    n_a <- nrow(IKY)
    left <- lapply(inv_sqrt_list, function(Mi) Mi %*% IKY)
    P <- matrix(0.0, nrow = V * n_a, ncol = V * n_a)
    for (v in seq_len(V)) {
        rows <- ((v - 1L) * n_a + 1L):(v * n_a)
        for (u in seq_len(V)) {
            if (u == v) next
            cols <- ((u - 1L) * n_a + 1L):(u * n_a)
            P[rows, cols] <- left[[v]] %*% inv_sqrt_list[[u]]
        }
    }
    P
}

# Step 5: top-K strictly positive eigenpairs of P. Negative eigenvalues
# correspond to anti-correlated directions and are never kept. A full eigen()
# is used: P is small (V n_a), and iterative solvers can struggle with its
# block structure.
.top_eigen <- function(P, K) {
    ev <- eigen(P, symmetric = TRUE)
    pos <- which(ev$values > 1e-10)
    if (length(pos) == 0L) stop("No positive eigenvalues in the dual matrix.")
    K_use <- min(K, length(pos))
    if (K_use < K) {
        warning(sprintf(
            "Only %d positive eigenvalues; K reduced %d -> %d.",
            length(pos), K, K_use
        ), call. = FALSE)
    }
    keep <- pos[seq_len(K_use)]
    list(values = ev$values[keep], vectors = ev$vectors[, keep, drop = FALSE])
}

# Step 6: anchor scores Z_v = M_v^{-1/2} phi_v, M-orthonormalized so that
# Z_v^T M_v Z_v = I.
.anchor_scores <- function(phi_v, inv_sqrt_v, chol_v) {
    K <- ncol(phi_v)
    Z_raw <- inv_sqrt_v %*% phi_v
    ZMZ <- crossprod(chol_v %*% Z_raw)
    ZMZ <- (ZMZ + t(ZMZ)) / 2
    R <- tryCatch(chol(ZMZ), error = function(e) NULL)
    if (!is.null(R)) {
        return(t(backsolve(R, t(Z_raw))))
    }
    ev <- eigen(ZMZ, symmetric = TRUE)
    d_inv <- 1 / sqrt(pmax(ev$values[seq_len(K)], 1e-14))
    Z_raw %*% ev$vectors[, seq_len(K), drop = FALSE] %*% diag(d_inv, K)
}

# Step 7 (optional): indices of the screen_top features with the largest
# one-way ANOVA F-statistic against Y_a, or all features if screen_top is
# Inf or >= p_v.
.anova_screen <- function(Xa, Y_a, screen_top) {
    p_v <- ncol(Xa)
    if (!is.finite(screen_top) || screen_top >= p_v) {
        return(seq_len(p_v))
    }
    y_f <- as.factor(Y_a)
    n_a <- nrow(Xa)
    grand_mean <- colMeans(Xa)
    ss_between <- numeric(p_v)
    ss_within <- numeric(p_v)
    for (lev in levels(y_f)) {
        idx <- which(y_f == lev)
        if (length(idx) < 2L) next
        Xc <- Xa[idx, , drop = FALSE]
        mean_c <- colMeans(Xc)
        ss_between <- ss_between + length(idx) * (mean_c - grand_mean)^2
        ss_within <- ss_within + colSums(sweep(Xc, 2, mean_c, "-")^2)
    }
    df_b <- nlevels(y_f) - 1L
    df_w <- n_a - nlevels(y_f)
    f <- (ss_between / max(df_b, 1L)) / (ss_within / max(df_w, 1L) + 1e-12)
    f[!is.finite(f)] <- 0
    order(f, decreasing = TRUE)[seq_len(screen_top)]
}

# Step 8: min-norm projection matrix W_v = X_s^T (X_s X_s^T + ridge)^{-1} Z_v
# on the screened columns X_s; screened-out features get zero rows, so W_v is
# always p_v x K.
.recover_W <- function(Xa, Za, Y_a, ridge_w, screen_top) {
    keep <- .anova_screen(Xa, Y_a, screen_top)
    Xs <- Xa[, keep, drop = FALSE]
    Ka <- tcrossprod(Xs)
    Ka_reg <- Ka + diag(ridge_w * max(diag(Ka)), nrow(Ka))
    W_s <- t(Xs) %*% tryCatch(
        solve(Ka_reg, Za),
        error = function(e) {
            sv <- svd(Ka_reg)
            sv$v %*% diag(1 / pmax(sv$d, 1e-10 * max(sv$d)), length(sv$d)) %*%
                t(sv$u) %*% Za
        }
    )
    W <- matrix(0, ncol(Xa), ncol(Za))
    W[keep, ] <- W_s
    W
}
