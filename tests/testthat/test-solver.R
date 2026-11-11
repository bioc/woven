# Tests for woven_mcca_dual() against small cases computed by hand.
#
# With lambdas = 0 the metric is M_v = (1 + 1e-8) I (the 1e-8 is the solver's
# numerical nudge), so every off-diagonal block of the dual matrix is
# IKY / (1 + 1e-8), with IKY = I + gamma_y * Ytilde Ytilde^T / n_a.
#
# For two balanced classes, the centered one-hot matrix Ytilde has columns
# c and -c with c_i = +1/2 (class 1) or -1/2 (class 2), so
# Ytilde Ytilde^T / n_a = 2 c c^T / n_a has one non-zero eigenvalue,
# 2 * (n_a / 4) / n_a = 1/2, with eigenvector u = c / ||c||.
# Hence the top eigenvalue of IKY is 1 + gamma_y / 2.
#
# The dual matrix is (J_V - I_V) kron IKY / (1 + 1e-8), whose top eigenvalue
# is (V - 1) * (1 + gamma_y / 2) / (1 + 1e-8), with eigenvector
# (u, ..., u) / sqrt(V). The M-orthonormal anchor scores are then
# Z_v = +/- u / sqrt(1 + 1e-8) for every modality, i.e. +/- 1 / sqrt(n_a)
# with the sign given by the class.

hand_fixture <- function(V, n_a = 20L, p = c(30L, 25L, 40L), seed = 1L) {
    set.seed(seed)
    X_list <- lapply(p[seq_len(V)], function(pv) matrix(rnorm(n_a * pv), n_a, pv))
    list(X_list = X_list, Y = rep(1:2, each = n_a / 2), n_a = n_a)
}

test_that("V = 2: top eigenvalue and anchor scores match the hand computation", {
    d <- hand_fixture(V = 2L)
    gamma <- 2
    fit <- suppressWarnings(woven_mcca_dual(d$X_list,
        anchor_idx = seq_len(d$n_a), Y = d$Y, K = 1L,
        lambdas = 0, gamma_y = gamma, verbose = FALSE
    ))
    expect_equal(fit$singular_values, (1 + gamma / 2) / (1 + 1e-8),
        tolerance = 1e-10
    )
    c_vec <- ifelse(d$Y == 1L, 0.5, -0.5)
    u <- c_vec / sqrt(sum(c_vec^2))
    for (v in 1:2) {
        z <- fit$Za_list[[v]][, 1]
        s <- sign(sum(z * u))
        expect_equal(z, s * u / sqrt(1 + 1e-8), tolerance = 1e-8)
    }
    # both modalities share the same eigenvector block, hence the same sign
    expect_equal(fit$Za_list[[1]], fit$Za_list[[2]], tolerance = 1e-8)
})

test_that("V = 3: top eigenvalue is (V - 1)(1 + gamma/2)", {
    d <- hand_fixture(V = 3L)
    gamma <- 2
    fit <- suppressWarnings(woven_mcca_dual(d$X_list,
        anchor_idx = seq_len(d$n_a), Y = d$Y, K = 1L,
        lambdas = 0, gamma_y = gamma, verbose = FALSE
    ))
    expect_equal(fit$singular_values, 2 * (1 + gamma / 2) / (1 + 1e-8),
        tolerance = 1e-10
    )
    expect_length(fit$W_list, 3L)
    expect_equal(vapply(fit$W_list, nrow, integer(1)), c(30L, 25L, 40L))
    for (v in 2:3) {
        expect_equal(fit$Za_list[[v]], fit$Za_list[[1]], tolerance = 1e-8)
    }
})

test_that("V = 3 with lambdas > 0: scores are M-orthonormal and eigenvalues sorted", {
    d <- hand_fixture(V = 3L, n_a = 30L)
    fit <- suppressWarnings(woven_mcca_dual(d$X_list,
        anchor_idx = seq_len(d$n_a), Y = rep(1:3, 10), K = 3L,
        lambdas = c(0.1, 1, 5), gamma_y = 1, k_nn = 5L, verbose = FALSE
    ))
    expect_equal(fit$K, 3L)
    expect_false(is.unsorted(rev(fit$singular_values)))
    for (v in 1:3) {
        Zv <- fit$Za_list[[v]]
        expect_equal(crossprod(Zv, fit$M_list[[v]] %*% Zv), diag(3),
            tolerance = 1e-4, ignore_attr = TRUE
        )
    }
})

test_that("projection matrices interpolate the anchor scores when p_v > n_a", {
    d <- hand_fixture(V = 2L)
    fit <- suppressWarnings(woven_mcca_dual(d$X_list,
        anchor_idx = seq_len(d$n_a), Y = d$Y, K = 2L,
        lambdas = 0.1, k_nn = 5L, verbose = FALSE
    ))
    for (v in 1:2) {
        expect_equal(fit$Xa_list[[v]] %*% fit$W_list[[v]], fit$Za_list[[v]],
            tolerance = 1e-5
        )
    }
})

test_that("screen_top zeroes the rows of non-selected features and keeps signal", {
    set.seed(2)
    n_a <- 40L
    Y <- rep(1:2, each = n_a / 2)
    X1 <- matrix(rnorm(n_a * 50), n_a, 50)
    X1[Y == 2, 1:3] <- X1[Y == 2, 1:3] + 4 # features 1-3 carry the signal
    X2 <- matrix(rnorm(n_a * 20), n_a, 20)
    fit <- suppressWarnings(woven_mcca_dual(list(X1, X2),
        anchor_idx = seq_len(n_a), Y = Y, K = 1L, screen_top = 5L,
        k_nn = 5L, verbose = FALSE
    ))
    nz <- which(rowSums(abs(fit$W_list[[1]])) > 0)
    expect_length(nz, 5L)
    expect_true(all(1:3 %in% nz))
    expect_equal(dim(fit$W_list[[1]]), c(50L, 1L))
    expect_setequal(woven:::.anova_screen(X1, Y, 3L), 1:3)
    expect_equal(woven:::.anova_screen(X1, Y, Inf), seq_len(50))
})

test_that("larger ridge_w shrinks the projection matrices", {
    d <- hand_fixture(V = 2L)
    norm_w <- function(r) {
        f <- suppressWarnings(woven_mcca_dual(d$X_list,
            anchor_idx = seq_len(d$n_a), Y = d$Y, K = 1L,
            ridge_w = r, k_nn = 5L, verbose = FALSE
        ))
        sum(f$W_list[[1]]^2)
    }
    expect_gt(norm_w(1e-8), norm_w(1))
    expect_gt(norm_w(1), norm_w(10))
})

test_that("K is reduced with a warning when there are too few positive eigenvalues", {
    # For V = 2 the dual matrix always has n_a positive eigenvalues, so the
    # reduction is exercised on the eigen step directly.
    expect_warning(
        ev <- woven:::.top_eigen(diag(c(3, 1, -1, -2)), K = 3L),
        "Only 2 positive eigenvalues; K reduced 3 -> 2"
    )
    expect_equal(ev$values, c(3, 1))
    expect_equal(dim(ev$vectors), c(4L, 2L))
})

test_that("woven_mcca_dual() validates its inputs", {
    d <- hand_fixture(V = 2L)
    expect_error(
        woven_mcca_dual(d$X_list[1], anchor_idx = 1:20, Y = d$Y, K = 1L),
        "at least 2 modalities"
    )
    expect_error(
        woven_mcca_dual(d$X_list,
            anchor_idx = 1:20, Y = d$Y[-1], K = 1L,
            verbose = FALSE
        ),
        "length\\(Y\\)"
    )
    expect_error(
        woven_mcca_dual(d$X_list,
            anchor_idx = 1:20, Y = d$Y, K = 1L,
            lambdas = c(1, 2, 3), verbose = FALSE
        ),
        "lambdas"
    )
    expect_error(
        woven_mcca_dual(d$X_list, anchor_idx = c(1:19, 25L), Y = d$Y, K = 1L),
        "outside"
    )
})

test_that("woven_mcca_dual() prints progress when verbose", {
    d <- hand_fixture(V = 2L)
    expect_message(
        suppressWarnings(woven_mcca_dual(d$X_list,
            anchor_idx = seq_len(20L), Y = d$Y, K = 1L, k_nn = 5L
        )),
        "MCCA-dual"
    )
})

test_that("internal step helpers return the expected structures", {
    Y_a <- rep(1:2, each = 4)
    IKY <- woven:::.label_kernel(Y_a, gamma_y = 2)
    expect_equal(dim(IKY), c(8L, 8L))
    expect_equal(max(eigen(IKY, symmetric = TRUE)$values), 2, tolerance = 1e-12)

    P <- woven:::.block_dual_matrix(list(diag(8), diag(8)), IKY)
    expect_equal(P[1:8, 9:16], IKY, ignore_attr = TRUE)
    expect_equal(P[1:8, 1:8], matrix(0, 8, 8))
    expect_true(isSymmetric(P))

    roots <- woven:::.metric_roots(list(diag(8) * 0), lambdas = 1, n_a = 8L)
    expect_equal(roots$sqrt[[1]] %*% roots$sqrt[[1]], roots$M[[1]],
        tolerance = 1e-12
    )
    expect_equal(roots$inv_sqrt[[1]] %*% roots$M[[1]] %*% roots$inv_sqrt[[1]],
        diag(8),
        tolerance = 1e-10
    )
    expect_error(woven:::.top_eigen(-diag(3), 1L), "No positive eigenvalues")
})
