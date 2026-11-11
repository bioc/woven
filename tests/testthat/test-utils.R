# Internal helpers in utils.R and laplacian.R.

test_that(".block_missing_rows() flags only all-NA rows", {
    X <- matrix(1:12, 4, 3)
    X[2, ] <- NA
    X[3, 1] <- NA
    expect_identical(woven:::.block_missing_rows(X), c(FALSE, TRUE, FALSE, FALSE))
})

test_that("na_impute_median() drops block-missing rows and all-NA columns", {
    X <- matrix(c(1, 2, NA, 4, 5, 6, 7, 8, NA, NA, NA, NA), 4, 3)
    X[4, ] <- NA
    out <- woven:::na_impute_median(X)
    expect_equal(dim(out), c(3L, 2L))
    expect_equal(out[3, 1], 1.5) # median of 1, 2
    expect_false(anyNA(out))
})

test_that("mat_sqrt() returns square root and inverse square root", {
    set.seed(1)
    A <- crossprod(matrix(rnorm(30), 10, 3)) + diag(3)
    r <- woven:::mat_sqrt(A)
    expect_equal(r$sqrt %*% r$sqrt, A, tolerance = 1e-10)
    expect_equal(r$inv_sqrt %*% A %*% r$inv_sqrt, diag(3), tolerance = 1e-10)
    # a singular matrix is floored, not inverted to Inf
    r0 <- woven:::mat_sqrt(diag(c(1, 0)))
    expect_true(all(is.finite(r0$inv_sqrt)))
})

test_that(".scale_fit() / .scale_apply() round-trip and keep NA rows", {
    set.seed(2)
    X <- matrix(rnorm(40, mean = 5, sd = 3), 10, 4)
    X[, 4] <- 7 # constant column keeps scale 1
    X[3, ] <- NA
    sf <- woven:::.scale_fit(list(a = X))
    expect_equal(unname(colMeans(sf$X$a[, 1:3], na.rm = TRUE)), rep(0, 3), tolerance = 1e-12)
    expect_equal(sf$scale[[1]][4], 1)
    expect_true(all(is.na(sf$X$a[3, ])))
    again <- woven:::.scale_apply(list(a = X), sf$center, sf$scale)
    expect_equal(again, sf$X)
})

test_that("build_laplacian() errors with no observed rows", {
    X <- matrix(NA_real_, 5, 3)
    expect_error(woven:::build_laplacian(X), "No observed rows")
})

test_that("build_laplacian() caps k at n_obs - 1 and uses a fixed sigma", {
    set.seed(3)
    X <- matrix(rnorm(5 * 4), 5, 4)
    L <- woven:::build_laplacian(X, k = 50L, sigma = 2)
    expect_equal(dim(L), c(5L, 5L))
    expect_true(all(abs(Matrix::rowSums(L)) < 1e-10))
    # identical rows: zero distances fall back to sigma = 1
    L1 <- woven:::build_laplacian(matrix(1, 4, 2), k = 2L)
    expect_true(all(is.finite(as.matrix(L1))))
})

test_that("build_laplacian() accepts a precomputed kernel for observed rows", {
    set.seed(4)
    X <- matrix(rnorm(6 * 3), 6, 3)
    X[2, ] <- NA
    Kmat <- exp(-as.matrix(dist(X[-2, ]))^2)
    L <- woven:::build_laplacian(X, kernel_matrix = Kmat)
    expect_equal(dim(L), c(6L, 6L))
    expect_equal(as.numeric(L[2, ]), rep(0, 6))
    expect_true(all(abs(Matrix::rowSums(L)) < 1e-10))
    expect_error(
        woven:::build_laplacian(X, kernel_matrix = Kmat[1:3, 1:3]),
        "kernel_matrix must be 5 x 5"
    )
})

test_that("woven_precompute() rejects non-list input and records settings", {
    expect_error(woven_precompute(matrix(1, 3, 3)), "list of matrices")
    set.seed(5)
    X <- list(a = matrix(rnorm(60), 20, 3), b = matrix(rnorm(40), 20, 2))
    pc <- woven_precompute(X, k_nn = 4L)
    expect_named(pc, c("a", "b"))
    expect_identical(attr(pc, "woven_k_nn"), 4L)
    expect_identical(attr(pc, "woven_n"), 20L)
})
