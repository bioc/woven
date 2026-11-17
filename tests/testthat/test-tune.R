# woven_tune(): cross-validated selection of gamma_y, ridge_w, screen_top.

test_that("the fast fold path equals refitting woven() per setting", {
    data(woven_example, package = "woven")
    X <- woven_example$X_missing
    Yf <- factor(woven_example$Y)
    set.seed(1)
    inp <- woven:::.woven_check_input(X, Yf, NULL, 3L, FALSE)
    folds <- woven:::.stratified_folds(inp$anchor_idx, Yf[inp$anchor_idx], 3L)
    grid <- expand.grid(
        gamma_y = c(1, 5), ridge_w = c(1e-8, 1), screen_top = c(Inf, 8),
        lambda = 0.1
    )
    fast <- woven:::.tune_fold(inp$X_list, Yf, folds[[2]], grid, 3L, 10L, TRUE)
    te <- folds[[2]]
    tr <- setdiff(seq_len(90), te)
    slow <- t(vapply(seq_len(nrow(grid)), function(i) {
        g <- grid[i, ]
        f <- woven(lapply(X, function(M) M[tr, , drop = FALSE]), Yf[tr],
            K = 3L, gamma_y = g$gamma_y, ridge_w = g$ridge_w,
            screen_top = g$screen_top, verbose = FALSE, plot = FALSE
        )
        Za <- latent_scores(f)[anchor_indices(f), ]
        m <- woven:::.lda_fit(Za, as.integer(Yf[tr])[anchor_indices(f)])
        woven:::.lda_predict(m, woven_scores(f, lapply(X, function(M) {
            M[te, , drop = FALSE]
        })))
    }, integer(length(te))))
    expect_identical(unname(fast), unname(slow))
})

test_that("woven_tune() returns the best setting, results and a refit", {
    data(woven_example, package = "woven")
    set.seed(2)
    tune <- woven_tune(woven_example$X_missing, woven_example$Y,
        K = 3L, gamma_y = c(1, 5), ridge_w = c(1e-8, 1),
        screen_top = c(Inf, 10, 100), verbose = FALSE, plot = FALSE
    )
    expect_s4_class(tune, "wovenTune")
    # screen_top = 100 exceeds every modality (<= 25 features): same as Inf
    expect_equal(nrow(tune$results), 8L)
    expect_false(is.unsorted(tune$results$ber))
    expect_equal(tune$best$ber, tune$results$ber[1])
    expect_s4_class(tune$fit, "woven")
    expect_equal(tune$fit@ridge_w, tune$best$ridge_w)
    expect_equal(tune$fit@gamma_y, tune$best$gamma_y)
    expect_equal(tune$chance, 2 / 3)
    expect_length(tune$folds, 3L)
    expect_setequal(unlist(tune$folds), anchor_indices(tune$fit))
    expect_output(show(tune), "Held-out BER")
    expect_error(tune$nope, "not part of")
    skip_if_not_installed("ggplot2")
    expect_s3_class(plot(tune), "ggplot")
})

test_that("woven_tune() finds real signal and handles repeats and small classes", {
    set.seed(3)
    n <- 90
    y <- rep(c("a", "b", "c"), each = 30)
    shift <- as.integer(factor(y))
    X <- list(
        A = matrix(rnorm(n * 40), n, 40),
        B = matrix(rnorm(n * 30), n, 30)
    )
    X$A[, 1:3] <- X$A[, 1:3] + 1.5 * shift
    X$B[, 1:3] <- X$B[, 1:3] - 1.5 * shift
    X$A[sample(n, 15), ] <- NA
    tune <- woven_tune(X, y, K = 2L, gamma_y = 1, ridge_w = c(1e-8, 1),
        screen_top = c(Inf, 5), n_repeats = 2L, verbose = FALSE,
        refit = FALSE, plot = FALSE
    )
    expect_lt(tune$best$ber, 0.2)
    expect_null(tune$fit)
    expect_length(tune$folds, 6L)

    y2 <- c(rep("a", 44), rep("b", 44), "c", "c")
    expect_message(
        woven_tune(lapply(X, function(M) {
            M[is.na(M)] <- 0
            M
        }), y2, K = 2L, gamma_y = 1, ridge_w = 1, screen_top = Inf,
            refit = FALSE, plot = FALSE
        ),
        "Using 2 folds"
    )
})
