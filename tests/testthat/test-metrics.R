test_that("woven_silhouette returns scalar in [-1, 1]", {
    set.seed(1)
    Z <- matrix(rnorm(100 * 5), 100, 5)
    labels <- rep(1:4, 25)
    s <- woven_silhouette(Z, labels)
    expect_true(is.finite(s))
    expect_true(s >= -1 && s <= 1)
})

test_that("woven_silhouette returns NA with one class", {
    Z <- matrix(rnorm(40 * 3), 40, 3)
    expect_true(is.na(woven_silhouette(Z, rep(1, 40))))
})

test_that("woven_davies_bouldin is non-negative", {
    set.seed(2)
    Z <- matrix(rnorm(80 * 4), 80, 4)
    labels <- rep(1:4, 20)
    db <- woven_davies_bouldin(Z, labels)
    expect_true(is.finite(db) && db >= 0)
})

test_that("woven_nmi is in [0, 1]", {
    set.seed(3)
    Z <- matrix(rnorm(60 * 3), 60, 3)
    labels <- rep(1:3, 20)
    nmi <- woven_nmi(Z, labels)
    expect_true(is.finite(nmi) && nmi >= 0 && nmi <= 1)
})

test_that("woven_nmi is 1 for perfectly separated clusters", {
    # Clusters perfectly separated — k-means should recover them exactly
    set.seed(4)
    Z <- rbind(
        matrix(rnorm(30 * 2, mean = 0), 30, 2),
        matrix(rnorm(30 * 2, mean = 10), 30, 2),
        matrix(rnorm(30 * 2, mean = 20), 30, 2)
    )
    labels <- rep(1:3, each = 30)
    nmi <- woven_nmi(Z, labels, n_cl = 3L)
    expect_gt(nmi, 0.95)
})

test_that("woven_rv returns value in [0, 1]", {
    set.seed(5)
    Z <- matrix(rnorm(50 * 4), 50, 4)
    Z_true <- matrix(rnorm(50 * 5), 50, 5)
    rv <- woven_rv(Z, Z_true)
    expect_true(is.finite(rv) && rv >= 0 && rv <= 1)
})

test_that("woven_rv is 1 when matrices are identical", {
    set.seed(6)
    Z <- matrix(rnorm(40 * 3), 40, 3)
    rv <- woven_rv(Z, Z)
    expect_gt(rv, 0.999)
})

test_that("woven_rv is higher for correlated than uncorrelated matrices", {
    set.seed(61)
    Z <- matrix(rnorm(50 * 4), 50, 4)
    Z_corr <- Z + matrix(rnorm(50 * 4, sd = 0.1), 50, 4) # nearly identical
    Z_rand <- matrix(rnorm(50 * 5), 50, 5) # random, uncorrelated
    expect_gt(woven_rv(Z, Z_corr), woven_rv(Z, Z_rand))
})

test_that("woven_ess_retention is n_used/n_total", {
    expect_equal(woven_ess_retention(145L, 300L), 145 / 300)
    expect_equal(woven_ess_retention(300L, 300L), 1.0)
})

test_that("woven_metrics returns named numeric vector from fit", {
    set.seed(7)
    n <- 60
    K <- 3L
    groups <- rep(c("A", "B", "C"), 20)
    X1 <- matrix(rnorm(n * 30), n, 30)
    X2 <- matrix(rnorm(n * 20), n, 20)

    fit <- woven(list(X1, X2), Y = groups, K = K, verbose = FALSE)
    m <- woven_metrics(fit, groups)

    expect_named(m, c("Silhouette", "Davies-Bouldin", "NMI", "ESS"))
    expect_true(all(is.finite(m)))
    expect_equal(m["ESS"], c(ESS = 1.0))
})

test_that("woven_davies_bouldin is lower for better-separated clusters", {
    set.seed(8)
    lab <- rep(1:3, each = 20)
    Z_far <- matrix(rnorm(120), 60, 2) + cbind(lab * 10, 0)
    Z_near <- matrix(rnorm(120), 60, 2) + cbind(lab * 1, 0)
    expect_lt(woven_davies_bouldin(Z_far, lab), woven_davies_bouldin(Z_near, lab))
    expect_true(is.na(woven_davies_bouldin(Z_far, rep(1, 60))))
    expect_true(is.na(woven_nmi(Z_far, rep(1, 60))))
})

test_that("woven_effect_bias() recovers a simulated treatment effect", {
    set.seed(9)
    n <- 2000
    Z <- matrix(rnorm(n * 2), n, 2)
    labels <- rep(c("g1", "g2"), each = n / 2)
    treatment <- rep(0:1, n / 2)
    true_eff <- c(1, 3)
    outcome <- true_eff[as.integer(factor(labels))] * treatment +
        0.5 * Z[, 1] + rnorm(n, sd = 0.3)
    expect_lt(woven_effect_bias(Z, outcome, treatment, labels, true_eff), 0.05)
    # a wrong "true" effect shows up as a large bias
    expect_gt(woven_effect_bias(Z, outcome, treatment, labels, c(2, 6)), 0.4)
    # subgroups that are too small are skipped
    expect_true(is.nan(woven_effect_bias(
        Z[1:4, ], outcome[1:4],
        treatment[1:4], labels[1:4], true_eff
    )))
})

test_that("woven_ess_retention() rejects impossible counts", {
    expect_error(woven_ess_retention(5, 4))
    expect_error(woven_ess_retention(1, 0))
})

test_that("woven_all_metrics() matches direct calls, including the Nystrom error", {
    data(woven_example, package = "woven")
    fit <- woven(woven_example$X_missing, Y = woven_example$Y, K = 3L, verbose = FALSE)
    Z <- latent_scores(fit)
    set.seed(4)
    truth <- matrix(rnorm(90 * 2), 90, 2)
    trt <- rep(0:1, 45)
    out_y <- rnorm(90)

    set.seed(21)
    all_m <- woven:::woven_all_metrics(Z, woven_example$Y,
        n_total = 90L,
        Z_true = truth, outcome = out_y, treatment = trt,
        true_effects = c(1, 1, 1), fit = fit,
        X_list = woven_example$X_missing, n_loo = 6L
    )
    expect_named(all_m, c(
        "silhouette", "davies_bouldin", "nmi", "ess_retention",
        "rv_coefficient", "effect_bias", "nystrom_error"
    ))

    # same RNG stream as inside woven_all_metrics(): nmi draws first, then
    # woven_nystrom_error() samples its held-out anchors
    set.seed(21)
    nmi <- woven_nmi(Z, woven_example$Y)
    ny <- woven_nystrom_error(fit, woven_example$X_missing, n_loo = 6L)
    expect_identical(all_m$nmi, nmi)
    expect_identical(all_m$nystrom_error, ny)
    expect_identical(all_m$silhouette, woven_silhouette(Z, woven_example$Y))
    expect_identical(all_m$rv_coefficient, woven_rv(Z, truth))
    expect_identical(
        all_m$effect_bias,
        woven_effect_bias(Z, out_y, trt, woven_example$Y, c(1, 1, 1))
    )
    expect_true(is.finite(all_m$nystrom_error))
})

test_that("woven_nystrom_error() validates n_loo and X_list", {
    data(woven_example, package = "woven")
    fit <- woven(woven_example$X_complete, Y = woven_example$Y, K = 3L, verbose = FALSE)
    expect_error(woven_nystrom_error(fit, woven_example$X_complete, n_loo = 0L), "n_loo")
    expect_error(woven_nystrom_error(fit, woven_example$X_complete, n_loo = 91L), "n_loo")
    expect_error(woven_nystrom_error(fit, woven_example$X_complete, n_loo = 2.5), "n_loo")
    expect_error(woven_nystrom_error(fit, woven_example$X_complete[1:2]), "X_list")
    set.seed(1)
    expect_true(is.finite(woven_nystrom_error(fit, woven_example$X_complete, n_loo = 90L)))
})

test_that("woven_nystrom_error() is small for a stable fit", {
    # strong shared signal: dropping one anchor barely moves its score
    data(woven_example, package = "woven")
    fit <- woven(woven_example$X_complete, Y = woven_example$Y, K = 2L, verbose = FALSE)
    set.seed(2)
    err <- woven_nystrom_error(fit, woven_example$X_complete, n_loo = 10L)
    spread <- sqrt(sum(apply(latent_scores(fit), 2, stats::var)))
    expect_lt(err, spread)
})
