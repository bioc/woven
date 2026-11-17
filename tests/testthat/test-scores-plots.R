# woven_scores(), woven_predict() and the woven_plot_*() functions.

sp_fit <- function() {
    data(woven_example, package = "woven")
    woven(woven_example$X_missing, Y = woven_example$Y, K = 3L, verbose = FALSE)
}

test_that("woven_scores() returns an n_new x K numeric matrix with names", {
    fit <- sp_fit()
    Xn <- lapply(woven_example$X_missing, function(X) X[1:7, ])
    Z <- woven_scores(fit, Xn)
    expect_true(is.matrix(Z) && is.numeric(Z))
    expect_equal(dim(Z), c(7L, 3L))
    expect_identical(colnames(Z), paste0("Dim", 1:3))
    expect_identical(rownames(Z), rownames(woven_example$X_missing$RNA)[1:7])
})

test_that("woven_scores() on the training data reproduces latent_scores()", {
    fit <- sp_fit()
    expect_equal(woven_scores(fit, woven_example$X_missing), latent_scores(fit),
        tolerance = 1e-12
    )
})

test_that("woven_scores() gives NA rows to subjects with no observed modality", {
    fit <- sp_fit()
    Xn <- lapply(woven_example$X_complete, function(X) X[1:4, ])
    Xn <- lapply(Xn, function(X) {
        X[2, ] <- NA
        X
    })
    Z <- woven_scores(fit, Xn)
    expect_true(all(is.na(Z[2, ])))
    expect_false(anyNA(Z[-2, ]))
})

test_that("a subject's score from one modality is that modality's projection", {
    fit <- sp_fit()
    Xn <- lapply(woven_example$X_complete, function(X) X[1, , drop = FALSE])
    Xn$Methylation[] <- NA
    Xn$Proteomics[] <- NA
    xs <- (Xn$RNA - fit@scale_center[[1]]) / fit@scale_scale[[1]]
    expect_equal(unname(woven_scores(fit, Xn)),
        unname(xs %*% projection_matrices(fit)$RNA),
        tolerance = 1e-12
    )
})

test_that("woven_scores() rejects mismatched inputs", {
    fit <- sp_fit()
    Xn <- lapply(woven_example$X_missing, function(X) X[1:3, ])
    expect_error(woven_scores(fit, unname(Xn[1:2])), "list of 3 matrices")
    expect_error(woven_scores(fit, list(Lipids = Xn$RNA)), "Unknown modality")
    Xn$Proteomics <- Xn$Proteomics[, 1:5]
    expect_error(woven_scores(fit, Xn), "Proteomics\\) has 5 features")
})

test_that("woven_predict() returns probabilities for both methods", {
    fit <- sp_fit()
    Xn <- lapply(woven_example$X_missing, function(X) X[1:10, ])
    for (m in c("centroid", "knn")) {
        p <- woven_predict(fit, Xn, method = m, k_pred = 3L)
        expect_s3_class(p, "data.frame")
        expect_equal(nrow(p), 10L)
        expect_named(p, c("predicted_class", "confidence", "p_AD", "p_CN", "p_MCI"))
        expect_equal(rowSums(p[, 3:5]), rep(1, 10),
            tolerance = 1e-12,
            ignore_attr = TRUE
        )
        expect_true(all(p$predicted_class %in% c("AD", "CN", "MCI")))
    }
    expect_error(woven_predict(fit, Xn, method = "svm"), "should be one of")
})

test_that("woven_predict() classifies the training subjects well above chance", {
    fit <- sp_fit()
    p <- woven_predict(fit, woven_example$X_missing)
    expect_gt(mean(p$predicted_class == woven_example$Y), 0.6) # chance = 1/3
})

test_that("woven_plot_*() return ggplot objects", {
    skip_if_not_installed("ggplot2")
    fit <- sp_fit()
    expect_s3_class(woven_plot_vip(fit), "ggplot")
    expect_s3_class(woven_plot_vip(fit, modality = "Proteomics", n_top = 5L), "ggplot")
    expect_s3_class(woven_plot_loadings(fit, dim = 2L), "ggplot")
    expect_s3_class(woven_plot_loadings(fit, dim = 1L, modality = c("RNA", "Proteomics")), "ggplot")
    expect_s3_class(woven_plot_variance(fit), "ggplot")
    expect_s3_class(plot(fit, labels = woven_example$Y, dims = c(2L, 3L)), "ggplot")
})

test_that("woven_plot_vip() data rank features by VIP", {
    skip_if_not_installed("ggplot2")
    fit <- sp_fit()
    p <- woven_plot_vip(fit, modality = 1L, n_top = 4L)
    expect_equal(nrow(p$data), 4L)
    expect_false(is.unsorted(rev(p$data$vip_score)))
})

test_that("woven_plot_*() error on a bad dim or modality", {
    skip_if_not_installed("ggplot2")
    fit <- sp_fit()
    expect_error(woven_plot_vip(fit, modality = 4L), "integer in 1..3")
    expect_error(woven_plot_vip(fit, modality = "Lipids"), "Modality 'Lipids' not found")
    expect_error(woven_plot_vip(fit, modality = NA), "single modality")
    expect_error(woven_plot_vip(fit, feature_names = letters[1:3]), "single modality")
    expect_error(
        woven_plot_vip(fit, modality = 1L, feature_names = letters[1:3]),
        "feature_names length"
    )
    expect_error(woven_plot_loadings(fit, dim = 4L), "dim must be a single integer in 1..3")
    expect_error(woven_plot_loadings(fit, dim = 0L), "dim must be")
    expect_error(woven_plot_loadings(fit, modality = "Lipids"), "not found")
    expect_error(
        woven_plot_loadings(fit, feature_names = list("a", "b", "c")),
        "feature_names\\[\\[1\\]\\] length"
    )
    expect_error(plot(fit, dims = c(1L, 9L)), "dims must be two integers")
})

test_that("woven_plot_loadings() accepts a character vector of feature names", {
    skip_if_not_installed("ggplot2")
    fit <- sp_fit()
    nm <- paste0("g", 1:25)
    p <- woven_plot_loadings(fit, modality = "RNA", feature_names = nm, n_top = 5L)
    expect_true(all(as.character(p$data$feature) %in% nm))
})

test_that("show() lists modalities and classes", {
    fit <- sp_fit()
    expect_output(show(fit), "Classes    : AD, CN, MCI")
    expect_output(show(fit), "RNA, Methylation, Proteomics")
})
