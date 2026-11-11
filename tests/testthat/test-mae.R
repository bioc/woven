test_that("woven() accepts a MultiAssayExperiment and matches the equivalent list fit", {
    skip_if_not_installed("MultiAssayExperiment")
    suppressWarnings(suppressPackageStartupMessages(library(MultiAssayExperiment)))

    set.seed(1)
    n <- 60
    p1 <- 20
    p2 <- 15
    K <- 2L
    Y <- rep(1:2, each = n / 2)
    X1 <- matrix(rnorm(n * p1), n, p1)
    X2 <- matrix(rnorm(n * p2), n, p2)
    miss <- matrix(runif(n * 2) < 0.3, n, 2)
    for (i in which(rowSums(miss) == 2)) miss[i, sample(2, 1)] <- FALSE
    X1[miss[, 1], ] <- NA
    X2[miss[, 2], ] <- NA

    fit_list <- woven(list(X1, X2), Y = Y, K = K, verbose = FALSE)

    ids <- paste0("S", seq_len(n))
    rownames(X1) <- ids
    rownames(X2) <- ids
    present1 <- ids[!is.na(X1[, 1])]
    present2 <- ids[!is.na(X2[, 1])]
    exps <- ExperimentList(list(
        Mod1 = t(X1[present1, , drop = FALSE]),
        Mod2 = t(X2[present2, , drop = FALSE])
    ))
    cd <- DataFrame(row.names = ids, group = Y)
    mae <- MultiAssayExperiment(exps, colData = cd)

    # Y as a colData column name (character length 1)
    fit_mae <- woven(mae, Y = "group", K = K, verbose = FALSE)
    expect_equal(dim(fit_mae@Z), dim(fit_list@Z))
    expect_equal(mean(!is.na(fit_mae@Z[, 1])), mean(!is.na(fit_list@Z[, 1])))
    for (k in seq_len(K)) {
        d_same <- max(abs(fit_list@Z[, k] - fit_mae@Z[, k]), na.rm = TRUE)
        d_flip <- max(abs(fit_list@Z[, k] + fit_mae@Z[, k]), na.rm = TRUE)
        expect_lt(min(d_same, d_flip), 1e-8)
    }

    # Y as an explicit vector (same subject order as colData(mae))
    fit_mae2 <- woven(mae, Y = Y, K = K, verbose = FALSE)
    expect_equal(dim(fit_mae2@Z), dim(fit_list@Z))
})

test_that(".mae_to_xlist() fills all-NA rows for subjects absent from an experiment", {
    skip_if_not_installed("MultiAssayExperiment")
    suppressWarnings(suppressPackageStartupMessages(library(MultiAssayExperiment)))

    exps <- ExperimentList(list(
        A = matrix(1:6, nrow = 3, dimnames = list(c("f1", "f2", "f3"), c("s1", "s2"))),
        B = matrix(7:15, nrow = 3, dimnames = list(c("f1", "f2", "f3"), c("s2", "s3", "s4")))
    ))
    cd <- DataFrame(row.names = c("s1", "s2", "s3", "s4"), group = c(1, 1, 2, 2))
    mae <- MultiAssayExperiment(exps, colData = cd)

    xl <- .mae_to_xlist(mae)
    expect_named(xl, c("A", "B"))
    expect_equal(dim(xl$A), c(4L, 3L))
    expect_equal(dim(xl$B), c(4L, 3L))
    expect_true(all(is.na(xl$A["s3", ])))
    expect_true(all(is.na(xl$A["s4", ])))
    expect_true(all(is.na(xl$B["s1", ])))
    expect_false(any(is.na(xl$A["s1", ])))
    expect_false(any(is.na(xl$B["s4", ])))
})

test_that(".mae_to_xlist() stops on replicate samples and names them", {
    skip_if_not_installed("MultiAssayExperiment")
    suppressWarnings(suppressPackageStartupMessages(library(MultiAssayExperiment)))
    A <- matrix(c(10, 11, 20, 21, 30, 31),
        nrow = 2,
        dimnames = list(c("f1", "f2"), c("s1", "s2", "s3"))
    )
    B <- matrix(1:6, nrow = 2, dimnames = list(c("g1", "g2"), c("t1", "t2", "t3")))
    smap <- DataFrame(
        assay = factor(c("A", "A", "A", "B", "B", "B")),
        primary = c("P1", "P1", "P2", "P1", "P2", "P3"),
        colname = c("s1", "s2", "s3", "t1", "t2", "t3")
    )
    cd <- DataFrame(row.names = c("P1", "P2", "P3"), group = c(1, 2, 1))
    mae <- MultiAssayExperiment(ExperimentList(list(A = A, B = B)),
        colData = cd, sampleMap = smap
    )
    expect_error(woven:::.mae_to_xlist(mae), "'A' has more than one sample for subject\\(s\\) P1")
    merged <- mergeReplicates(mae)
    xl <- woven:::.mae_to_xlist(merged)
    expect_equal(unname(xl$A["P1", ]), c(15, 16))
    expect_true(all(is.na(xl$A["P3", ])))
})

test_that("woven() row names follow colData for a MultiAssayExperiment", {
    skip_if_not_installed("MultiAssayExperiment")
    suppressWarnings(suppressPackageStartupMessages(library(MultiAssayExperiment)))
    data(woven_example, package = "woven")
    ids <- sprintf("S%03d", 1:90)
    exps <- lapply(woven_example$X_missing, function(X) {
        keep <- rowSums(is.na(X)) < ncol(X)
        t(X[keep, , drop = FALSE])
    })
    mae <- MultiAssayExperiment(ExperimentList(exps),
        colData = DataFrame(row.names = ids, dx = woven_example$Y)
    )
    fit <- woven(mae, Y = "dx", K = 3L, verbose = FALSE)
    expect_identical(rownames(latent_scores(fit)), rownames(colData(mae)))
    fit_list <- woven(woven_example$X_missing, Y = woven_example$Y, K = 3L, verbose = FALSE)
    expect_equal(unname(latent_scores(fit)), unname(latent_scores(fit_list)),
        tolerance = 1e-10
    )
})

test_that("woven_scores()/woven_predict() accept a MultiAssayExperiment", {
    skip_if_not_installed("MultiAssayExperiment")
    suppressWarnings(suppressPackageStartupMessages(library(MultiAssayExperiment)))
    data(woven_example, package = "woven")
    ids <- sprintf("S%03d", 1:90)
    exps <- lapply(woven_example$X_missing, function(X) {
        keep <- rowSums(is.na(X)) < ncol(X)
        t(X[keep, , drop = FALSE])
    })
    mae <- MultiAssayExperiment(ExperimentList(exps),
        colData = DataFrame(row.names = ids, dx = woven_example$Y)
    )
    fit <- woven(mae[, 1:60], Y = "dx", K = 3L, verbose = FALSE)
    # experiments given in a different order and features shuffled
    test <- mae[, 61:90, rev(names(mae))]
    test[["RNA"]] <- test[["RNA"]][rev(rownames(test[["RNA"]])), ]
    Z_mae <- woven_scores(fit, test)
    X_test <- lapply(woven_example$X_missing, function(X) X[61:90, ])
    expect_equal(unname(Z_mae), unname(woven_scores(fit, X_test)), tolerance = 1e-12)
    expect_identical(rownames(Z_mae), ids[61:90])
    p <- woven_predict(fit, test)
    expect_identical(rownames(p), ids[61:90])
    only_rna <- suppressWarnings(suppressMessages(mae[, , "RNA"]))
    expect_error(woven_scores(fit, only_rna), "lacks experiment")
})
