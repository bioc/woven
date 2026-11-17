# Input alignment, progress output, stored labels, top features and the
# plot(type = ) entry point.

ux_data <- function() {
    data(woven_example, package = "woven")
    woven_example
}

test_that("separate tables matched by row name give the same fit as padded input", {
    ex <- ux_data()
    ids <- rownames(ex$X_missing$RNA)
    fit_pad <- woven(ex$X_missing, Y = ex$Y, K = 3L, verbose = FALSE)

    set.seed(3)
    tabs <- lapply(ex$X_missing, function(M) {
        M <- M[!is.na(M[, 1]), , drop = FALSE]
        M[sample(nrow(M)), , drop = FALSE]
    })
    tabs$RNA <- as.data.frame(tabs$RNA)
    y_named <- setNames(ex$Y, ids)[sample(length(ids))]
    fit_tab <- woven(tabs, Y = y_named, K = 3L, verbose = FALSE)

    expect_setequal(rownames(latent_scores(fit_tab)), ids)
    # Reordering rows can flip the sign of a latent dimension; match signs.
    Za <- latent_scores(fit_tab)[ids, ]
    Zb <- latent_scores(fit_pad)[ids, ]
    sgn <- sign(colSums(Za * Zb))
    expect_equal(sweep(Za, 2L, sgn, "*"), Zb, tolerance = 1e-6)
    expect_identical(unname(as.character(fit_tab@Y[ids])), ex$Y)
})

test_that("tables covering different subjects are padded with missing blocks", {
    set.seed(4)
    a <- matrix(rnorm(30 * 6), 30, 6, dimnames = list(paste0("S", 1:30), NULL))
    b <- matrix(rnorm(25 * 5), 25, 5, dimnames = list(paste0("S", 11:35), NULL))
    y <- setNames(rep(c("x", "y"), length.out = 35), paste0("S", 1:35))
    fit <- woven(list(A = a, B = b), Y = y, K = 2L, verbose = FALSE)
    expect_equal(fit@n, 35L)
    expect_equal(length(anchor_indices(fit)), 20L)
    expect_equal(sum(fit@observed[, "A"]), 30L)
    expect_equal(sum(fit@observed[, "B"]), 25L)
    expect_false(anyNA(latent_scores(fit)))
    expect_identical(rownames(projection_matrices(fit)$A), paste0("A_", 1:6))
})

test_that("bad tabular input gives clear errors", {
    set.seed(5)
    a <- data.frame(id = letters[1:20], x = rnorm(20))
    b <- matrix(rnorm(40), 20, 2)
    expect_error(woven(list(a, b), Y = rep(1:2, 10)), "non-numeric column")
    d1 <- matrix(rnorm(40), 20, 2, dimnames = list(rep(c("p", "q"), 10), NULL))
    d2 <- matrix(rnorm(40), 20, 2, dimnames = list(paste0("s", 1:20), NULL))
    expect_error(woven(list(d1, d2), Y = rep(1:2, 10)), "duplicated row names")
    expect_error(woven(data.frame(x = 1:3), Y = 1:3), "list of matrices")
})

test_that("woven() reports each step and a closing summary", {
    ex <- ux_data()
    msgs <- character(0)
    withCallingHandlers(
        woven(ex$X_missing, Y = ex$Y, K = 3L, plot = FALSE),
        message = function(m) {
            msgs <<- c(msgs, conditionMessage(m))
            invokeRestart("muffleMessage")
        }
    )
    all_msgs <- paste(msgs, collapse = "")
    expect_match(all_msgs, "WOVEN \\| 90 subjects \\| 3 modalities")
    expect_match(all_msgs, "Anchors \\(observed in every modality\\): 60 of 90")
    expect_match(all_msgs, "\\[1/3\\] Similarity graph 3/3: Proteomics")
    expect_match(all_msgs, "\\[2/3\\]")
    expect_match(all_msgs, "complete-case method would keep 60")
    expect_match(all_msgs, "Top features:")
    expect_silent(woven(ex$X_missing, Y = ex$Y, K = 3L, verbose = FALSE, plot = FALSE))
})

test_that("woven_top_features() ranks by VIP with the highest class", {
    ex <- ux_data()
    fit <- woven(ex$X_missing, Y = ex$Y, K = 3L, verbose = FALSE)
    tf <- woven_top_features(fit, n = 4)
    expect_named(tf, c("modality", "rank", "feature", "vip", "highest_in"))
    expect_equal(nrow(tf), 12L)
    for (m in fit@mod_names) {
        expect_false(is.unsorted(rev(tf$vip[tf$modality == m])))
    }
    expect_true(all(tf$highest_in %in% c("AD", "CN", "MCI")))
    one <- woven_top_features(fit, n = Inf, modality = "Proteomics")
    expect_equal(nrow(one), 15L)
    vip <- woven:::.vip_scores(fit, 3L)
    expect_equal(mean(vip^2), 1, tolerance = 1e-10)
    cm <- fit@class_means$RNA
    top <- tf$feature[tf$modality == "RNA"][1]
    expect_identical(tf$highest_in[1], colnames(cm)[which.max(cm[top, ])])
    expect_error(woven_top_features(fit, n = 0), "n must be")
})

test_that("stored labels are the defaults for metrics and plots", {
    ex <- ux_data()
    fit <- woven(ex$X_missing, Y = ex$Y, K = 3L, verbose = FALSE)
    set.seed(1)
    m1 <- woven_metrics(fit)
    set.seed(1)
    m2 <- woven_metrics(fit, ex$Y)
    expect_equal(m1, m2)
    skip_if_not_installed("ggplot2")
    p <- plot(fit)
    expect_setequal(levels(p$data$group), c("AD", "CN", "MCI"))
    expect_error(plot(fit, labels = 1:3), "one entry per subject")
})

test_that("plot(type = ) dispatches to every plot", {
    skip_if_not_installed("ggplot2")
    ex <- ux_data()
    fit <- woven(ex$X_missing, Y = ex$Y, K = 3L, verbose = FALSE)
    for (ty in c("scores", "features", "loadings", "coverage", "variance")) {
        expect_s3_class(plot(fit, type = ty), "ggplot")
    }
    pf <- plot(fit, type = "features", n_top = 5, modality = c("RNA", "Proteomics"))
    expect_equal(nrow(pf$data), 10L)
    expect_s3_class(woven_plot_vip(fit, color_by = "modality"), "ggplot")
    pc <- woven_plot_coverage(fit)
    expect_equal(sum(pc$data$n), 90L)
    expect_match(pc$labels$subtitle, "Complete-case keeps 60 of 90")
    expect_error(plot(fit, type = "heatmap"), "should be one of")
})

test_that("new data can be separate tables and may omit modalities", {
    ex <- ux_data()
    fit <- woven(ex$X_missing, Y = ex$Y, K = 3L, verbose = FALSE)
    Xc <- lapply(ex$X_complete, function(X) X[1:6, ])
    full <- woven_scores(fit, Xc)
    # RNA only, by name: same as padding the other modalities with NA rows
    pad <- Xc
    pad$Methylation[] <- NA
    pad$Proteomics[] <- NA
    expect_equal(woven_scores(fit, list(RNA = Xc$RNA)), woven_scores(fit, pad))
    # blocks in a different order, as data frames, rows shuffled
    shuffled <- list(
        Proteomics = as.data.frame(Xc$Proteomics[6:1, ]),
        RNA = Xc$RNA, Methylation = Xc$Methylation
    )
    expect_equal(woven_scores(fit, shuffled)[rownames(full), ], full)
    p <- woven_predict(fit, list(RNA = Xc$RNA))
    expect_equal(nrow(p), 6L)
})

test_that("plots adapt to K = 1, many classes, continuous labels", {
    skip_if_not_installed("ggplot2")
    set.seed(8)
    n <- 120
    y <- rep(paste0("G", 1:10), length.out = n)
    X <- list(
        A = matrix(rnorm(n * 12), n, 12) + as.integer(factor(y)),
        B = matrix(rnorm(n * 9), n, 9)
    )
    X$B[1:20, ] <- NA
    fit1 <- woven(X, Y = y, K = 1L, verbose = FALSE)
    expect_s3_class(plot(fit1), "ggplot")
    expect_s3_class(plot(fit1, type = "variance"), "ggplot")
    fit3 <- woven(X, Y = y, K = 3L, verbose = FALSE)
    p <- plot(fit3)
    expect_s3_class(p, "ggplot")
    expect_length(unique(ggplot2::ggplot_build(p)$data[[1]]$colour), 10L)
    expect_s3_class(plot(fit3, labels = runif(n)), "ggplot")
    expect_s3_class(plot(fit3, highlight_anchors = FALSE), "ggplot")
})

test_that("coverage plot lumps rare patterns", {
    skip_if_not_installed("ggplot2")
    set.seed(9)
    n <- 150
    X <- lapply(1:5, function(v) matrix(rnorm(n * 6), n, 6))
    names(X) <- paste0("M", 1:5)
    m <- matrix(runif(n * 5) < 0.3, n)
    for (i in which(rowSums(m) == 5)) m[i, 1] <- FALSE
    for (v in 1:5) X[[v]][m[, v], ] <- NA
    fit <- woven(X, Y = rep(1:3, 50), K = 2L, verbose = FALSE)
    pc <- woven_plot_coverage(fit, max_patterns = 6L)
    expect_lte(nlevels(pc$data$pattern), 6L)
    expect_equal(sum(pc$data$n), n)
    expect_true(any(grepl("^Other patterns", levels(pc$data$pattern))))
})

test_that("features screened out are not plotted", {
    skip_if_not_installed("ggplot2")
    data(woven_example, package = "woven")
    fit <- woven(woven_example$X_missing, Y = woven_example$Y, K = 3L,
        screen_top = 5, verbose = FALSE
    )
    p <- woven_plot_vip(fit, n_top = 20)
    expect_equal(nrow(p$data), 15L) # 5 kept per modality
    expect_true(all(p$data$vip_score > 0))
})
