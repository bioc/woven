# Input validation, precomp checks and reproducibility of woven().

val_fixture <- function(n = 40L, seed = 5L) {
    set.seed(seed)
    X1 <- matrix(rnorm(n * 12), n, 12)
    X2 <- matrix(rnorm(n * 10), n, 10)
    rownames(X1) <- rownames(X2) <- sprintf("P%02d", seq_len(n))
    X1[c(1, 2), ] <- NA # subjects P01, P02 miss modality 1
    X2[c(3), ] <- NA # subject P03 misses modality 2
    list(X = list(RNA = X1, Prot = X2), Y = rep(1:2, length.out = n), n = n)
}

test_that("anchor_idx including block-missing subjects is rejected and named", {
    d <- val_fixture()
    err <- expect_error(
        woven(d$X, Y = d$Y, anchor_idx = 1:20, K = 2L, verbose = FALSE),
        "missing a whole modality block"
    )
    expect_match(conditionMessage(err), "RNA: P01, P02")
    expect_match(conditionMessage(err), "Prot: P03")
})

test_that("anchor_idx out of range, duplicated or non-integer is rejected", {
    d <- val_fixture()
    expect_error(
        woven(d$X, Y = d$Y, anchor_idx = c(4:20, 41L), K = 2L, verbose = FALSE),
        "outside 1..40"
    )
    expect_error(
        woven(d$X, Y = d$Y, anchor_idx = c(4:20, 4L), K = 2L, verbose = FALSE),
        "duplicated subjects: P04"
    )
    expect_error(
        woven(d$X, Y = d$Y, anchor_idx = c(4.5, 5:20), K = 2L, verbose = FALSE),
        "integer row indices"
    )
    expect_error(
        woven(d$X, Y = d$Y, anchor_idx = c(NA, 5:20), K = 2L, verbose = FALSE),
        "integer row indices"
    )
})

test_that("unnamed inputs report subject indices instead of names", {
    d <- val_fixture()
    X <- lapply(d$X, unname)
    names(X) <- NULL
    err <- expect_error(
        woven(X, Y = d$Y, anchor_idx = 1:20, K = 2L, verbose = FALSE)
    )
    expect_match(conditionMessage(err), "Modality 1: 1, 2")
})

test_that("a valid user-supplied anchor_idx equals automatic detection", {
    d <- val_fixture()
    fit_auto <- woven(d$X, Y = d$Y, K = 2L, verbose = FALSE)
    fit_user <- woven(d$X, Y = d$Y, anchor_idx = 4:40, K = 2L, verbose = FALSE)
    expect_identical(anchor_indices(fit_auto), 4:40)
    expect_identical(latent_scores(fit_auto), latent_scores(fit_user))
})

test_that("woven() errors when no subject is observed in every modality", {
    d <- val_fixture()
    d$X$Prot[4:40, ] <- NA
    d$X$RNA[3, ] <- d$X$RNA[4, ] # P03 now only has RNA; nobody has both
    expect_error(
        woven(d$X, Y = d$Y, K = 2L, verbose = FALSE),
        "No fully-observed subjects"
    )
})

test_that("precomp built with a different scale is an error", {
    d <- val_fixture()
    pc_raw <- woven_precompute(d$X, scale = FALSE)
    pc_std <- woven_precompute(d$X, scale = TRUE)
    expect_identical(attr(pc_raw, "woven_scale"), FALSE)
    expect_identical(attr(pc_std, "woven_scale"), TRUE)
    expect_error(
        woven(d$X, Y = d$Y, K = 2L, precomp = pc_raw, scale = TRUE, verbose = FALSE),
        "scale = FALSE but woven\\(\\) was called with scale = TRUE"
    )
    expect_error(
        woven(d$X, Y = d$Y, K = 2L, precomp = pc_std, scale = FALSE, verbose = FALSE),
        "scale = TRUE but woven\\(\\) was called with scale = FALSE"
    )
    # matching settings give the same fit as building the graphs internally
    f_pc <- woven(d$X, Y = d$Y, K = 2L, precomp = pc_raw, scale = FALSE, verbose = FALSE)
    f_in <- woven(d$X, Y = d$Y, K = 2L, scale = FALSE, verbose = FALSE)
    expect_equal(latent_scores(f_pc), latent_scores(f_in), tolerance = 1e-10)
})

test_that("precomp with the wrong number of subjects is an error", {
    d <- val_fixture()
    pc <- woven_precompute(lapply(d$X, function(X) X[1:30, ]))
    expect_error(
        woven(d$X, Y = d$Y, K = 2L, precomp = pc, verbose = FALSE),
        "40 x 40"
    )
})

test_that("a hand-built precomp without a scale record warns", {
    d <- val_fixture()
    pc <- woven_precompute(d$X)
    attr(pc, "woven_scale") <- NULL
    expect_warning(
        woven(d$X, Y = d$Y, K = 2L, precomp = pc, verbose = FALSE),
        "does not record how it was scaled"
    )
})

test_that("k_nn differing from precomp warns and the fit records precomp's k_nn", {
    d <- val_fixture()
    pc <- woven_precompute(d$X, k_nn = 7L)
    expect_warning(
        fit <- woven(d$X, Y = d$Y, K = 2L, precomp = pc, k_nn = 10L, verbose = FALSE),
        "k_nn = 10 is ignored"
    )
    expect_identical(fit@k_nn, 7L)
})

test_that("K larger than a modality's feature count warns", {
    d <- val_fixture()
    d$X$Prot <- d$X$Prot[, 1:2]
    expect_warning(
        woven(d$X, Y = d$Y, K = 3L, verbose = FALSE),
        "exceeds the number of features in Prot \\(p=2\\)"
    )
})

test_that("other argument errors are informative", {
    d <- val_fixture()
    expect_error(woven(d$X$RNA, Y = d$Y), "list of matrices")
    expect_error(woven(d$X[1], Y = d$Y), "at least 2 modalities")
    expect_error(woven(list(d$X$RNA, "a"), Y = d$Y), "not a matrix or data frame")
    expect_error(
        woven(list(unname(d$X$RNA), unname(d$X$Prot)[1:30, ]), Y = d$Y),
        "same number of rows"
    )
    expect_error(woven(d$X, Y = d$Y[-1]), "length\\(Y\\)")
    expect_error(woven(d$X, Y = replace(d$Y, 5, NA)), "NA values")
    expect_error(woven(d$X, Y = d$Y, K = 0L), "K must be >= 1")
    expect_error(woven(d$X, Y = d$Y, K = 40L), "less than n")
})

test_that("Y_labels holds class names; numeric Y gives numbers as strings", {
    d <- val_fixture()
    f_num <- woven(d$X, Y = d$Y, K = 2L, verbose = FALSE)
    expect_identical(f_num@Y_labels, c("1", "2"))
    f_chr <- woven(d$X, Y = c("b", "a")[d$Y], K = 2L, verbose = FALSE)
    expect_identical(f_chr@Y_labels, c("a", "b"))
    pred <- woven_predict(f_chr, lapply(d$X, function(X) X[4:8, ]))
    expect_true(all(pred$predicted_class %in% c("a", "b")))
})

test_that("woven() is deterministic, and seeded metrics are reproducible", {
    data(woven_example, package = "woven")
    f1 <- woven(woven_example$X_missing, Y = woven_example$Y, K = 3L, verbose = FALSE)
    f2 <- woven(woven_example$X_missing, Y = woven_example$Y, K = 3L, verbose = FALSE)
    expect_identical(f1, f2)

    set.seed(11)
    m1 <- woven_metrics(f1, woven_example$Y)
    set.seed(11)
    m2 <- woven_metrics(f1, woven_example$Y)
    expect_identical(m1, m2)

    set.seed(3)
    e1 <- woven_nystrom_error(f1, woven_example$X_missing, n_loo = 4L)
    set.seed(3)
    e2 <- woven_nystrom_error(f1, woven_example$X_missing, n_loo = 4L)
    expect_identical(e1, e2)
})

test_that("verbose mode reports anchors and scoring", {
    d <- val_fixture()
    expect_message(
        woven(d$X, Y = d$Y, K = 2L, plot = FALSE),
        "Anchors \\(observed in every modality\\): 37 of 40"
    )
})
