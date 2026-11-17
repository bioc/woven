# input.R -- turn what users actually have into woven's internal format
#
# Internally woven() works on a list of V numeric matrices with the same rows
# (subjects) in the same order, where a subject missing a modality has an
# all-NA row. Users rarely have data in that shape: each assay usually comes
# as its own table covering its own subset of subjects. The helpers here
# accept that form directly (matrices or data frames, matched by row names)
# and build the aligned list, so no manual padding is needed.

# Convert one modality to a numeric matrix. Data frames are accepted when
# every column is numeric; anything else stops with a message naming the
# modality and the offending columns.
.as_block_matrix <- function(X, label) {
    if (is.data.frame(X)) {
        num <- vapply(X, is.numeric, logical(1L))
        if (!all(num)) {
            fmt <- paste0(
                "Modality %s has non-numeric column(s): %s. Keep only ",
                "numeric features (move IDs to row names, labels to Y)."
            )
            stop(sprintf(
                fmt, label,
                paste(utils::head(names(X)[!num], 5L), collapse = ", ")
            ), call. = FALSE)
        }
        X <- as.matrix(X)
    }
    if (!is.matrix(X)) {
        fmt <- paste0(
            "Modality %s is not a matrix or data frame (got %s). ",
            "Convert with as.matrix()."
        )
        stop(sprintf(fmt, label, class(X)[1L]), call. = FALSE)
    }
    if (!is.numeric(X)) {
        stop(sprintf("Modality %s must be numeric.", label), call. = FALSE)
    }
    storage.mode(X) <- "double"
    X
}

# Align modalities on subject IDs. If every modality has row names and they
# are not already identical, the result covers the union of subjects (in
# order of first appearance) and a subject absent from a modality gets an
# all-NA row there. Matrices without row names must already share rows.
# Returns the aligned list and whether alignment was needed.
.align_blocks <- function(X_list) {
    rn <- lapply(X_list, rownames)
    has_rn <- !vapply(rn, is.null, logical(1L))
    same <- all(has_rn) &&
        all(vapply(rn, identical, logical(1L), rn[[1L]]))
    if (!all(has_rn) || same) {
        return(list(X_list = X_list, aligned = FALSE))
    }
    for (v in seq_along(X_list)) {
        dup <- unique(rn[[v]][duplicated(rn[[v]])])
        if (length(dup) > 0L) {
            fmt <- paste0(
                "Modality %s has duplicated row names (subject IDs): %s. ",
                "Each subject needs one row per modality."
            )
            lab <- if (is.null(names(X_list)) || !nzchar(names(X_list)[v])) {
                sprintf("[[%d]]", v)
            } else {
                sprintf("'%s'", names(X_list)[v])
            }
            stop(sprintf(
                fmt, lab,
                paste(utils::head(dup, 5L), collapse = ", ")
            ), call. = FALSE)
        }
    }
    ids <- unique(unlist(rn, use.names = FALSE))
    out <- lapply(X_list, function(X) {
        M <- matrix(NA_real_,
            nrow = length(ids), ncol = ncol(X),
            dimnames = list(ids, colnames(X))
        )
        M[rownames(X), ] <- X
        M
    })
    names(out) <- names(X_list)
    list(X_list = out, aligned = TRUE)
}

# Give every modality a name and every feature a name, so plots and tables
# are readable and features from different modalities never collide.
.name_blocks <- function(X_list) {
    V <- length(X_list)
    nm <- names(X_list)
    if (is.null(nm)) nm <- rep("", V)
    nm[is.na(nm) | nm == ""] <- paste0("Modality ", seq_len(V))[
        is.na(nm) | nm == ""
    ]
    names(X_list) <- make.unique(nm)
    for (v in seq_len(V)) {
        cn <- colnames(X_list[[v]])
        if (is.null(cn) || all(cn == "")) {
            colnames(X_list[[v]]) <- paste0(
                names(X_list)[v], "_", seq_len(ncol(X_list[[v]]))
            )
        }
    }
    X_list
}

# Match labels to subjects. A named Y whose names cover the subject IDs is
# reordered by name (so its order need not match the data); otherwise Y is
# taken in row order and must have one entry per subject.
.align_labels <- function(Y, ids, n) {
    if (!is.null(names(Y)) && !is.null(ids) && all(ids %in% names(Y))) {
        return(Y[ids])
    }
    if (length(Y) != n) {
        fmt <- paste0(
            "length(Y) = %d but the data have %d subjects. Give one label ",
            "per subject, in row order, or a named vector whose names are ",
            "the subject IDs (row names)."
        )
        stop(sprintf(fmt, length(Y), n), call. = FALSE)
    }
    Y
}
