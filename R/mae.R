# mae.R -- MultiAssayExperiment input support for woven()
#
# MultiAssayExperiment (Ramos et al. 2017) is the Bioconductor container for
# the setting woven() targets: several assays measured on overlapping but
# non-identical sets of patients. This file converts a MultiAssayExperiment
# into the named list of subject x feature matrices (with all-NA rows for
# subjects missing an assay) that the rest of the package uses, so no other
# code needs to know about MultiAssayExperiment.

#' Convert a MultiAssayExperiment to woven's internal list format
#'
#' Aligns every experiment in \code{mae} to the patient roster
#' \code{rownames(colData(mae))} using \code{sampleMap(mae)}, transposes each
#' assay to subject x feature, and fills an all-NA row for any subject that
#' has no sample in that experiment.
#'
#' The roster is the subjects in \code{colData(mae)}. MultiAssayExperiment
#' drops \code{colData} rows that have no sample in any experiment, so every
#' subject in the result has at least one observed experiment.
#'
#' Each subject may have at most one sample per experiment. If an experiment
#' has replicate samples for a subject (several columns mapped to the same
#' primary), this function stops and names them; combine replicates first,
#' for example with \code{MultiAssayExperiment::mergeReplicates()} (see
#' \code{MultiAssayExperiment::replicated()} to find them).
#'
#' @param mae a MultiAssayExperiment object
#' @return named list of V numeric matrices, each n x p_v, with rows in the
#'   order of \code{rownames(colData(mae))}
#' @noRd
.mae_to_xlist <- function(mae) {
    # SummarizedExperiment is a dependency of MultiAssayExperiment, so it is
    # installed whenever MultiAssayExperiment is; both are checked explicitly.
    for (pkg in c("MultiAssayExperiment", "SummarizedExperiment")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            fmt <- paste0(
                "Package '%s' is required to pass a MultiAssayExperiment ",
                "to woven(). Install it with BiocManager::install('%s')."
            )
            stop(sprintf(
                fmt,
                pkg, pkg
            ), call. = FALSE)
        }
    }
    primary <- rownames(MultiAssayExperiment::colData(mae))
    exps <- MultiAssayExperiment::experiments(mae)
    smap <- as.data.frame(MultiAssayExperiment::sampleMap(mae))

    X_list <- lapply(names(exps), function(nm) {
        ex <- exps[[nm]]
        m <- if (methods::is(ex, "SummarizedExperiment")) {
            SummarizedExperiment::assay(ex)
        } else {
            as.matrix(ex)
        }
        map_i <- smap[smap$assay == nm, , drop = FALSE]
        map_i <- map_i[map_i$colname %in% colnames(m), , drop = FALSE]
        dup <- unique(map_i$primary[duplicated(map_i$primary)])
        if (length(dup) > 0L) {
            fmt <- paste0(
                "Experiment '%s' has more than one sample for subject(s) ",
                "%s. woven() needs one sample per subject per experiment: ",
                "combine replicates first, e.g. with ",
                "MultiAssayExperiment::mergeReplicates()."
            )
            stop(sprintf(
                fmt,
                nm, paste(utils::head(dup, 10L), collapse = ", ")
            ), call. = FALSE)
        }
        mt <- t(m) # features x samples -> samples x features
        rownames(mt) <- map_i$primary[match(rownames(mt), map_i$colname)]
        out <- matrix(NA_real_,
            nrow = length(primary), ncol = ncol(mt),
            dimnames = list(primary, colnames(mt))
        )
        common <- intersect(primary, rownames(mt))
        out[common, ] <- mt[common, , drop = FALSE]
        out
    })
    names(X_list) <- names(exps)
    X_list
}
