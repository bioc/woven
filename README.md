# woven

**Weighted Omics View Embedding via Nystrom**: supervised multi-omics
integration when some subjects are missing whole data types.

## The problem

In many cohorts, not every subject has every omics assay. Methods differ in
how they treat those subjects:

- **DIABLO** (mixOmics) is supervised but requires every subject to have
  every data type, so subjects missing one are removed.
- **MOFA2** accepts samples with whole views missing (coded as `NA`), but
  it is unsupervised: its factors are not steered towards an outcome.
- **missRows** handles missing individuals by multiple imputation within
  multiple factor analysis; it is unsupervised and imputes the missing
  blocks.

If the subjects with missing data differ from the rest, removing them biases
the analysis.

## What woven does

woven fits a shared, supervised latent space on the **anchor** subjects
(observed in every data type) and then scores every other subject from the
data types it has, using the learned projection matrices. No subject with at
least one data type is dropped, and no missing block is imputed.

- Supervised: class labels enter the cross-covariance, so latent dimensions
  are both shared across data types and class-discriminative.
- Closed form: one eigendecomposition; no iterations or random restarts.
- Tuned by cross-validation: `woven_tune()` chooses the ridge, feature
  screening and supervision settings by held-out balanced error rate. The
  defaults leave ridge and screening off, so tuning usually improves
  accuracy on high-dimensional data.
- Interpretable: one linear projection matrix per data type.
- Bioconductor input: takes a `MultiAssayExperiment` directly, or separate
  tables matched by subject ID.
- Readable output: progress and a coverage summary while fitting, one
  `plot()` entry point, and a top-features table.

## Installation

```r
# Bioconductor
BiocManager::install("woven")

# Development version
remotes::install_github("NathanBresette/woven")
```

## Quick start

```r
library(woven)
library(MultiAssayExperiment)
data(miniACC) # TCGA adrenocortical carcinoma subset

acc <- miniACC[, !is.na(miniACC$C1A.C1B),
               c("RNASeq2GeneNorm", "miRNASeqGene", "RPPAArray")]
for (e in c("RNASeq2GeneNorm", "miRNASeqGene"))
    acc[[e]] <- log2(assay(acc[[e]]) + 1)

fit <- woven(acc, Y = "C1A.C1B", K = 2L) # 78 tumours, 33 without RPPA
# or let cross-validation choose the settings, then use tune$fit:
# tune <- woven_tune(acc, Y = "C1A.C1B", K = 2L); fit <- tune$fit
plot(fit)                        # latent space, coloured by subtype
plot(fit, type = "features")     # top features per data type
plot(fit, type = "coverage")     # who has which data type
woven_top_features(fit)          # the same ranking as a table
head(latent_scores(fit))         # one row per subject
```

Without a MultiAssayExperiment, pass one table per data type (matrices or
data frames, subjects in rows). Tables with row names (subject IDs) may cover
different subjects in any order; woven matches them and treats a subject
absent from a table as missing that data type:

```r
fit <- woven(list(RNA = rna_table, Protein = protein_table),
             Y = labels_named_by_subject_id)
woven_predict(fit, list(RNA = new_rna)) # classify subjects with RNA only
```

## Documentation

```r
vignette("woven-intro", package = "woven")
```

See `?woven` for the main entry point, and `?woven_scores`,
`?woven_predict`, `?woven_top_features` and `?woven_metrics` for the
supporting functions.

## Citation

Bresette N, Lin A-L, Cheng J. Weighted Omics View Embedding via Nystrom
(WOVEN): supervised multi-omics integration for block-missing data. *In
preparation*, 2026.

## License

MIT
