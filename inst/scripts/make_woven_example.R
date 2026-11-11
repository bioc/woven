#!/usr/bin/env Rscript
# Generates the woven_example dataset shipped in data/woven_example.rda.
# All values are simulated (no real subjects); see ?woven_example.
#
# Run from the package source root:
#   Rscript inst/scripts/make_woven_example.R
# Output: data/woven_example.rda (identical to the shipped file; the script is
# fully seeded with set.seed(42) and set.seed(7)).

set.seed(42)
n <- 90L # 3 groups x 30 subjects
p1 <- 25L # RNA-seq (reduced)
p2 <- 20L # Methylation
p3 <- 15L # Proteomics
K_true <- 3L

# Group labels
Y <- rep(c("CN", "MCI", "AD"), each = 30L)

# Shared latent signal per group
Z_true <- matrix(0, n, K_true)
for (g in seq_len(3L)) {
    idx <- ((g - 1L) * 30L + 1L):(g * 30L)
    mu <- c(if (g == 1L) c(2, 0, 0) else if (g == 2L) c(0, 2, 0) else c(0, 0, 2))
    Z_true[idx, ] <- matrix(rnorm(30L * K_true, mean = mu, sd = 0.4), 30L, K_true)
}

make_view <- function(p) {
    W <- matrix(rnorm(K_true * p, sd = 0.5), K_true, p)
    X <- Z_true %*% W + matrix(rnorm(n * p, sd = 0.6), n, p)
    X <- scale(X)
    rownames(X) <- paste0("S", sprintf("%03d", seq_len(n)))
    colnames(X) <- paste0(deparse(substitute(p)), "_f", seq_len(p))
    X
}

RNA <- make_view(p1)
colnames(RNA) <- paste0("gene_", seq_len(p1))
Methyl <- make_view(p2)
colnames(Methyl) <- paste0("cpg_", seq_len(p2))
Prot <- make_view(p3)
colnames(Prot) <- paste0("prot_", seq_len(p3))

# Induce ~33% block-missing: each subject missing at most one modality
set.seed(7)
missing_subj <- sample(seq_len(n), size = 30L)
miss_mod <- sample(1:3, size = 30L, replace = TRUE)
RNA_miss <- RNA
Methyl_miss <- Methyl
Prot_miss <- Prot
for (i in seq_along(missing_subj)) {
    s <- missing_subj[i]
    if (miss_mod[i] == 1L) RNA_miss[s, ] <- NA_real_
    if (miss_mod[i] == 2L) Methyl_miss[s, ] <- NA_real_
    if (miss_mod[i] == 3L) Prot_miss[s, ] <- NA_real_
}


woven_example <- list(
    X_complete = list(RNA = RNA, Methylation = Methyl, Proteomics = Prot),
    X_missing  = list(RNA = RNA_miss, Methylation = Methyl_miss, Proteomics = Prot_miss),
    Y          = Y
)

save(woven_example,
    file = file.path("data", "woven_example.rda"),
    compress = "bzip2"
)
cat("woven_example saved to data/woven_example.rda\n")
