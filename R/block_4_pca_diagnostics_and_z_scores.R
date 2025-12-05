##############################################
#### 4) PCA diagnostics + z-scores ###########
##############################################

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(tidyr)
})

if (!exists("out_dir")) out_dir <- "qc_outputs"
qc_path <- file.path(out_dir, "03_log10_QC_matrix.csv")
stopifnot(file.exists(qc_path))

# ---- tiny bootstrap imputer (PCA only)
bootstrap_impute <- function(M) {
  M <- as.matrix(M); storage.mode(M) <- "double"
  for (j in seq_len(ncol(M))) {
    miss <- is.na(M[, j])
    if (any(miss)) {
      pool <- M[!miss, j]
      if (length(pool)) M[miss, j] <- sample(pool, sum(miss), replace = TRUE)
    }
  }
  M
}

# ---- load matrix
qc_full <- readr::read_csv(qc_path, show_col_types = FALSE)
stopifnot("sample_id" %in% names(qc_full))

# metadata columns (keep whatever exists)
meta_cols <- intersect(c("sample_id", ".__subset", "contributor_code", "visit", "sample_type"),
                       names(qc_full))
feat_cols <- setdiff(names(qc_full), meta_cols)  # proteins only

# ---- PCA on imputed matrix (all subsets together)
M <- qc_full %>% select(all_of(feat_cols))
row_ids <- qc_full$sample_id
M_imp <- bootstrap_impute(M)
rownames(M_imp) <- row_ids

use_irlba <- requireNamespace("irlba", quietly = TRUE) &&
  (nrow(M_imp) > 5000 || ncol(M_imp) > 3000)

if (use_irlba) {
  message("Using irlba::prcomp_irlba (top 10 PCs).")
  pc <- irlba::prcomp_irlba(M_imp, n = 10, center = TRUE, scale. = TRUE)
} else {
  message("Using stats::prcomp.")
  pc <- prcomp(M_imp, center = TRUE, scale. = TRUE)
}

# Scree + save
pct <- 100 * pc$sdev^2 / sum(pc$sdev^2)
scree_df <- data.frame(PC = seq_along(pct), VarExp = pct)
ggplot(scree_df[1:min(20, nrow(scree_df)),], aes(PC, VarExp)) +
  geom_col() + geom_line() + geom_point() +
  labs(title = "Scree (top PCs)", y = "% variance explained") +
  theme_minimal()
ggsave(file.path(out_dir, "PCA_scree.png"), width = 7, height = 5, dpi = 120)

# PC scores (first 3) + subset/color
scores <- as.data.frame(pc$x[, 1:min(3, ncol(pc$x)), drop = FALSE])
scores$sample_id <- rownames(M_imp)

# join subset/sample_type only if present
join_cols <- intersect(c("sample_id", ".__subset", "sample_type"), names(qc_full))
scores <- scores %>% left_join(qc_full %>% select(all_of(join_cols)), by = "sample_id")

readr::write_csv(scores, file.path(out_dir, "04_pca_scores.csv"))

p <- ggplot(scores, aes(PC1, PC2, color = .__subset %||% "all")) +
  geom_point(size = 1.2, alpha = .85) +
  theme_minimal() +
  labs(title = "PCA: PC1 vs PC2", color = "subset")
ggsave(file.path(out_dir, "PCA_PC1_PC2.png"), p, width = 7, height = 5, dpi = 120)

# ---- Z-scores within each subset (on log10 scale)
# IMPORTANT: do NOT return .__subset inside group_modify’s data; dplyr adds it back
meta_cols_nosubset <- setdiff(meta_cols, ".__subset")

zs_by_subset <- qc_full %>%
  {
    if (".__subset" %in% names(.)) group_by(., .__subset) else mutate(., .__subset := "all") %>% group_by(.__subset)
  } %>%
  group_modify(~{
    X  <- as.matrix(select(.x, all_of(feat_cols)))
    Xz <- scale(X, center = TRUE, scale = TRUE)
    meta_keep <- select(.x, all_of(meta_cols_nosubset))
    bind_cols(meta_keep, as.data.frame(Xz), .name_repair = "minimal")
  }) %>%
  ungroup() %>%
  relocate(sample_id, .__subset,
           any_of(setdiff(meta_cols, c("sample_id", ".__subset"))))

readr::write_csv(zs_by_subset, file.path(out_dir, "04_zscores_by_subset.csv"))
