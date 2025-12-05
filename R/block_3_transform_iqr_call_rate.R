##############################################
#### 3) Transform + IQR + call-rate QC #######
##############################################

suppressPackageStartupMessages({ library(dplyr); library(readr) })

# ---- config ----
if (!exists("out_dir")) out_dir <- "qc_outputs_filtered_by_ad_sex_apoe"
in_path <- file.path(out_dir, "02_expr_annotated_filtered.csv")
stopifnot(file.exists(in_path))


# ---- helpers ----
call_rate_cols <- function(mat) 100 * (1 - colMeans(is.na(mat)))
call_rate_rows <- function(mat) 100 * (1 - rowMeans(is.na(mat)))

log10_safely <- function(M) {
  M <- suppressWarnings(log10(M))
  M[!is.finite(M)] <- NA_real_
  M
}

iqr_mask <- function(M) {
  # Column-wise IQR outlier masking: replace values outside [Q1-1.5*IQR, Q3+1.5*IQR] with NA
  q1  <- apply(M, 2, quantile, probs = 0.25, na.rm = TRUE)
  q3  <- apply(M, 2, quantile, probs = 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  for (j in seq_len(ncol(M))) {
    y <- M[, j]
    if (all(is.na(y))) next
    lb <- q1[j] - 1.5 * iqr[j]
    ub <- q3[j] + 1.5 * iqr[j]
    idx <- !is.na(y) & (y < lb | y > ub)
    if (any(idx)) y[idx] <- NA_real_
    M[, j] <- y
  }
  M
}

# ---- load annotated expression (SOMAmer IDs retained) ----
expr_filt <- readr::read_csv(in_path, show_col_types = FALSE)

# Meta columns present from Block 2
meta_cols <- intersect(c("sample_id","contributor_code","visit","sample_type"), names(expr_filt))
stopifnot("sample_id" %in% meta_cols)

# Ensure analyte columns are numeric
analyte_cols <- setdiff(names(expr_filt), meta_cols)
expr_filt <- expr_filt %>%
  mutate(across(all_of(analyte_cols), ~ suppressWarnings(as.numeric(.))))

# Split by sample_type if available; otherwise single group "all"
subset_col <- if ("sample_type" %in% names(expr_filt)) "sample_type" else NULL
split_list <- if (!is.null(subset_col)) split(expr_filt, expr_filt[[subset_col]]) else list(all = expr_filt)

do_qc <- function(df, tag) {
  stopifnot("sample_id" %in% names(df))
  # Separate meta and analytes
  Y <- df %>% select(-any_of(meta_cols)) %>% as.matrix()
  storage.mode(Y) <- "double"
  rownames(Y) <- df$sample_id
  
  # 1) log10 transform
  Y <- log10_safely(Y)
  
  # 2) IQR outlier masking
  Y <- iqr_mask(Y)
  
  # 3) Call-rate pass 1 (65%)
  cr_c1 <- call_rate_cols(Y)
  cr_r1 <- call_rate_rows(Y)
  keep_c1 <- names(cr_c1)[cr_c1 >= 65]
  keep_r1 <- rownames(Y)[cr_r1 >= 65]
  Y1 <- Y[keep_r1, keep_c1, drop = FALSE]
  
  # 4) Call-rate pass 2 (85%)
  cr_c2 <- call_rate_cols(Y1)
  cr_r2 <- call_rate_rows(Y1)
  keep_c2 <- names(cr_c2)[cr_c2 >= 85]
  keep_r2 <- rownames(Y1)[cr_r2 >= 85]
  Y2 <- Y1[keep_r2, keep_c2, drop = FALSE]
  
  # Save diagnostics
  readr::write_csv(
    data.frame(analyte = names(cr_c1), call_rate = unname(cr_c1)),
    file.path(out_dir, sprintf("CR_%s_analytes_pass1_65.csv", tag))
  )
  readr::write_csv(
    data.frame(sample_id = rownames(Y), call_rate = unname(cr_r1)),
    file.path(out_dir, sprintf("CR_%s_samples_pass1_65.csv", tag))
  )
  readr::write_csv(
    data.frame(analyte = names(cr_c2), call_rate = unname(cr_c2)),
    file.path(out_dir, sprintf("CR_%s_analytes_pass2_85.csv", tag))
  )
  readr::write_csv(
    data.frame(sample_id = rownames(Y1), call_rate = unname(cr_r2)),
    file.path(out_dir, sprintf("CR_%s_samples_pass2_85.csv", tag))
  )
  
  # Carry meta for kept samples
  meta_keep <- df %>%
    select(any_of(meta_cols)) %>%
    filter(sample_id %in% rownames(Y2))
  
  list(M_log10_qc = Y2, meta = meta_keep)
}

# Run QC per subset
qc_list <- lapply(names(split_list), function(nm) do_qc(split_list[[nm]], nm))
names(qc_list) <- names(split_list)

# Union of all features across subsets
all_feats <- Reduce(union, lapply(qc_list, function(x) colnames(x$M_log10_qc)))

# Build a combined data frame
qc_full <- dplyr::bind_rows(lapply(names(qc_list), function(nm) {
  M <- qc_list[[nm]]$M_log10_qc
  # add missing features as NA to align columns
  add <- setdiff(all_feats, colnames(M))
  if (length(add)) {
    M <- cbind(M, matrix(NA_real_, nrow = nrow(M), ncol = length(add),
                         dimnames = list(rownames(M), add)))
  }
  # order columns lexicographically for reproducibility
  M <- M[, sort(colnames(M)), drop = FALSE]
  
  df <- as.data.frame(M)
  df$sample_id <- rownames(M)
  df$.__subset <- nm
  
  # attach meta aligned by sample_id
  df <- qc_list[[nm]]$meta %>%
    right_join(df, by = "sample_id") %>%
    relocate(sample_id, .__subset, any_of(setdiff(meta_cols, "sample_id")))
  df
}))

# Write unified QC matrix
out_qc <- file.path(out_dir, "03_log10_QC_matrix.csv")
readr::write_csv(qc_full, out_qc)
message("Saved QC matrix: ", out_qc)
message("QC dims: ", nrow(qc_full), " samples x ", length(setdiff(names(qc_full), meta_cols)) , " analytes (post-QC)")
#QC dims: 1474 samples x 7284 analytes (post-QC)