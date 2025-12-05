##############################################
#### 5B) AD-female APOE dataset (log2 RFU) ###
##############################################
suppressPackageStartupMessages({ library(dplyr); library(readr); library(stringr) })

out_dir  <- get0("out_dir", ifnotfound = "qc_outputs")
qc_path  <- file.path(out_dir, "03_log10_QC_matrix.csv")   # from Block 3
pheno_path <- "clinicalV1_filtered_ad_apoe.csv"
stopifnot(file.exists(qc_path), file.exists(pheno_path))

# Load QC matrix (log10 scale) and phenotype
qc <- readr::read_csv(qc_path, show_col_types = FALSE)
ph <- readr::read_csv(pheno_path, show_col_types = FALSE)

# Identify meta vs analytes in QC
meta_cols_qc <- intersect(c("sample_id",".__subset","contributor_code","visit","sample_type"), names(qc))
stopifnot("sample_id" %in% meta_cols_qc)
analyte_cols <- setdiff(names(qc), meta_cols_qc)

# Keep only SOMAmer analytes
is_somamer <- function(x) grepl("^seq[._][0-9]+[._][0-9]+$", tolower(x))
analyte_cols <- analyte_cols[is_somamer(analyte_cols)]

# Convert log10 to log2 so effects are in standard units
log10_to_log2 <- 1 / log10(2)  # ~3.321928
qc[analyte_cols] <- lapply(qc[analyte_cols], function(v) as.numeric(v) * log10_to_log2)

# ---- phenotype prep & filters ----
norm_apoe <- function(x){
  x <- tolower(as.character(x))
  x <- gsub("[^0-9e/]", "", x)
  x <- gsub("e", "", x)
  x <- gsub("(\\d)[\\s]*[-_:.,][\\s]*(\\d)", "\\1/\\2", x)
  x
}
has_allele <- function(gt, a) grepl(paste0("(^|/)", a, "(/|$)"), gt)

# Standardize IDs/columns
names(ph) <- tolower(gsub("[^A-Za-z0-9]+","_", names(ph)))
if (!"sample_id" %in% names(ph)) {
  cand <- names(ph)[grepl("^sample(_id)?$|somalogic|^sid$|^id$", names(ph))]
  stopifnot(length(cand) > 0); ph <- ph %>% rename(sample_id = !!cand[1])
}
ad_col  <- if ("ad" %in% names(ph)) "ad" else stop("Phenotype missing 'AD' column.")
sex_col <- if ("sex" %in% names(ph)) "sex" else stop("Phenotype missing 'Sex' column.")
mmse_col <- if ("cognitive_test_score" %in% names(ph)) "cognitive_test_score" else stop("Phenotype missing MMSE column (cognitive_test_score).")
apoe_cands <- names(ph)[grepl("apoe", names(ph), ignore.case = TRUE)]
stopifnot(length(apoe_cands) > 0)
apoe_col <- apoe_cands[1]

# Merge & filter: AD females with MMSE < 24; group as E3E3 vs E4_plus (exclude E2)
dat <- qc %>%
  left_join(ph %>% select(sample_id, all_of(c(ad_col, sex_col, apoe_col, mmse_col))), by = "sample_id") %>%
  mutate(
    ad_flag  = .data[[ad_col]] == 1,
    sex_f    = .data[[sex_col]] == 2,
    mmse_ok  = suppressWarnings(as.numeric(.data[[mmse_col]])) < 24,
    apoe_std = norm_apoe(.data[[apoe_col]]),
    has_e2   = has_allele(apoe_std, "2"),
    has_e4   = has_allele(apoe_std, "4"),
    has_e3   = has_allele(apoe_std, "3"),
    apoe_group = dplyr::case_when(
      !has_e2 & has_e4 ~ "E4_plus",
      !has_e2 & !has_e4 & has_e3 ~ "E3E3",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(ad_flag, sex_f, mmse_ok, !is.na(apoe_group)) %>%
  relocate(sample_id, apoe_group, .__subset, contributor_code, sample_type, visit)

# Save Block 5B dataset (log2 RFU)
outf <- file.path(out_dir, "05B_AD_female_APOE_log2.csv")
readr::write_csv(dat, outf)
message("Saved: ", outf)

##############################################
#### 6) limma DE on log2 + Volcano ##########
##############################################
suppressPackageStartupMessages({ library(dplyr); library(readr); library(limma); library(ggplot2) })
have_EV <- requireNamespace("EnhancedVolcano", quietly = TRUE)

# Reload the 5B dataset (ensures 'dat' is present even in a fresh session)
dat <- readr::read_csv(file.path(out_dir, "05B_AD_female_APOE_log2.csv"), show_col_types = FALSE)

# Build expression matrix X (rows=features = SOMAmer, cols=samples) from *this* dataset
is_somamer <- function(x) grepl("^seq[._][0-9]+[._][0-9]+$", tolower(x))
feat_cols <- names(dat)[is_somamer(names(dat))]
stopifnot(length(feat_cols) > 0)

X <- t(as.matrix(sapply(dat[feat_cols], function(v) suppressWarnings(as.numeric(v)))))
rownames(X) <- feat_cols
colnames(X) <- dat$sample_id

# Drop zero-variance features
keep_feat <- apply(X, 1, function(z) isTRUE(var(z, na.rm = TRUE) > 0))
X <- X[keep_feat, , drop = FALSE]

# ---- Design: E3/3 (control) vs E4 carriers (case) with pruning ----
grp <- factor(dat$apoe_group, levels = c("E3E3","E4_plus"))

covars <- data.frame(sample_id = dat$sample_id, stringsAsFactors = FALSE)
if ("contributor_code" %in% names(dat)) covars$contributor <- factor(dat$contributor_code)
if (".__subset" %in% names(dat))        covars$subset      <- factor(dat$.__subset)
if ("PC2" %in% names(dat))              covars$PC2         <- dat$PC2 else if ("pc2" %in% names(dat)) covars$PC2 <- dat$pc2
if ("PC3" %in% names(dat))              covars$PC3         <- dat$PC3 else if ("pc3" %in% names(dat)) covars$PC3 <- dat$pc3
if ("age_at_visit" %in% names(dat))     covars$age         <- dat$age_at_visit else
  if ("age" %in% names(dat))            covars$age         <- dat$age

design_df <- cbind(sample_id = dat$sample_id, grp = grp, covars[setdiff(names(covars), "sample_id")])
keep_samp <- complete.cases(design_df)
if (!all(keep_samp)) message("Dropping ", sum(!keep_samp), " samples due to missing covariates.")
design_df <- design_df[keep_samp, , drop = FALSE]
X <- X[, design_df$sample_id, drop = FALSE]

# ---- Robust covariate pruning (vapply) ----
dfmod <- design_df
char_cols <- setdiff(names(dfmod)[vapply(dfmod, is.character, logical(1))], c("sample_id"))
if (length(char_cols)) dfmod[char_cols] <- lapply(dfmod[char_cols], factor)

all_na   <- names(dfmod)[vapply(dfmod, function(x) all(is.na(x)), logical(1))]
fac_cols <- names(dfmod)[vapply(dfmod, is.factor,  logical(1))]
one_lvl  <- setdiff(fac_cols[vapply(dfmod[fac_cols], function(x) nlevels(x) < 2, logical(1))], "grp")
num_cols <- names(dfmod)[vapply(dfmod, is.numeric, logical(1))]
zero_var <- character(0)
if (length(num_cols)) {
  zero_var_flags <- vapply(dfmod[num_cols], function(x) {
    x <- suppressWarnings(as.numeric(x))
    v <- stats::var(x, na.rm = TRUE)
    isTRUE(v == 0) || !any(is.finite(x))
  }, logical(1))
  zero_var <- names(zero_var_flags)[zero_var_flags]
}
drop_cols <- setdiff(unique(c(all_na, one_lvl, zero_var)), c("grp","sample_id"))
if (length(drop_cols)) {
  message("Pruning covariates: ", paste(drop_cols, collapse = ", "))
  dfmod <- dfmod[, setdiff(names(dfmod), drop_cols), drop = FALSE]
}

dfmod$grp <- droplevels(dfmod$grp)
stopifnot(nlevels(dfmod$grp) >= 2)

design <- model.matrix(~ 0 + grp + ., data = dfmod[, setdiff(names(dfmod), "sample_id"), drop = FALSE])
colnames(design)[colnames(design) == "grpE3E3"]    <- "E3E3"
colnames(design)[colnames(design) == "grpE4_plus"] <- "E4_plus"

# Align X again to dfmod sample order (just to be extra safe)
X <- X[, dfmod$sample_id, drop = FALSE]

# ---- Fit limma on log2 intensities -------------------

# (Optional) drop features with heavy missingness; tweak threshold if you like
keep_nona <- rowMeans(!is.na(X)) >= 0.80   # keep features with >=80% observed values
X <- X[keep_nona, , drop = FALSE]

# Sample-level precision weights (helps when some samples are noisier)
aw <- limma::arrayWeights(X, design)  # length = ncol(X)

# Expand to an observation-weight matrix for lmFit
W <- matrix(aw, nrow = nrow(X), ncol = ncol(X), byrow = TRUE)

# Fit: contrast = E3E3 − E4_plus  (positive log2FC = higher in E3/3)
fit  <- limma::lmFit(X, design, weights = W)
cont <- makeContrasts(E4plus_vs_E3E3 = E4_plus - E3E3, levels = design)
fit2 <- limma::eBayes(limma::contrasts.fit(fit, cont), trend = TRUE, robust = TRUE)

# Results
res <- limma::topTable(fit2, coef = "E3E3_vs_E4plus", number = Inf, sort.by = "P")
res$feature_id <- rownames(res)
res$FDR <- p.adjust(res$P.Value, "BH")

# Save + diagnostic
out_dir <- get0("out_dir", ifnotfound = "qc_outputs")
readr::write_csv(res, file.path(out_dir, "06B_limma_log2_E3E3_vs_E4plus.csv"))
png(file.path(out_dir, "06B_limma_plotSA.png"), width = 800, height = 600)
limma::plotSA(fit2, main = "Moderated s2 vs mean (trend=TRUE, robust=TRUE)")
dev.off()


############################################################
## Exact-match annotation + volcano (labels = target)     ##
############################################################
# Requires 'analyte_info' in memory with: column_name, target, entrez_gene_symbol, uni_prot, etc.
stopifnot(exists("analyte_info"))
req_ann <- c("column_name","target","entrez_gene_symbol","uni_prot")
stopifnot(all(req_ann %in% names(analyte_info)))

ann_slice <- analyte_info %>%
  transmute(column_name = as.character(column_name),
            target, entrez_gene_symbol, uni_prot)

res2 <- readr::read_csv(res_path, show_col_types = FALSE)
res_annot <- res2 %>%
  mutate(feature_id = as.character(feature_id)) %>%
  left_join(ann_slice, by = c("feature_id" = "column_name")) %>%
  mutate(label = dplyr::coalesce(
    ifelse(!is.na(target) & target != "", target, NA_character_),
    ifelse(!is.na(entrez_gene_symbol) & entrez_gene_symbol != "", entrez_gene_symbol, NA_character_),
    feature_id
  ))

readr::write_csv(res_annot, file.path(out_dir, "06B_limma_log2_E3E3_vs_E4plus_annotated.csv"))

# Volcano with target labels (EnhancedVolcano if available; else ggplot)
if (have_EV) {
  library(EnhancedVolcano)
  sig <- res_annot %>% filter(is.finite(FDR), FDR < 0.05) %>% arrange(FDR, dplyr::desc(abs(logFC)))
  p <- EnhancedVolcano::EnhancedVolcano(
    res_annot,
    lab = res_annot$label, x = "logFC", y = "FDR",
    pCutoff = 0.05, FCcutoff = 0,
    selectLab = head(unique(res_annot$label[match(sig$feature_id, res_annot$feature_id)]), 20),
    drawConnectors = TRUE, boxedLabels = TRUE, widthConnectors = 0.4, colAlpha = 0.6, labSize = 3,
    title = "E3/3 vs E4 carriers — log2 scale",
    subtitle = "+log2FC = higher in E3/3 controls",
    ylab = expression(-log[10]("FDR"))
  )
  ggsave(file.path(out_dir, "09B_volcano_log2_target.png"), p, width = 10, height = 7, dpi = 150)
} else {
  res_annot$neglog10FDR <- -log10(pmax(res_annot$FDR, 1e-300))
  sig20 <- res_annot %>% filter(FDR < 0.05) %>% arrange(FDR, dplyr::desc(abs(logFC))) %>% head(20)
  if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
  library(ggrepel)
  p <- ggplot(res_annot, aes(logFC, neglog10FDR)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text_repel(data = sig20, aes(label = label), size = 3, max.overlaps = 30) +
    labs(title = "E3/3 vs E4 carriers — log2 scale",
         subtitle = "+log2FC = higher in E3/3 controls",
         x = "log2 fold-change (E3/3 − E4 carriers)", y = expression(-log[10]("FDR"))) +
    theme_bw()
  ggsave(file.path(out_dir, "09B_volcano_log2_target.png"), p, width = 10, height = 7, dpi = 150)
}

# ===== Summaries, effect-size filter, per-gene collapse, enrichment helpers ===
suppressPackageStartupMessages({ library(dplyr); library(readr); library(ggplot2) })

out_dir <- get0("out_dir", ifnotfound = "qc_outputs")
res_ann <- readr::read_csv(file.path(out_dir, "06B_limma_log2_E3E3_vs_E4plus_annotated.csv"),
                           show_col_types = FALSE)

# 1) Direction counts and top tables (FDR<0.05)
sig <- res_ann %>% filter(is.finite(FDR), FDR < 0.05, is.finite(logFC))
c(E3E3_higher = sum(sig$logFC > 0), E4carriers_higher = sum(sig$logFC < 0))

top_E3E3 <- sig %>% arrange(desc(logFC), FDR) %>% head(30)
top_E4   <- sig %>% arrange(logFC, FDR)       %>% head(30)
readr::write_csv(top_E3E3, file.path(out_dir, "08_top30_higher_in_E3E3_log2.csv"))
readr::write_csv(top_E4,   file.path(out_dir, "08_top30_higher_in_E4carriers_log2.csv"))

# 2) Add an effect-size floor (|log2FC| >= 0.2) for reporting
sig_strict <- res_ann %>% filter(FDR < 0.05, abs(logFC) >= 0.2)
readr::write_csv(sig_strict, file.path(out_dir, "08_sig_FDRlt0.05_lfc0.2.csv"))

# 3) Per-gene consolidation (one row per entrez_gene_symbol)
collapsed <- res_ann %>%
  filter(!is.na(entrez_gene_symbol) & entrez_gene_symbol != "") %>%
  group_by(entrez_gene_symbol) %>%
  slice_max(order_by = abs(t), n = 1, with_ties = FALSE) %>%
  ungroup()
readr::write_csv(collapsed, file.path(out_dir, "08_collapsed_genelevel.csv"))

# Check for discordant-sign targets with multiple aptamers
discordant <- res_ann %>%
  filter(!is.na(entrez_gene_symbol) & entrez_gene_symbol != "") %>%
  group_by(entrez_gene_symbol) %>%
  summarize(n = n(),
            n_pos = sum(logFC > 0, na.rm=TRUE),
            n_neg = sum(logFC < 0, na.rm=TRUE),
            mixed = n_pos > 0 & n_neg > 0,
            .groups="drop") %>%
  filter(mixed)
readr::write_csv(discordant, file.path(out_dir, "08_discordant_aptamers.csv"))

# 4) Quick group means (so tables show magnitude on the same scale used in DE)
dat_log2 <- readr::read_csv(file.path(out_dir, "05B_AD_female_APOE_log2.csv"), show_col_types = FALSE)
is_somamer <- function(x) grepl("^seq[._][0-9]+[._][0-9]+$", tolower(x))
feat_cols <- names(dat_log2)[is_somamer(names(dat_log2))]

means <- dat_log2 %>%
  select(apoe_group, all_of(feat_cols)) %>%
  group_by(apoe_group) %>%
  summarise(across(all_of(feat_cols), ~ mean(as.numeric(.), na.rm = TRUE)), .groups="drop") %>%
  tidyr::pivot_longer(-apoe_group, names_to="feature_id", values_to="group_mean") %>%
  tidyr::pivot_wider(names_from=apoe_group, values_from=group_mean, names_prefix="mean_")

report <- res_ann %>%
  left_join(means, by="feature_id") %>%
  mutate(direction = ifelse(logFC > 0, "higher_in_E3E3", "higher_in_E4carriers"),
         mean_diff = mean_E3E3 - mean_E4_plus) %>%
  arrange(FDR, desc(abs(logFC)))
readr::write_csv(report, file.path(out_dir, "08_report_log2_with_group_means.csv"))

# 5) Optional: minimum effect testing with limma::treat (lfc = 0.2)
# (run this only if you still have 'fit'/'design' and can refit quickly)
# library(limma)
# fit_treat <- limma::treat(limma::contrasts.fit(limma::lmFit(X, design), makeContrasts(E3E3 - E4_plus, levels=design)), lfc = 0.2)
# res_treat <- limma::topTable(limma::eBayes(fit_treat), number = Inf)
# readr::write_csv(res_treat, file.path(out_dir, "06B_limma_log2_treat_lfc0.2.csv"))
