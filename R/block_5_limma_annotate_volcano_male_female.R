# scripts/05B_MvF_AD_APOE4no24_limma.R
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr)
  library(limma); library(ggplot2)
})

out_dir    <- get0("out_dir", ifnotfound = "qc_outputs")
qc_path    <- file.path(out_dir, "03_log10_QC_matrix.csv")
pheno_path <- "clinicalV1_filtered_ad_apoe.csv"
stopifnot(file.exists(qc_path), file.exists(pheno_path))

# ---------- USER: exact column names & codes ----------
sample_id_col_user <- "sample_id"
ad_col_user        <- "ad"
sex_col_user       <- "sex"
apoe_col_user      <- "apoe"                    # set your APOE col
mmse_col_user      <- "cognitive_test_score"    # per your data
pd_col_user        <- "pd"
ftd_col_user       <- "ftd"
als_col_user       <- "als"
sex_male_code      <- 1
sex_female_code    <- 2
# ------------------------------------------------------

# ---------- helpers ----------
clean_names_local <- function(x){
  x <- gsub("[^A-Za-z0-9]+","_",x); x <- gsub("_+","_",x)
  x <- gsub("^_|_$","",x); tolower(x)
}
.norm_name <- function(nm) clean_names_local(nm)
.clean_id  <- function(x) tolower(trimws(as.character(x)))
is_somamer <- function(x) grepl("^seq[._][0-9]+[._][0-9]+$", tolower(x))
norm_apoe <- function(x){
  x <- tolower(as.character(x))
  x <- gsub("[^0-9/]", "", x)
  x <- gsub("(\\d)[\\s]*[-_:.,][\\s]*(\\d)", "\\1/\\2", x)
  x
}
has_allele <- function(gt, a) grepl(paste0("(^|/)", a, "(/|$)"), gt)

# ---------- load & prep QC (log10 → log2) ----------
qc <- readr::read_csv(qc_path, show_col_types = FALSE)
names(qc) <- clean_names_local(names(qc))
stopifnot("sample_id" %in% names(qc))

meta_cols_qc <- intersect(c("sample_id",".__subset","contributor_code","visit","sample_type"), names(qc))
analyte_cols <- setdiff(names(qc), meta_cols_qc)
analyte_cols <- analyte_cols[is_somamer(analyte_cols)]
stopifnot(length(analyte_cols) > 0)

log10_to_log2 <- 1 / log10(2)
qc[analyte_cols] <- lapply(qc[analyte_cols], function(v) as.numeric(v) * log10_to_log2)

# ---------- phenotype ----------
ph <- readr::read_csv(pheno_path, show_col_types = FALSE)
names(ph) <- clean_names_local(names(ph))

sample_id_col <- .norm_name(sample_id_col_user)
ad_col        <- .norm_name(ad_col_user)
sex_col       <- .norm_name(sex_col_user)
apoe_col      <- .norm_name(apoe_col_user)
mmse_col      <- .norm_name(mmse_col_user)
pd_col        <- .norm_name(pd_col_user)
ftd_col       <- .norm_name(ftd_col_user)
als_col       <- .norm_name(als_col_user)

need <- c(sample_id_col, ad_col, sex_col, apoe_col, mmse_col, pd_col, ftd_col, als_col)
miss <- setdiff(need, names(ph))
if (length(miss)) stop("Missing phenotype columns: ", paste(miss, collapse = ", "))

ph[[sample_id_col]] <- .clean_id(ph[[sample_id_col]])
qc[["sample_id"]]   <- .clean_id(qc[["sample_id"]])

# ---------- merge & filter ----------
dat <- qc %>%
  left_join(ph %>% select(all_of(need)), by = c("sample_id" = sample_id_col)) %>%
  mutate(
    ad_flag  = .data[[ad_col]] == 1,
    sex_m    = suppressWarnings(.data[[sex_col]] == sex_male_code),
    sex_f    = suppressWarnings(.data[[sex_col]] == sex_female_code),
    mmse_val = suppressWarnings(as.numeric(.data[[mmse_col]])),
    apoe_std = norm_apoe(.data[[apoe_col]]),
    has_e2   = has_allele(apoe_std, "2"),
    has_e4   = has_allele(apoe_std, "4"),
    pd_val   = suppressWarnings(as.numeric(.data[[pd_col]])),
    ftd_val  = suppressWarnings(as.numeric(.data[[ftd_col]])),
    als_val  = suppressWarnings(as.numeric(.data[[als_col]]))
  ) %>%
  filter(
    ad_flag,
    (sex_m | sex_f),
    !is.na(mmse_val) & mmse_val < 24,
    !is.na(apoe_std) & has_e4 & !has_e2,    # exclude 2/4
    !is.na(pd_val)  & pd_val  == 0,
    !is.na(ftd_val) & ftd_val == 0,
    !is.na(als_val) & als_val == 0
  ) %>%
  mutate(
    group = ifelse(sex_m, "Male_E4", "Female_E4"),
    group = factor(group, levels = c("Female_E4","Male_E4"))
  ) %>%
  relocate(sample_id, group,
           any_of(c(".__subset","contributor_code","sample_type","visit",
                    mmse_col, pd_col, ftd_col, als_col)))

outf5 <- file.path(out_dir, "05B_AD_APOE4no24_MvF_log2.csv")
readr::write_csv(dat, outf5)
message("Saved: ", outf5)

# ==================== LIMMA ====================
have_EV <- requireNamespace("EnhancedVolcano", quietly = TRUE)
dat <- readr::read_csv(outf5, show_col_types = FALSE)
feat_cols <- names(dat)[is_somamer(names(dat))]
stopifnot(length(feat_cols) > 0)

X <- t(as.matrix(sapply(dat[feat_cols], function(v) suppressWarnings(as.numeric(v)))))
rownames(X) <- feat_cols
colnames(X) <- dat$sample_id

# drop 0-var
keep_feat <- apply(X, 1, function(z) isTRUE(var(z, na.rm = TRUE) > 0))
X <- X[keep_feat, , drop = FALSE]

# group & covariates
grp <- factor(dat$group, levels = c("Female_E4","Male_E4"))

covars <- data.frame(sample_id = dat$sample_id, stringsAsFactors = FALSE)
if ("contributor_code" %in% names(dat)) covars$contributor <- factor(dat$contributor_code)
if (".__subset" %in% names(dat))        covars$subset      <- factor(dat$.__subset)
if ("PC2" %in% names(dat))              covars$PC2         <- dat$PC2 else if ("pc2" %in% names(dat)) covars$PC2 <- dat$pc2
if ("PC3" %in% names(dat))              covars$PC3         <- dat$PC3 else if ("pc3" %in% names(dat)) covars$PC3 <- dat$pc3
if ("age_at_visit" %in% names(dat))     covars$age         <- dat$age_at_visit else if ("age" %in% names(dat)) covars$age <- dat$age

design_df <- cbind(sample_id = dat$sample_id, grp = grp, covars[setdiff(names(covars), "sample_id")])
keep_samp <- complete.cases(design_df)
if (!all(keep_samp)) message("Dropping ", sum(!keep_samp), " samples due to missing covariates.")
design_df <- design_df[keep_samp, , drop = FALSE]
X <- X[, design_df$sample_id, drop = FALSE]

# prune bad covariates
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
    x <- suppressWarnings(as.numeric(x)); v <- stats::var(x, na.rm = TRUE)
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
colnames(design)[colnames(design) == "grpFemale_E4"] <- "Female_E4"
colnames(design)[colnames(design) == "grpMale_E4"]   <- "Male_E4"

# align again
X <- X[, dfmod$sample_id, drop = FALSE]

# missingness filter
keep_nona <- rowMeans(!is.na(X)) >= 0.80
X <- X[keep_nona, , drop = FALSE]

# weights & fit
aw <- limma::arrayWeights(X, design)
W  <- matrix(aw, nrow = nrow(X), ncol = ncol(X), byrow = TRUE)

fit  <- limma::lmFit(X, design, weights = W)
cont <- limma::makeContrasts(Male_minus_Female = Male_E4 - Female_E4, levels = design)
fit2 <- limma::eBayes(limma::contrasts.fit(fit, cont), trend = TRUE, robust = TRUE)

res <- limma::topTable(fit2, coef = "Male_minus_Female", number = Inf, sort.by = "P")
res$feature_id <- rownames(res)
res$FDR <- p.adjust(res$P.Value, "BH")

res_path <- file.path(out_dir, "06B_limma_log2_Male_minus_Female.csv")
readr::write_csv(res, res_path)
png(file.path(out_dir, "06B_limma_plotSA_MvF.png"), width = 800, height = 600)
limma::plotSA(fit2, main = "Moderated s2 vs mean (trend=TRUE, robust=TRUE)")
dev.off()

# ----------------- annotation + volcano -----------------
analyte_info <- read_csv("SomalogicAnalyteInfoV1_anonymized.csv")
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

readr::write_csv(res_annot, file.path(out_dir, "06B_limma_log2_Male_minus_Female_annotated.csv"))

if (requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  library(EnhancedVolcano)
  sig <- res_annot %>% filter(is.finite(FDR), FDR < 0.05) %>% arrange(FDR, dplyr::desc(abs(logFC)))
  p <- EnhancedVolcano::EnhancedVolcano(
    res_annot,
    lab = res_annot$label, x = "logFC", y = "FDR",
    pCutoff = 0.05, FCcutoff = 0,
    selectLab = head(unique(res_annot$label[match(sig$feature_id, res_annot$feature_id)]), 20),
    drawConnectors = TRUE, boxedLabels = TRUE, widthConnectors = 0.4, colAlpha = 0.6, labSize = 3,
    title = "Male (ε4, AD, MMSE<24) − Female (ε4, AD, MMSE<24)",
    subtitle = "+log2FC = higher in Males",
    ylab = expression(-log[10]("FDR"))
  )
  ggsave(file.path(out_dir, "09B_volcano_log2_Male_minus_Female.png"), p, width = 10, height = 7, dpi = 150)
} else {
  res_annot$neglog10FDR <- -log10(pmax(res_annot$FDR, 1e-300))
  if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
  library(ggrepel)
  sig20 <- res_annot %>% filter(FDR < 0.05) %>% arrange(FDR, dplyr::desc(abs(logFC))) %>% head(20)
  p <- ggplot(res_annot, aes(logFC, neglog10FDR)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text_repel(data = sig20, aes(label = label), size = 3, max.overlaps = 30) +
    labs(title = "Male (ε4, AD, MMSE<24) − Female (ε4, AD, MMSE<24)",
         subtitle = "+log2FC = higher in Males",
         x = "log2 fold-change (Male − Female)", y = expression(-log[10]("FDR"))) +
    theme_bw()
  ggsave(file.path(out_dir, "09B_volcano_log2_Male_minus_Female.png"), p, width = 10, height = 7, dpi = 150)
}

# ----------------- reports -----------------
res_ann <- readr::read_csv(file.path(out_dir, "06B_limma_log2_Male_minus_Female_annotated.csv"), show_col_types = FALSE)

sig <- res_ann %>% filter(is.finite(FDR), FDR < 0.05, is.finite(logFC))
c(Male_higher = sum(sig$logFC > 0), Female_higher = sum(sig$logFC < 0))

top_Male   <- sig %>% arrange(desc(logFC), FDR) %>% head(30)
top_Female <- sig %>% arrange(logFC, FDR)       %>% head(30)
readr::write_csv(top_Male,   file.path(out_dir, "08_top30_higher_in_Male_log2.csv"))
readr::write_csv(top_Female, file.path(out_dir, "08_top30_higher_in_Female_log2.csv"))

feat_cols <- names(dat)[is_somamer(names(dat))]
means <- dat %>%
  select(group, all_of(feat_cols)) %>%
  group_by(group) %>%
  summarise(across(all_of(feat_cols), ~ mean(as.numeric(.), na.rm = TRUE)), .groups="drop") %>%
  tidyr::pivot_longer(-group, names_to="feature_id", values_to="group_mean") %>%
  tidyr::pivot_wider(names_from=group, values_from=group_mean, names_prefix="mean_")

report <- res_ann %>%
  left_join(means, by="feature_id") %>%
  mutate(direction = ifelse(logFC > 0, "higher_in_Male", "higher_in_Female"),
         mean_diff = mean_Male_E4 - mean_Female_E4) %>%
  arrange(FDR, desc(abs(logFC)))
readr::write_csv(report, file.path(out_dir, "08_report_log2_with_group_means_MvF.csv"))
