##############################################
#### 4) AD-female, MMSE<24, APOE groups ######
##############################################

suppressPackageStartupMessages({ library(dplyr); library(readr); library(stringr) })
if (!exists("out_dir")) out_dir <- "qc_outputs"

zs_path  <- file.path(out_dir, "04_zscores_by_subset.csv")
pca_path <- file.path(out_dir, "04_pca_scores.csv")
pheno_path <- "clinicalV1_filtered_ad_apoe.csv"
stopifnot(file.exists(zs_path), file.exists(pheno_path))

clean_names_local <- function(x){
  x <- gsub("[^A-Za-z0-9]+","_",x); x <- gsub("_+","_",x)
  x <- gsub("^_|_$","",x); tolower(x)
}
.clean_id <- function(x) tolower(trimws(as.character(x)))
norm_apoe <- function(x){
  x <- tolower(as.character(x))
  x <- gsub("[^0-9/]", "", x)      # keep digits and slash
  x <- gsub("(\\d)[\\s]*[-_:.,][\\s]*(\\d)", "\\1/\\2", x)
  x
}
has_allele <- function(gt, a) grepl(paste0("(^|/)", a, "(/|$)"), gt)

# Load
zs <- readr::read_csv(zs_path, show_col_types = FALSE)
ph <- readr::read_csv(pheno_path, show_col_types = FALSE)
names(zs) <- clean_names_local(names(zs))
names(ph) <- clean_names_local(names(ph))

# sample_id
if (!"sample_id" %in% names(ph)) {
  cand <- names(ph)[grepl("^sample(_id)?$|somalogic|^sid$|^id$", names(ph))]
  stopifnot(length(cand) > 0); ph <- ph %>% rename(sample_id = !!cand[1])
}
ph$sample_id <- .clean_id(ph$sample_id)
zs$sample_id <- .clean_id(zs$sample_id)

# required phenotype columns (by your dictionary)
stopifnot("ad"  %in% names(ph), "sex" %in% names(ph))
# APOE column
apoe_col <- names(ph)[grepl("apoe", names(ph), ignore.case = TRUE)][1]
stopifnot(length(apoe_col) == 1)

# MMSE column 
mmse_col <- if ("mmse" %in% names(ph)) "mmse" else if ("cognitive_test_score" %in% names(ph)) "cognitive_test_score" else if ("moca" %in% names(ph)) "moca" else NA_character_

# Merge & filter: AD females, MMSE<24, APOE groups (exclude any E2)
dat <- zs %>%
  left_join(ph %>% select(sample_id, ad, sex, all_of(c(apoe_col, mmse_col) %>% na.omit())), by = "sample_id") %>%
  mutate(
    ad_flag  = ad == 1,
    sex_f    = sex == 2,
    mmse_val = if (!is.na(mmse_col)) suppressWarnings(as.numeric(.data[[mmse_col]])) else NA_real_,
    apoe_std = norm_apoe(.data[[apoe_col]]),
    has_e2   = has_allele(apoe_std, "2"),
    has_e4   = has_allele(apoe_std, "4"),
    has_e3   = has_allele(apoe_std, "3")
  ) %>%
  filter(ad_flag, sex_f, !is.na(apoe_std), !has_e2) %>%
  { if (!is.na(mmse_col)) filter(., !is.na(mmse_val) & mmse_val < 24) else . } %>%
  mutate(
    apoe_group = case_when(
      has_e4 ~ "E4_plus",          # includes 3/4 or 4/4 (no E2)
      has_e3 & !has_e4 ~ "E3E3",   # 3/3 (no E2)
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(apoe_group))

# bring PCs if present
if (file.exists(pca_path)) {
  pcs <- readr::read_csv(pca_path, show_col_types = FALSE) %>% select(sample_id, starts_with("PC"))
  names(pcs) <- clean_names_local(names(pcs))
  dat <- dat %>% left_join(pcs, by = "sample_id")
}

# keep meta columns up front
front_meta <- c("sample_id","apoe_group",".__subset","contributor_code","sample_type",
                if (!is.na(mmse_col)) mmse_col else NULL)
front_meta <- intersect(front_meta, names(dat))
dat <- dat %>% relocate(all_of(front_meta))
readr::write_csv(dat, file.path(out_dir, "05_AD_female_APOE_dataset.csv"))


