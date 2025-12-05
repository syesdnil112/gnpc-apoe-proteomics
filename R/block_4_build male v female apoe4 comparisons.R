# scripts/05_AD_APOE4_no24_MvF_dataset.R
# Purpose: AD (==1), MMSE<24, APOE ε4 carriers (exclude 2/4), pd==0 & ftd==0 & als==0
#          Compare Male vs Female. 

suppressPackageStartupMessages({ library(dplyr); library(readr); library(stringr) })
if (!exists("out_dir")) out_dir <- "qc_outputs"

# ------------------------ USER INPUTS (edit these) ------------------------
# Provide the exact column names as they appear in CSVs.
sample_id_col_user <- "sample_id"
ad_col_user        <- "ad"
sex_col_user       <- "sex"
apoe_col_user      <- "apoe"                  
mmse_col_user      <- "cognitive_test_score" 
pd_col_user        <- "pd"
ftd_col_user       <- "ftd"
als_col_user       <- "als"

# Sex codes in your data (integers/strings after coercion)
sex_male_code   <- 1
sex_female_code <- 2
# -------------------------------------------------------------------------

zs_path    <- file.path(out_dir, "04_zscores_by_subset.csv")
pca_path   <- file.path(out_dir, "04_pca_scores.csv")
pheno_path <- "clinicalV1_filtered_ad_apoe.csv"
stopifnot(file.exists(zs_path), file.exists(pheno_path))

# --- Helpers (name cleaning consistent with prior scripts)
clean_names_local <- function(x){
  x <- gsub("[^A-Za-z0-9]+","_",x); x <- gsub("_+","_",x)
  x <- gsub("^_|_$","",x); tolower(x)
}
.clean_id <- function(x) tolower(trimws(as.character(x)))
norm_apoe <- function(x){
  x <- tolower(as.character(x))
  x <- gsub("[^0-9/]", "", x)                  # keep digits and slash
  x <- gsub("(\\d)[\\s]*[-_:.,][\\s]*(\\d)", "\\1/\\2", x)
  x
}
has_allele <- function(gt, a) grepl(paste0("(^|/)", a, "(/|$)"), gt)

# Normalize user-provided names to match cleaned columns
norm_name <- function(nm) clean_names_local(nm)

# --- Load
zs <- readr::read_csv(zs_path, show_col_types = FALSE)
ph <- readr::read_csv(pheno_path, show_col_types = FALSE)
names(zs) <- clean_names_local(names(zs))
names(ph) <- clean_names_local(names(ph))

# Map user names → cleaned names
sample_id_col <- norm_name(sample_id_col_user)
ad_col        <- norm_name(ad_col_user)
sex_col       <- norm_name(sex_col_user)
apoe_col      <- norm_name(apoe_col_user)
mmse_col      <- norm_name(mmse_col_user)
pd_col        <- norm_name(pd_col_user)
ftd_col       <- norm_name(ftd_col_user)
als_col       <- norm_name(als_col_user)

# Validate presence 
needed <- c(sample_id_col, ad_col, sex_col, apoe_col, mmse_col, pd_col, ftd_col, als_col)
missing <- setdiff(needed, names(ph))
if (length(missing)) stop("Missing phenotype columns: ", paste(missing, collapse = ", "))

# Ensure sample_id exists in zs or discover a compatible one
if (!sample_id_col %in% names(zs)) {
  stop("Expression table is missing '", sample_id_col, "'. Provide the correct sample ID column name.")
}

# Harmonize IDs
ph[[sample_id_col]] <- .clean_id(ph[[sample_id_col]])
zs[[sample_id_col]] <- .clean_id(zs[[sample_id_col]])

# Merge & compute
dat <- zs %>%
  left_join(
    ph %>% select(all_of(c(sample_id_col, ad_col, sex_col, apoe_col, mmse_col, pd_col, ftd_col, als_col))),
    by = sample_id_col
  ) %>%
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
  # Filters: AD, sex in {M,F}, MMSE<24, APOE ε4 carriers excluding any ε2 (removes 2/4),
  #          and no comorbid neurodegenerative disease
  filter(
    ad_flag,
    (sex_m | sex_f),
    !is.na(mmse_val) & mmse_val < 24,
    !is.na(apoe_std) & has_e4 & !has_e2,
    !is.na(pd_val)  & pd_val  == 0,
    !is.na(ftd_val) & ftd_val == 0,
    !is.na(als_val) & als_val == 0
  ) %>%
  mutate(
    group = ifelse(sex_m, "Male_E4", "Female_E4"),
    group = factor(group, levels = c("Female_E4","Male_E4"))
  )

# Optional PCs
if (file.exists(pca_path)) {
  pcs <- readr::read_csv(pca_path, show_col_types = FALSE)
  names(pcs) <- clean_names_local(names(pcs))
  if (!(sample_id_col %in% names(pcs))) {
    stop("PCA file missing '", sample_id_col, "'. Align PCA 'sample_id' column name.")
  }
  pcs <- pcs %>% select(all_of(c(sample_id_col, grep("^pc\\d+$|^pc", names(pcs), value = TRUE))))
  dat <- dat %>% left_join(pcs, by = sample_id_col)
}

# Meta to front
front_meta <- c(sample_id_col, "group", mmse_col, pd_col, ftd_col, als_col, ".__subset", "contributor_code", "sample_type")
front_meta <- intersect(front_meta, names(dat))
dat <- dat %>% relocate(all_of(front_meta))

# Write
out_file <- file.path(out_dir, "05_AD_APOE4_no24_MvF_dataset.csv")
readr::write_csv(dat, out_file)
message("Wrote: ", out_file, " | N=", nrow(dat),
        " (Female_E4=", sum(dat$group=="Female_E4"), ", Male_E4=", sum(dat$group=="Male_E4"), ")")
