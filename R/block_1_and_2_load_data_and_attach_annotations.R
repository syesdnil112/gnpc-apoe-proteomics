##############################################
######## 1) Load plates & merge columns ######
##############################################

data_dir <- "unmodified_plates"  # Somalogic01V1_anonymized_filtered_ad.csv ... Somalogic12V1_anonymized_filtered_ad.csv
plate_files <- file.path(data_dir, sprintf("Somalogic%02dV1_anonymized_filtered_ad.csv", 1:12))
stopifnot(all(file.exists(plate_files)))

suppressPackageStartupMessages({ library(data.table); library(dplyr); library(readr) })

.clean_id <- function(x) tolower(trimws(as.character(x)))

clean_names_local <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)   # non-alnum -> "_"
  x <- gsub("_+", "_", x)              # collapse multiple "_"
  x <- gsub("^_|_$", "", x)            # trim leading/trailing "_"
  tolower(x)
}

read_plate <- function(f) {
  dt <- data.table::fread(f, check.names = FALSE, showProgress = FALSE)
  
  # ensure first column is sample_id
  if (!"sample_id" %in% names(dt)) names(dt)[1] <- "sample_id"
  
  # standardize column names (e.g., "contributor code" -> "contributor_code")
  names(dt) <- clean_names_local(names(dt))
  
  dt$sample_id <- .clean_id(dt$sample_id)
  
  # Identify metadata & analyte columns
  meta_cols <- intersect(c("contributor_code","visit","sample_type"), names(dt))
  analyte_cols <- setdiff(names(dt), c("sample_id", meta_cols))
  
  # Make analytes numeric, then convert sentinel -1 -> NA **only** on analytes
  dt <- dt %>%
    mutate(across(all_of(analyte_cols), ~ suppressWarnings(as.numeric(.)))) %>%
    mutate(across(all_of(analyte_cols), ~ dplyr::na_if(., -1)))
  
  list(
    meta = dt %>% select(sample_id, all_of(meta_cols)),
    expr = dt %>% select(sample_id, all_of(analyte_cols))
  )
}

plates <- lapply(plate_files, read_plate)

# Common sample IDs across all plates
common_ids <- Reduce(intersect, lapply(plates, function(x) x$expr$sample_id))
if (!length(common_ids)) stop("No overlapping sample IDs across Somalogic plates.")

# Align & column-bind analytes
expr_list <- lapply(plates, function(x) x$expr %>% filter(sample_id %in% common_ids) %>% arrange(sample_id))
merged <- expr_list[[1]]
for (k in 2:length(expr_list)) {
  merged <- merged %>% left_join(expr_list[[k]], by = "sample_id")
}

# Attach **one** set of metadata (from the first plate) aligned to common_ids
meta_df <- plates[[1]]$meta %>% filter(sample_id %in% common_ids) %>% arrange(sample_id)
merged <- meta_df %>% left_join(merged, by = "sample_id")

# De-duplicate sample rows and analyte columns (safeguards)
merged <- merged %>% distinct(sample_id, .keep_all = TRUE)
merged <- merged %>% select(sample_id, contributor_code, visit, sample_type,
                            any_of(unique(setdiff(names(merged),
                                                  c("sample_id","contributor_code","visit","sample_type")))))

# where to save outputs
out_dir <- "qc_outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
readr::write_csv(merged, file.path(out_dir, "00_merged_raw.csv"))
dim(merged)

##############################################
#### 2) Attach annotations & keep proteins ####
##############################################

if (!exists("out_dir")) out_dir <- "qc_outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
merged_path <- file.path(out_dir, "00_merged_raw.csv")
anno_path   <- "SomalogicAnalyteInfoV1_anonymized.csv"

stopifnot(file.exists(merged_path), file.exists(anno_path))

suppressPackageStartupMessages({ library(data.table); library(dplyr); library(readr) })

clean_names_local <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x); x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x); tolower(x)
}

# 2.1 Load merged matrix (meta + analytes) and clean names
merged <- readr::read_csv(merged_path, show_col_types = FALSE)
names(merged) <- clean_names_local(names(merged))
meta_cols    <- intersect(c("sample_id","contributor_code","visit","sample_type"), names(merged))
analyte_cols <- setdiff(names(merged), meta_cols)

# 2.2 Load & standardize annotations
analyte_info<- data.table::fread(anno_path, check.names = FALSE, showProgress = FALSE) %>% as.data.frame()
names(analyte_info) <- clean_names_local(names(analyte_info))
stopifnot("column_name" %in% names(analyte_info))   # SOMAmer/SomaScan key

# 2.3 Filter to analytes present in the data and human proteins
anno_present <- analyte_info %>%
  filter(column_name %in% analyte_cols) %>%
  { if ("organism" %in% names(.)) filter(., tolower(organism) == "human") else . } %>%
  { if ("type" %in% names(.)) filter(., tolower(type) == "protein") else . }

kept_ids <- unique(anno_present$column_name)
message("Analytes kept after annotation filtering: ", length(kept_ids)) 
#7289 analytes kept after annotation filtering

# IMPORTANT: keep expression columns as SOMAmer IDs (column_name). Do NOT rename here.
expr_mat <- merged %>% select(all_of(meta_cols), all_of(kept_ids))

# Save matrix + the filtered annotation table
readr::write_csv(expr_mat,   file.path(out_dir, "01_expr_annotated_filtered.csv"))
readr::write_csv(anno_present, file.path(out_dir, "01_analyte_info_postfilter.csv"))

