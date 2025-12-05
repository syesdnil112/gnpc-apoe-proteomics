# scripts/07_reactome_brain_MvF.R
# Brain-related Reactome enrichment for Male(ε4, AD, MMSE<24) vs Female(ε4, AD, MMSE<24)

suppressPackageStartupMessages({ library(dplyr); library(readr); library(ggplot2); library(limma) })
out_dir <- get0("out_dir", ifnotfound = "qc_outputs")

# ---- Inputs ----
res_path      <- file.path(out_dir, "06B_limma_log2_Male_minus_Female_annotated.csv")
analyte_path  <- "SomalogicAnalyteInfoV1_anonymized.csv"
reactome_gmt  <- list.files(file.path(out_dir, "gene_sets"), pattern = "Reactome.*\\.gmt$", full.names = TRUE)[1]

stopifnot(file.exists(res_path))
stopifnot(file.exists(analyte_path))
stopifnot(!is.na(reactome_gmt) && file.exists(reactome_gmt))

# ---- Load results (Male − Female) ----
res <- readr::read_csv(res_path, show_col_types = FALSE)

# ---- Read analyte info for ID mapping ----
an <- readr::read_csv(analyte_path, show_col_types = FALSE)
# Require these exact columns; map if different before running.
req_ann <- c("column_name","target","entrez_gene_symbol","uni_prot")
miss <- setdiff(req_ann, names(an))
if (length(miss)) stop("Analyte info missing columns: ", paste(miss, collapse = ", "))

sym_map <- an %>%
  transmute(
    uni = as.character(uni_prot),
    sym = as.character(entrez_gene_symbol)
  ) %>%
  filter(!is.na(uni), uni != "", !is.na(sym), sym != "") %>%
  distinct()

# ---- Read GMT ----
read_gmt <- function(path){
  con <- file(path,"r"); on.exit(close(con))
  sets <- list(); nm <- character(); i <- 0L
  while(length(ln <- readLines(con, n=1L, warn=FALSE)) > 0L){
    p <- strsplit(ln, "\t")[[1]]
    if (length(p) >= 3) { i <- i + 1L; nm[i] <- p[1]; sets[[i]] <- unique(p[-c(1,2)]) }
  }
  names(sets) <- nm; sets
}
rx_sets_raw <- read_gmt(reactome_gmt)

# ---- Detect ID type and map to SYMBOLs if UniProt-like ----
is_uniprot_like <- function(x) grepl("^[A-NR-Z][0-9][A-Z0-9]{3}[0-9](-\\d+)?$", x)
id_tokens <- unlist(rx_sets_raw, use.names = FALSE)
frac_uniprot <- mean(is_uniprot_like(head(id_tokens, 5000)), na.rm = TRUE)

to_symbols <- function(vec){
  if (isTRUE(frac_uniprot > 0.5)) {
    df <- tibble(uni = vec) %>% left_join(sym_map, by = "uni")
    out <- df$sym
    out[is.na(out) | out == ""] <- df$uni[is.na(out) | out == ""]
    unique(out)
  } else {
    unique(vec)
  }
}
rx_sets <- lapply(rx_sets_raw, to_symbols)

# ---- Build ranked gene-level statistics (t) ----
gene_df <- res %>%
  filter(!is.na(entrez_gene_symbol), entrez_gene_symbol != "", is.finite(t)) %>%
  group_by(entrez_gene_symbol) %>%
  slice_max(order_by = abs(t), n = 1, with_ties = FALSE) %>%
  ungroup()

stat <- gene_df$t
names(stat) <- gene_df$entrez_gene_symbol
universe <- names(stat)

# ---- Brain-like set filter ----
brain_keywords <- c(
  "nervous system","neuro","neuron","synap","axon","dendrit","glia","astrocy","microglia",
  "oligodendro","myelin","vesicle","neurotrans","dopamine","glutamat","gaba","cholin",
  "seroton","noradrenerg","neurotroph","brain","hippocamp","cortex","cerebell","tau","amyloid"
)
rx_name <- names(rx_sets)
brain_idx <- grepl(paste(brain_keywords, collapse="|"), rx_name, ignore.case = TRUE)
rx_sets_brain <- rx_sets[brain_idx]
message("Brain-like Reactome pathways kept: ", length(rx_sets_brain), " / ", length(rx_sets))

# ---- Build index on universe; drop small sets (<10) ----
index <- lapply(rx_sets_brain, function(gs) which(universe %in% gs))
keep <- vapply(index, function(ix) length(ix) >= 10, logical(1))
index <- index[keep]
rx_sets_brain <- rx_sets_brain[keep]

# ---- cameraPR (pre-ranked t) ----
cam_b <- limma::cameraPR(statistic = stat, index = index, inter.gene.cor = 0.01, use.ranks = FALSE)
cam_b$Set <- names(index)
cam_b$Direction <- vapply(index, function(ix) if (mean(stat[ix], na.rm = TRUE) > 0) "Male_up" else "Female_up", character(1))
cam_b <- cam_b %>% relocate(Set) %>% arrange(FDR, PValue)
readr::write_csv(cam_b, file.path(out_dir, "07_reactome_cameraPR_brain_only_MvF.csv"))

# ---- Plot top 20 ----
cam_b$negLog10FDR <- -log10(pmax(cam_b$FDR, 1e-300))
g_cam_b <- ggplot(head(cam_b, 20),
                  aes(negLog10FDR, reorder(Set, negLog10FDR), size = NGenes, color = Direction)) +
  geom_point() +
  labs(title = "Reactome — brain-related pathways (cameraPR) M − F",
       subtitle = "+ = Male_up (t > 0)",
       x = expression(-log[10]("FDR")), y = NULL) +
  theme_bw()
ggsave(file.path(out_dir, "07_reactome_cameraPR_brain_only_MvF_top20.png"), g_cam_b, width = 15, height = 6, dpi = 150)

# ---- ORA for up-genes each side (brain-only sets) ----
sig_genes <- res %>% filter(is.finite(FDR), FDR < 0.05, !is.na(entrez_gene_symbol), entrez_gene_symbol != "")
sig_set <- unique(sig_genes$entrez_gene_symbol)

sig_up_Male   <- names(stat)[stat > 0 & names(stat) %in% sig_set]
sig_up_Female <- names(stat)[stat < 0 & names(stat) %in% sig_set]

do_ora <- function(sig, uni, sets){
  resl <- lapply(names(sets), function(s){
    gs <- sets[[s]]
    in_set <- uni %in% gs
    A <- sum(in_set & uni %in% sig); B <- sum(in_set & !(uni %in% sig))
    C <- sum(!in_set & uni %in% sig); D <- sum(!in_set & !(uni %in% sig))
    ft <- fisher.test(matrix(c(A,B,C,D), 2), alternative = "greater")
    c(Set = s, Overlap = A, SetSize = A + B, PValue = ft$p.value, OddsRatio = unname(ft$estimate))
  })
  df <- as.data.frame(do.call(rbind, resl), stringsAsFactors = FALSE)
  df$Overlap  <- as.numeric(df$Overlap)
  df$SetSize  <- as.numeric(df$SetSize)
  df$PValue   <- as.numeric(df$PValue)
  df$OddsRatio<- as.numeric(df$OddsRatio)
  df$FDR <- p.adjust(df$PValue, "BH")
  df %>% arrange(FDR, PValue)
}

ora_Male   <- do_ora(sig_up_Male,   universe, rx_sets_brain)
ora_Female <- do_ora(sig_up_Female, universe, rx_sets_brain)
readr::write_csv(ora_Male,   file.path(out_dir, "07_reactome_ORA_brain_only_up_in_Male.csv"))
readr::write_csv(ora_Female, file.path(out_dir, "07_reactome_ORA_brain_only_up_in_Female.csv"))
