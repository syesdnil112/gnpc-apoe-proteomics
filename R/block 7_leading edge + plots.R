############################################################
## Drill into brain Reactome results: leading edge + plots
############################################################
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(tidyr); library(stringr); library(limma)
})

out_dir <- get0("out_dir", ifnotfound = "qc_outputs")

# --- Load annotated DE and log2 dataset (for sample-level plots) -------------
res_ann_path <- file.path(out_dir, "06B_limma_log2_E3E3_vs_E4plus_annotated.csv")
dat_log2_path <- file.path(out_dir, "05B_AD_female_APOE_log2.csv")
stopifnot(file.exists(res_ann_path), file.exists(dat_log2_path))
res <- readr::read_csv(res_ann_path, show_col_types = FALSE)
dat_log2 <- readr::read_csv(dat_log2_path, show_col_types = FALSE)

# --- Recreate gene-level ranked stats used by cameraPR -----------------------
gene_df <- res %>%
  filter(!is.na(entrez_gene_symbol), entrez_gene_symbol != "", is.finite(t)) %>%
  group_by(entrez_gene_symbol) %>%               # 1 aptamer per gene
  slice_max(order_by = abs(t), n = 1, with_ties = FALSE) %>%
  ungroup()

stat <- gene_df$t
names(stat) <- gene_df$entrez_gene_symbol
universe <- names(stat)

# --- Recreate brain-filtered Reactome sets + cameraPR if not present ---------
# Expect "rx_sets_brain" and "cam_b" from the earlier step; rebuild if absent.
need_build <- !exists("rx_sets_brain") || !exists("cam_b") || is.null(cam_b) || !nrow(cam_b)

if (need_build) {
  # Locate the Reactome GMT you downloaded earlier
  reactome_gmt <- list.files(file.path(out_dir, "gene_sets"), pattern = "Reactome.*\\.gmt$", full.names = TRUE)[1]
  stopifnot(file.exists(reactome_gmt))
  
  # Read GMT
  read_gmt <- function(path){
    con <- file(path,"r"); on.exit(close(con)); sets <- list(); nm <- character(); i <- 0L
    while(length(ln <- readLines(con, n=1L, warn=FALSE))>0L){
      p <- strsplit(ln, "\t")[[1]]; if(length(p) >= 3){ i <- i+1L; nm[i] <- p[1]; sets[[i]] <- unique(p[-c(1,2)]) }
    }
    names(sets) <- nm; sets
  }
  rx_sets_raw <- read_gmt(reactome_gmt)
  
  # Map UniProt -> SYMBOL if needed, using your analyte_info
  stopifnot(exists("analyte_info"))
  sym_map <- analyte_info %>% transmute(uni = as.character(uni_prot),
                                        sym = as.character(entrez_gene_symbol)) %>%
    filter(!is.na(uni), uni!="", !is.na(sym), sym!="") %>% distinct()
  is_uniprot_like <- function(x) grepl("^[A-NR-Z][0-9][A-Z0-9]{3}[0-9](-\\d+)?$", x)
  id_tokens <- unlist(rx_sets_raw, use.names = FALSE)
  frac_uniprot <- mean(is_uniprot_like(head(id_tokens, 5000)), na.rm = TRUE)
  
  to_symbols <- function(vec){
    if (frac_uniprot > 0.5) {
      df <- tibble(uni = vec) %>% left_join(sym_map, by = "uni")
      out <- df$sym; out[is.na(out) | out==""] <- df$uni[is.na(out) | out==""]
      unique(out)
    } else unique(vec)
  }
  rx_sets <- lapply(rx_sets_raw, to_symbols)
  
  # Brain name filter (same keywords as before)
  brain_keywords <- c(
    "nervous system","neuro","neuron","synap","axon","dendrit","glia","astrocy","microglia",
    "oligodendro","myelin","vesicle","neurotrans","dopamine","glutamat","gaba","cholin",
    "seroton","noradrenerg","neurotroph","brain","hippocamp","cortex","cerebell","tau","amyloid"
  )
  rx_name <- names(rx_sets)
  brain_idx <- grepl(paste(brain_keywords, collapse="|"), rx_name, ignore.case = TRUE)
  rx_sets_brain <- rx_sets[brain_idx]
  
  # Build index on universe and keep sets with >=10 genes in universe
  index <- lapply(rx_sets_brain, function(gs) which(universe %in% gs))
  keep <- vapply(index, function(ix) length(ix) >= 10, logical(1))
  index <- index[keep]; rx_sets_brain <- rx_sets_brain[keep]
  
  # cameraPR
  cam_b <- limma::cameraPR(statistic = stat, index = index, inter.gene.cor = 0.01, use.ranks = FALSE)
  cam_b$Set <- names(index)
  cam_b$Direction <- vapply(index, function(ix) if (mean(stat[ix], na.rm=TRUE) > 0) "E3E3_up" else "E4carriers_up", character(1))
  cam_b <- cam_b %>% relocate(Set) %>% arrange(FDR, PValue)
}

# ---- Leading-edge genes for the top Reactome brain pathways (FIXED) ---------
suppressPackageStartupMessages({ library(dplyr); library(readr) })

stopifnot(exists("cam_b"), exists("rx_sets_brain"), exists("stat"), exists("universe"))

top_n    <- min(10, nrow(cam_b))
top_sets <- cam_b$Set[seq_len(top_n)]

# Rebuild index in case environment changed
index <- lapply(rx_sets_brain, function(gs) which(universe %in% gs))

leading_edge <- lapply(top_sets, function(s) {
  ix <- index[[s]]
  if (is.null(ix) || length(ix) == 0) return(NULL)
  genes  <- universe[ix]
  dirpos <- mean(stat[ix], na.rm = TRUE) > 0
  ord    <- if (dirpos) order(-stat[genes]) else order(stat[genes])
  data.frame(
    Set       = s,
    SYMBOL    = genes[ord],
    t         = unname(stat[genes][ord]),
    Direction = if (dirpos) "E3E3_up" else "E4carriers_up",
    stringsAsFactors = FALSE
  )
})

# Base-R safe bind + head-per-set
leading_edge <- leading_edge[!vapply(leading_edge, is.null, logical(1))]
if (length(leading_edge) == 0) stop("No leading-edge sets available.")

le_tbl <- do.call(rbind, lapply(leading_edge, function(df) df[seq_len(min(30, nrow(df))), , drop = FALSE]))

out_dir <- get0("out_dir", ifnotfound = "qc_outputs")
readr::write_csv(le_tbl, file.path(out_dir, "08_reactome_brain_leading_edge_top10sets.csv"))
message("Wrote leading-edge table: ", file.path(out_dir, "08_reactome_brain_leading_edge_top10sets.csv"))


# --- 2) Sample-level plots for leading-edge genes (violin/box, faceted) -----
# Map each SYMBOL -> the aptamer (feature_id) we used in the gene collapse
map_gene_to_feature <- res %>%
  filter(entrez_gene_symbol %in% unique(le_tbl$SYMBOL)) %>%
  group_by(entrez_gene_symbol) %>%
  slice_max(order_by = abs(t), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(entrez_gene_symbol, feature_id, target)

# Pull expression for those features from the log2 dataset
feat_cols <- unique(map_gene_to_feature$feature_id)
long_expr <- dat_log2 %>%
  select(sample_id, apoe_group, all_of(feat_cols)) %>%
  pivot_longer(-c(sample_id, apoe_group), names_to = "feature_id", values_to = "log2_RFU") %>%
  left_join(map_gene_to_feature, by = "feature_id")

# Faceted violins for the top 12 leading-edge genes across top sets
genes_to_plot <- unique(le_tbl$SYMBOL)[1:min(12, length(unique(le_tbl$SYMBOL)))]
plot_df <- long_expr %>% filter(entrez_gene_symbol %in% genes_to_plot)

ggplot(plot_df, aes(apoe_group, log2_RFU, fill = apoe_group)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.size = 0.4, alpha = 0.9) +
  facet_wrap(~ entrez_gene_symbol, scales = "free_y") +
  labs(title = "Leading-edge genes from brain Reactome sets",
       subtitle = "log2 RFU by group; E3/3 vs E4 carriers",
       x = NULL, y = "log2 RFU") +
  theme_bw() + theme(legend.position = "none") ->
  g_violin
ggsave(file.path(out_dir, "08_leading_edge_violin_top12.png"), g_violin, width = 12, height = 8, dpi = 150)

# --- 3) Complementary ORA on the same brain-only sets ------------------------
sig_up_E3 <- gene_df %>% filter(FDR < 0.05, logFC > 0) %>% pull(entrez_gene_symbol) %>% unique()
sig_up_E4 <- gene_df %>% filter(FDR < 0.05, logFC < 0) %>% pull(entrez_gene_symbol) %>% unique()

do_ora <- function(sig, uni, sets){
  resl <- lapply(names(sets), function(s){
    gs <- sets[[s]]; in_set <- uni %in% gs
    A <- sum(in_set & uni %in% sig); B <- sum(in_set & !(uni %in% sig))
    C <- sum(!in_set & uni %in% sig); D <- sum(!in_set & !(uni %in% sig))
    ft <- fisher.test(matrix(c(A,B,C,D), 2), alternative = "greater")
    c(Set=s, Overlap=A, SetSize=A+B, PValue=ft$p.value, OddsRatio=unname(ft$estimate))
  })
  df <- as.data.frame(do.call(rbind, resl), stringsAsFactors = FALSE)
  df$Overlap <- as.numeric(df$Overlap); df$SetSize <- as.numeric(df$SetSize)
  df$PValue <- as.numeric(df$PValue);   df$OddsRatio <- as.numeric(df$OddsRatio)
  df$FDR <- p.adjust(df$PValue, "BH"); df %>% arrange(FDR, PValue)
}
ora_E3_b <- do_ora(sig_up_E3, universe, rx_sets_brain)
ora_E4_b <- do_ora(sig_up_E4, universe, rx_sets_brain)
readr::write_csv(ora_E3_b, file.path(out_dir, "08_reactome_ORA_brain_only_up_in_E3E3.csv"))
readr::write_csv(ora_E4_b, file.path(out_dir, "08_reactome_ORA_brain_only_up_in_E4carriers.csv"))

message("Wrote:",
        "\n • 08_reactome_brain_leading_edge_top10sets.csv",
        "\n • 08_leading_edge_violin_top12.png",
        "\n • 08_reactome_ORA_brain_only_up_in_E3E3.csv",
        "\n • 08_reactome_ORA_brain_only_up_in_E4carriers.csv")


# 1) Summarize the 12 plotted genes with DE stats + group means
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr) })

out_dir <- get0("out_dir", ifnotfound = "qc_outputs")
res_ann <- readr::read_csv(file.path(out_dir, "06B_limma_log2_E3E3_vs_E4plus_annotated.csv"),
                           show_col_types = FALSE)
dat_log2 <- readr::read_csv(file.path(out_dir, "05B_AD_female_APOE_log2.csv"),
                            show_col_types = FALSE)

genes_shown <- unique(res_ann$entrez_gene_symbol[match(
  c("CASK","DLG3","DLG4","ERBB4","GRIA4","HRAS","KRAS","LIN7A","MAPT","NEFL","PRKAR1A","RPS6KA6"),
  res_ann$entrez_gene_symbol
)])
genes_shown <- genes_shown[!is.na(genes_shown)]

# one aptamer per gene (the one used by enrichment)
per_gene <- res_ann %>%
  filter(entrez_gene_symbol %in% genes_shown) %>%
  group_by(entrez_gene_symbol) %>%
  slice_max(order_by = abs(t), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(entrez_gene_symbol, feature_id, target, logFC, AveExpr, t, P.Value, FDR)

# group means on the same log2 scale
feat_cols <- per_gene$feature_id
means <- dat_log2 %>%
  select(apoe_group, all_of(feat_cols)) %>%
  group_by(apoe_group) %>%
  summarise(across(all_of(feat_cols), ~ mean(as.numeric(.), na.rm = TRUE)), .groups="drop") %>%
  pivot_longer(-apoe_group, names_to="feature_id", values_to="group_mean") %>%
  pivot_wider(names_from = apoe_group, values_from = group_mean,
              names_prefix = "mean_")

tbl <- per_gene %>%
  left_join(means, by="feature_id") %>%
  mutate(direction = ifelse(logFC > 0, "higher_in_E3E3", "higher_in_E4"),
         abs_logFC = abs(logFC)) %>%
  arrange(FDR, desc(abs_logFC))

readr::write_csv(tbl, file.path(out_dir, "09_leading_edge_summary.csv"))
tbl %>% select(entrez_gene_symbol, target, logFC, FDR, mean_E3E3, mean_E4_plus, direction)
