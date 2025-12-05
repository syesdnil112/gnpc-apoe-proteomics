# ==== Brain-only Reactome enrichment (name-filtered) =========================
suppressPackageStartupMessages({ library(dplyr); library(readr); library(ggplot2); library(limma) })
out_dir <- get0("out_dir", ifnotfound = "qc_outputs")

# Expect these from the previous Reactome code; rebuild if absent:
if (!exists("rx_sets")) {
  # Recreate from your saved annotated DE and downloaded GMT
  res <- readr::read_csv(file.path(out_dir, "06B_limma_log2_E3E3_vs_E4plus_annotated.csv"), show_col_types = FALSE)
  reactome_gmt <- list.files(file.path(out_dir, "gene_sets"), pattern = "Reactome.*\\.gmt$", full.names = TRUE)[1]
  stopifnot(file.exists(reactome_gmt))
  read_gmt <- function(path){
    con <- file(path,"r"); on.exit(close(con)); sets <- list(); nm <- character(); i <- 0L
    while(length(ln <- readLines(con, n=1L, warn=FALSE))>0L){
      p <- strsplit(ln, "\t")[[1]]; if(length(p) >= 3){ i <- i+1L; nm[i] <- p[1]; sets[[i]] <- unique(p[-c(1,2)]) }
    }
    names(sets) <- nm; sets
  }
  rx_sets_raw <- read_gmt(reactome_gmt)
  
  # Map to SYMBOLs if GMT is UniProt (via analyte_info)
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
  
  # Build ranked stats once
  gene_df <- res %>% filter(!is.na(entrez_gene_symbol), entrez_gene_symbol!="", is.finite(t)) %>%
    group_by(entrez_gene_symbol) %>% slice_max(order_by = abs(t), n = 1, with_ties = FALSE) %>% ungroup()
  stat <- gene_df$t; names(stat) <- gene_df$entrez_gene_symbol
  universe <- names(stat)
}

# ---- filter set names to "brain-like" pathways
brain_keywords <- c(
  "nervous system","neuro","neuron","synap","axon","dendrit","glia","astrocy","microglia",
  "oligodendro","myelin","vesicle","neurotrans","dopamine","glutamat","gaba","cholin",
  "seroton","noradrenerg","neurotroph","brain","hippocamp","cortex","cerebell","tau","amyloid"
)
rx_name <- names(rx_sets)
brain_idx <- grepl(paste(brain_keywords, collapse="|"), rx_name, ignore.case = TRUE)
rx_sets_brain <- rx_sets[brain_idx]
message("Brain-like Reactome pathways kept: ", length(rx_sets_brain), " / ", length(rx_sets))

# Build index on your universe and drop too-small sets
index <- lapply(rx_sets_brain, function(gs) which(universe %in% gs))
keep <- vapply(index, function(ix) length(ix) >= 10, logical(1))
index <- index[keep]; rx_sets_brain <- rx_sets_brain[keep]

# cameraPR (pre-ranked t)
cam_b <- limma::cameraPR(statistic = stat, index = index, inter.gene.cor = 0.01, use.ranks = FALSE)
cam_b$Set <- names(index)
cam_b$Direction <- vapply(index, function(ix) if (mean(stat[ix], na.rm=TRUE) > 0) "E3E3_up" else "E4carriers_up", character(1))
cam_b <- cam_b %>% relocate(Set) %>% arrange(FDR, PValue)
readr::write_csv(cam_b, file.path(out_dir, "07_reactome_cameraPR_brain_only.csv"))

# Plot
cam_b$negLog10FDR <- -log10(pmax(cam_b$FDR, 1e-300))
ggplot(head(cam_b, 20),
       aes(negLog10FDR, reorder(Set, negLog10FDR), size = NGenes, color = Direction)) +
  geom_point() + labs(title = "Reactome — brain-related pathways (cameraPR)",
                      x = expression(-log[10]("FDR")), y = NULL) + theme_bw() ->
  g_cam_b
ggsave(file.path(out_dir, "07_reactome_cameraPR_brain_only_top20.png"), g_cam_b, width = 15, height = 6, dpi = 150)

# ORA for up-genes each side (using same brain-only sets)
sig_up_E3 <- names(stat)[stat > 0 & names(stat) %in% (res$entrez_gene_symbol[res$FDR < 0.05])]
sig_up_E4 <- names(stat)[stat < 0 & names(stat) %in% (res$entrez_gene_symbol[res$FDR < 0.05])]
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
readr::write_csv(ora_E3_b, file.path(out_dir, "07_reactome_ORA_brain_only_up_in_E3E3.csv"))
readr::write_csv(ora_E4_b, file.path(out_dir, "07_reactome_ORA_brain_only_up_in_E4carriers.csv"))
