# =============================================================================
# 06_DEG_SNP_association.R
#
# Purpose: Test associations between top 25 DEGs and 15 SNPs using
#          linear models adjusted for confounders.
#          Also runs unnamed gene lookup via Ensembl/biomaRt.
#
# Note:    Top 25 DEGs are predominantly non-coding RNA and pseudogenes.
#          Heatmap removed as uninformative. Results saved as CSV only.
#
# Input:   rds_objects/logCPM.rds        (from 02_transcriptomics.R)
#          rds_objects/gene_clinical.rds  (from 02_transcriptomics.R)
#          rds_objects/deg_results.rds    (from 02_transcriptomics.R)
#
# Output:  results/top25_DEG_SNP_associations.csv
#          results/top25_DEG_SNP_significant_hits.csv
#          results/unnamed_genes_lookup.csv
# =============================================================================

library(limma)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
library(tidyverse)
library(ggrepel)
library(patchwork)
library(RColorBrewer)

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

cat("Libraries loaded\n")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ── Load objects ──────────────────────────────────────────────────────────────

logCPM        <- readRDS("rds_objects/logCPM.rds")
gene_clinical <- readRDS("rds_objects/gene_clinical.rds")
deg_results   <- readRDS("rds_objects/deg_results.rds")

cat("logCPM loaded:",        nrow(logCPM), "genes x", ncol(logCPM), "samples\n")
cat("gene_clinical loaded:", nrow(gene_clinical), "samples\n")
cat("deg_results loaded:",   nrow(deg_results), "genes\n")

snp_cols <- c("GSTM1","GSTT1","GSTP1","NQO1","COMT","MTHFR",
              "CAT1","Glu298Asp","XRCC1","BCMO1","SLC23A1",
              "ZBED3","APOE","HNF1A","TCF7L2")

# ── F1. Extract top 25 DEGs ───────────────────────────────────────────────────

top25 <- deg_results %>%
  rownames_to_column("ENSEMBL") %>%
  arrange(P.Value) %>%
  head(25) %>%
  dplyr::select(ENSEMBL, logFC, P.Value, adj.P.Val)

symbols_25 <- bitr(top25$ENSEMBL, fromType = "ENSEMBL",
                   toType = "SYMBOL", OrgDb = org.Hs.eg.db)

top25 <- top25 %>%
  left_join(symbols_25, by = "ENSEMBL") %>%
  mutate(gene_label = ifelse(is.na(SYMBOL), ENSEMBL, SYMBOL))

cat("\nTop 25 DEGs:\n")
print(top25 %>% dplyr::select(gene_label, logFC, P.Value))

# ── F2. Gene x SNP association loop ──────────────────────────────────────────

snps_found <- snp_cols[snp_cols %in% colnames(gene_clinical)]
cat("SNPs available:", length(snps_found), "/ 15\n")
cat("Running", nrow(top25), "x", length(snps_found), "tests...\n")

association_results <- map_dfr(top25$ENSEMBL, function(gene_id) {
  if (!gene_id %in% rownames(logCPM)) return(NULL)
  expr_vals <- logCPM[gene_id, gene_clinical$gene_col]
  
  map_dfr(snps_found, function(snp) {
    snp_vals <- gene_clinical[[snp]]
    if (all(is.na(snp_vals))) return(NULL)
    
    df_test <- data.frame(
      expr   = as.numeric(expr_vals),
      snp    = as.numeric(as.character(snp_vals)),
      Age_z  = gene_clinical$Age_z,
      Sex    = gene_clinical$Sex,
      BMI_z  = gene_clinical$BMI_z,
      Season = gene_clinical$Season
    ) %>% drop_na()
    
    if (nrow(df_test) < 20) return(NULL)
    
    model <- tryCatch(
      lm(expr ~ snp + Age_z + Sex + BMI_z + Season, data = df_test),
      error = function(e) NULL
    )
    if (is.null(model)) return(NULL)
    coef_s  <- summary(model)$coefficients
    if (!"snp" %in% rownames(coef_s)) return(NULL)
    
    data.frame(
      ENSEMBL = gene_id,
      SNP     = snp,
      beta    = round(coef_s["snp", "Estimate"],   4),
      se      = round(coef_s["snp", "Std. Error"], 4),
      p_value = round(coef_s["snp", "Pr(>|t|)"],  4)
    )
  })
})

association_results <- association_results %>%
  left_join(top25 %>% dplyr::select(ENSEMBL, gene_label, logFC),
            by = "ENSEMBL") %>%
  mutate(
    FDR = round(p.adjust(p_value, method = "BH"), 4),
    sig = case_when(
      FDR     < 0.05 ~ "FDR<0.05",
      p_value < 0.05 ~ "p<0.05",
      p_value < 0.1  ~ "trend",
      TRUE           ~ "NS"
    )
  )

cat("Total tests:", nrow(association_results), "\n")
cat("FDR<0.05:",    sum(association_results$FDR     < 0.05, na.rm = TRUE), "\n")
cat("p<0.05:",      sum(association_results$p_value < 0.05, na.rm = TRUE), "\n")

write.csv(association_results,
          "results/top25_DEG_SNP_associations.csv", row.names = FALSE)

# ── F3. Significant hits ──────────────────────────────────────────────────────
# Heatmap removed — top 25 DEGs are predominantly non-coding RNA and
# pseudogenes, making the heatmap uninformative for thesis presentation.
# Results saved as CSV for reference only.

# Significant hits
sig_hits <- association_results %>%
  dplyr::filter(p_value < 0.05) %>%
  arrange(p_value) %>%
  dplyr::select(gene_label, SNP, beta, se, p_value, FDR, sig)

cat("\nSignificant hits (p<0.05):\n")
if (nrow(sig_hits) > 0) print(sig_hits, row.names = FALSE) else cat("None\n")

write.csv(sig_hits, "results/top25_DEG_SNP_significant_hits.csv",
          row.names = FALSE)
cat("Saved: results/top25_DEG_SNP_significant_hits.csv\n")

# ── G. Unnamed gene lookup ────────────────────────────────────────────────────
# Included here since it operates on top25 from this script

tryCatch({
  mart <- useEnsembl(biomart  = "ensembl",
                     dataset  = "hsapiens_gene_ensembl",
                     mirror   = "useast")
  
  unnamed_ids <- top25 %>%
    dplyr::filter(str_starts(gene_label, "ENSG")) %>%
    pull(ENSEMBL)
  
  if (length(unnamed_ids) > 0) {
    lookup <- getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol",
                     "description", "gene_biotype"),
      filters    = "ensembl_gene_id",
      values     = unnamed_ids,
      mart       = mart
    )
    cat("\nUnnamed gene lookup:\n")
    print(lookup)
    write.csv(lookup, "results/unnamed_genes_lookup.csv", row.names = FALSE)
    cat("Saved: results/unnamed_genes_lookup.csv\n")
  } else {
    cat("No unnamed genes in top 25\n")
  }
}, error = function(e) {
  cat("Ensembl server unavailable — look up manually at ensembl.org\n")
  cat("IDs to look up:\n")
  print(top25 %>%
          dplyr::filter(str_starts(gene_label, "ENSG")) %>%
          pull(ENSEMBL))
})

cat("\n06_DEG_SNP_association.R complete\n")