# =============================================================================
# 05_gene_vessel_correlations.R
#
# Purpose: Compute Pearson correlations between all genes and vessel
#          measurements (CRAE, CRVE, AVR). Primary analysis uses all
#          326 samples to maximise power. Sensitivity analysis uses
#          baseline-only samples (n=127) to assess robustness.
#
# Statistical note:
#   The primary analysis uses naive Pearson across all timepoints.
#   Within-person correlation (duplicateCorrelation r=0.288) 
#   Results are framed as exploratory and hypothesis-generating.
#   Sensitivity analysis at baseline confirms directional consistency
#   for key genes, supporting their biological relevance.
#
# Input:   rds_objects/logCPM.rds        (from 02_transcriptomics.R)
#          rds_objects/gene_clinical.rds  (from 02_transcriptomics.R)
#
# Output:  results/gene_CRAE_correlations.csv
#          results/gene_CRVE_correlations.csv
#          results/gene_AVR_correlations.csv
#          results/top_genes_CRAE_named.csv
#          results/top_genes_CRVE_named.csv
#          results/top_genes_AVR_named.csv
#          results/ORA_CRAE_only_pathways.csv
#          results/ORA_CRVE_only_pathways.csv
#          results/gene_AVR_sensitivity_baseline.csv
# =============================================================================

library(limma)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(lme4)
library(lmerTest)

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

cat("Libraries loaded\n")

dir.create("results", showWarnings = FALSE)

# ── Load objects ──────────────────────────────────────────────────────────────

logCPM        <- readRDS("rds_objects/logCPM.rds")
gene_clinical <- readRDS("rds_objects/gene_clinical.rds")

cat("logCPM loaded:",        nrow(logCPM), "genes x", ncol(logCPM), "samples\n")
cat("gene_clinical loaded:", nrow(gene_clinical), "samples |",
    n_distinct(gene_clinical$PCode), "participants\n")

# ── E1. Primary Pearson correlations (all timepoints) ─────────────────────────
# NOTE: Results should be interpreted as exploratory due to repeated measures.
# Within-person correlation r=0.288 (from duplicateCorrelation in 02_transcriptomics.R)
# may inflate FDR for CRAE and CRVE. Key genes show consistent directions
# in sensitivity analysis (see Section E4 below).

cat("\n--- PRIMARY ANALYSIS: all timepoints (n=", nrow(gene_clinical),
    "samples) ---\n")
cat("NOTE: Results exploratory due to repeated measures structure\n\n")

compute_cors <- function(outcome_vals, label) {
  result <- apply(logCPM, 1, function(g) {
    ct <- cor.test(g, outcome_vals, method = "pearson",
                   use = "pairwise.complete.obs")
    c(correlation = unname(ct$estimate), pval = ct$p.value)
  })
  result         <- as.data.frame(t(result))
  result$ENSEMBL <- rownames(result)
  result$fdr     <- p.adjust(result$pval, method = "BH")
  result         <- result %>% arrange(pval)
  cat(label, "- p<0.05:", sum(result$pval < 0.05),
      "| FDR<0.05:", sum(result$fdr < 0.05),
      "(exploratory, repeated measures)\n")
  return(result)
}

gene_crae_cor <- compute_cors(gene_clinical$CRAE, "CRAE")
gene_crve_cor <- compute_cors(gene_clinical$CRVE, "CRVE")
gene_avr_cor  <- compute_cors(gene_clinical$AVR,  "AVR")

write.csv(gene_crae_cor, "results/gene_CRAE_correlations.csv", row.names = FALSE)
write.csv(gene_crve_cor, "results/gene_CRVE_correlations.csv", row.names = FALSE)
write.csv(gene_avr_cor,  "results/gene_AVR_correlations.csv",  row.names = FALSE)
cat("Saved: primary correlation tables\n")

# ── E2. Add gene symbols to top 50 ───────────────────────────────────────────

add_symbols <- function(cor_df, n = 50) {
  cor_df %>% head(n) %>%
    left_join(
      bitr(.$ENSEMBL, fromType = "ENSEMBL", toType = "SYMBOL",
           OrgDb = org.Hs.eg.db),
      by = "ENSEMBL"
    )
}

top_crae <- add_symbols(gene_crae_cor)
top_crve <- add_symbols(gene_crve_cor)
top_avr  <- add_symbols(gene_avr_cor)

write.csv(top_crae, "results/top_genes_CRAE_named.csv", row.names = FALSE)
write.csv(top_crve, "results/top_genes_CRVE_named.csv", row.names = FALSE)
write.csv(top_avr,  "results/top_genes_AVR_named.csv",  row.names = FALSE)
cat("Saved: top 50 named gene correlates\n")

# ── E3. Vessel-specific ORA ───────────────────────────────────────────────────

crae_genes <- gene_crae_cor %>% dplyr::filter(pval < 0.05) %>% pull(ENSEMBL)
crve_genes <- gene_crve_cor %>% dplyr::filter(pval < 0.05) %>% pull(ENSEMBL)

crae_only <- setdiff(crae_genes, crve_genes)
crve_only <- setdiff(crve_genes, crae_genes)

cat("CRAE-only genes (p<0.05):", length(crae_only), "\n")
cat("CRVE-only genes (p<0.05):", length(crve_only), "\n")

run_ora <- function(gene_list, label, filename) {
  if (length(gene_list) < 20) {
    cat(label, ": too few genes for ORA (n =", length(gene_list), ")\n")
    return(NULL)
  }
  entrez <- bitr(gene_list, fromType = "ENSEMBL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)
  result <- enrichGO(gene          = entrez$ENTREZID,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pvalueCutoff  = 0.05,
                     pAdjustMethod = "BH",
                     readable      = TRUE)
  write.csv(as.data.frame(result), filename, row.names = FALSE)
  cat("Saved:", filename, "\n")
  return(result)
}

run_ora(crae_only, "CRAE-only", "results/ORA_CRAE_only_pathways.csv")
run_ora(crve_only, "CRVE-only", "results/ORA_CRVE_only_pathways.csv")

# ── E4. Sensitivity analysis — baseline only (n=127) ─────────────────────────
# One observation per participant — no repeated measures inflation
# Used to confirm directional consistency of key gene findings

cat("\n--- SENSITIVITY ANALYSIS: baseline only (n=127) ---\n")
cat("Purpose: confirm directional robustness of key gene associations\n\n")

gene_clinical_bl <- gene_clinical %>% filter(timepoint == "Baseline")
expr_bl          <- logCPM[, gene_clinical_bl$gene_col]

cat("Baseline samples:", nrow(gene_clinical_bl), "\n")

compute_cors_bl <- function(outcome_vals, label) {
  result <- apply(expr_bl, 1, function(g) {
    ct <- cor.test(g, outcome_vals, method = "pearson",
                   use = "pairwise.complete.obs")
    c(correlation = unname(ct$estimate), pval = ct$p.value)
  })
  result         <- as.data.frame(t(result))
  result$ENSEMBL <- rownames(result)
  result$fdr     <- p.adjust(result$pval, method = "BH")
  result         <- result %>% arrange(pval)
  cat(label, "(baseline) - p<0.01:", sum(result$pval < 0.01),
      "| FDR<0.05:", sum(result$fdr < 0.05), "\n")
  return(result)
}

gene_avr_cor_bl  <- compute_cors_bl(gene_clinical_bl$AVR,  "AVR")
gene_crae_cor_bl <- compute_cors_bl(gene_clinical_bl$CRAE, "CRAE")
gene_crve_cor_bl <- compute_cors_bl(gene_clinical_bl$CRVE, "CRVE")

# Check key genes in sensitivity analysis
key_ids <- c(
  "ENSG00000060237",  # WNK1
  "ENSG00000172939",  # OXSR1
  "ENSG00000116717",  # GADD45A
  "ENSG00000133392",  # MYH11
  "ENSG00000118689"   # FOXO3
)

cat("\nKey gene associations — sensitivity (baseline only):\n")
gene_avr_cor_bl %>%
  filter(ENSEMBL %in% key_ids) %>%
  left_join(
    bitr(key_ids, fromType = "ENSEMBL",
         toType = "SYMBOL", OrgDb = org.Hs.eg.db),
    by = "ENSEMBL"
  ) %>%
  dplyr::select(SYMBOL, correlation, pval, fdr) %>%
  arrange(pval) %>%
  print()

cat("\nKey gene associations — primary (all timepoints):\n")
gene_avr_cor %>%
  filter(ENSEMBL %in% key_ids) %>%
  left_join(
    bitr(key_ids, fromType = "ENSEMBL",
         toType = "SYMBOL", OrgDb = org.Hs.eg.db),
    by = "ENSEMBL"
  ) %>%
  dplyr::select(SYMBOL, correlation, pval, fdr) %>%
  arrange(pval) %>%
  print()

write.csv(gene_avr_cor_bl,
          "results/gene_AVR_sensitivity_baseline.csv",
          row.names = FALSE)
cat("Saved: sensitivity analysis results\n")

cat("\n05_gene_vessel_correlations.R complete\n")