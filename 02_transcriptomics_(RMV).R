# =============================================================================
# 02_transcriptomics.R
#
# Purpose: Load gene expression data, link to masterfile, run TMM
#          normalisation, voom transformation, iterative
#          duplicateCorrelation, DEG analysis and GSEA.
#
# Input:   rds_objects/master.rds  (from 01_data_preparation.R)
#          Gene expression xlsx    (file chooser dialog will open)
#          Sample key CSV          (file chooser dialog will open)
#            → MiBlend_SampleKey_Masterlist_220714_LongFormat_edit_230310.csv
#            → Required to recover DI values for SeqBatch 2/3 samples whose
#              filenames omit the _D component (e.g. P120_T3_S21 vs P47_T1_D0_S136)
#
# Saves:   rds_objects/logCPM.rds
#          rds_objects/gene_clinical.rds
#          rds_objects/deg_results.rds
#
# Output:  results/DEG_results_AllDIs.csv
#          results/GSEA_KEGG_results.csv
#          results/GSEA_GO_results.csv
#          figures/plot_volcano.png
#          figures/plot_GSEA_KEGG_dotplot.png
# =============================================================================

library(limma)
library(edgeR)
library(statmod)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(fgsea)
library(biomaRt)
library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(ggrepel)
library(patchwork)
library(RColorBrewer)

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

cat("Libraries loaded\n")

dir.create("rds_objects", showWarnings = FALSE)
dir.create("results",     showWarnings = FALSE)
dir.create("figures",     showWarnings = FALSE)

# ── Load master ───────────────────────────────────────────────────────────────

master <- readRDS("rds_objects/master.rds")
cat("master.rds loaded:", nrow(master), "rows\n")

snp_cols <- c("GSTM1","GSTT1","GSTP1","NQO1","COMT","MTHFR",
              "CAT1","Glu298Asp","XRCC1","BCMO1","SLC23A1",
              "ZBED3","APOE","HNF1A","TCF7L2")

# ── B1. Load gene expression ──────────────────────────────────────────────────

cat("\nSelect gene expression file (.xlsx)...\n")
expr_raw <- read_excel(file.choose(), sheet=1, skip=1)

expr_df           <- as.data.frame(expr_raw)
rownames(expr_df) <- expr_df[, 1]
expr_df           <- expr_df[, -1]
expr_df           <- expr_df %>%
  mutate(across(everything(), ~ as.numeric(as.character(.))))

cat("Genes:", nrow(expr_df), "| Samples:", ncol(expr_df), "\n")

# ── B2. Link gene file to masterfile ─────────────────────────────────────────
#
# The gene expression file contains samples from two sequencing batches:
#   SeqBatch 1: filenames include DI  e.g. P47_T1_D0_S136
#   SeqBatch 2/3: filenames omit DI   e.g. P120_T3_S21
#
# Without the DI field, inner_join() on PCode + PERIOD + DI fails for
# SeqBatch 2/3 samples, dropping ~96 participants. The sample key file
# (provided by Simone van Breda) maps every SampleSeqName to its correct
# DI value regardless of batch, allowing all participants to be recovered.

snp_data <- master %>%
  distinct(PCode, .keep_all = TRUE) %>%
  mutate(PCode_chr = as.character(PCode)) %>%
  dplyr::select(PCode_chr, all_of(snp_cols))

# ── B2a. Load sample key and build DI lookup ──────────────────────────────────
cat("\nSelect sample key CSV (MiBlend_SampleKey_Masterlist)...\n")
sample_key <- read.csv(file.choose(), stringsAsFactors = FALSE)

# Clean DI column: remove "D" prefix → numeric (D0 → 0, D1 → 1 etc.)
sample_key <- sample_key %>%
  mutate(
    DI_numeric     = as.numeric(str_remove(DI, "D")),
    PERIOD_numeric = as.integer(TD)
  )

# Lookup table: SampleSeqName → DI and PERIOD
di_lookup <- sample_key %>%
  dplyr::select(SampleSeqName, DI_numeric, PERIOD_numeric) %>%
  dplyr::rename(
    gene_col        = SampleSeqName,
    DI_from_key     = DI_numeric,
    PERIOD_from_key = PERIOD_numeric
  )

cat("Sample key loaded:", nrow(di_lookup), "entries\n")

# ── B2b. Build sample_map with DI recovered for all batches ──────────────────
sample_map <- tibble(gene_col = colnames(expr_df)) %>%
  mutate(
    PCode  = as.integer(str_extract(gene_col, "(?<=P)\\d+")),
    PERIOD = as.integer(str_extract(gene_col, "(?<=_T)\\d+")),
    DI     = as.numeric(str_extract(gene_col, "(?<=_D)\\d+"))
  ) %>%
  left_join(di_lookup, by = "gene_col") %>%
  mutate(
    # Use DI from filename where available; fill from sample key otherwise
    DI     = if_else(is.na(DI),     DI_from_key,     DI),
    PERIOD = if_else(is.na(PERIOD), PERIOD_from_key, PERIOD)
  ) %>%
  dplyr::select(gene_col, PCode, PERIOD, DI)

# Rows with NA PCode are header/metadata artefacts from xlsx import — drop
sample_map <- sample_map %>% filter(!is.na(PCode))

cat("sample_map: ", nrow(sample_map), "samples |",
    n_distinct(sample_map$PCode), "unique participants\n")
cat("Samples with missing DI:", sum(is.na(sample_map$DI)), "\n")

# ── B2c. Join to masterfile ───────────────────────────────────────────────────
master_num <- master %>%
  mutate(PCode_n  = as.integer(as.character(PCode)),
         PERIOD_n = as.integer(as.character(PERIOD)))

gene_clinical <- sample_map %>%
  inner_join(
    master_num %>%
      dplyr::select(PCode_n, PERIOD_n, DI, CRAE, CRVE, AVR,
                    timepoint, responder_group, overall_score,
                    group_detox, group_methylation,
                    group_antioxidant, group_carotenoid,
                    group_metabolic, Age_z, Sex, BMI_z, WBC_z, Season),
    by = c("PCode" = "PCode_n", "PERIOD" = "PERIOD_n", "DI" = "DI")
  ) %>%
  mutate(
    PCode_chr = as.character(PCode),
    PCode     = droplevels(as.factor(PCode)),
    Sex       = droplevels(as.factor(Sex)),
    Season    = droplevels(as.factor(Season)),
    subset    = droplevels(as.factor(case_when(
      timepoint == "Baseline" ~ paste0("Baseline_D", DI),
      TRUE                    ~ paste0("PostTest_D", DI)
    )))
  ) %>%
  left_join(snp_data, by = "PCode_chr") %>%
  dplyr::select(-PCode_chr)

cat("Before NA removal:", nrow(gene_clinical), "samples |",
    n_distinct(gene_clinical$PCode), "participants\n")

complete_rows <- complete.cases(
  gene_clinical[, c("subset", "Sex", "Age_z", "Season", "WBC_z")]
)
gene_clinical <- gene_clinical[complete_rows, ]
expr_clean    <- as.matrix(expr_df[, gene_clinical$gene_col])

cat("After NA removal:", nrow(gene_clinical), "samples |",
    n_distinct(gene_clinical$PCode), "participants\n")
cat("SNPs in gene_clinical:", sum(snp_cols %in% colnames(gene_clinical)), "/ 15\n")
stopifnot(all(colnames(expr_clean) == gene_clinical$gene_col))
cat("Alignment confirmed\n")

# ── B3. TMM normalisation + voom ─────────────────────────────────────────────
# Using the updated version with NA checks and zero-variance gene filter

dge <- DGEList(counts = expr_clean)
dge <- calcNormFactors(dge, method = "TMM")

design_deg <- model.matrix(
  ~ 0 + subset + Sex + Age_z + Season + WBC_z,
  data = gene_clinical
)

# Remove samples with NA in design matrix
keep_samples <- complete.cases(design_deg)
if (sum(!keep_samples) > 0) {
  cat("Removing", sum(!keep_samples), "samples with NA in design matrix\n")
  gene_clinical <- gene_clinical[keep_samples, ]
  expr_clean    <- expr_clean[, keep_samples]
  dge           <- DGEList(counts = expr_clean)
  dge           <- calcNormFactors(dge, method = "TMM")
  design_deg    <- model.matrix(
    ~ 0 + subset + Sex + Age_z + Season + WBC_z,
    data = gene_clinical
  )
}

cat("Samples after NA removal:", ncol(dge$counts), "\n")
cat("Design matrix rows:", nrow(design_deg), "\n")
cat("Match:", ncol(dge$counts) == nrow(design_deg), "\n")

# Remove zero-variance genes
keep_genes <- rowSums(dge$counts > 0) >= 3
dge        <- dge[keep_genes, ]
cat("Genes after filtering:", nrow(dge$counts), "\n")

# voom round 1
v1 <- voom(dge, design = design_deg, plot = FALSE)

# duplicateCorrelation round 1
cat("duplicateCorrelation round 1 (~3-5 min)...\n")
dc1 <- duplicateCorrelation(v1, design = design_deg,
                            block = gene_clinical$PCode)
cat("Correlation r1:", round(dc1$consensus.correlation, 3), "\n")

if (is.na(dc1$consensus.correlation) ||
    abs(dc1$consensus.correlation) >= 1) {
  cat("WARNING: Invalid correlation in round 1, using 0\n")
  dc1$consensus.correlation <- 0
}

# voom round 2
v2  <- voom(dge, design = design_deg, block = gene_clinical$PCode,
            correlation = dc1$consensus.correlation, plot = FALSE)
dc2 <- duplicateCorrelation(v2, design = design_deg,
                            block = gene_clinical$PCode)
cat("Correlation r2:", round(dc2$consensus.correlation, 3), "\n")

if (is.na(dc2$consensus.correlation) ||
    abs(dc2$consensus.correlation) >= 1) {
  cat("WARNING: Invalid correlation in round 2, using round 1 value\n")
  dc2$consensus.correlation <- dc1$consensus.correlation
}

logCPM <- v2$E
cat("logCPM stored:", nrow(logCPM), "genes x", ncol(logCPM), "samples\n")

# ── B4. Limma DEG analysis ────────────────────────────────────────────────────

fit <- lmFit(v2, design = design_deg, block = gene_clinical$PCode,
             correlation = dc2$consensus.correlation)

post_cols <- grep("PostTest", colnames(design_deg), value = TRUE)
base_cols <- grep("Baseline", colnames(design_deg), value = TRUE)

alldis_str <- paste0(
  "(", paste(post_cols, collapse=" + "), ")/", length(post_cols),
  " - (",
  paste(base_cols, collapse=" + "), ")/", length(base_cols)
)

contr_matrix           <- makeContrasts(contrasts = alldis_str,
                                        levels = colnames(design_deg))
colnames(contr_matrix) <- "AllDIs"
fit2                   <- contrasts.fit(fit, contr_matrix)
fit_ebayes             <- eBayes(fit2)

deg_results <- topTable(fit_ebayes,
                        coef    = "AllDIs",
                        n       = Inf,
                        sort.by = "P")

cat("\nTop 10 DEGs:\n")
print(head(deg_results[, c("logFC","AveExpr","P.Value","adj.P.Val")], 10))

write.csv(deg_results, "results/DEG_results_AllDIs.csv")
cat("Saved: results/DEG_results_AllDIs.csv\n")

# ── B5. Volcano plot ──────────────────────────────────────────────────────────

pathway_genes <- c(
  RHOA  = "ENSG00000067560", RAC1  = "ENSG00000136238",
  CDC42 = "ENSG00000070831", ROCK1 = "ENSG00000067900",
  ROCK2 = "ENSG00000134318", PAK1  = "ENSG00000149269",
  NOS3  = "ENSG00000164867", NFE2L2 = "ENSG00000116044"
)

deg_plot <- deg_results %>%
  rownames_to_column("ENSEMBL") %>%
  mutate(
    is_key   = ENSEMBL %in% pathway_genes,
    Gene     = ifelse(is_key,
                      names(pathway_genes)[match(ENSEMBL, pathway_genes)],
                      NA),
    category = case_when(
      is_key                               ~ "Key pathway gene",
      adj.P.Val < 0.05 & abs(logFC) > 0.5 ~ "Significant DEG",
      TRUE                                  ~ "Not significant"
    ) %>% factor(levels = c("Not significant","Significant DEG",
                            "Key pathway gene"))
  )

ggplot(deg_plot, aes(x = logFC, y = -log10(P.Value), colour = category)) +
  geom_point(data = dplyr::filter(deg_plot, category == "Not significant"),
             alpha = 0.2, size = 0.6) +
  geom_point(data = dplyr::filter(deg_plot, category == "Significant DEG"),
             size = 1.5) +
  geom_point(data = dplyr::filter(deg_plot, is_key), size = 3.5) +
  geom_label_repel(data = dplyr::filter(deg_plot, is_key),
                   aes(label = Gene), size = 3.5, box.padding = 0.6) +
  scale_colour_manual(
    values = c("Not significant"  = "grey80",
               "Significant DEG"  = "#F39C12",
               "Key pathway gene" = "#E74C3C"),
    name = NULL) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "red", alpha = 0.6) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
             colour = "steelblue", alpha = 0.6) +
  labs(title    = "DEG: Post-Intervention vs Baseline (AllDIs pooled)",
       subtitle = "voom + duplicateCorrelation",
       x = "Log2 Fold-Change", y = "-log10(P-value)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave("figures/plot_volcano.png", width = 9, height = 6, dpi = 300)
cat("Saved: figures/plot_volcano.png\n")

# ── B6. GSEA ─────────────────────────────────────────────────────────────────

all_entrez <- bitr(rownames(deg_results), fromType = "ENSEMBL",
                   toType = "ENTREZID", OrgDb = org.Hs.eg.db)

gene_ranks <- deg_results %>%
  rownames_to_column("ENSEMBL") %>%
  inner_join(all_entrez, by = "ENSEMBL") %>%
  arrange(desc(logFC)) %>%
  distinct(ENTREZID, .keep_all = TRUE)

ranked_vec <- setNames(gene_ranks$logFC, gene_ranks$ENTREZID)

gsea_kegg <- gseKEGG(geneList = ranked_vec, organism = "hsa",
                     pvalueCutoff = 0.05, verbose = FALSE)
gsea_go   <- gseGO(geneList = ranked_vec, OrgDb = org.Hs.eg.db,
                   ont = "BP", pvalueCutoff = 0.05, verbose = FALSE)

cat("GSEA KEGG significant:", nrow(as.data.frame(gsea_kegg)), "\n")
cat("GSEA GO significant:",   nrow(as.data.frame(gsea_go)),   "\n")

write.csv(as.data.frame(gsea_kegg),
          "results/GSEA_KEGG_results.csv", row.names = FALSE)
write.csv(as.data.frame(gsea_go),
          "results/GSEA_GO_results.csv",   row.names = FALSE)

if (nrow(as.data.frame(gsea_kegg)) > 0) {
  dotplot(gsea_kegg, showCategory = 20,
          title = "GSEA KEGG Pathways") + theme_minimal()
  ggsave("figures/plot_GSEA_KEGG_dotplot.png", width = 11, height = 8, dpi = 300)
}
cat("Saved: GSEA outputs\n")

# ── Save objects ──────────────────────────────────────────────────────────────

saveRDS(logCPM,        "rds_objects/logCPM.rds")
saveRDS(gene_clinical, "rds_objects/gene_clinical.rds")
saveRDS(deg_results,   "rds_objects/deg_results.rds")

cat("\nSaved:\n")
cat("  rds_objects/logCPM.rds\n")
cat("  rds_objects/gene_clinical.rds\n")
cat("  rds_objects/deg_results.rds\n")
cat("\n02_transcriptomics.R complete\n")