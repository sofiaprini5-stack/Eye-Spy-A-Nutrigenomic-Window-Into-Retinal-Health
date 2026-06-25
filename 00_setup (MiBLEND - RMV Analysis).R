# =============================================================================
# 00_setup.R
# MiBLEND Thesis — Sofia M. Prini
# Maastricht University | 2025-2026
#
# Load all libraries, define file paths, and set shared options.
# Source this file at the top of every other script with:
# source("00_setup.R")
#
# Run this block ONCE if packages are not yet installed:
# install.packages(c("tidyverse","readxl","lme4","lmerTest",
#                    "ggplot2","ggrepel","patchwork","RColorBrewer","statmod"))
# install.packages("BiocManager")
# BiocManager::install(c("limma","edgeR","clusterProfiler",
#                        "org.Hs.eg.db","enrichplot","fgsea","biomaRt"))
# =============================================================================

# ── Libraries ─────────────────────────────────────────────────────────────────

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

# Resolve function conflicts explicitly
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

cat("Libraries loaded\n")

# ── File paths ────────────────────────────────────────────────────────────────
# Update these paths if your files are in a different location

PATHS <- list(
  # Input files
  masterfile  = file.choose,   # called as PATHS$masterfile() when needed
  expr_file   = file.choose,   # called as PATHS$expr_file() when needed
  
  # Saved objects directory — all .rds files go here
  rds_dir     = "rds_objects/",
  
  # Output directories
  results_dir = "results/",
  figures_dir = "figures/"
)

# Create output directories if they don't exist
dir.create(PATHS$rds_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(PATHS$results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(PATHS$figures_dir, showWarnings = FALSE, recursive = TRUE)

# ── SNP column names ──────────────────────────────────────────────────────────
# Used across multiple scripts — defined once here

SNP_COLS <- c("GSTM1", "GSTT1", "GSTP1", "NQO1", "COMT", "MTHFR",
              "CAT1", "Glu298Asp", "XRCC1", "BCMO1", "SLC23A1",
              "ZBED3", "APOE", "HNF1A", "TCF7L2")

BINARY_SNPS   <- c("GSTM1", "GSTT1")
ADDITIVE_SNPS <- setdiff(SNP_COLS, BINARY_SNPS)

# ── Plot theme ────────────────────────────────────────────────────────────────
# Shared ggplot2 theme applied to all figures

theme_miblend <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(
      legend.position  = "bottom",
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold"),
      plot.caption     = element_text(size = 8, colour = "grey50")
    )
}

# ── Colour palette ────────────────────────────────────────────────────────────
# Consistent colours used across all figures

COLOURS <- list(
  responder      = "#2471A3",   # blue  — Responder group
  poor_responder = "#C0392B",   # red   — Poor Responder group
  significant    = "#2171b5",   # blue  — significant result
  trend          = "#fd8d3c",   # orange — trend (p<0.1)
  ns             = "grey60"     # grey  — not significant
)

cat("Setup complete\n")
cat("Directories ready:", PATHS$rds_dir, PATHS$results_dir, PATHS$figures_dir, "\n")
