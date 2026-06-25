# =============================================================================
# 04_vessel_stratified.R
#
# Purpose: Reverse engineering analysis — classify participants by baseline
#          vessel calibre, run chi-square screening across 15 SNPs,
#          followed by logistic regression for nominally significant SNPs.
#          Generate enhanced forest plots with allele counts.
#
# Input:   rds_objects/master.rds  (from 01_data_preparation.R)
#
# Output:  results/fdr_corrected_chisquare_results.csv
#          results/logistic_regression_CRAE_results.csv
#          results/logistic_regression_AVR_results.csv
#          results/top25_DEG_SNP_significant_hits.csv
#          figures/plot_forest_v3_CRAE_allcounts.png
#          figures/plot_forest_v3_AVR_allcounts.png
# =============================================================================

library(tidyverse)
library(readxl)
library(ggrepel)
library(patchwork)
library(RColorBrewer)

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

cat("Libraries loaded\n")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ── Load master ───────────────────────────────────────────────────────────────

master <- readRDS("rds_objects/master.rds")
cat("master.rds loaded:", nrow(master), "rows\n")

snp_cols      <- c("GSTM1","GSTT1","GSTP1","NQO1","COMT","MTHFR",
                   "CAT1","Glu298Asp","XRCC1","BCMO1","SLC23A1",
                   "ZBED3","APOE","HNF1A","TCF7L2")
binary_snps   <- c("GSTM1","GSTT1")
additive_snps <- setdiff(snp_cols, binary_snps)

# ── D1. Baseline data + SNP numeric coding ────────────────────────────────────

baseline <- master %>%
  filter(PERIOD == 1) %>%
  distinct(PCode, .keep_all = TRUE) %>%
  mutate(
    Age_z = as.numeric(scale(Age)),
    BMI_z = as.numeric(scale(BMI))
  )

for (snp in binary_snps) {
  baseline[[paste0(snp, "_num")]] <- as.numeric(as.character(baseline[[snp]]))
}
for (snp in additive_snps) {
  baseline[[paste0(snp, "_num")]] <- case_when(
    baseline[[snp]] == 2 ~ 0,
    baseline[[snp]] == 3 ~ 1,
    baseline[[snp]] == 4 ~ 2
  )
}

# ── D2. Median splits ─────────────────────────────────────────────────────────

crae_median <- median(baseline$CRAE, na.rm = TRUE)
avr_median  <- median(baseline$AVR,  na.rm = TRUE)

baseline <- baseline %>%
  mutate(
    CRAE_wide = as.integer(CRAE > crae_median),
    AVR_high  = as.integer(AVR  > avr_median)
  )

cat("CRAE median:", round(crae_median, 1), "um\n")
cat("AVR median:",  round(avr_median,  3), "\n")

# ── D3. Chi-square + FDR ──────────────────────────────────────────────────────

all_snps <- c(binary_snps, additive_snps)

run_chisq <- function(data, group_col, snp) {
  ct <- table(data[[group_col]], data[[snp]])
  if (ncol(ct) < 2 || nrow(ct) < 2) return(NA)
  test <- tryCatch(
    chisq.test(ct),
    warning = function(w) chisq.test(ct, simulate.p.value = TRUE)
  )
  return(test$p.value)
}

crae_pvals <- sapply(all_snps, function(s) run_chisq(baseline, "CRAE_wide", s))
avr_pvals  <- sapply(all_snps, function(s) run_chisq(baseline, "AVR_high",  s))

crae_fdr <- p.adjust(crae_pvals, method = "BH")
avr_fdr  <- p.adjust(avr_pvals,  method = "BH")

fdr_results <- data.frame(
  SNP      = all_snps,
  CRAE_p   = round(crae_pvals, 4),
  CRAE_FDR = round(crae_fdr,   4),
  CRAE_sig = case_when(
    crae_fdr   < 0.05 ~ "FDR<0.05",
    crae_pvals < 0.05 ~ "nominal",
    crae_pvals < 0.1  ~ "trend",
    TRUE              ~ "NS"),
  AVR_p    = round(avr_pvals,  4),
  AVR_FDR  = round(avr_fdr,   4),
  AVR_sig  = case_when(
    avr_fdr   < 0.05 ~ "FDR<0.05",
    avr_pvals < 0.05 ~ "nominal",
    avr_pvals < 0.1  ~ "trend",
    TRUE             ~ "NS")
)

cat("\nChi-square + FDR results:\n")
print(fdr_results, row.names = FALSE)
write.csv(fdr_results, "results/fdr_corrected_chisquare_results.csv",
          row.names = FALSE)
cat("Saved: results/fdr_corrected_chisquare_results.csv\n")

# ── D4. Logistic regression ───────────────────────────────────────────────────

covariates <- c("Age_z", "Sex", "BMI_z", "Season")

run_logistic <- function(data, outcome_col, snp) {
  snp_num     <- paste0(snp, "_num")
  formula_str <- paste(outcome_col, "~", snp_num, "+",
                       paste(covariates, collapse = " + "))
  model <- tryCatch(
    glm(as.formula(formula_str), data = data, family = binomial),
    error = function(e) NULL
  )
  if (is.null(model)) return(NULL)
  coef_s  <- summary(model)$coefficients
  snp_row <- coef_s[snp_num, ]
  data.frame(
    SNP      = snp,
    OR       = round(exp(snp_row["Estimate"]), 3),
    CI_lower = round(exp(snp_row["Estimate"] - 1.96 * snp_row["Std. Error"]), 3),
    CI_upper = round(exp(snp_row["Estimate"] + 1.96 * snp_row["Std. Error"]), 3),
    p_value  = round(snp_row["Pr(>|z|)"], 4)
  )
}

crae_log <- do.call(rbind, lapply(all_snps, function(s)
  run_logistic(baseline, "CRAE_wide", s))) %>%
  mutate(FDR = round(p.adjust(p_value, method = "BH"), 4),
         sig = case_when(
           FDR     < 0.05 ~ "Significant",
           p_value < 0.1  ~ "NS (p<0.1)",
           TRUE           ~ "NS"))

avr_log <- do.call(rbind, lapply(all_snps, function(s)
  run_logistic(baseline, "AVR_high", s))) %>%
  mutate(FDR = round(p.adjust(p_value, method = "BH"), 4),
         sig = case_when(
           FDR     < 0.05 ~ "Significant",
           p_value < 0.1  ~ "NS (p<0.1)",
           TRUE           ~ "NS"))

write.csv(crae_log, "results/logistic_regression_CRAE_results.csv",
          row.names = FALSE)
write.csv(avr_log,  "results/logistic_regression_AVR_results.csv",
          row.names = FALSE)
cat("Saved: logistic regression results\n")

# ── D5. Enhanced forest plots with allele counts ──────────────────────────────

allele_labels <- list(
  GSTM1     = c("0"="Null",        "1"="Present"),
  GSTT1     = c("0"="Null",        "1"="Present"),
  GSTP1     = c("2"="Ile/Ile",     "3"="Ile/Val",     "4"="Val/Val"),
  NQO1      = c("2"="CC",          "3"="CT",           "4"="TT"),
  COMT      = c("2"="Met/Met",     "3"="Val/Met",      "4"="Val/Val"),
  MTHFR     = c("2"="CC",          "3"="CT",           "4"="TT"),
  CAT1      = c("2"="CC",          "3"="CT",           "4"="TT"),
  Glu298Asp = c("2"="Glu/Glu",     "3"="Glu/Asp",     "4"="Asp/Asp"),
  XRCC1     = c("2"="Arg/Arg",     "3"="Arg/Gln",     "4"="Gln/Gln"),
  BCMO1     = c("2"="Wild-type",   "3"="Heterozygous", "4"="Variant"),
  SLC23A1   = c("2"="Wild-type",   "3"="Heterozygous", "4"="Variant"),
  ZBED3     = c("2"="Wild-type",   "3"="Heterozygous", "4"="Variant"),
  APOE      = c("2"="e2/e2",       "3"="e3/e3",        "4"="e4/e4"),
  HNF1A     = c("2"="Wild-type",   "3"="Heterozygous", "4"="Variant"),
  TCF7L2    = c("2"="CC",          "3"="CT",           "4"="TT")
)

count_df <- map_dfr(all_snps, function(snp) {
  labels <- allele_labels[[snp]]
  baseline %>%
    count(.data[[snp]]) %>%
    mutate(
      SNP    = snp,
      allele = recode(as.character(.data[[snp]]), !!!labels),
      role   = case_when(
        as.character(.data[[snp]]) %in% c("2","1") ~ "Protective",
        as.character(.data[[snp]]) == "3"          ~ "Heterozygous",
        as.character(.data[[snp]]) %in% c("4","0") ~ "Risk",
        TRUE ~ "Unknown"
      )
    ) %>%
    dplyr::select(SNP, allele, n, role)
})

make_combined_plot <- function(log_results, count_data,
                               outcome_label, outcome_note) {
  snp_order <- rev(log_results$SNP)
  
  count_plot <- count_data %>%
    dplyr::filter(SNP %in% log_results$SNP) %>%
    mutate(
      SNP  = factor(SNP,  levels = snp_order),
      role = factor(role, levels = c("Protective","Heterozygous","Risk"))
    ) %>%
    ggplot(aes(x = n, y = SNP, fill = role)) +
    geom_col(width = 0.6, colour = "white", linewidth = 0.3) +
    geom_text(aes(label = paste0(allele, "\nn=", n)),
              position = position_stack(vjust = 0.5),
              size = 2.8, colour = "white", fontface = "bold",
              lineheight = 0.9) +
    scale_fill_manual(
      values = c("Protective"   = "#2471A3",
                 "Heterozygous" = "#F39C12",
                 "Risk"         = "#C0392B"),
      name = "Allele group") +
    labs(x = "N participants", y = NULL,
         title = "Genotype distribution") +
    theme_minimal(base_size = 11) +
    theme(axis.text.y     = element_text(size = 10, face = "bold"),
          legend.position = "bottom",
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank())
  
  forest_plot <- log_results %>%
    mutate(SNP = factor(SNP, levels = snp_order)) %>%
    ggplot(aes(x = OR, y = SNP, colour = sig)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               colour = "grey50", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper),
                   height = 0.25, linewidth = 0.8) +
    geom_point(size = 3.5) +
    scale_colour_manual(
      values = c("Significant" = "#2171b5",
                 "NS (p<0.1)"  = "#fd8d3c",
                 "NS"          = "grey60"),
      name = NULL) +
    coord_cartesian(xlim = c(0, 5)) +
    labs(x = "Odds Ratio (95% CI)", y = NULL,
         title    = paste0("OR: ", outcome_label),
         subtitle = outcome_note) +
    theme_minimal(base_size = 11) +
    theme(axis.text.y      = element_blank(),
          axis.ticks.y     = element_blank(),
          legend.position  = "bottom",
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank())
  
  count_plot + forest_plot +
    plot_layout(widths = c(1.8, 1)) +
    plot_annotation(
      title   = paste0("SNP Association with ", outcome_label,
                       " — All 15 SNPs"),
      caption = "Adjusted for Age, Sex, BMI, Season | BH-FDR applied"
    )
}

p_crae_forest <- make_combined_plot(
  crae_log, count_df,
  "Wide vs Narrow Arterioles (CRAE)",
  "OR > 1 = protective allele -> wider arterioles")

p_avr_forest <- make_combined_plot(
  avr_log, count_df,
  "High vs Low AVR",
  "OR > 1 = each additional risk allele associated with higher AVR odds")

ggsave("figures/plot_forest_v3_CRAE_allcounts.png", p_crae_forest,
       width = 14, height = 10, dpi = 200)
ggsave("figures/plot_forest_v3_AVR_allcounts.png",  p_avr_forest,
       width = 14, height = 10, dpi = 200)
cat("Saved: enhanced forest plots\n")

cat("\n04_vessel_stratified.R complete\n")
