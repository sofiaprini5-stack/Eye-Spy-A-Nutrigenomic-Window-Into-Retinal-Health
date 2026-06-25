# =============================================================================
# 03_clinical_models.R
# MiBLEND Thesis — Sofia M. Prini
# Maastricht University | 2025-2026
#
# Purpose: Run linear mixed models (LMMs) for CRAE, CRVE and AVR.
#          Generate trajectory plots and genetic group forest plot.
#
# Input:   rds_objects/master.rds  (from 01_data_preparation.R)
#
# Output:  results/results_model*.txt  (all LMM summaries)
#          figures/plot_CRAE_trajectories.png
#          figures/plot_CRVE_trajectories.png
#          figures/plot_AVR_trajectories.png
#          figures/plot_forest_genetic_groups.png
#          figures/plot_Glu298Asp_CRAE_trajectories.png
#          figures/plot_Glu298Asp_CRVE_trajectories.png
#          figures/plot_COMT_CRAE_trajectories.png
#          figures/plot_COMT_CRVE_trajectories.png
# =============================================================================

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

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ── Load master ───────────────────────────────────────────────────────────────

master <- readRDS("rds_objects/master.rds")
cat("master.rds loaded:", nrow(master), "rows\n")

# ── C1. Trajectory plots ──────────────────────────────────────────────────────

p_crae <- master %>%
  dplyr::filter(!is.na(responder_group), !is.na(timepoint_label)) %>%
  ggplot(aes(x = timepoint_label, y = CRAE, group = PCode,
             colour = responder_group)) +
  geom_line(alpha = 0.2, linewidth = 0.4) +
  stat_summary(aes(group = responder_group, fill = responder_group),
               fun.data = mean_se, geom = "ribbon", alpha = 0.2, colour = NA) +
  stat_summary(aes(group = responder_group), fun = mean,
               geom = "line", linewidth = 2) +
  stat_summary(aes(group = responder_group), fun = mean,
               geom = "point", size = 3.5) +
  scale_colour_manual(
    values = c("Poor_Responder" = "#C0392B", "Responder" = "#2471A3"),
    name   = "Genetic Group",
    labels = c("Poor Responder", "Responder")) +
  scale_fill_manual(
    values = c("Poor_Responder" = "#C0392B", "Responder" = "#2471A3"),
    guide  = "none") +
  labs(title    = "Retinal Arteriolar (CRAE) Trajectories by Genetic Group",
       subtitle = "Shaded band = mean ± SE | Thin lines = individuals",
       x = "Timepoint", y = "CRAE (µm)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

p_crve <- master %>%
  dplyr::filter(!is.na(responder_group), !is.na(timepoint_label)) %>%
  ggplot(aes(x = timepoint_label, y = CRVE, group = PCode,
             colour = responder_group)) +
  geom_line(alpha = 0.2, linewidth = 0.4) +
  stat_summary(aes(group = responder_group, fill = responder_group),
               fun.data = mean_se, geom = "ribbon", alpha = 0.2, colour = NA) +
  stat_summary(aes(group = responder_group), fun = mean,
               geom = "line", linewidth = 2) +
  stat_summary(aes(group = responder_group), fun = mean,
               geom = "point", size = 3.5) +
  scale_colour_manual(
    values = c("Poor_Responder" = "#C0392B", "Responder" = "#2471A3"),
    name   = "Genetic Group",
    labels = c("Poor Responder", "Responder")) +
  scale_fill_manual(
    values = c("Poor_Responder" = "#C0392B", "Responder" = "#2471A3"),
    guide  = "none") +
  labs(title    = "Retinal Venular (CRVE) Trajectories by Genetic Group",
       subtitle = "Shaded band = mean ± SE | Thin lines = individuals",
       x = "Timepoint", y = "CRVE (µm)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave("figures/plot_CRAE_trajectories.png", p_crae, width = 9, height = 6, dpi = 300)
ggsave("figures/plot_CRVE_trajectories.png", p_crve, width = 9, height = 6, dpi = 300)
cat("Saved: trajectory plots\n")

# ── C2. Linear mixed models ───────────────────────────────────────────────────

# CRAE models
model1 <- lmer(CRAE ~ PERIOD + Age_z + Sex + BMI_z + WBC_z + Season +
                 (1|PCode), data = master, REML = TRUE)

model2 <- lmer(CRAE ~ PERIOD + overall_score + Age_z + Sex + BMI_z +
                 WBC_z + Season + (1|PCode), data = master, REML = TRUE)

model3 <- lmer(CRAE ~ PERIOD * responder_group + Age_z + Sex + BMI_z +
                 WBC_z + Season + (1|PCode),
               data = master %>% dplyr::filter(!is.na(responder_group)),
               REML = TRUE)

model4 <- lmer(CRAE ~ PERIOD + group_detox + group_methylation +
                 group_antioxidant + group_carotenoid + group_metabolic +
                 Age_z + Sex + BMI_z + WBC_z + Season + (1|PCode),
               data = master, REML = TRUE)

# CRVE models
model_CRVE_A <- lmer(CRVE ~ PERIOD + Age_z + Sex + BMI_z + WBC_z +
                       Season + (1|PCode), data = master, REML = TRUE)

model_CRVE_B <- lmer(CRVE ~ PERIOD * responder_group + Age_z + Sex +
                       BMI_z + WBC_z + Season + (1|PCode),
                     data = master %>% dplyr::filter(!is.na(responder_group)),
                     REML = TRUE)

model_CRVE_C <- lmer(CRVE ~ PERIOD + group_detox + group_methylation +
                       group_antioxidant + group_carotenoid +
                       group_metabolic + Age_z + Sex + BMI_z + WBC_z +
                       Season + (1|PCode), data = master, REML = TRUE)

# AVR model
model_AVR <- lmer(AVR ~ PERIOD * responder_group + Age_z + Sex + BMI_z +
                    WBC_z + Season + (1|PCode),
                  data = master %>% dplyr::filter(!is.na(responder_group)),
                  REML = TRUE)

# Save all model summaries
for (nm in c("model1","model2","model3","model4",
             "model_CRVE_A","model_CRVE_B","model_CRVE_C","model_AVR")) {
  sink(paste0("results/results_", nm, ".txt"))
  print(summary(get(nm)))
  sink()
}
cat("Saved: all LMM results\n")

# ── C3. Forest plot — genetic groups ──────────────────────────────────────────

coef_forest <- as.data.frame(coef(summary(model4))) %>%
  rownames_to_column("term") %>%
  dplyr::rename(estimate = Estimate, se = `Std. Error`, pval = `Pr(>|t|)`) %>%
  dplyr::filter(str_detect(term, "group_")) %>%
  mutate(
    label = recode(term,
                   "group_detox"       = "Phase II Detox\n(GSTM1, GSTT1, GSTP1, NQO1)",
                   "group_methylation" = "Methylation\n(COMT, MTHFR)",
                   "group_antioxidant" = "Antioxidant + Vascular\n(CAT1, eNOS/Glu298Asp, XRCC1)",
                   "group_carotenoid"  = "Carotenoid Uptake\n(BCMO1, SLC23A1, ZBED3)",
                   "group_metabolic"   = "Metabolic/CVD Risk\n(APOE, HNF1A, TCF7L2)"),
    significant = pval < 0.05,
    ci_low      = estimate - 1.96 * se,
    ci_high     = estimate + 1.96 * se
  )

ggplot(coef_forest, aes(x = estimate, y = reorder(label, estimate),
                        colour = significant)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.25) +
  geom_point(size = 4) +
  scale_colour_manual(
    values = c("FALSE" = "#E07070", "TRUE" = "#4A90D9"),
    labels = c("p >= 0.05", "p < 0.05"), name = NULL) +
  labs(title    = "Effect of Genetic Groups on CRAE",
       subtitle = "Positive = wider arterioles | Bars = 95% CI",
       x = "Effect on CRAE (µm)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 10))

ggsave("figures/plot_forest_genetic_groups.png", width = 9, height = 6, dpi = 300)
cat("Saved: figures/plot_forest_genetic_groups.png\n")

# ── C4. AVR trajectory plot ───────────────────────────────────────────────────

p_avr <- master %>%
  dplyr::filter(!is.na(responder_group), !is.na(timepoint_label)) %>%
  ggplot(aes(x = timepoint_label, y = AVR, group = PCode,
             colour = responder_group)) +
  geom_line(alpha = 0.2, linewidth = 0.4) +
  stat_summary(aes(group = responder_group, fill = responder_group),
               fun.data = mean_se, geom = "ribbon", alpha = 0.2, colour = NA) +
  stat_summary(aes(group = responder_group), fun = mean,
               geom = "line", linewidth = 2) +
  stat_summary(aes(group = responder_group), fun = mean,
               geom = "point", size = 3.5) +
  scale_colour_manual(
    values = c("Poor_Responder" = "#C0392B", "Responder" = "#2471A3"),
    name   = "Genetic Group",
    labels = c("Poor Responder", "Responder")) +
  scale_fill_manual(
    values = c("Poor_Responder" = "#C0392B", "Responder" = "#2471A3"),
    guide  = "none") +
  labs(title    = "Arteriole-to-Venule Ratio (AVR) Trajectories by Genetic Group",
       subtitle = "Shaded band = mean ± SE | Thin lines = individuals",
       x = "Timepoint", y = "AVR") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave("figures/plot_AVR_trajectories.png", p_avr, width = 9, height = 6, dpi = 300)
cat("Saved: figures/plot_AVR_trajectories.png\n")

# ── C5. Glu298Asp trajectory plots ───────────────────────────────────────────

master_glu <- master %>%
  dplyr::filter(!is.na(Glu298Asp), !is.na(timepoint_label)) %>%
  mutate(
    Glu298Asp_label = factor(
      case_when(
        Glu298Asp == 2 ~ "Glu/Glu (protective)",
        Glu298Asp == 3 ~ "Glu/Asp (heterozygous)",
        Glu298Asp == 4 ~ "Asp/Asp (risk)"
      ),
      levels = c("Glu/Glu (protective)", "Glu/Asp (heterozygous)", "Asp/Asp (risk)")
    )
  )

glu_colours <- c("Glu/Glu (protective)"   = "#2471A3",
                 "Glu/Asp (heterozygous)" = "#F39C12",
                 "Asp/Asp (risk)"         = "#C0392B")

# Glu298Asp CRAE
p_glu_crae <- master_glu %>%
  ggplot(aes(x = timepoint_label, y = CRAE,
             group = PCode, colour = Glu298Asp_label)) +
  geom_line(alpha = 0.2, linewidth = 0.4) +
  stat_summary(aes(group = Glu298Asp_label, fill = Glu298Asp_label),
               fun.data = mean_se, geom = "ribbon", alpha = 0.15, colour = NA) +
  stat_summary(aes(group = Glu298Asp_label), fun = mean,
               geom = "line", linewidth = 2) +
  stat_summary(aes(group = Glu298Asp_label), fun = mean,
               geom = "point", size = 3.5) +
  scale_colour_manual(values = glu_colours, name = "Genotype") +
  scale_fill_manual(values = glu_colours, guide = "none") +
  labs(title    = "CRAE Trajectories by Glu298Asp Genotype",
       subtitle = "Shaded band = mean ± SE | Thin lines = individuals",
       x = "Timepoint", y = "CRAE (µm)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

# Glu298Asp CRVE
p_glu_crve <- master_glu %>%
  ggplot(aes(x = timepoint_label, y = CRVE,
             group = PCode, colour = Glu298Asp_label)) +
  geom_line(alpha = 0.2, linewidth = 0.4) +
  stat_summary(aes(group = Glu298Asp_label, fill = Glu298Asp_label),
               fun.data = mean_se, geom = "ribbon", alpha = 0.15, colour = NA) +
  stat_summary(aes(group = Glu298Asp_label), fun = mean,
               geom = "line", linewidth = 2) +
  stat_summary(aes(group = Glu298Asp_label), fun = mean,
               geom = "point", size = 3.5) +
  scale_colour_manual(values = glu_colours, name = "Genotype") +
  scale_fill_manual(values = glu_colours, guide = "none") +
  labs(title    = "CRVE Trajectories by Glu298Asp Genotype",
       subtitle = "Shaded band = mean ± SE | Thin lines = individuals",
       x = "Timepoint", y = "CRVE (µm)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave("figures/plot_Glu298Asp_CRAE_trajectories.png",
       p_glu_crae, width = 9, height = 6, dpi = 300)
ggsave("figures/plot_Glu298Asp_CRVE_trajectories.png",
       p_glu_crve, width = 9, height = 6, dpi = 300)
cat("Saved: Glu298Asp trajectory plots\n")

# ── C6. COMT trajectory plots ─────────────────────────────────────────────────

master_comt <- master %>%
  dplyr::filter(!is.na(COMT), !is.na(timepoint_label)) %>%
  mutate(
    COMT_label = factor(
      case_when(
        COMT == 2 ~ "Met/Met (protective)",
        COMT == 3 ~ "Val/Met (intermediate)",
        COMT == 4 ~ "Val/Val (risk)"
      ),
      levels = c("Met/Met (protective)", "Val/Met (intermediate)", "Val/Val (risk)")
    )
  )

# Compute n per group for labels
comt_n <- master_comt %>%
  dplyr::filter(timepoint_label == "Baseline") %>%
  count(COMT_label) %>%
  mutate(label = paste0("n=", n))

comt_colours <- c("Met/Met (protective)"   = "#2471A3",
                  "Val/Met (intermediate)" = "#F39C12",
                  "Val/Val (risk)"         = "#C0392B")

make_comt_plot <- function(outcome, ylab, title) {
  # Get baseline means for label placement
  baseline_means <- master_comt %>%
    dplyr::filter(timepoint_label == "Baseline") %>%
    group_by(COMT_label) %>%
    summarise(mean_val = mean(.data[[outcome]], na.rm = TRUE)) %>%
    left_join(comt_n, by = "COMT_label")
  
  p <- master_comt %>%
    ggplot(aes(x = timepoint_label, y = .data[[outcome]],
               group = PCode, colour = COMT_label)) +
    geom_line(alpha = 0.2, linewidth = 0.4) +
    stat_summary(aes(group = COMT_label, fill = COMT_label),
                 fun.data = mean_se, geom = "ribbon", alpha = 0.15, colour = NA) +
    stat_summary(aes(group = COMT_label), fun = mean,
                 geom = "line", linewidth = 2) +
    stat_summary(aes(group = COMT_label), fun = mean,
                 geom = "point", size = 3.5) +
    geom_text(data = baseline_means,
              aes(x = "Baseline", y = mean_val, label = label,
                  colour = COMT_label),
              hjust = 1.3, size = 3.5, inherit.aes = FALSE) +
    scale_colour_manual(values = comt_colours, name = "Genotype") +
    scale_fill_manual(values = comt_colours, guide = "none") +
    labs(title    = title,
         subtitle = "Blue=protective | Orange=heterozygous | Red=risk",
         x = "Timepoint", y = ylab) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  return(p)
}

p_comt_crae <- make_comt_plot("CRAE", "CRAE (µm)", "COMT and CRAE Trajectories")
p_comt_crve <- make_comt_plot("CRVE", "CRVE (µm)", "COMT and CRVE Trajectories")

ggsave("figures/plot_COMT_CRAE_trajectories.png",
       p_comt_crae, width = 9, height = 6, dpi = 300)
ggsave("figures/plot_COMT_CRVE_trajectories.png",
       p_comt_crve, width = 9, height = 6, dpi = 300)
cat("Saved: COMT trajectory plots\n")

cat("\n03_clinical_models.R complete\n")


# Summary table: mean vessel calibre per timepoint per responder group
summary_by_group <- master %>%
  filter(!is.na(responder_group)) %>%
  group_by(PERIOD, responder_group) %>%
  summarise(
    n = n(),
    mean_CRAE = round(mean(CRAE, na.rm = TRUE), 2),
    sd_CRAE = round(sd(CRAE, na.rm = TRUE), 2),
    mean_CRVE = round(mean(CRVE, na.rm = TRUE), 2),
    sd_CRVE = round(sd(CRVE, na.rm = TRUE), 2),
    mean_AVR = round(mean(AVR, na.rm = TRUE), 3),
    sd_AVR = round(sd(AVR, na.rm = TRUE), 3),
    .groups = "drop"
  )

print(summary_by_group)

write.csv(summary_by_group,
          "summary_vessel_by_group_timepoint.csv",
          row.names = FALSE)