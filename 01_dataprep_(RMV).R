# =============================================================================
# 01_data_preparation.R
#
# Load and clean the masterfile, compute genetic group scores,
# classify participants as Responders / Poor Responders.
#
# Input:   Masterfile xlsx (file chooser dialog will open)
# Saves:   rds_objects/master.rds
# Output:  results/genetic_score_summary.csv
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

# Create output directories
dir.create("rds_objects", showWarnings = FALSE)
dir.create("results",     showWarnings = FALSE)
dir.create("figures",     showWarnings = FALSE)

# ── A1. Load masterfile ───────────────────────────────────────────────────────

master <- read_excel(file.choose())

cat("Masterfile:", nrow(master), "rows x", ncol(master), "cols\n")

# ── A2. Clean and type-convert ────────────────────────────────────────────────

master <- master %>%
  mutate(across(everything(), ~ na_if(as.character(.), " "))) %>%
  type_convert() %>%
  mutate(
    PCode     = as.factor(PCode),
    Sex       = as.factor(Sex),
    Season    = as.factor(Season),
    PERIOD    = relevel(as.factor(PERIOD), ref = "1"),
    across(c(DI, CRAE, CRVE, AVR, Age, BMI, WBC_Count,
             Lutein, Lycopene, bCARO, aCARO,
             VitCcorr, TEAC), as.numeric)
  ) %>%
  mutate(
    Age_z = as.numeric(scale(Age)),
    BMI_z = as.numeric(scale(BMI)),
    WBC_z = as.numeric(scale(WBC_Count)),
    timepoint = if_else(PERIOD == 1 & DI == 0,
                        "Baseline", "Post_Intervention") %>%
      factor(levels = c("Baseline", "Post_Intervention")),
    timepoint_label = case_when(
      timepoint == "Baseline" ~ "Baseline",
      PERIOD == "2"           ~ "Post 1",
      PERIOD == "3"           ~ "Post 2",
      TRUE ~ paste("Post", as.character(PERIOD))
    ) %>% factor(levels = c("Baseline", "Post 1", "Post 2"))
  )

cat("Timepoint split:\n")
print(table(master$timepoint))

# ── A3. Genetic scoring ───────────────────────────────────────────────────────

master <- master %>%
  mutate(
    # Group 1: Phase II Detoxification
    GSTM1_s = case_when(GSTM1==1~1, GSTM1==0~0, TRUE~NA_real_),
    GSTT1_s = case_when(GSTT1==1~1, GSTT1==0~0, TRUE~NA_real_),
    GSTP1_s = case_when(GSTP1==2~1, GSTP1==3~0.5, GSTP1==4~0, TRUE~NA_real_),
    NQO1_s  = case_when(NQO1==2~1,  NQO1==3~0.5,  NQO1==4~0,  TRUE~NA_real_),
    group_detox = rowMeans(cbind(GSTM1_s,GSTT1_s,GSTP1_s,NQO1_s), na.rm=TRUE),
    
    # Group 2: Methylation
    COMT_s  = case_when(COMT==2~1,  COMT==3~0.5,  COMT==4~0,  TRUE~NA_real_),
    MTHFR_s = case_when(MTHFR==2~1, MTHFR==3~0.5, MTHFR==4~0, TRUE~NA_real_),
    group_methylation = rowMeans(cbind(COMT_s,MTHFR_s), na.rm=TRUE),
    
    # Group 3: Antioxidant + Vascular
    CAT1_s      = case_when(CAT1==2~1,      CAT1==3~0.5,      CAT1==4~0,      TRUE~NA_real_),
    Glu298Asp_s = case_when(Glu298Asp==2~1, Glu298Asp==3~0.5, Glu298Asp==4~0, TRUE~NA_real_),
    XRCC1_s     = case_when(XRCC1==2~1,     XRCC1==3~0.5,     XRCC1==4~0,     TRUE~NA_real_),
    group_antioxidant = rowMeans(cbind(CAT1_s,Glu298Asp_s,XRCC1_s), na.rm=TRUE),
    
    # Group 4: Carotenoid
    BCMO1_s   = case_when(BCMO1==2~1,   BCMO1==3~0.5,   BCMO1==4~0,   TRUE~NA_real_),
    SLC23A1_s = case_when(SLC23A1==2~1, SLC23A1==3~0.5, SLC23A1==4~0, TRUE~NA_real_),
    ZBED3_s   = case_when(ZBED3==2~1,   ZBED3==3~0.5,   ZBED3==4~0,   TRUE~NA_real_),
    group_carotenoid = rowMeans(cbind(BCMO1_s,SLC23A1_s,ZBED3_s), na.rm=TRUE),
    
    # Group 5: Metabolic/CVD
    APOE_s   = case_when(APOE==2~1,   APOE==3~0.5,   APOE==4~0,   TRUE~NA_real_),
    HNF1A_s  = case_when(HNF1A==2~1,  HNF1A==3~0.5,  HNF1A==4~0,  TRUE~NA_real_),
    TCF7L2_s = case_when(TCF7L2==2~1, TCF7L2==3~0.5, TCF7L2==4~0, TRUE~NA_real_),
    group_metabolic = rowMeans(cbind(APOE_s,HNF1A_s,TCF7L2_s), na.rm=TRUE),
    
    overall_score = rowMeans(cbind(group_detox, group_methylation,
                                   group_antioxidant, group_carotenoid,
                                   group_metabolic), na.rm=TRUE)
  )

# ── A4. Responder classification ──────────────────────────────────────────────

median_score <- master %>%
  distinct(PCode, overall_score) %>%
  pull(overall_score) %>%
  median(na.rm=TRUE)

master <- master %>%
  mutate(responder_group = case_when(
    overall_score >  median_score ~ "Responder",
    overall_score <= median_score ~ "Poor_Responder",
    TRUE ~ NA_character_
  ) %>% factor(levels = c("Poor_Responder","Responder")))

cat("Median score:", round(median_score, 3), "\n")
print(master %>% distinct(PCode, responder_group) %>% count(responder_group))

# ── A5. Save ──────────────────────────────────────────────────────────────────

saveRDS(master, file = "rds_objects/master.rds")
cat("Saved: rds_objects/master.rds\n")

write.csv(
  master %>% distinct(PCode, overall_score, responder_group,
                      group_detox, group_methylation, group_antioxidant,
                      group_carotenoid, group_metabolic),
  "results/genetic_score_summary.csv",
  row.names = FALSE
)
cat("Saved: results/genetic_score_summary.csv\n")

cat("\n01_data_preparation.R complete\n")
cat("Participants:", nlevels(master$PCode), "\n")
cat("Observations:", nrow(master), "\n")