# =============================================================================
# 07_BOPs_(RMV).R
# MiBLEND Thesis — Sofia M. Prini
# Maastricht University | 2025-2026
#
# Purpose: Complete Cytoscape script for all 5 BOPs and physiological map
#          Single merged script — replaces ALL_BOPs_complete.R
#
# BOPs:
#   BOP1 — WP1995: Effects of Nitric Oxide (NOS3/Glu298Asp)
#   BOP2 — WP5603: Caffeine in Myocytes (MYH11/NOTCH3/Ca2+)
#   BOP3 — WP1742: TP53 Network (GADD45A/FOXO3/MYC)
#   BOP4 — Custom: WNK1-OXSR1 kinase cascade
#   BOP5 — WP231 + WP3932: Immune/Inflammatory (JAK2/DAPP1/MAP3K8/TAB3)
#
# Verified analysis values (post BMI correction):
#   CRAE Post1: +2.97um, p 8.06e-05
#   CRAE Post2: +4.31um, p 1.19e-07
#   AVR  Post1: +0.024,  p 2.28e-05
#   AVR  Post2: +0.022,  p 1.45e-04
#   Glu298Asp AVR: OR 2.14, p 0.012, adj.p 0.177 (nominal only)
#   NQO1 CRAE: OR 0.165, p 0.017, adj.p 0.100 (nominal only)
#
# Colour scheme (consistent across all BOPs):
#   Red    (#C0392B) = FDR-significant AVR correlate (this study)
#   Gold   (#F1C40F) = literature-supported node
#   Teal   (#148F77) = downstream physiological outcome
#   Pink   (#FADBD8) = GSEA pathway finding
#   White             = biological context, not directly tested
#   Blue   (#2471A3) = stress/inflammation BOP nodes
#
# IMPORTANT: No equals signs in addAnnotationText() strings
#   Cytoscape parses key=value as JSON causing duplicate key errors
#   Use colons or parentheses: r -0.350, adj.p 4.2e-6
#
# Run order:
#   1. Open Cytoscape 3.10.4
#   2. Source this script — all networks imported automatically
#
# Output: /Users/sofi/Desktop/BOP_figures/
#   BOP1_NOS_pathway_final.png
#   BOP2_myocytes_final.png
#   BOP3_TP53_final.png
#   BOP4_WNK1_OXSR1_final.png
#   BOP5a_TNFalpha_final.png
#   BOP5b_PI3K_Akt_final.png
#   physiological_map_final.png
# =============================================================================

library(RCy3)
library(tidyverse)

cytoscapePing()
cat("Connected to Cytoscape\n")
cat("Available networks:\n")
print(getNetworkList())

out_dir <- "/Users/sofi/Desktop/BOP_figures/"
dir.create(out_dir, showWarnings = FALSE)

# Network names — update if getNetworkList() shows different names
BOP1_name <- "Effects of nitric oxide - Homo sapiens"
BOP2_name <- "Caffeine in myocytes - Homo sapiens"
BOP3_name <- "TP53 network - Homo sapiens"
BOP4_name <- "WNK1-OXSR1 signalling - Custom BOP4"
BOP5_TNF  <- "TNF-alpha signaling - Homo sapiens"
BOP5_PI3K <- "Focal adhesion: PI3K-Akt-mTOR-signaling - Homo sapiens"

# Helper: delete network if exists
delete_if_exists <- function(name) {
  nets <- getNetworkList()
  exists <- if (is.data.frame(nets)) name %in% nets$name else
    name %in% as.character(nets)
  if (exists) { deleteNetwork(name); cat("Deleted:", name, "\n") }
}

# =============================================================================
# IMPORT ALL WIKIPATHWAYS NETWORKS
# =============================================================================

cat("\nImporting all WikiPathways networks...\n")

cat("Importing BOP1: WP1995 - Effects of Nitric Oxide...\n")
commandsRun("wikipathways import-as-pathway id=WP1995")
Sys.sleep(4)

cat("Importing BOP2: WP5603 - Caffeine in Myocytes...\n")
commandsRun("wikipathways import-as-pathway id=WP5603")
Sys.sleep(4)

cat("Importing BOP3: WP1742 - TP53 Network...\n")
commandsRun("wikipathways import-as-pathway id=WP1742")
Sys.sleep(4)

cat("Importing BOP5a: WP231 - TNF-alpha Signalling...\n")
commandsRun("wikipathways import-as-pathway id=WP231")
Sys.sleep(4)

cat("Importing BOP5b: WP3932 - PI3K-Akt-mTOR Signalling...\n")
commandsRun("wikipathways import-as-pathway id=WP3932")
Sys.sleep(4)

cat("All WikiPathways networks imported\n")
cat("Available networks:\n")
print(getNetworkList())

# =============================================================================
# BOP 1 — WP1995: Effects of Nitric Oxide
# Glu298Asp: OR 2.14, p 0.012, adj.p 0.177 (nominal only)
# =============================================================================

setCurrentNetwork(BOP1_name)

setNodeColorBypass("NOS3", "#C0392B")
setNodeBorderColorBypass("NOS3", "#922B21")
setNodeBorderWidthBypass("NOS3", 4)
setNodeLabelColorBypass("NOS3", "#FFFFFF")

setNodeColorBypass("Nitric oxide", "#F1C40F")
setNodeBorderColorBypass("Nitric oxide", "#B7950B")
setNodeBorderWidthBypass("Nitric oxide", 3)
setNodeColorBypass("L-Arginine", "#F1C40F")
setNodeBorderColorBypass("L-Arginine", "#B7950B")
setNodeColorBypass("Citrulline", "#F1C40F")
setNodeBorderColorBypass("Citrulline", "#B7950B")

nos3_pos <- getNodePosition("NOS3", network = BOP1_name)
no_pos   <- getNodePosition("Nitric oxide", network = BOP1_name)

addAnnotationText(
  text     = "Glu298Asp (rs1799983)\nOR 2.14, p 0.012, adj.p 0.177 (nominal)\nAsp allele nominally associated\nwith higher AVR group\n(unexpected direction — see Discussion)",
  x.pos    = nos3_pos$x_location - 250,
  y.pos    = nos3_pos$y_location - 100,
  fontSize = 13, color = "#922B21",
  network  = BOP1_name)

addAnnotationText(
  text     = "NO to cGMP signalling\nto vasodilation to CRAE increase\nCRAE Post1 +2.97um (p 8.06e-05)\nCRAE Post2 +4.31um (p 1.19e-07)\n(Forstermann & Munzel 2006)",
  x.pos    = no_pos$x_location + 200,
  y.pos    = no_pos$y_location + 20,
  fontSize = 13, color = "#B7950B",
  network  = BOP1_name)

addAnnotationText(
  text     = "Cross-BOP connections:\nTNFa (BOP5) uncouples eNOS\nleading to reduced NO\nNO modulates WNK1 activity (BOP4)\n\nLEGEND\nRed  : Observed in this study\nGold : Literature supported\nWhite: Biological context",
  x.pos    = nos3_pos$x_location - 400,
  y.pos    = nos3_pos$y_location + 250,
  fontSize = 12, color = "#2C3E50",
  network  = BOP1_name)

exportImage(paste0(out_dir, "BOP1_NOS_pathway.png"),
            type = "PNG", resolution = 300, zoom = 250)
message("BOP1 exported.")

# =============================================================================
# BOP 2 — WP5603: Caffeine in Myocytes
# MYH11: r +0.213, adj.p 0.016 (FDR-significant, exploratory)
# Ca2+: KEGG adj.p 0.0015, GO adj.p 0.015
# Bridge: WNK->OSR1->TRPV4->Ca2+ (Garrud 2024)
# =============================================================================

setCurrentNetwork(BOP2_name)

setNodeColorBypass("MYH11", "#C0392B")
setNodeBorderColorBypass("MYH11", "#922B21")
setNodeBorderWidthBypass("MYH11", 4)
setNodeLabelColorBypass("MYH11", "#FFFFFF")

setNodeColorBypass("Ca\u00b2\u207a", "#C0392B")
setNodeBorderColorBypass("Ca\u00b2\u207a", "#922B21")
setNodeBorderWidthBypass("Ca\u00b2\u207a", 4)
setNodeLabelColorBypass("Ca\u00b2\u207a", "#FFFFFF")

for (g in c("RYR1","RYR2","RYR3")) {
  setNodeColorBypass(g, "#F1C40F")
  setNodeBorderColorBypass(g, "#B7950B")
  setNodeBorderWidthBypass(g, 3)
}
setNodeBorderColorBypass("Caffeine", "#B7950B")
setNodeBorderWidthBypass("Caffeine", 3)

ca_pos  <- tryCatch(getNodePosition("Ca\u00b2\u207a", network = BOP2_name),
                    error = function(e) list(x_location=195, y_location=252))
myh_pos <- tryCatch(getNodePosition("MYH11", network = BOP2_name),
                    error = function(e) list(x_location=781, y_location=582))

addAnnotationText(
  text     = "Ca2+ signalling KEGG (adj.p 0.0015)\nGO calcium adhesion (adj.p 0.015)\nsmooth muscle relaxation\narteriolar vasodilation to CRAE increase",
  x.pos    = ca_pos$x_location + 120,
  y.pos    = ca_pos$y_location - 140,
  fontSize = 13, color = "#922B21",
  network  = BOP2_name)

addAnnotationText(
  text     = "MYH11: AVR r +0.213, adj.p 0.016\nSmooth muscle myosin heavy chain\nDifferentiated contractile SMC phenotype\nRay et al. Sci Rep 2020\n\nNOTCH3: r +0.317 (baseline sensitivity)\nSelectively expressed in retinal mural cells\nRequired for mural cell recruitment\nLiu et al. Circ Res 2010",
  x.pos    = myh_pos$x_location - 60,
  y.pos    = myh_pos$y_location - 100,
  fontSize = 13, color = "#922B21",
  network  = BOP2_name)

addAnnotationText(
  text     = "BOP4 bridge (Garrud et al. PNAS 2024):\nWNK to OSR1 to TRPV4 to Ca2+\nto vasodilation in ECs\n\nBOP5 bridge:\nJAK2 to PI3K to Akt to eNOS\nto NO to reduced Ca2+ contraction\n\nLEGEND\nRed  : Observed in this study\nGold : Literature supported\nWhite: Biological context\nTeal : Cross-BOP connection",
  x.pos    = ca_pos$x_location - 300,
  y.pos    = ca_pos$y_location + 200,
  fontSize = 12, color = "#2C3E50",
  network  = BOP2_name)

exportImage(paste0(out_dir, "BOP2_myocytes.png"),
            type = "PNG", resolution = 300, zoom = 250)
message("BOP2 exported.")

# =============================================================================
# BOP 3 — WP1742: TP53 Network
# GADD45A: r -0.350, adj.p 4.2e-6 (TOP-RANKED AVR correlate #1)
# FOXO3:   r -0.209, adj.p 0.019 (exploratory)
# MXI1:    r -0.225, adj.p 0.011 (exploratory)
# HIPK1:   r -0.282, adj.p 0.001 (top 50)
# CDKN2D:  r -0.270, adj.p 0.002 (top 50)
# C5:      r -0.288, adj.p 0.001 (top 50, SASP connection to BOP5)
# =============================================================================

setCurrentNetwork(BOP3_name)

setNodeColorBypass("GADD45A", "#C0392B")
setNodeBorderColorBypass("GADD45A", "#922B21")
setNodeBorderWidthBypass("GADD45A", 4)
setNodeLabelColorBypass("GADD45A", "#FFFFFF")

setNodeColorBypass("MYC", "#C0392B")
setNodeBorderColorBypass("MYC", "#922B21")
setNodeBorderWidthBypass("MYC", 4)
setNodeLabelColorBypass("MYC", "#FFFFFF")

tryCatch({
  setNodeColorBypass("TP53", "#F1C40F")
  setNodeBorderColorBypass("TP53", "#B7950B")
  setNodeBorderWidthBypass("TP53", 3)
}, error = function(e) cat("TP53 node not found\n"))

gadd_pos <- getNodePosition("GADD45A", network = BOP3_name)
myc_pos  <- getNodePosition("MYC",     network = BOP3_name)

addAnnotationText(
  text     = "GADD45A: AVR r -0.350, adj.p 4.2e-6\nTOP-RANKED AVR correlate (#1 of top 50)\nFOXO3a to GADD45A\n(Tran et al. Science 2002)\nHIPK1 (r -0.282) and CDKN2D (r -0.270)\nalso in top 50 named correlates",
  x.pos    = gadd_pos$x_location - 220,
  y.pos    = gadd_pos$y_location + 80,
  fontSize = 13, color = "#922B21",
  network  = BOP3_name)

addAnnotationText(
  text     = "FOXO3a to MXI1 to represses MYC\nFOXO3: r -0.209, adj.p 0.019\nMXI1:  r -0.225, adj.p 0.011\n(both exploratory — outside top 50)\nDelpuech et al. MCB 2007",
  x.pos    = myc_pos$x_location + 120,
  y.pos    = myc_pos$y_location - 100,
  fontSize = 13, color = "#922B21",
  network  = BOP3_name)

addAnnotationText(
  text     = "Cross-BOP connection:\nGADD45A drives stress-induced senescence\nSASP connects BOP3 to BOP5\nC5 (r -0.288) and CDKN2D (r -0.270)\nappear in both cell cycle (WP179)\nand SASP (WP3391) pathways\n\nLEGEND\nRed  : Observed in this study\nGold : Literature supported\nWhite: Biological context\n\nExploratory: FDR-significant\nbut outside top 50 named correlates",
  x.pos    = gadd_pos$x_location - 380,
  y.pos    = gadd_pos$y_location + 280,
  fontSize = 12, color = "#2C3E50",
  network  = BOP3_name)

exportImage(paste0(out_dir, "BOP3_TP53.png"),
            type = "PNG", resolution = 300, zoom = 250)
message("BOP3 exported.")

# =============================================================================
# BOP 4 — CUSTOM: WNK1-OXSR1 Signalling
# WNK1:  r -0.216, adj.p 0.015 (exploratory)
# OXSR1: r -0.234, adj.p 0.009 (exploratory)
# Chloride transport GO adj.p 0.020
# Bridge: Garrud 2024 WNK->OSR1->TRPV4->Ca2+ links BOP4 to BOP2
# =============================================================================

delete_if_exists(BOP4_name)

bop4_nodes <- data.frame(
  id    = c("WNK1","OXSR1","STK39","SLC12A2","SLC12A7",
            "Hyperosmotic_stress","Chloride_sensing",
            "Ion_homeostasis","Blood_pressure","AVR_outcome"),
  label = c("WNK1 (r -0.216, exploratory)",
            "OXSR1 (r -0.234, exploratory)",
            "STK39 (SPAK)",
            "SLC12A2 (NKCC1)",
            "SLC12A7 (KCC4)",
            "Hyperosmotic stress",
            "Chloride sensing\nGO adj.p 0.020",
            "Ion homeostasis & cell volume",
            "Blood pressure regulation",
            "AVR (this study)"),
  stringsAsFactors = FALSE)

bop4_edges <- data.frame(
  source      = c("Hyperosmotic_stress","Chloride_sensing",
                  "WNK1","WNK1",
                  "OXSR1","OXSR1","STK39","STK39",
                  "SLC12A2","SLC12A7",
                  "Ion_homeostasis","Blood_pressure"),
  target      = c("WNK1","WNK1",
                  "OXSR1","STK39",
                  "SLC12A2","SLC12A7","SLC12A2","SLC12A7",
                  "Ion_homeostasis","Ion_homeostasis",
                  "Blood_pressure","AVR_outcome"),
  interaction = c("activates","activates",
                  "phosphorylates","phosphorylates",
                  "activates","inhibits","activates","inhibits",
                  "regulates","regulates",
                  "determines","reflected in"),
  stringsAsFactors = FALSE)

createNetworkFromDataFrames(
  nodes = bop4_nodes, edges = bop4_edges,
  title = BOP4_name, collection = "MiBLEND BOPs")

for (g in c("STK39","SLC12A2","SLC12A7")) {
  setNodeColorBypass(g, "#F1C40F", network = BOP4_name)
  setNodeBorderColorBypass(g, "#B7950B", network = BOP4_name)
  setNodeBorderWidthBypass(g, 3, network = BOP4_name)
}

setNodeColorBypass("Chloride_sensing", "#FADBD8", network = BOP4_name)
setNodeBorderColorBypass("Chloride_sensing", "#C0392B", network = BOP4_name)
setNodeBorderWidthBypass("Chloride_sensing", 3, network = BOP4_name)

for (g in c("Ion_homeostasis","Blood_pressure","AVR_outcome")) {
  setNodeColorBypass(g, "#148F77", network = BOP4_name)
  setNodeBorderColorBypass(g, "#0E6655", network = BOP4_name)
  setNodeLabelColorBypass(g, "#FFFFFF", network = BOP4_name)
}

setNodePropertyBypass("Hyperosmotic_stress", list("#FFFFFF"),
                      "NODE_FILL_COLOR", network = BOP4_name)
setNodePropertyBypass("Hyperosmotic_stress", list("#808080"),
                      "NODE_BORDER_PAINT", network = BOP4_name)

layoutNetwork("hierarchical", network = BOP4_name)

# Sizes — individual calls to avoid vector bypass errors
for (g in c("WNK1","OXSR1")) {
  setNodeWidthBypass(g, 200, network = BOP4_name)
  setNodeHeightBypass(g, 60, network = BOP4_name)
}
for (g in c("STK39","SLC12A2","SLC12A7")) {
  setNodeWidthBypass(g, 170, network = BOP4_name)
  setNodeHeightBypass(g, 55, network = BOP4_name)
}
for (g in c("Ion_homeostasis","Blood_pressure","AVR_outcome",
            "Hyperosmotic_stress","Chloride_sensing")) {
  setNodeWidthBypass(g, 220, network = BOP4_name)
  setNodeHeightBypass(g, 60, network = BOP4_name)
}

# Labels
setNodeLabelBypass("WNK1",               "WNK1 (r -0.216, exploratory)", network = BOP4_name)
setNodeLabelBypass("OXSR1",              "OXSR1 (r -0.234, exploratory)", network = BOP4_name)
setNodeLabelBypass("STK39",              "STK39 (SPAK)", network = BOP4_name)
setNodeLabelBypass("SLC12A2",            "SLC12A2 (NKCC1)", network = BOP4_name)
setNodeLabelBypass("SLC12A7",            "SLC12A7 (KCC4)", network = BOP4_name)
setNodeLabelBypass("Hyperosmotic_stress","Hyperosmotic stress", network = BOP4_name)
setNodeLabelBypass("Chloride_sensing",   "Chloride sensing\nGO adj.p 0.020", network = BOP4_name)
setNodeLabelBypass("Ion_homeostasis",    "Ion homeostasis & cell volume", network = BOP4_name)
setNodeLabelBypass("Blood_pressure",     "Blood pressure regulation", network = BOP4_name)
setNodeLabelBypass("AVR_outcome",        "AVR (this study)", network = BOP4_name)

# Colours — individual calls to avoid vector bypass errors
for (g in c("WNK1","OXSR1")) {
  setNodeColorBypass(g, "#C0392B", network = BOP4_name)
  setNodeBorderColorBypass(g, "#922B21", network = BOP4_name)
  setNodeBorderWidthBypass(g, 4, network = BOP4_name)
  setNodeLabelColorBypass(g, "#FFFFFF", network = BOP4_name)
}

for (g in c("STK39","SLC12A2","SLC12A7")) {
  setNodeColorBypass(g, "#F1C40F", network = BOP4_name)
  setNodeBorderColorBypass(g, "#B7950B", network = BOP4_name)
  setNodeBorderWidthBypass(g, 3, network = BOP4_name)
}

setNodeColorBypass("Chloride_sensing", "#FADBD8", network = BOP4_name)
setNodeBorderColorBypass("Chloride_sensing", "#C0392B", network = BOP4_name)
setNodeBorderWidthBypass("Chloride_sensing", 3, network = BOP4_name)

for (g in c("Ion_homeostasis","Blood_pressure","AVR_outcome")) {
  setNodeColorBypass(g, "#148F77", network = BOP4_name)
  setNodeBorderColorBypass(g, "#0E6655", network = BOP4_name)
  setNodeLabelColorBypass(g, "#FFFFFF", network = BOP4_name)
}

setNodePropertyBypass("Hyperosmotic_stress", list("#FFFFFF"),
                      "NODE_FILL_COLOR", network = BOP4_name)
setNodePropertyBypass("Hyperosmotic_stress", list("#808080"),
                      "NODE_BORDER_PAINT", network = BOP4_name)

bop4_all_edges <- getAllEdges(network = BOP4_name)
setEdgeColorBypass(bop4_all_edges, rep("#000000", length(bop4_all_edges)))
inh <- bop4_all_edges[grep("inhibits", bop4_all_edges)]
if (length(inh) > 0) {
  setEdgeColorBypass(inh, rep("#C0392B", length(inh)))
  setEdgeLineStyleBypass(inh, rep("LONG_DASH", length(inh)))
}

layoutNetwork("hierarchical", network = BOP4_name)

addAnnotationText(
  text     = "LEGEND\nRed  : FDR-significant AVR correlate\nGold : Literature supported\nTeal : Downstream outcome\nWhite: Upstream context\nPink : GSEA finding\n\nSolid : activates or regulates\nDashed: inhibits\n\nExploratory: FDR-significant but\noutside top 50 named correlates",
  x.pos = -450, y.pos = 200,
  fontSize = 12, color = "#2C3E50",
  network = BOP4_name)

addAnnotationText(
  text     = "Cross-BOP bridge (Garrud et al. PNAS 2024):\nWNK to OSR1 to TRPV4 to Ca2+ to vasodilation\nLinks BOP4 to BOP2\n\nSources:\nRichardson & Alessi J Cell Sci 2008\nShekarabi et al. Cell Metab 2017\nBergaya et al. Hypertension 2011\nGarrud et al. PNAS 2024",
  x.pos = -450, y.pos = 500,
  fontSize = 11, color = "#717D7E",
  network = BOP4_name)

exportImage(paste0(out_dir, "BOP4_WNK1_OXSR1.png"),
            type = "PNG", resolution = 300, zoom = 250)
message("BOP4 exported.")

# =============================================================================
# BOP 5a — WP231: TNF-alpha Signalling
# TAB3:   r -0.280, adj.p 0.001
# MAP3K8: r -0.267, adj.p 0.002
# Cross-BOP: TNFa uncouples NOS3 (BOP1)
# =============================================================================

tryCatch({
  setCurrentNetwork(BOP5_TNF)
  
  for (g in c("TAB3","MAP3K8")) {
    tryCatch({
      setNodeColorBypass(g, "#C0392B")
      setNodeBorderColorBypass(g, "#922B21")
      setNodeBorderWidthBypass(g, 4)
      setNodeLabelColorBypass(g, "#FFFFFF")
    }, error = function(e) cat("Node not found:", g, "\n"))
  }
  
  for (g in c("TNF","TNFR1","TNFR2","NFKB1","RELA","TRAF2","RIPK1","CASP8")) {
    tryCatch({
      setNodeColorBypass(g, "#F1C40F")
      setNodeBorderColorBypass(g, "#B7950B")
      setNodeBorderWidthBypass(g, 3)
    }, error = function(e) NULL)
  }
  
  addAnnotationText(
    text     = "BOP5a: TNF-alpha signalling (WP231)\nTAB3:   AVR r -0.280, adj.p 0.001\nMAP3K8: AVR r -0.267, adj.p 0.002\nDAPP1:  AVR r -0.314 (ranked #2 overall)\nPI3K downstream adaptor in leukocytes\nRegulates microvascular permeability via ROS\nHao et al. Front Immunol 2020\nAll FDR-significant top-50 AVR correlates",
    x.pos = 50, y.pos = 50,
    fontSize = 13, color = "#922B21",
    network = BOP5_TNF)
  
  addAnnotationText(
    text     = "Cross-BOP connection to BOP1:\nTNFa activates NF-kB which suppresses NOS3\nleading to eNOS uncoupling and reduced NO\nSource: Forstermann & Munzel 2006\n\nCross-BOP connection to BOP5b:\nTNFa activates JAK2 signalling\n\nLEGEND\nRed  : FDR-significant AVR correlate\nGold : Literature supported\nWhite: Biological context",
    x.pos = 50, y.pos = 350,
    fontSize = 12, color = "#2C3E50",
    network = BOP5_TNF)
  
  exportImage(paste0(out_dir, "BOP5a_TNFalpha.png"),
              type = "PNG", resolution = 300, zoom = 250)
  message("BOP5a TNF-alpha exported.")
  
}, error = function(e) cat("BOP5a error:", e$message,
                           "\nCheck network name with getNetworkList()\n"))

# =============================================================================
# BOP 5b — WP3932: PI3K-Akt-mTOR Signalling
# JAK2:  r -0.253, adj.p 0.004
# RAB8A: r -0.270, adj.p 0.002
# GNG5:  r -0.269, adj.p 0.002
# RAB10: r -0.259, adj.p 0.003
# Cross-BOP: JAK2->PI3K->Akt->eNOS->NO bridges BOP5 to BOP1 and BOP2
# =============================================================================

tryCatch({
  setCurrentNetwork(BOP5_PI3K)
  
  for (g in c("JAK2","RAB8A","GNG5","RAB10")) {
    tryCatch({
      setNodeColorBypass(g, "#C0392B")
      setNodeBorderColorBypass(g, "#922B21")
      setNodeBorderWidthBypass(g, 4)
      setNodeLabelColorBypass(g, "#FFFFFF")
    }, error = function(e) cat("Node not found:", g, "\n"))
  }
  
  for (g in c("PIK3CA","PIK3CB","AKT1","AKT2","AKT3","MTOR","NOS3","PTEN")) {
    tryCatch({
      setNodeColorBypass(g, "#F1C40F")
      setNodeBorderColorBypass(g, "#B7950B")
      setNodeBorderWidthBypass(g, 3)
    }, error = function(e) NULL)
  }
  
  addAnnotationText(
    text     = "BOP5b: PI3K-Akt-mTOR signalling (WP3932)\nJAK2:  AVR r -0.253, adj.p 0.004\nRAB8A: AVR r -0.270, adj.p 0.002\nGNG5:  AVR r -0.269, adj.p 0.002\nRAB10: AVR r -0.259, adj.p 0.003\nAll FDR-significant top-50 AVR correlates",
    x.pos = 50, y.pos = 50,
    fontSize = 13, color = "#922B21",
    network = BOP5_PI3K)
  
  addAnnotationText(
    text     = "Cross-BOP connection to BOP1 and BOP2:\nJAK2 to PI3K to Akt to NOS3 phosphorylation\nto increased NO to vasodilation\nto reduced Ca2+ smooth muscle contraction\n\nLEGEND\nRed  : FDR-significant AVR correlate\nGold : Literature supported (PI3K pathway)\nWhite: Biological context",
    x.pos = 50, y.pos = 400,
    fontSize = 12, color = "#148F77",
    network = BOP5_PI3K)
  
  exportImage(paste0(out_dir, "BOP5b_PI3K_Akt.png"),
              type = "PNG", resolution = 300, zoom = 250)
  message("BOP5b PI3K-Akt exported.")
  
}, error = function(e) cat("BOP5b error:", e$message,
                           "\nCheck network name with getNetworkList()\n"))

# =============================================================================
# PHYSIOLOGICAL MAP — all 5 BOPs with inter-BOP gene connections
# =============================================================================

map_name <- "MiBLEND Physiological Map"
delete_if_exists(map_name)

map_nodes <- data.frame(
  id    = c(
    "Intervention",
    "BOP1","BOP2","BOP3","BOP4","BOP5",
    "NOS3","WNK1","OXSR1","Ca2plus","MYH11",
    "GADD45A","FOXO3","C5","CDKN2D",
    "JAK2","MAP3K8","TAB3","TNFa",
    "TRPV4","Chloride",
    "CRAE_out","AVR_out"),
  label = c(
    "Phytochemical\nintervention",
    "BOP1\nNOS3/eNOS\nWP1995",
    "BOP2\nCa2+/MYH11\nWP5603",
    "BOP3\nGADD45A/p53\nWP1742",
    "BOP4\nWNK1-OXSR1\nCustom",
    "BOP5\nImmune/Inflam.\nWP231/WP3932",
    "NOS3\n(Glu298Asp OR 2.14)",
    "WNK1\n(r -0.216)",
    "OXSR1\n(r -0.234)",
    "Ca2+\n(KEGG adj.p 0.0015)",
    "MYH11\n(r +0.213)",
    "GADD45A\n(r -0.350 #1)",
    "FOXO3\n(r -0.209)",
    "C5\n(r -0.288)",
    "CDKN2D\n(r -0.270)",
    "JAK2\n(r -0.253)",
    "MAP3K8\n(r -0.267)",
    "TAB3\n(r -0.280)",
    "TNF-alpha\n(literature)",
    "TRPV4\n(Garrud 2024)",
    "Chloride\n(GO adj.p 0.020)",
    "CRAE\nPost1 +2.97um, Post2 +4.31um",
    "AVR\nPost1 +0.024, Post2 +0.022"),
  stringsAsFactors = FALSE)

map_edges <- data.frame(
  source = c(
    "Intervention","Intervention","Intervention","Intervention","Intervention",
    "BOP1","BOP4","BOP4","BOP2","BOP2",
    "BOP3","BOP3","BOP5","BOP5","BOP5","BOP5",
    "NOS3","WNK1","OXSR1","TRPV4",
    "GADD45A","C5",
    "JAK2","TNFa",
    "Chloride",
    "Ca2plus","MYH11","NOS3","GADD45A","CRAE_out"),
  target = c(
    "BOP1","BOP2","BOP3","BOP4","BOP5",
    "NOS3","WNK1","OXSR1","Ca2plus","MYH11",
    "GADD45A","FOXO3","JAK2","MAP3K8","TAB3","TNFa",
    "WNK1","TRPV4","TRPV4","Ca2plus",
    "C5","CDKN2D",
    "Ca2plus","NOS3",
    "WNK1",
    "CRAE_out","CRAE_out","AVR_out","AVR_out","AVR_out"),
  interaction = c(
    "activates","activates","activates","activates","activates",
    "contains","contains","contains","contains","contains",
    "contains","contains","contains","contains","contains","contains",
    "NO_suppresses_WNK","WNK_activates","OXSR1_activates","TRPV4_releases",
    "SASP_releases","SASP_releases",
    "JAK2_PI3K_eNOS","TNFa_uncouples_eNOS",
    "chloride_activates",
    "relaxation","contraction","nominally_associated","stress_marker","increases"),
  stringsAsFactors = FALSE)

createNetworkFromDataFrames(
  nodes = map_nodes, edges = map_edges,
  title = map_name, collection = "MiBLEND BOPs")

layoutNetwork("force-directed", network = map_name)

# BOP nodes — vascular tone red
for (n in c("BOP1","BOP4","BOP2")) {
  setNodeColorBypass(n, "#C0392B", network = map_name)
  setNodeBorderColorBypass(n, "#922B21", network = map_name)
  setNodeBorderWidthBypass(n, 4, network = map_name)
  setNodeLabelColorBypass(n, "#FFFFFF", network = map_name)
  setNodeWidthBypass(n, 160, network = map_name)
  setNodeHeightBypass(n, 80, network = map_name)
}

# BOP nodes — stress/inflammation blue
for (n in c("BOP3","BOP5")) {
  setNodeColorBypass(n, "#2471A3", network = map_name)
  setNodeBorderColorBypass(n, "#1A5276", network = map_name)
  setNodeBorderWidthBypass(n, 4, network = map_name)
  setNodeLabelColorBypass(n, "#FFFFFF", network = map_name)
  setNodeWidthBypass(n, 160, network = map_name)
  setNodeHeightBypass(n, 80, network = map_name)
}

# Intervention — teal
setNodePropertyBypass("Intervention", list("#D5F5E3"),
                      "NODE_FILL_COLOR", network = map_name)
setNodeBorderColorBypass("Intervention", "#148F77", network = map_name)
setNodeBorderWidthBypass("Intervention", 3, network = map_name)
setNodeWidthBypass("Intervention", 160, network = map_name)
setNodeHeightBypass("Intervention", 70, network = map_name)

# Gene nodes — red
for (n in c("NOS3","WNK1","OXSR1","Ca2plus","MYH11",
            "GADD45A","FOXO3","C5","CDKN2D",
            "JAK2","MAP3K8","TAB3")) {
  setNodeColorBypass(n, "#C0392B", network = map_name)
  setNodeBorderColorBypass(n, "#922B21", network = map_name)
  setNodeBorderWidthBypass(n, 2, network = map_name)
  setNodeLabelColorBypass(n, "#FFFFFF", network = map_name)
  setNodeWidthBypass(n, 140, network = map_name)
  setNodeHeightBypass(n, 55, network = map_name)
}

# Gold — literature mediators
for (n in c("TNFa","TRPV4","Chloride")) {
  setNodeColorBypass(n, "#F1C40F", network = map_name)
  setNodeBorderColorBypass(n, "#B7950B", network = map_name)
  setNodeBorderWidthBypass(n, 2, network = map_name)
  setNodeWidthBypass(n, 140, network = map_name)
  setNodeHeightBypass(n, 55, network = map_name)
}

# Outcomes — teal
for (n in c("CRAE_out","AVR_out")) {
  setNodeColorBypass(n, "#148F77", network = map_name)
  setNodeBorderColorBypass(n, "#0E6655", network = map_name)
  setNodeBorderWidthBypass(n, 3, network = map_name)
  setNodeLabelColorBypass(n, "#FFFFFF", network = map_name)
  setNodeWidthBypass(n, 160, network = map_name)
  setNodeHeightBypass(n, 65, network = map_name)
}

# Edge styling
map_all_edges <- getAllEdges(network = map_name)
setEdgeColorBypass(map_all_edges, rep("#555555", length(map_all_edges)),
                   network = map_name)

garrud <- map_all_edges[grepl("TRPV4", map_all_edges)]
if (length(garrud) > 0) {
  setEdgeColorBypass(garrud, rep("#148F77", length(garrud)), network = map_name)
  setEdgeLineStyleBypass(garrud, rep("LONG_DASH", length(garrud)), network = map_name)
}

tfa <- map_all_edges[grepl("uncouples", map_all_edges)]
if (length(tfa) > 0) {
  setEdgeColorBypass(tfa, rep("#C0392B", length(tfa)), network = map_name)
  setEdgeLineStyleBypass(tfa, rep("LONG_DASH", length(tfa)), network = map_name)
}

jak <- map_all_edges[grepl("JAK2_PI3K", map_all_edges)]
if (length(jak) > 0) {
  setEdgeColorBypass(jak, rep("#2471A3", length(jak)), network = map_name)
  setEdgeLineStyleBypass(jak, rep("LONG_DASH", length(jak)), network = map_name)
}

sasp <- map_all_edges[grepl("SASP", map_all_edges)]
if (length(sasp) > 0) {
  setEdgeColorBypass(sasp, rep("#8E44AD", length(sasp)), network = map_name)
  setEdgeLineStyleBypass(sasp, rep("LONG_DASH", length(sasp)), network = map_name)
}

layoutNetwork("force-directed", network = map_name)

addAnnotationText(
  text     = "LEGEND\nRed (large)   : Vascular tone BOP\nBlue (large)  : Stress/inflammation BOP\nRed (small)   : FDR-significant AVR correlate gene\nGold          : Literature-supported mediator\nTeal          : Clinical outcome or intervention\n\nSolid grey    : BOP membership or primary pathway\nRed dashed    : TNFa inhibits NOS3 (BOP5->BOP1)\nTeal dashed   : Garrud 2024 bridge (OXSR1->TRPV4->Ca2+)\nBlue dashed   : JAK2->PI3K->eNOS (BOP5->BOP1/BOP2)\nPurple dashed : SASP senescence (BOP3->BOP5)",
  x.pos = 50, y.pos = 50,
  fontSize = 12, color = "#2C3E50",
  network = map_name)

addAnnotationText(
  text     = "Tip: Manually arrange BOP nodes into two-axis structure\nLeft axis  (red):  Intervention -> BOP1 -> BOP4 -> BOP2 -> outcomes\nRight axis (blue): Intervention -> BOP3 -> BOP5 -> BOP2\nGene nodes positioned between connected BOPs",
  x.pos = 50, y.pos = 700,
  fontSize = 11, color = "#717D7E",
  network = map_name)

exportImage(paste0(out_dir, "physiological_map.png"),
            type = "PNG", resolution = 300, zoom = 200)
message("Physiological map exported.")

# =============================================================================
# CROP WHITESPACE — all BOP figures
# =============================================================================

crop_whitespace <- function(in_path, out_path, fuzz = 15, padding = 40) {
  img          <- png::readPNG(in_path)
  r_ch         <- img[,,1] * 255
  g_ch         <- img[,,2] * 255
  b_ch         <- img[,,3] * 255
  near_white   <- (r_ch >= (255-fuzz)) & (g_ch >= (255-fuzz)) & (b_ch >= (255-fuzz))
  content_rows <- which(!apply(near_white, 1, all))
  content_cols <- which(!apply(near_white, 2, all))
  if (length(content_rows) == 0 || length(content_cols) == 0) {
    file.copy(in_path, out_path, overwrite = TRUE)
    return(invisible(NULL))
  }
  r1 <- max(1, min(content_rows) - padding)
  r2 <- min(nrow(img), max(content_rows) + padding)
  c1 <- max(1, min(content_cols) - padding)
  c2 <- min(ncol(img), max(content_cols) + padding)
  img_cropped <- img[r1:r2, c1:c2, 1:dim(img)[3]]
  png::writePNG(img_cropped, out_path)
  message(basename(in_path), ": cropped")
}

fnames <- c("BOP1_NOS_pathway.png","BOP2_myocytes.png",
            "BOP3_TP53.png","BOP4_WNK1_OXSR1.png",
            "BOP5a_TNFalpha.png","BOP5b_PI3K_Akt.png",
            "physiological_map.png")

for (fname in fnames) {
  in_p  <- paste0(out_dir, fname)
  out_p <- paste0(out_dir, sub(".png", "_final.png", fname))
  if (file.exists(in_p)) crop_whitespace(in_p, out_p)
}

cat("\nALL_BOPs_complete.R done\n")
cat("Figures saved in:", out_dir, "\n")