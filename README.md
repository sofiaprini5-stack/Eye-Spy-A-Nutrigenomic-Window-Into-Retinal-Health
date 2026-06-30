# Eye Spy: Opening a Nutrigenomic Window onto Retinal Vascular Health

This repository contains the analysis pipeline for a bachelor's thesis investigating the effects of a phytochemical-rich dietary intervention on retinal microvascular caliber, genetic moderation of intervention response, and transcriptomic pathway signatures, using data from the MiBLEND study (Maastricht University).

Author: Sofia M. Prini
Programme: BSc Biomedical Sciences, Faculty of Health, Medicine and Life Sciences (FHML), Maastricht University
Supervisors: Danyel Jennen & Simone van Breda (Translational Genomics)

Overview

This study is a secondary analysis of the MiBLEND randomised crossover dietary intervention, examining whether phytochemical-rich smoothie blends improve retinal arteriolar and venular caliber (CRAE, CRVE, AVR), whether genetic profile across 15 functionally grouped SNPs moderates this response, and what transcriptomic pathway signatures underlie the observed vascular changes. Five Beneficial Outcome Pathways (BOPs) were constructed to contextualise key molecular findings within established vascular biology, each structured around a Molecular Initiating Event (MIE), Key Event(s) (KE), and a vascular Outcome.

Analysis Pipeline

The pipeline runs sequentially from data preparation through pathway-level visualisation: linking and cleaning the source datasets, scoring genetic variants into functional groups, fitting linear mixed models for vessel caliber outcomes, running transcriptomic differential expression and gene set enrichment analysis, computing gene-vessel correlations, building protein-protein interaction networks, and finally constructing the BOP figures. Full descriptions of individual steps are provided in the Methods section of the thesis.

Requirements

Analyses were performed in R. Key packages include `lme4` and `lmerTest` for linear mixed models, `edgeR` and `limma` for transcriptomic normalisation and differential expression, `clusterProfiler` for GSEA and ORA, and `RCy3` with Cytoscape for PPI network and WikiPathways analysis. A full package and version list is provided in the thesis appendix.

Data Availability

Raw participant-level data are not included in this repository due to participant privacy and data sharing agreements associated with the original MiBLEND study. Source data structure and variable definitions are described in the Methods section of the thesis.

Reference

Original study: DeBenedictis JN, Murrell C, Hauser D, van Herwijnen M, Elen B, de Kok TM, van Breda SG. Effects of Different Combinations of Phytochemical-Rich Fruits and Vegetables on Chronic Disease Risk Markers and Gene Expression Changes: Insights from the MiBLEND Study, a Randomized Trial. Antioxidants (Basel). 2024;13(8):915.
