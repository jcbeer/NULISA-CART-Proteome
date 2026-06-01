################################################################################
# NULISA-CART-Proteome
# s00_setup.R — Project configuration, libraries, and shared utilities
#
# Dynamic proteomic signatures distinguish CAR T-cell expansion and acute 
# toxicities in large B cell lymphoma patients
#
# Blaum EM, Beer JC, Mian A, Pramanik D, Kirkpatrick R, Ariyapala IS, 
# Hao Q, Feng W, Ma XJ, Luo Y, Hill BT, Gupta N
#
# This script is sourced at the top of every analysis script.
# It loads all required packages, sets up directories, defines shared 
# constants (timepoint labels, color palettes), and provides common 
# utility functions used across multiple analyses.
################################################################################

# =============================================================================
# 1. PACKAGES
# =============================================================================
# This project uses renv for reproducible package management.
# To install all dependencies at the correct versions, run:
#   renv::restore()
#
# If setting up outside of renv, the following manual installation steps apply:
#
#   # 1. Install devtools
#   if (!requireNamespace("devtools", quietly = TRUE))
#     install.packages("devtools")
#
#   # 2. Install BiocManager for Bioconductor packages
#   if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#
#   # 3. Install Bioconductor packages
#   BiocManager::install(c("ComplexHeatmap", "clusterProfiler", 
#                           "org.Hs.eg.db", "ReactomePA"))
#
#   # 4. Install ggalt from CRAN snapshot (required by NULISAseqR)
#   install.packages("ggalt", repos = "http://packagemanager.posit.co/cran/2025-08-02")
#
#   # 5. Install PCAtools (Alamar fork, required by NULISAseqR)
#   devtools::install_github("Alamar-Biosciences/PCAtools")
#
#   # 6. Install NULISAseqR (Alamar Biosciences)
#   devtools::install_github("Alamar-Biosciences/NULISAseqR")
#
#   # 7. Install remaining CRAN packages
#   install.packages(c("tidyverse", "readxl", "data.table", "survival", "survminer",
#                       "ggsurvfit", "msigdbr", "circlize", "RColorBrewer", "rcartocolor",
#                       "ggrepel", "cowplot", "gridExtra", "scales", "openxlsx"))
#
#   For full installation details, see:
#   https://nulisaseqr.alamarbio.com/user-guide/index.html#installation

# Core data manipulation
library(tidyverse)      # includes dplyr, tidyr, ggplot2, tibble, readr, purrr, stringr
library(readxl)          # Excel file reading (not auto-loaded by tidyverse)
library(data.table)

# NULISA-specific
library(NULISAseqR)     # v1.5.0 — LME modeling for NULISAseq data
library(XML)             # undeclared NULISAseqR dependency
library(fields)          # undeclared NULISAseqR dependency
library(PCAtools)        # Alamar fork, undeclared NULISAseqR dependency

# Statistical modeling
library(survival)       # v3.8-3 — Cox PH models, Surv objects
library(survminer)      # v0.5.1 — survival curve visualization
library(ggsurvfit)       # survfit2 and ggsurvfit for KM curves

# Pathway enrichment
library(clusterProfiler) # v4.16.0 — ORA and GSEA
library(org.Hs.eg.db)    # v3.21.0 — human gene annotation
library(ReactomePA)      # v1.52.0 — Reactome pathway analysis
library(msigdbr)         # v25.1.1 — MSigDB gene sets (Hallmark)

# Visualization
library(ComplexHeatmap)  # v2.24.1 — heatmaps with annotations
library(circlize)        # colorRamp2 for heatmap color scales
library(RColorBrewer)
library(rcartocolor)     # additional color palettes
library(ggrepel)         # non-overlapping text labels
library(cowplot)          # multi-panel figure composition
library(gridExtra)       # arranging multiple plots
library(grid)            # textGrob and other grid utilities
library(scales)          # percent formatting, color scales

# Output
library(openxlsx)        # Excel file creation for supplementary tables

# Project paths
library(here)            # robust path handling relative to project root


# =============================================================================
# 2. DIRECTORY STRUCTURE
# =============================================================================
# Base project directory — automatically detected from project root.
# All other paths are relative to this.
projectDir <- here::here()

# Data directories
dataDir         <- file.path(projectDir, "data")
dataDir_private <- file.path(projectDir, "data", "private")  # .gitignored — patient-level data
dataDir_public  <- file.path(projectDir, "data", "public")   # aggregate outputs (tracked)

# Figure directories
figureDir          <- file.path(projectDir, "figures")
figureDir_main     <- file.path(figureDir, "main")
figureDir_extended <- file.path(figureDir, "extended_data")

# Supplementary tables output
suppDir <- file.path(projectDir, "supplementary_tables")
suppDir_working <- file.path(suppDir, "working")  # intermediate tables (not referenced in manuscript)

# Create directories if they don't exist
dirs_to_create <- c(
  dataDir, dataDir_private, dataDir_public,
  figureDir, figureDir_main, figureDir_extended,
  suppDir, suppDir_working
)
for (d in dirs_to_create) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

cat("Project directory:", projectDir, "\n")
cat("Private data:     ", dataDir_private, "\n")
cat("Public data:      ", dataDir_public, "\n")
cat("Figures:          ", figureDir, "\n")
cat("Supp tables:      ", suppDir, "\n")


# =============================================================================
# 3. CONSTANTS
# =============================================================================

# --- Timepoint definitions ---
# Factor levels for the categorical time variable (in chronological order)
TIMEPOINT_LEVELS <- c("DA", "DAY_0", "DAY_1_2", "DAY_3_5", "DAY_6_9", "DAY_11_16")

# Display labels for figures
TIMEPOINT_LABELS <- c("Day of Apheresis", "Day 0", "Day 1-2", 
                       "Day 3-5", "Day 6-9", "Day 11-16")

# Short labels for compact figures (heatmaps, trajectories)
TIMEPOINT_LABELS_SHORT <- c("DA", "0", "1-2", "3-5", "6-9", "11-16")

# --- Cohort size ---
N_PATIENTS    <- 80
N_PROTEINS    <- 204
N_GENES       <- 201  # unique genes (heterodimers share proteins)
N_EXPANS_PTS  <- 71   # patients with expansion data

# --- Binary group definitions ---
# CRS: grade 2-4 vs 0-1 (clinically significant threshold)
# ICANS: grade 1-4 vs 0 (any CNS involvement is clinically relevant)
# Expansion: above vs below median at day 6-9

# --- Statistical thresholds ---
LME_SIG_THRESHOLD    <- 0.05  # unadjusted p for LRT filtering
FDR_THRESHOLD        <- 0.05  # BH-corrected FDR for volcano plots
PATHWAY_P_THRESHOLD  <- 0.1   # relaxed threshold for pathway enrichment
COX_SIG_THRESHOLD    <- 0.05  # unadjusted p for Cox models

# --- Clustering rules ---
# Based on biological plausibility and parsimony:
#   n > 100 significant proteins  -> 6 clusters
#   100 >= n >= 10                -> 4 clusters
#   n < 10                        -> 2 clusters
get_n_clusters <- function(n_sig_proteins) {
  if (n_sig_proteins > 100) return(6)
  if (n_sig_proteins >= 10) return(4)
  return(2)
}

# --- Excluded timepoints ---
# D-5 and DA-5 are excluded from all LME analyses
EXCLUDED_TIMEPOINTS <- c("D-5", "DA-5")




# =============================================================================
# 4. SHARED UTILITY FUNCTIONS
# =============================================================================

# ---------------------------------------------------------------------------
# save_pdf: wrapper for cairo_pdf to handle Unicode symbols (mu, ge, etc.)
# ---------------------------------------------------------------------------
save_pdf <- function(filename, plot_expr, width = 10, height = 8, dir = figureDir) {
  filepath <- file.path(dir, filename)
  cairo_pdf(filepath, width = width, height = height)
  eval(plot_expr)
  dev.off()
  cat("Saved:", filepath, "\n")
}


# ---------------------------------------------------------------------------
# prep_lme_data: standard data preparation for LME modeling
# Filters out excluded timepoints and sets zeros to NA
# ---------------------------------------------------------------------------
prep_lme_data <- function(cc_data, samples) {
  # Remove excluded timepoints (D-5 and DA-5)
  samples_filtered <- samples[!(samples$day %in% EXCLUDED_TIMEPOINTS), ]
  cc_data_filtered <- cc_data[, samples_filtered$SampleName]
  
  # Set zero NPQ values to NA (0.07% of values)
  cc_data_filtered[cc_data_filtered == 0] <- NA
  
  cat(sprintf("Prepared data: %d proteins x %d samples (%d patients)\n",
              nrow(cc_data_filtered), ncol(cc_data_filtered), 
              length(unique(samples_filtered$patientID))))
  cat(sprintf("Zero values set to NA: %d (%.2f%%)\n",
              sum(cc_data_filtered == 0, na.rm = TRUE),
              sum(is.na(cc_data_filtered)) / length(unlist(cc_data_filtered)) * 100))
  
  return(list(data = cc_data_filtered, samples = samples_filtered))
}


# ---------------------------------------------------------------------------
# write_supp_table: write a named list of data frames to a multi-sheet Excel
# ---------------------------------------------------------------------------
write_supp_table <- function(sheet_list, filename, dir = suppDir) {
  filepath <- file.path(dir, filename)
  wb <- createWorkbook()
  
  for (sheet_name in names(sheet_list)) {
    # Excel sheet names limited to 31 characters
    safe_name <- substr(sheet_name, 1, 31)
    addWorksheet(wb, safe_name)
    writeData(wb, safe_name, sheet_list[[sheet_name]])
  }
  
  saveWorkbook(wb, filepath, overwrite = TRUE)
  cat("Saved supplementary table:", filepath, "\n")
}


# =============================================================================
# 5. SESSION INFO LOGGING
# =============================================================================
# Print key package versions for reproducibility
cat("\n", strrep("=", 60), "\n")
cat("NULISA-CART-Proteome: s00_setup.R loaded\n")
cat(strrep("=", 60), "\n")
cat("R version:          ", as.character(getRversion()), "\n")
cat("NULISAseqR:         ", as.character(packageVersion("NULISAseqR")), "\n")
cat("ComplexHeatmap:     ", as.character(packageVersion("ComplexHeatmap")), "\n")
cat("clusterProfiler:    ", as.character(packageVersion("clusterProfiler")), "\n")
cat("survival:           ", as.character(packageVersion("survival")), "\n")
cat("survminer:          ", as.character(packageVersion("survminer")), "\n")
cat("ggplot2:            ", as.character(packageVersion("ggplot2")), "\n")
cat(strrep("=", 60), "\n\n")
