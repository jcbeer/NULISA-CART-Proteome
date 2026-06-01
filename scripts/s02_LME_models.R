################################################################################
# NULISA-CART-Proteome
# s02_LME_models.R
#
# Fits four sets of linear mixed-effects (LME) models:
#   1. Time effect only (n=80): NPQ ~ age + sex + day_cat + (1|patientID)
#   2. CRS × time (n=80): NPQ ~ age + sex + day_cat * CRS_grade_2_4
#   3. ICANS × time (n=80): NPQ ~ age + sex + day_cat * ICANS_grade_1_4
#   4. Expansion × time (n=71): NPQ ~ age + sex + day_cat * expansion_above_median
#
# Generates:
#   - Volcano plots for each model set (Extended Data Figs 2A, 3A, 5A, 9A)
#   - Coefficient matrices for clustering (saved for s03)
#   - Supplementary Table 2: LME model coefficients (all 4 model sets)
#
# Inputs:
#   - data/private/cleaned_data.RData
#
# Outputs:
#   - data/public/LME_model_results.RData (modelStats, LRTstats, coefficients)
#   - data/public/LME_coefs_for_clustering.RData (coefficient matrices for s03)
#   - supplementary_tables/SuppTable2_LME_model_coefficients.xlsx
################################################################################

source(file.path("scripts", "s00_setup.R"))

# Load cleaned data
load(file.path(dataDir_private, "cleaned_data.RData"))


# =============================================================================
# VOLCANO PLOT FUNCTION
# =============================================================================
# Uses alamarColorPalette for point colors to match manuscript figures.
# Nominal p < 0.05: light colors (circles)
# FDR-corrected p < 0.05: dark colors (diamonds)

plot_volcano <- function(coefs,
                        p_vals,
                        p_vals_FDR = NULL,
                        target_labels,
                        title = NULL,
                        xlabel = expression('log'[2]*'(fold change)'),
                        ylabel = expression('-log'[10]*'(p-value)'),
                        xlimits = NULL,
                        ylimits = NULL,
                        sig_threshold = 0.05,
                        sig_label = 'p = 0.05',
                        sig_line_color = 'black',
                        sig_label_color = 'black',
                        label_all_targets = FALSE,
                        target_labels_off = FALSE,
                        target_label_colors = NULL,
                        target_point_colors = NULL,
                        target_point_shapes = 20,
                        target_label_size = 2,
                        target_label_segment_color = 'grey',
                        segment_size = 0.5,
                        max.overlaps = Inf,
                        force = 1,
                        force_pull = 1,
                        plot_title_font_size = 14,
                        axis_label_font_size = 12,
                        tick_label_font_size = 12,
                        plot_aspect_ratio = 1,
                        log_y = TRUE,
                        add_axis_lines = TRUE) {
  
  # Blank out non-significant labels (unless all requested)
  if (label_all_targets == FALSE) {
    target_labels[p_vals > sig_threshold] <- ''
  }
  if (target_labels_off == TRUE) {
    target_labels <- rep('', length(target_labels))
  }
  
  # Create point colors using Alamar palette
  red_blue <- alamarColorPalette(n = 2, palette = 2, nReps = 5, tint = 'light')
  if (is.null(target_point_colors)) {
    target_point_colors <- rep('grey', length(target_labels))
    target_point_colors[p_vals < sig_threshold & coefs > 0] <- red_blue[[2]][1]
    target_point_colors[p_vals < sig_threshold & coefs < 0] <- red_blue[[1]][1]
    if (!is.null(p_vals_FDR)) {
      target_point_colors[p_vals < sig_threshold & coefs > 0] <- red_blue[[2]][3]
      target_point_colors[p_vals < sig_threshold & coefs < 0] <- red_blue[[1]][3]
      target_point_colors[p_vals_FDR < sig_threshold & coefs > 0] <- red_blue[[2]][1]
      target_point_colors[p_vals_FDR < sig_threshold & coefs < 0] <- red_blue[[1]][1]
      target_point_shapes <- rep(19, length(target_point_colors))
      target_point_shapes[p_vals_FDR < sig_threshold] <- 18
    }
  }
  point_colors <- unique(target_point_colors)
  names(point_colors) <- point_colors
  
  # Create label colors
  if (is.null(target_label_colors)) {
    target_label_colors <- target_point_colors
  }
  label_colors <- unique(target_label_colors)
  names(label_colors) <- label_colors
  
  # Set axis limits
  if (is.null(xlimits)) {
    xlimits <- c(min(coefs), max(coefs))
    xmin <- min(coefs)
  } else {
    xmin <- xlimits[1]
  }
  if (is.null(ylimits)) {
    if (log_y == TRUE) ymax <- max(-log10(p_vals))
    if (log_y == FALSE) ymax <- max(p_vals)
    ylimits <- c(0, ymax)
  }
  
  # Axis line settings
  if (add_axis_lines == TRUE) {
    axis_line_color <- 'grey'
    axis_line_width <- 1
  } else {
    axis_line_color <- 'white'
    axis_line_width <- 0
  }
  
  # Compute y values
  if (log_y == TRUE) y_vals <- -log10(p_vals)
  if (log_y == FALSE) y_vals <- p_vals
  
  # Build plot data
  plot_data <- data.frame(
    coefs = coefs,
    minus_log10_p_vals = y_vals,
    target_labels = target_labels,
    target_point_colors = target_point_colors,
    target_label_colors = target_label_colors,
    target_point_shapes = target_point_shapes,
    stringsAsFactors = FALSE
  )
  
  volcano_plot <- ggplot(plot_data, aes(x = coefs, y = minus_log10_p_vals, label = target_labels)) +
    geom_vline(xintercept = 0, color = 'grey40', linewidth = 0.8) +
    geom_hline(yintercept = 0, color = 'grey40', linewidth = 0.8) +
    geom_point(aes(color = target_point_colors, shape = target_point_shapes), size = 2) +
    scale_color_manual(values = point_colors) +
    scale_shape_identity() +
    scale_x_continuous(limits = xlimits) +
    scale_y_continuous(limits = ylimits) +
    labs(x = xlabel, y = ylabel, title = title) +
    theme_light() +
    theme(
      axis.line = element_line(color = axis_line_color, linewidth = axis_line_width),
      plot.title = element_text(hjust = 0.5, face = 'bold', size = plot_title_font_size),
      axis.title = element_text(size = axis_label_font_size),
      axis.text = element_text(size = tick_label_font_size),
      aspect.ratio = plot_aspect_ratio,
      plot.margin = margin(4, 4, 4, 4),
      legend.position = 'none'
    ) +
    geom_hline(yintercept = -log10(sig_threshold), color = sig_line_color) +
    annotate('text', x = xmin, y = -log10(sig_threshold), label = sig_label,
             color = sig_label_color, size = 2, hjust = 0.25, vjust = -0.5,
             fontface = 'italic') +
    geom_text_repel(
      size = target_label_size,
      max.overlaps = max.overlaps,
      force = force,
      force_pull = force_pull,
      segment.color = target_label_segment_color,
      segment.size = segment_size,
      color = target_label_colors
    ) +
    coord_cartesian(clip = 'off')
  
  return(volcano_plot)
}


# =============================================================================
# 1. DATA PREPARATION
# =============================================================================
cat("=== Preparing data for LME models ===\n")

# Use the prep_lme_data helper from s00 (filters excluded timepoints, zeros → NA)
prepped <- prep_lme_data(NPQ_data, sample_metadata)
NPQ_lme   <- prepped$data
samples_lme <- prepped$samples

# Ensure grouping variables are factors for interaction models
samples_lme$CRS_grade_2_4 <- factor(samples_lme$CRS_grade_2_4)
samples_lme$ICANS_grade_1_4 <- factor(samples_lme$ICANS_grade_1_4)

# Verify day_cat factor levels and reference level
cat("day_cat levels:", levels(samples_lme$day_cat), "\n")
cat("day_cat reference level:", levels(samples_lme$day_cat)[1], "\n")


# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================

# Extract interaction coefficients from lmerNULISAseq output into a clean matrix
# for downstream clustering. Returns a data frame with target + one column per
# timepoint, plus the LRT p-value.
extract_interaction_coefs <- function(lme_result, group_var_coef_prefix,
                                      coef_labels = c('DA', 'D0', 'D1_2', 'D3_5', 'D6_9', 'D11_16'),
                                      output_prefix = "time") {
  
  # Build column names for the interaction coefficients
  # DA coefficient = main effect of group (at reference level)
  # Other timepoints = day_cat{LEVEL}.{group_var} interaction terms
  da_col <- paste0(group_var_coef_prefix, "_coef")
  interaction_cols <- paste0("day_catDAY_", 
                              c("0", "1_2", "3_5", "6_9", "11_16"),
                              ".", group_var_coef_prefix, "_coef")
  all_cols <- c("target", da_col, interaction_cols)
  
  # Extract coefficients
  coefs <- lme_result$modelStats[, all_cols]
  colnames(coefs) <- c("target", paste0(output_prefix, "_", coef_labels))
  
  # Merge LRT p-values
  coefs <- merge(coefs, lme_result$LRTstats[, c("target", "Chisq_test_pval")], all = TRUE)
  
  return(coefs)
}

# Count proteins with at least one significant coefficient across the 
# timepoint-specific terms. This is the number reported in the Results
# (e.g., "121 proteins associated with severe CRS").
count_sig_proteins <- function(lme_result, pval_col_pattern, sig_threshold = 0.05,
                                pval_suffix = "_pval_unadj$", label = "unadjusted") {
  pval_cols <- grep(pval_col_pattern, names(lme_result$modelStats), value = TRUE)
  # Keep only columns matching the suffix (unadjusted or FDR)
  pval_cols <- pval_cols[grepl(pval_suffix, pval_cols)]
  
  pval_matrix <- lme_result$modelStats[, pval_cols, drop = FALSE]
  any_sig <- apply(pval_matrix, 1, function(row) any(row < sig_threshold, na.rm = TRUE))
  
  n_sig <- sum(any_sig)
  sig_targets <- lme_result$modelStats$target[any_sig]
  
  cat(sprintf("  Proteins with >= 1 sig coefficient (%s p < %.2f): %d of %d\n",
              label, sig_threshold, n_sig, nrow(lme_result$modelStats)))
  
  return(list(n = n_sig, targets = sig_targets))
}

# Generate a 3x2 grid of volcano plots for one model set (one per timepoint)
generate_volcano_grid <- function(lme_result, coef_prefix, pval_prefix,
                                   fdr_prefix = NULL, main_title,
                                   xlimits = NULL, ylimits = NULL) {
  
  # Timepoint suffixes for coefficient/pval column names
  tp_suffixes <- c("", "day_catDAY_0.", "day_catDAY_1_2.", 
                    "day_catDAY_3_5.", "day_catDAY_6_9.", "day_catDAY_11_16.")
  tp_labels <- TIMEPOINT_LABELS
  
  plots <- list()
  for (i in seq_along(tp_suffixes)) {
    coef_col <- paste0(tp_suffixes[i], coef_prefix, "_coef")
    pval_col <- paste0(tp_suffixes[i], coef_prefix, "_pval_unadj")
    fdr_col  <- if (!is.null(fdr_prefix)) paste0(tp_suffixes[i], coef_prefix, "_pval_FDR") else NULL
    
    fdr_vals <- if (!is.null(fdr_col)) lme_result$modelStats[[fdr_col]] else NULL
    
    # ylabel on leftmost plots in 3x2 grid: DA (i=1) and Day 3-5 (i=4)
    yl <- if (i %in% c(1, 4)) expression('-log'[10]*'(p-value)') else ''
    # xlabel on bottom row only (positions 4, 5, 6)
    xl <- if (i %in% 4:6) expression('log'[2]*'(fold change)') else ''
    
    plots[[i]] <- plot_volcano(
      coefs = lme_result$modelStats[[coef_col]],
      p_vals = lme_result$modelStats[[pval_col]],
      p_vals_FDR = fdr_vals,
      target_labels = lme_result$modelStats$target,
      title = tp_labels[i],
      xlimits = xlimits,
      ylimits = ylimits,
      sig_label = '',
      xlabel = xl,
      ylabel = yl,
      max.overlaps = 10
    ) +
      annotate('text', x = 2.2, 
               y = -log10(0.05), label = 'p = 0.05',
               color = 'black', size = 3, hjust = 0.25, vjust = -0.5,
               fontface = 'italic')
  }
  
  # Arrange in 3x2 grid
  grid.arrange(
    plots[[1]], plots[[2]], plots[[3]],
    plots[[4]], plots[[5]], plots[[6]],
    ncol = 3,
    top = textGrob(main_title, gp = gpar(col = "black", fontsize = 20, fontface = "bold"))
  )
  
  # Individual full-page plots with full axis labels
  for (i in seq_along(plots)) {
    print(plots[[i]] + 
            labs(x = expression('log'[2]*'(fold change)'),
                 y = expression('-log'[10]*'(p-value)')))
  }
  
  invisible(plots)
}


# =============================================================================
# 3. MODEL SET 1: Time effect only (all 80 patients)
# =============================================================================
cat("\n=== Model Set 1: Time effect (n=80) ===\n")

lme_time <- lmerNULISAseq(
  data = NPQ_lme,
  sampleInfo = samples_lme,
  sampleName_var = 'SampleName',
  modelFormula_fixed = 'age + sex + day_cat',
  reduced_modelFormula_fixed = 'age + sex',
  modelFormula_random = '(1|patientID)',
  return_model_fits = TRUE
)

# Extract time coefficients: intercept + day effects
# For time model, we want the mean NPQ at each timepoint
time_coef_cols <- c('target', 'intercept_coef', 'day_catDAY_0_coef',
                     'day_catDAY_1_2_coef', 'day_catDAY_3_5_coef',
                     'day_catDAY_6_9_coef', 'day_catDAY_11_16_coef')
time_coefs <- lme_time$modelStats[, time_coef_cols]
colnames(time_coefs) <- c('target', paste0('time_', c('DA', 'D0', 'D1_2', 'D3_5', 'D6_9', 'D11_16')))

# Convert to mean NPQ at each timepoint: add intercept (DA) to all other timepoints
time_coefs[, 3:ncol(time_coefs)] <- time_coefs[, 3:ncol(time_coefs)] + time_coefs$time_DA

# Merge LRT p-values
time_coefs <- merge(time_coefs, lme_time$LRTstats[, c('target', 'Chisq_test_pval')], all = TRUE)

cat(sprintf("  Significant proteins (LRT p < 0.05): %d of %d\n",
            sum(time_coefs$Chisq_test_pval < 0.05, na.rm = TRUE), nrow(time_coefs)))

# Count proteins with at least one significant timepoint coefficient
time_sig <- count_sig_proteins(lme_time, "day_cat.*_pval_unadj")
time_sig_fdr <- count_sig_proteins(lme_time, "day_cat.*_pval_FDR", 
                                    pval_suffix = "_pval_FDR$", label = "FDR")

# Volcano plots (Extended Data Fig 2A) — 5 panels (no DA panel for time-only)
pdf(file.path(figureDir_extended, 'ExtDataFig2A_time_volcano.pdf'),
    width = 15, height = 10)

tp_coef_names <- c('day_catDAY_0', 'day_catDAY_1_2', 'day_catDAY_3_5',
                    'day_catDAY_6_9', 'day_catDAY_11_16')
tp_labels_no_DA <- TIMEPOINT_LABELS[2:6]

time_volcanoes <- list()
for (i in seq_along(tp_coef_names)) {
  coef_col <- paste0(tp_coef_names[i], '_coef')
  pval_col <- paste0(tp_coef_names[i], '_pval_unadj')
  fdr_col  <- paste0(tp_coef_names[i], '_pval_FDR')
  
  # ylabel on leftmost plots: Day 0 (i=1, top-left real plot) and Day 3-5 (i=3, bottom-left)
  yl <- if (i %in% c(1, 3)) expression('-log'[10]*'(p-value)') else ''
  xl <- if (i %in% 3:5) expression('log'[2]*'(fold change)') else ''
  
  time_volcanoes[[i]] <- plot_volcano(
    coefs = lme_time$modelStats[[coef_col]],
    p_vals = lme_time$modelStats[[pval_col]],
    p_vals_FDR = lme_time$modelStats[[fdr_col]],
    target_labels = lme_time$modelStats$target,
    title = tp_labels_no_DA[i],
    ylimits = c(0, 80), xlimits = c(-4, 4),
    sig_label = '',
    xlabel = xl,
    ylabel = yl,
    max.overlaps = 10
  ) +
    annotate('text', x = 2.2, y = -log10(0.05), label = 'p = 0.05',
             color = 'black', size = 3, hjust = 0.25, vjust = -0.5,
             fontface = 'italic')
}

grid.arrange(
  nullGrob(), time_volcanoes[[1]],
  time_volcanoes[[2]], time_volcanoes[[3]],
  time_volcanoes[[4]], time_volcanoes[[5]],
  ncol = 3,
  top = textGrob(expression(bold('Main effect of time (relative to Day of Apheresis)')),
                 gp = gpar(col = "black", fontsize = 20))
)

# Individual full-page plots for each timepoint (with full axis labels)
for (i in seq_along(time_volcanoes)) {
  print(time_volcanoes[[i]] + 
          labs(x = expression('log'[2]*'(fold change)'),
               y = expression('-log'[10]*'(p-value)')))
}

dev.off()
cat("Saved: ExtDataFig2A_time_volcano.pdf\n")


# =============================================================================
# 4. MODEL SET 2: CRS × time interaction (all 80 patients)
# =============================================================================
cat("\n=== Model Set 2: CRS × time (n=80) ===\n")

# Filter to complete cases for CRS
samples_CRS <- samples_lme[!is.na(samples_lme$CRS_grade_2_4), ]
NPQ_CRS <- NPQ_lme[, samples_CRS$SampleName]

lme_CRS <- lmerNULISAseq(
  data = NPQ_CRS,
  sampleInfo = samples_CRS,
  sampleName_var = 'SampleName',
  modelFormula_fixed = 'age + sex + day_cat * CRS_grade_2_4',
  modelFormula_random = '(1|patientID)',
  reduced_modelFormula_fixed = 'age + sex + day_cat + CRS_grade_2_4',
  return_model_fits = TRUE
)

# Extract CRS interaction coefficients
CRS_coefs <- extract_interaction_coefs(
  lme_CRS, 
  group_var_coef_prefix = "CRS_grade_2_41",
  output_prefix = "CRS_time"
)

cat(sprintf("  Significant proteins (LRT p < 0.05): %d of %d\n",
            sum(CRS_coefs$Chisq_test_pval < 0.05, na.rm = TRUE), nrow(CRS_coefs)))

# Count proteins with at least one significant CRS interaction coefficient
CRS_sig <- count_sig_proteins(lme_CRS, "CRS_grade_2_41")
CRS_sig_fdr <- count_sig_proteins(lme_CRS, "CRS_grade_2_41", 
                                   pval_suffix = "_pval_FDR$", label = "FDR")

# Volcano plots (Extended Data Fig 3A)
pdf(file.path(figureDir_extended, 'ExtDataFig3A_CRS_volcano.pdf'),
    width = 15, height = 10)
generate_volcano_grid(
  lme_CRS, 
  coef_prefix = "CRS_grade_2_41",
  main_title = 'Maximum CRS (Grade 2-4 vs 0-1) × time interaction',
  xlimits = c(-1.5, 2.5), ylimits = c(0, 10.3)
)
dev.off()
cat("Saved: ExtDataFig3A_CRS_volcano.pdf\n")


# =============================================================================
# 5. MODEL SET 3: ICANS × time interaction (all 80 patients)
# =============================================================================
cat("\n=== Model Set 3: ICANS × time (n=80) ===\n")

# Filter to complete cases for ICANS
samples_ICANS <- samples_lme[!is.na(samples_lme$ICANS_grade_1_4), ]
NPQ_ICANS <- NPQ_lme[, samples_ICANS$SampleName]

lme_ICANS <- lmerNULISAseq(
  data = NPQ_ICANS,
  sampleInfo = samples_ICANS,
  sampleName_var = 'SampleName',
  modelFormula_fixed = 'age + sex + day_cat * ICANS_grade_1_4',
  modelFormula_random = '(1|patientID)',
  reduced_modelFormula_fixed = 'age + sex + day_cat + ICANS_grade_1_4',
  return_model_fits = TRUE
)

# Extract ICANS interaction coefficients
ICANS_coefs <- extract_interaction_coefs(
  lme_ICANS,
  group_var_coef_prefix = "ICANS_grade_1_41",
  output_prefix = "ICANS_time"
)

cat(sprintf("  Significant proteins (LRT p < 0.05): %d of %d\n",
            sum(ICANS_coefs$Chisq_test_pval < 0.05, na.rm = TRUE), nrow(ICANS_coefs)))

# Count proteins with at least one significant ICANS interaction coefficient
ICANS_sig <- count_sig_proteins(lme_ICANS, "ICANS_grade_1_41")
ICANS_sig_fdr <- count_sig_proteins(lme_ICANS, "ICANS_grade_1_41",
                                     pval_suffix = "_pval_FDR$", label = "FDR")

# Volcano plots (Extended Data Fig 5A)
pdf(file.path(figureDir_extended, 'ExtDataFig5A_ICANS_volcano.pdf'),
    width = 15, height = 10)
generate_volcano_grid(
  lme_ICANS,
  coef_prefix = "ICANS_grade_1_41",
  main_title = 'Maximum ICANS (Grade 1-4 vs 0) × time interaction',
  xlimits = c(-3, 2.6), ylimits = c(0, 14.5)
)
dev.off()
cat("Saved: ExtDataFig5A_ICANS_volcano.pdf\n")


# =============================================================================
# 6. MODEL SET 4: Expansion × time interaction (n=71 patients)
# =============================================================================
cat("\n=== Model Set 4: Expansion × time (n=71) ===\n")

# Filter to patients with expansion data
samples_expans <- samples_lme[!is.na(samples_lme$expansion_above_median), ]
NPQ_expans <- NPQ_lme[, samples_expans$SampleName]

cat(sprintf("  Patients with expansion data: %d\n",
            length(unique(samples_expans$patientID))))

lme_expans <- lmerNULISAseq(
  data = NPQ_expans,
  sampleInfo = samples_expans,
  sampleName_var = 'SampleName',
  modelFormula_fixed = 'age + sex + day_cat * expansion_above_median',
  modelFormula_random = '(1|patientID)',
  reduced_modelFormula_fixed = 'age + sex + day_cat + expansion_above_median',
  return_model_fits = TRUE
)

# Extract expansion interaction coefficients
expans_coefs <- extract_interaction_coefs(
  lme_expans,
  group_var_coef_prefix = "expansion_above_median",
  output_prefix = "expans_time"
)

cat(sprintf("  Significant proteins (LRT p < 0.05): %d of %d\n",
            sum(expans_coefs$Chisq_test_pval < 0.05, na.rm = TRUE), nrow(expans_coefs)))

# Count proteins with at least one significant expansion interaction coefficient
expans_sig <- count_sig_proteins(lme_expans, "expansion_above_median")
expans_sig_fdr <- count_sig_proteins(lme_expans, "expansion_above_median",
                                      pval_suffix = "_pval_FDR$", label = "FDR")

# Volcano plots (Extended Data Fig 9A)
pdf(file.path(figureDir_extended, 'ExtDataFig9A_expansion_volcano.pdf'),
    width = 15, height = 10)
generate_volcano_grid(
  lme_expans,
  coef_prefix = "expansion_above_median",
  main_title = 'CAR T-cell expansion (above vs below median) × time interaction',
  xlimits = c(-3.5, 3.5), ylimits = c(0, 14)
)
dev.off()
cat("Saved: ExtDataFig9A_expansion_volcano.pdf\n")


# =============================================================================
# 7. SAVE COEFFICIENTS FOR CLUSTERING (s03)
# =============================================================================
cat("\n=== Saving coefficients for clustering ===\n")

# Package all coefficient matrices and significant protein lists for s03
LME_coefs <- list(
  time_coefs    = time_coefs,
  CRS_coefs     = CRS_coefs,
  ICANS_coefs   = ICANS_coefs,
  expans_coefs  = expans_coefs
)

LME_sig_proteins <- list(
  time    = time_sig$targets,
  CRS     = CRS_sig$targets,
  ICANS   = ICANS_sig$targets,
  expans  = expans_sig$targets
)

LME_LRT_stats <- list(
  time    = lme_time$LRTstats,
  CRS     = lme_CRS$LRTstats,
  ICANS   = lme_ICANS$LRTstats,
  expans  = lme_expans$LRTstats
)

save(LME_coefs, LME_sig_proteins, LME_LRT_stats,
     file = file.path(dataDir_public, "LME_coefs_for_clustering.RData"))
cat("Saved: LME_coefs_for_clustering.RData\n")

# Also save full model results for downstream scripts that need modelStats
LME_results <- list(
  time  = list(modelStats = lme_time$modelStats, LRTstats = lme_time$LRTstats),
  CRS   = list(modelStats = lme_CRS$modelStats, LRTstats = lme_CRS$LRTstats),
  ICANS = list(modelStats = lme_ICANS$modelStats, LRTstats = lme_ICANS$LRTstats),
  expans = list(modelStats = lme_expans$modelStats, LRTstats = lme_expans$LRTstats)
)

save(LME_results,
     file = file.path(dataDir_public, "LME_model_results.RData"))
cat("Saved: LME_model_results.RData\n")


# =============================================================================
# 8. SUPPLEMENTARY TABLE 2: LME model coefficients and LRT statistics
# =============================================================================
cat("\n=== Generating Supplementary Table 2 ===\n")

# Merge LRT stats into modelStats for each model set
merge_lrt <- function(lme_result) {
  merge(lme_result$modelStats, lme_result$LRTstats, by = "target", all = TRUE)
}

write_supp_table(
  list(
    "Overall_Time"         = merge_lrt(lme_time),
    "CRS_Interaction"      = merge_lrt(lme_CRS),
    "ICANS_Interaction"    = merge_lrt(lme_ICANS),
    "Expansion_Interaction" = merge_lrt(lme_expans)
  ),
  "SuppTable2_LME_model_coefficients.xlsx"
)


cat("\n", strrep("=", 60), "\n")
cat("s02_LME_models.R complete.\n")
cat(strrep("=", 60), "\n")
cat("  Significant proteins (LRT / unadj / FDR, p < 0.05):\n")
cat(sprintf("  Model 1 (Time):      %d / %d / %d of %d proteins\n",
            sum(time_coefs$Chisq_test_pval < 0.05, na.rm = TRUE), time_sig$n, time_sig_fdr$n, nrow(time_coefs)))
cat(sprintf("  Model 2 (CRS):       %d / %d / %d of %d proteins\n",
            sum(CRS_coefs$Chisq_test_pval < 0.05, na.rm = TRUE), CRS_sig$n, CRS_sig_fdr$n, nrow(CRS_coefs)))
cat(sprintf("  Model 3 (ICANS):     %d / %d / %d of %d proteins\n",
            sum(ICANS_coefs$Chisq_test_pval < 0.05, na.rm = TRUE), ICANS_sig$n, ICANS_sig_fdr$n, nrow(ICANS_coefs)))
cat(sprintf("  Model 4 (Expansion): %d / %d / %d of %d proteins\n",
            sum(expans_coefs$Chisq_test_pval < 0.05, na.rm = TRUE), expans_sig$n, expans_sig_fdr$n, nrow(expans_coefs)))
