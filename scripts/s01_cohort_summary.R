################################################################################
# NULISA-CART-Proteome
# s01_cohort_summary.R
#
# Descriptive cohort summary and visualization.
#
# Generates:
#   - Table 1: Baseline and post-treatment clinical characteristics
#   - Fig 1: CRS/ICANS swimmer plots (kinetics of toxicity and treatment)
#   - Extended Data Fig 1: Sample collection summary
#   - CRS/ICANS distribution bar plots (Ext Data Fig 4a)
#   - Fig 2A: Full-cohort NPQ heatmap (all 204 proteins x all timepoints)
#   - Supplementary Table 1: Protein panel annotation
#
# Note: Clinical association tests (toxicity x response, CRS x ICANS,
# cross-tab heatmaps) are in s06_clinical_associations.R.
#
# Inputs:
#   - data/private/cleaned_data.RData
#   - data/public/swimmer_plot_data.xlsx
#   - data/public/sample_collection.xlsx
################################################################################

source(file.path("scripts", "s00_setup.R"))

# Load cleaned data
load(file.path(dataDir_private, "cleaned_data.RData"))

pm <- patient_metadata  # shorthand


# =============================================================================
# 1. TABLE 1: Baseline and post-treatment clinical characteristics
# =============================================================================
cat("=== Table 1 ===\n")

# Age
mean(pm$age); sd(pm$age); median(pm$age); min(pm$age); max(pm$age)

# Sex
table(pm$sex)
table(pm$sex) / sum(table(pm$sex))

# Prior therapies
table(pm$n_prior_therapies, useNA = 'always')
median(pm$n_prior_therapies, na.rm = TRUE)
min(pm$n_prior_therapies, na.rm = TRUE)
max(pm$n_prior_therapies, na.rm = TRUE)

# Prior stem cell transplantation
table(pm$prior_auto_sct, useNA = 'always')
table(pm$prior_auto_sct, useNA = 'always') / sum(table(pm$prior_auto_sct, useNA = 'always'))

# Disease type
table(pm$disease_type, useNA = 'always')
round(table(pm$disease_type, useNA = 'always') / sum(table(pm$disease_type, useNA = 'always')) * 100, 1)

# CAR-T product
table(pm$cart_product, useNA = 'always')
round(table(pm$cart_product, useNA = 'always') / sum(table(pm$cart_product, useNA = 'always')) * 100, 1)

# 3-month response
table(pm$response_3mo, useNA = 'always')
table(pm$response_3mo, useNA = 'always') / sum(table(pm$response_3mo, useNA = 'always')) * 100

# CRS
table(pm$max_CRS, useNA = 'always')
table(pm$CRS_grade_2_4, useNA = 'always')
round(table(pm$CRS_grade_2_4, useNA = 'always') / sum(table(pm$CRS_grade_2_4, useNA = 'always')) * 100, 1)
mean(pm$max_CRS); sd(pm$max_CRS); median(pm$max_CRS); min(pm$max_CRS); max(pm$max_CRS)

# Time to CRS >= 1
table(pm$time_to_CRS_any, useNA = 'always')
median(pm$time_to_CRS_any, na.rm = TRUE)
quantile(pm$time_to_CRS_any, probs = c(0.25, 0.75), na.rm = TRUE)

# ICANS
table(pm$max_ICANS, useNA = 'always')
table(pm$ICANS_grade_1_4, useNA = 'always')
round(table(pm$ICANS_grade_1_4, useNA = 'always') / sum(table(pm$ICANS_grade_1_4, useNA = 'always')) * 100, 1)
mean(pm$max_ICANS); sd(pm$max_ICANS); median(pm$max_ICANS); min(pm$max_ICANS); max(pm$max_ICANS)

# Time to ICANS >= 1
table(pm$time_to_ICANS_any, useNA = 'always')
median(pm$time_to_ICANS_any, na.rm = TRUE)
quantile(pm$time_to_ICANS_any, probs = c(0.25, 0.75), na.rm = TRUE)


# =============================================================================
# 2. FIG 1: CRS/ICANS SWIMMER PLOTS
# =============================================================================
cat("\n=== Fig 1: CRS/ICANS swimmer plots ===\n")

# --- Load swimmer plot data ---
# Data from Emily Blaum: daily CRS/ICANS grades and treatment administration
# for 80 patients over 14 days post-CAR T-cell infusion.
# Columns: ID, therapy_day, crs_grade, toci, siltux, icans_grade, steroids, anakinra
swimmer_data <- read_excel(file.path(dataDir_public, "swimmer_plot_data.xlsx"))

# Create treatment day columns (NA if treatment not given on that day)
swimmer_data <- swimmer_data %>%
  mutate(
    toci_day     = case_when(toci == 1 ~ therapy_day),
    siltux_day   = case_when(siltux == 1 ~ therapy_day),
    steroids_day = case_when(steroids == 1 ~ therapy_day),
    anakinra_day = case_when(anakinra == 1 ~ therapy_day)
  )

# Convert grades to character for discrete color mapping
swimmer_data$crs_grade   <- as.character(swimmer_data$crs_grade)
swimmer_data$icans_grade <- as.character(swimmer_data$icans_grade)

# Shared theme for swimmer plot panels
swimmer_theme <- theme_minimal(base_size = 8) +
  theme(
    axis.title   = element_text(size = 14),
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 3, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    plot.title   = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10)
  )

# Shared x-axis for days 0-14
day_breaks <- 0:14
day_labels <- as.character(day_breaks)

# --- Grade color palettes ---
# CRS: 0-4 + NA
crs_grade_colors <- c(
  "0"  = "grey40",
  "1"  = "deepskyblue2",
  "2"  = "darkmagenta",
  "3"  = "orange",
  "4"  = "chartreuse3",
  "NA" = "grey80"
)

# ICANS: 0-4 + NA
icans_grade_colors <- c(
  "0"  = "grey40",
  "1"  = "deepskyblue2",
  "2"  = "darkmagenta",
  "3"  = "orange",
  "4"  = "chartreuse3",
  "NA" = "grey80"
)

# --- Panel A: CRS swimmer plot ---
fig1a <- ggplot(swimmer_data, aes(y = ID, group = ID)) +
  swimmer_theme +
  geom_line(aes(x = therapy_day, color = crs_grade), linewidth = 0.9) +
  geom_point(aes(x = toci_day), shape = 15, size = 1.3, color = "black",
             na.rm = TRUE) +
  geom_point(aes(x = siltux_day), shape = 15, size = 1.3, color = "red",
             na.rm = TRUE) +
  scale_color_manual(
    name   = "CRS grade",
    values = crs_grade_colors,
    breaks = c("0", "1", "2", "3", "4", "NA"),
    labels = c("0", "1", "2", "3", "4", "NA"),
    na.translate = FALSE
  ) +
  scale_x_continuous(breaks = day_breaks, labels = day_labels) +
  labs(
    x = "Days post CAR T-cell infusion",
    y = "Sample ID"
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 2)))

# --- Panel B: ICANS swimmer plot ---
fig1b <- ggplot(swimmer_data, aes(y = ID, group = ID)) +
  swimmer_theme +
  geom_line(aes(x = therapy_day, color = icans_grade), linewidth = 0.9) +
  geom_point(aes(x = steroids_day), shape = 15, size = 1.3, color = "black",
             na.rm = TRUE) +
  geom_point(aes(x = anakinra_day), shape = 15, size = 1.3, color = "red",
             na.rm = TRUE) +
  scale_color_manual(
    name   = "ICANS grade",
    values = icans_grade_colors,
    breaks = c("0", "1", "2", "3", "4", "NA"),
    labels = c("0", "1", "2", "3", "4", "NA"),
    na.translate = FALSE
  ) +
  scale_x_continuous(breaks = day_breaks, labels = day_labels) +
  labs(
    x = "Days post CAR T-cell infusion",
    y = "Sample ID"
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 2)))

# --- Save combined figure ---
# Note: Treatment markers are shown as squares on the swimmer lanes.
# CRS panel: black = tocilizumab, red = siltuximab
# ICANS panel: black = steroids, red = anakinra
pdf(file.path(figureDir_main, "Fig1_CRS_ICANS_swimmer_plots.pdf"),
    width = 10, height = 16)
print(cowplot::plot_grid(fig1a, fig1b, ncol = 1, labels = c("a", "b"),
                         label_size = 18))
dev.off()
cat("Saved: Fig1_CRS_ICANS_swimmer_plots.pdf\n")


# =============================================================================
# 3. EXTENDED DATA FIG 1: SAMPLE COLLECTION SUMMARY
# =============================================================================
cat("\n=== Extended Data Fig 1: Sample collection ===\n")

# --- Load sample collection data ---
# Data from Emily Blaum: daily sample collection records for 80 patients.
# Columns: ID, therapy_day, sample_collect (1 = collected, 0 = not)
# therapy_day = -5 is a placeholder for the day of apheresis (actual timing
# varies between patients, approximately 21 days before infusion).
sample_collection <- read_excel(file.path(dataDir_public, "sample_collection.xlsx"))

# Create column with day only when sample was collected (NA otherwise)
sample_collection <- sample_collection %>%
  mutate(sample_day = case_when(sample_collect == 1 ~ therapy_day))

# --- Panel A: Sample collection dot plot ---
ext_fig1a <- ggplot(sample_collection, aes(y = ID, group = ID, x = therapy_day)) +
  theme_bw(base_size = 10) +
  geom_point(aes(x = sample_day), shape = 16, size = 1, color = "black",
             na.rm = TRUE) +
  scale_x_continuous(
    breaks = c(-5, 0:16),
    labels = c("DA", as.character(0:16))
  ) +
  labs(
    x = "Day of therapy",
    y = "Sample ID"
  ) +
  theme(
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.y  = element_text(size = 3),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.x  = element_text(size = 10)
  )

# --- Panel B: Sample size summary table ---
# Print sample counts per time interval
table(sample_metadata$day_cat)

nobs <- data.frame(
  `Time interval` = TIMEPOINT_LABELS,
  `Sample size (n)` = as.numeric(table(sample_metadata$day_cat)),
  check.names = FALSE
)
cat("Sample sizes per time interval:\n")
print(nobs)

# Create table as a grob for the figure
nobs_table <- gridExtra::tableGrob(
  nobs,
  rows  = NULL,
  theme = gridExtra::ttheme_minimal(base_size = 12)
)

# --- Save combined figure ---
pdf(file.path(figureDir_extended, "ExtDataFig1_sample_collection.pdf"),
    width = 14, height = 16)
print(cowplot::plot_grid(
  ext_fig1a, nobs_table,
  ncol = 1, labels = c("a", "b"), label_size = 18,
  rel_heights = c(4, 1)
))
dev.off()
cat("Saved: ExtDataFig1_sample_collection.pdf\n")


# =============================================================================
# 4. CRS/ICANS DISTRIBUTION BAR PLOTS (Ext Data Fig 4a)
# =============================================================================
cat("\n=== CRS/ICANS distributions ===\n")

# --- Alamar color palette ---
cols <- alamarColorPalette(11, nReps = 10)
maxICANS_col <- colorRamp2(c(0, 4), c(cols[[1]][10], cols[[1]][2]))
maxCRS_col   <- colorRamp2(c(0, 4), c(cols[[2]][10], cols[[2]][2]))

# --- maxCRS bar plot ---
maxCRS_plot <- function() {
  par(cex = 1.5, las = 1, mar = c(3, 3.5, 4.5, 0.5))
  tab <- table(pm$max_CRS)
  b <- barplot(tab, main = '', col = maxCRS_col(0:4))
  title(main = 'maxCRS', line = 2)
  par(xpd = TRUE)
  for (i in seq_along(tab)) {
    text(b[i, 1], tab[i], labels = tab[i], pos = 3)
  }
  mtext('Frequency', side = 2, line = 2.5, las = 0, cex = 1.5)
  par(xpd = FALSE)
  grid()
  barplot(tab, main = '', col = maxCRS_col(0:4),
          xaxt = 'n', yaxt = 'n', ylab = '', xlab = '', add = TRUE)
}

# --- maxICANS bar plot ---
maxICANS_plot <- function() {
  par(cex = 1.5, las = 1, mar = c(3, 3.5, 4.5, 0.5))
  tab <- table(pm$max_ICANS)
  b <- barplot(tab, main = '', col = maxICANS_col(0:4), ylim = c(0, 50))
  title(main = 'maxICANS', line = 2)
  par(xpd = TRUE)
  for (i in seq_along(tab)) {
    text(b[i, 1], tab[i], labels = tab[i], pos = 3)
  }
  mtext('Frequency', side = 2, line = 2.5, las = 0, cex = 1.5)
  par(xpd = FALSE)
  grid()
  barplot(tab, main = '', col = maxICANS_col(0:4), ylim = c(0, 50),
          xaxt = 'n', yaxt = 'n', ylab = '', xlab = '', add = TRUE)
}

# --- Save ---
pdf(file.path(figureDir_extended, 'ExtDataFig4a_maxCRS_ICANS_distributions.pdf'),
    width = 10, height = 5)
par(mfrow = c(1, 2))
maxCRS_plot()
maxICANS_plot()
dev.off()
cat("Saved: ExtDataFig4a_maxCRS_ICANS_distributions.pdf\n")


# =============================================================================
# 5. FIG 2A: Full-cohort NPQ heatmap
# =============================================================================
cat("\n=== Fig 2A: Full-cohort NPQ heatmap ===\n")

# Z-score each protein relative to Day of Apheresis mean and SD
DA_samples <- sample_metadata$SampleName[sample_metadata$day_cat == "DA"]
DA_means <- rowMeans(NPQ_data[, DA_samples], na.rm = TRUE)
DA_sds   <- apply(NPQ_data[, DA_samples], 1, sd, na.rm = TRUE)
NPQ_scaled <- t(scale(t(NPQ_data), center = DA_means, scale = DA_sds))

# Prepare heatmap sample metadata
heatmap_samples <- sample_metadata
heatmap_samples$sex_label <- factor(
  ifelse(heatmap_samples$sex == "M", "Male", "Female"),
  levels = c("Male", "Female")
)

# Day category display labels
day_cat_labels <- setNames(TIMEPOINT_LABELS, TIMEPOINT_LEVELS)
heatmap_samples$day_cat_label <- factor(
  day_cat_labels[as.character(heatmap_samples$day_cat)],
  levels = TIMEPOINT_LABELS
)

# Expansion annotation (Above/Below median / No data)
heatmap_samples$expansion <- as.character(heatmap_samples$expansion_above_median)
heatmap_samples$expansion[heatmap_samples$expansion == "0"] <- "Below median"
heatmap_samples$expansion[heatmap_samples$expansion == "1"] <- "Above median"
heatmap_samples$expansion[is.na(heatmap_samples$expansion_above_median)] <- "No data"
heatmap_samples$expansion <- factor(heatmap_samples$expansion,
                                     levels = c("Above median", "Below median", "No data"))

# Merge max CRS / ICANS from patient_metadata
heatmap_samples <- heatmap_samples %>%
  left_join(pm[, c("patientID", "max_CRS", "max_ICANS")], by = "patientID")

# Ensure column ordering matches NPQ_scaled
NPQ_scaled <- NPQ_scaled[, heatmap_samples$SampleName]

# --- Alamar color palettes for annotations ---
cols <- alamarColorPalette(11, nReps = 10)

age_colors <- colorRamp2(
  c(min(heatmap_samples$age, na.rm = TRUE), max(heatmap_samples$age, na.rm = TRUE)),
  c(cols[[8]][10], cols[[8]][2])
)
sex_colors <- c(cols[[6]][9], cols[[5]][5])
maxCRS_col <- colorRamp2(
  c(min(heatmap_samples$max_CRS, na.rm = TRUE), max(heatmap_samples$max_CRS, na.rm = TRUE)),
  c(cols[[2]][10], cols[[2]][2])
)
maxICANS_col <- colorRamp2(
  c(min(heatmap_samples$max_ICANS, na.rm = TRUE), max(heatmap_samples$max_ICANS, na.rm = TRUE)),
  c(cols[[1]][10], cols[[1]][2])
)
expans_col <- c(`Below median` = cols[[4]][9], `Above median` = cols[[4]][2], `No data` = 'grey')

# Column annotation
col_ha <- HeatmapAnnotation(
  `Age (years)` = heatmap_samples$age,
  Sex = heatmap_samples$sex_label,
  maxCRS = heatmap_samples$max_CRS,
  maxICANS = heatmap_samples$max_ICANS,
  `Expansion` = heatmap_samples$expansion,
  col = list(
    `Age (years)` = age_colors,
    Sex = c('Male' = sex_colors[2], 'Female' = sex_colors[1]),
    maxCRS = maxCRS_col,
    maxICANS = maxICANS_col,
    `Expansion` = expans_col
  ),
  show_legend = TRUE,
  annotation_legend_param = list(
    maxICANS = list(color_bar = 'discrete'),
    maxCRS = list(color_bar = 'discrete')
  ),
  annotation_name_side = "right"
)

# Heatmap color scale
col_fun <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))

# Draw heatmap
h <- Heatmap(
  NPQ_scaled,
  col = col_fun,
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_split = heatmap_samples$day_cat_label,
  cluster_column_slices = FALSE,
  clustering_method_rows = 'ward.D2',
  clustering_distance_rows = 'euclidean',
  clustering_method_columns = 'ward.D2',
  clustering_distance_columns = 'euclidean',
  row_names_gp = gpar(fontsize = 3),
  top_annotation = col_ha,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(title = 'NPQ z-score'),
  row_split = NULL,
  row_title = "Proteins",
  row_title_rot = 90,
  column_title_side = 'top'
)

# Save
pdf(file.path(figureDir_main, 'Fig2A_cohort_heatmap.pdf'),
    width = 20, height = 10)
draw(h, merge_legend = TRUE)
dev.off()
cat("Saved: Fig2A_cohort_heatmap.pdf\n")


# =============================================================================
# 6. SUPPLEMENTARY TABLE 1: Protein panel annotation
# =============================================================================
write_supp_table(
  list("Protein Panel Annotation" = target_metadata),
  "SuppTable1_protein_panel_annotation.xlsx"
)


cat("\n", strrep("=", 60), "\n")
cat("s01_cohort_summary.R complete.\n")
cat(strrep("=", 60), "\n")
