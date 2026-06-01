################################################################################
# NULISA-CART-Proteome
# s09_decouple_expansion_toxicity.R
#
# Decoupling analysis: tests whether vascular recovery trajectories distinguish
# patients who experience toxicity from those who do not, among high vs low
# CAR T-cell expanders. Uses cluster-level scores from the Effector Engine
# (expansion model) and Vascular Collapse (time model) clusters.
#
# Generates:
#   - Fig 6: 4-panel trajectory plot (Effector vs Vascular × High vs Low expansion)
#   - Supplementary Table 6: Decoupling analysis scores and Wilcoxon test results
#
# Inputs:
#   - data/private/cleaned_data.RData
#   - data/public/cluster_results.RData (from s03)
################################################################################

source(file.path("scripts", "s00_setup.R"))

# Load data
load(file.path(dataDir_private, "cleaned_data.RData"))
load(file.path(dataDir_public, "cluster_results.RData"))


# =============================================================================
# 1. IDENTIFY CLUSTER TARGETS
# =============================================================================
cat("=== Identifying cluster targets ===\n")

# Find Effector Engine targets (from expansion model)
expans_df <- cluster_results$expans$df
effector_cluster_num <- expans_df$cluster[expans_df$cluster_name == "The Effector Engine"][1]
effector_targets <- expans_df$target[expans_df$cluster == effector_cluster_num]
cat(sprintf("  Effector Engine: %d proteins (expansion cluster %d)\n",
            length(effector_targets), effector_cluster_num))

# Find Vascular Collapse targets (from time model)
time_df <- cluster_results$time$df
vascular_cluster_num <- time_df$cluster[time_df$cluster_name == "Vascular Collapse"][1]
vascular_targets <- time_df$target[time_df$cluster == vascular_cluster_num]
cat(sprintf("  Vascular Collapse: %d proteins (time cluster %d)\n",
            length(vascular_targets), vascular_cluster_num))


# =============================================================================
# 2. Z-SCORE USING DA BASELINE AND COMPUTE CLUSTER SCORES
# =============================================================================
cat("\n=== Computing cluster scores ===\n")

compute_cluster_score <- function(NPQ_data, sample_metadata, targets, score_name) {
  # Extract NPQ for cluster targets
  npq <- as.data.frame(t(NPQ_data[targets, ]))
  npq$SampleName <- rownames(npq)
  
  # Z-score each protein using DA baseline mean and SD
  DA_samples <- sample_metadata$SampleName[sample_metadata$day_cat == "DA"]
  
  z_cols <- c()
  for (target in targets) {
    da_vals <- npq[[target]][npq$SampleName %in% DA_samples]
    da_mean <- mean(da_vals, na.rm = TRUE)
    da_sd   <- sd(da_vals, na.rm = TRUE)
    z_col <- paste0(target, "_z")
    npq[[z_col]] <- (npq[[target]] - da_mean) / da_sd
    z_cols <- c(z_cols, z_col)
  }
  
  # Compute cluster score = mean z-score across proteins
  npq[[score_name]] <- rowMeans(npq[, z_cols], na.rm = TRUE)
  
  return(npq[, c("SampleName", score_name)])
}

effector_scores <- compute_cluster_score(NPQ_data, sample_metadata, 
                                          effector_targets, "effector_score")
vascular_scores <- compute_cluster_score(NPQ_data, sample_metadata,
                                          vascular_targets, "vascular_score")

# Merge scores with sample metadata
combined <- sample_metadata %>%
  left_join(effector_scores, by = "SampleName") %>%
  left_join(vascular_scores, by = "SampleName")

# Filter to patients with expansion data and create group variables
combined <- combined %>%
  filter(!is.na(expansion_above_median)) %>%
  mutate(
    expansion_group = factor(
      ifelse(expansion_above_median == 1, "High Expansion", "Low Expansion"),
      levels = c("High Expansion", "Low Expansion")
    ),
    any_severe_tox = ifelse(CRS_grade_2_4 == 1 | ICANS_grade_1_4 == 1,
                             "Severe Toxicity", "No/Mild Toxicity"),
    day_cat = factor(day_cat, levels = TIMEPOINT_LEVELS)
  )

# Sample sizes
cat("\n  Sample sizes by quadrant:\n")
combined %>%
  group_by(expansion_group, any_severe_tox) %>%
  summarise(n_patients = n_distinct(patientID), .groups = "drop") %>%
  print()


# =============================================================================
# 3. COMPUTE GROUP MEANS AND WILCOXON TESTS
# =============================================================================
cat("\n=== Running Wilcoxon rank-sum tests ===\n")

# Pivot to long format for both clusters
score_long <- combined %>%
  pivot_longer(cols = c(effector_score, vascular_score),
               names_to = "cluster", values_to = "score") %>%
  mutate(cluster = factor(
    ifelse(cluster == "effector_score", "Effector", "Vascular"),
    levels = c("Effector", "Vascular")
  ))

# Group means and SE
group_summary <- score_long %>%
  group_by(day_cat, expansion_group, any_severe_tox, cluster) %>%
  summarise(
    mean_val = mean(score, na.rm = TRUE),
    se_val = sd(score, na.rm = TRUE) / sqrt(n()),
    n_samples = n(),
    n_patients = n_distinct(patientID),
    .groups = "drop"
  )

# Wilcoxon tests: toxicity vs no toxicity at each timepoint × expansion × cluster
pval_results <- score_long %>%
  group_by(expansion_group, cluster, day_cat) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(score ~ any_severe_tox)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

# BH correction within each panel (stored for reference)
pval_results <- pval_results %>%
  group_by(expansion_group, cluster) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup()

cat("\n  Wilcoxon results:\n")
print(pval_results %>% arrange(expansion_group, cluster, day_cat), n = Inf)


# =============================================================================
# 4. FIGURE 6: 4-PANEL TRAJECTORY PLOT
# =============================================================================
cat("\n=== Generating Fig 6 ===\n")

# Use raw p-values for annotations (set to "p_adj" for BH-corrected)
ANNOTATION_P_COL <- "p_value"

# Color mapping with sample sizes in legend
color_values <- c(
  "High Expansion: No/Mild Toxicity" = "#2166ac",
  "High Expansion: Severe Toxicity"  = "#b2182b",
  "Low Expansion: No/Mild Toxicity"  = "#92c5de",
  "Low Expansion: Severe Toxicity"   = "#e8a0a0"
)

# Build color group variable
group_summary <- group_summary %>%
  mutate(color_group = paste0(expansion_group, ": ", any_severe_tox))

# Sample sizes for legend labels
n_per_group <- combined %>%
  group_by(expansion_group, any_severe_tox) %>%
  summarise(n_pat = n_distinct(patientID), .groups = "drop") %>%
  mutate(color_group = paste0(expansion_group, ": ", any_severe_tox))

label_map <- setNames(
  paste0(n_per_group$color_group, " (n = ", n_per_group$n_pat, ")"),
  n_per_group$color_group
)

color_values_labeled <- setNames(
  color_values[names(color_values)],
  label_map[names(color_values)]
)

group_summary <- group_summary %>%
  mutate(color_group_label = factor(label_map[color_group],
                                     levels = label_map[names(color_values)]))

# Significance annotations
y_positions <- group_summary %>%
  group_by(expansion_group, cluster, day_cat) %>%
  summarise(y_top = max(mean_val + se_val, na.rm = TRUE), .groups = "drop")

panel_ranges <- group_summary %>%
  group_by(expansion_group, cluster) %>%
  summarise(y_range = max(mean_val + se_val, na.rm = TRUE) - 
              min(mean_val - se_val, na.rm = TRUE), .groups = "drop")

pval_annotations <- pval_results %>%
  left_join(y_positions, by = c("expansion_group", "cluster", "day_cat")) %>%
  left_join(panel_ranges, by = c("expansion_group", "cluster")) %>%
  mutate(
    p_for_annot = .data[[ANNOTATION_P_COL]],
    sig_label = case_when(
      is.na(p_for_annot) ~ "",
      p_for_annot < 0.05 ~ "*",
      p_for_annot < 0.1  ~ "dot",
      TRUE ~ ""
    ),
    y_nudge = y_range * 0.12
  ) %>%
  filter(sig_label != "")

# Build plot
p <- ggplot(group_summary,
            aes(x = day_cat, y = mean_val,
                color = color_group_label, group = color_group_label)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
                width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_grid(cluster ~ expansion_group,
             scales = "free_y",
             labeller = labeller(
               cluster = c("Effector" = "The Effector Engine",
                            "Vascular" = "Vascular Collapse")
             )) +
  scale_color_manual(values = color_values_labeled, name = "") +
  scale_x_discrete(labels = TIMEPOINT_LABELS_SHORT) +
  labs(
    title = "Decoupling Expansion from Toxicity",
    x = NULL,
    y = "Mean Z-Score \u00b1 SE"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray90"),
    legend.position = "bottom",
    legend.text = element_text(size = 9)
  ) +
  guides(color = guide_legend(nrow = 2, byrow = FALSE))

# Add significance annotations
if (nrow(pval_annotations) > 0) {
  annot_star <- pval_annotations %>% filter(sig_label == "*")
  annot_dot  <- pval_annotations %>% filter(sig_label == "dot")
  
  if (nrow(annot_star) > 0) {
    p <- p + geom_text(data = annot_star,
                        aes(x = day_cat, y = y_top, label = sig_label),
                        inherit.aes = FALSE, vjust = -0.3, size = 8, color = "black")
  }
  if (nrow(annot_dot) > 0) {
    p <- p + geom_point(data = annot_dot,
                         aes(x = day_cat, y = y_top + y_nudge),
                         inherit.aes = FALSE, shape = 19, size = 1, color = "black")
  }
  
  p <- p + scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
  p_label <- ifelse(ANNOTATION_P_COL == "p_adj", "BH-adjusted", "unadjusted")
  p <- p + labs(caption = paste0("Wilcoxon rank-sum test (", p_label, 
                                  "):  * p < 0.05,  \u2022 0.05 \u2264 p < 0.1")) +
    theme(plot.caption = element_text(hjust = 0.5, size = 10, face = "italic"))
}

# Save
ggsave(file.path(figureDir_main, 'Fig6_decouple_expansion_toxicity.pdf'),
       p, width = 7, height = 6, device = cairo_pdf)
cat("Saved: Fig6_decouple_expansion_toxicity.pdf\n")


# =============================================================================
# 5. SAVE RESULTS
# =============================================================================
cat("\n=== Saving results ===\n")

decoupling_results <- list(
  combined_scores = combined,
  group_summary = group_summary,
  pval_results = pval_results,
  effector_targets = effector_targets,
  vascular_targets = vascular_targets
)

save(decoupling_results,
     file = file.path(dataDir_public, "decoupling_results.RData"))
cat("Saved: decoupling_results.RData\n")


# =============================================================================
# 6. SUPPLEMENTARY TABLE 6: Decoupling analysis
# =============================================================================
cat("\n=== Generating Supplementary Table 6 ===\n")

# Sheet 1: Wilcoxon test results
wilcoxon_sheet <- pval_results %>%
  arrange(expansion_group, cluster, day_cat)

# Sheet 2: Group summary statistics
summary_sheet <- group_summary %>%
  dplyr::select(day_cat, expansion_group, any_severe_tox, cluster,
                mean_val, se_val, n_samples, n_patients) %>%
  arrange(cluster, expansion_group, any_severe_tox, day_cat)

# Sheet 3: Cluster target lists
target_sheet <- data.frame(
  Cluster = c(rep("Effector Engine", length(effector_targets)),
              rep("Vascular Collapse", length(vascular_targets))),
  Target = c(sort(effector_targets), sort(vascular_targets))
)

write_supp_table(
  list(
    "Wilcoxon_Tests" = wilcoxon_sheet,
    "Group_Summary"  = summary_sheet,
    "Cluster_Targets" = target_sheet
  ),
  "SuppTable6_decoupling_analysis.xlsx"
)


cat("\n", strrep("=", 60), "\n")
cat("s09_decouple_expansion_toxicity.R complete.\n")
cat(strrep("=", 60), "\n")
