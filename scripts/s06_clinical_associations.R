################################################################################
# NULISA-CART-Proteome
# s06_clinical_associations.R
#
# Clinical association tests between toxicity, response, and expansion.
#
# Generates:
#   - Association tests: response × CRS, response × ICANS
#   - Association tests: CRS × ICANS (binary and full grade)
#   - CRS × ICANS cross-tab heatmaps: grade and time-to-onset (Ext Data Fig 4b)
#   - Time-to-CRS × time-to-ICANS Fisher's test
#   - Expansion × toxicity × response associations (Fig 5, Ext Data Fig 10)
#
# Inputs:
#   - data/private/cleaned_data.RData
################################################################################

source(file.path("scripts", "s00_setup.R"))

# Load cleaned data
load(file.path(dataDir_private, "cleaned_data.RData"))

pm <- patient_metadata  # shorthand


# =============================================================================
# 1. RESPONSE × TOXICITY
# =============================================================================
cat("=== Association of toxicity and response ===\n")

# CRS vs response
cat("\n--- Response × CRS (Grade 2-4 vs 0-1) ---\n")
addmargins(table(pm$response_3mo, pm$CRS_grade_2_4, useNA = 'always'))
addmargins(table(pm$response_3mo, pm$CRS_grade_2_4))
print(fisher.test(table(pm$response_3mo, pm$CRS_grade_2_4)))

# ICANS vs response
cat("\n--- Response × ICANS (Grade 1-4 vs 0) ---\n")
addmargins(table(pm$response_3mo, pm$ICANS_grade_1_4, useNA = 'always'))
addmargins(table(pm$response_3mo, pm$ICANS_grade_1_4))
print(fisher.test(table(pm$response_3mo, pm$ICANS_grade_1_4)))


# =============================================================================
# 2. CRS × ICANS ASSOCIATION
# =============================================================================
cat("\n=== Association of CRS and ICANS ===\n")

# Binary
cat("\n--- Binary: CRS (Grade 2-4) × ICANS (Grade 1-4) ---\n")
addmargins(table(pm$CRS_grade_2_4, pm$ICANS_grade_1_4, useNA = 'always'))
addmargins(table(pm$CRS_grade_2_4, pm$ICANS_grade_1_4))
print(fisher.test(table(pm$CRS_grade_2_4, pm$ICANS_grade_1_4)))

# Full grade
cat("\n--- Full grade: maxCRS × maxICANS ---\n")
print(fisher.test(table(pm$max_CRS, pm$max_ICANS)))


# =============================================================================
# 3. CRS × ICANS CROSS-TAB HEATMAPS (Ext Data Fig 4b)
# =============================================================================
cat("\n=== CRS × ICANS cross-tab heatmaps ===\n")

# --- Cross-tab: maxCRS grade vs maxICANS grade ---
cross_tab_df <- as.data.frame(table(pm$max_CRS, pm$max_ICANS, useNA = 'always'))
cross_tab_df$Percentage <- cross_tab_df$Freq / sum(cross_tab_df$Freq) * 100

cross_tab <- ggplot(cross_tab_df, aes(Var1, Var2, fill = Percentage)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), size = 2) +
  scale_fill_gradient(low = "white", high = "steelblue", name = "Percent") +
  labs(x = "Maximum CRS grade", y = "Maximum ICANS grade", fill = "Percent") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "bottom")

# --- Cross-tab: time-to-CRS vs time-to-ICANS ---
df_time <- data.frame(maxCRS = pm$time_to_CRS_any,
                      maxICANS = pm$time_to_ICANS_any)

tab_time <- as.data.frame(table(df_time$maxCRS, df_time$maxICANS, useNA = 'always'))
tab_time$Percentage <- tab_time$Freq / sum(tab_time$Freq) * 100

# Expand to full grid
full_grid <- expand.grid(Var1 = c(0:12, NA), Var2 = c(0:10, NA))
tab_time$Var1 <- as.integer(as.character(tab_time$Var1))
tab_time$Var2 <- as.integer(as.character(tab_time$Var2))
tab_filled <- full_grid %>%
  left_join(tab_time, by = c("Var1", "Var2")) %>%
  mutate(Freq = replace_na(Freq, 0),
         Percentage = replace_na(Percentage, 0))
tab_filled$Var1 <- factor(tab_filled$Var1)
tab_filled$Var2 <- factor(tab_filled$Var2)

cross_tab_time <- ggplot(tab_filled, aes(Var1, Var2, fill = Percentage)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), size = 2) +
  scale_fill_gradient(low = "white", high = "steelblue", name = "Percent") +
  labs(x = expression("Time to CRS " >= 1),
       y = expression("Time to ICANS " >= 1),
       fill = "Percent") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "bottom")

# --- Save ---
pdf(file.path(figureDir_extended, 'ExtDataFig4B_CRS_ICANS_crosstabs.pdf'),
    width = 12, height = 6)
grid.arrange(cross_tab, cross_tab_time, ncol = 2)
dev.off()
cat("Saved: ExtDataFig4B_CRS_ICANS_crosstabs.pdf\n")

# --- Time-to-onset Fisher's test ---
cat("\nFisher's exact test: time-to-CRS × time-to-ICANS\n")
time_CRS <- pm$time_to_CRS_any
time_ICANS <- pm$time_to_ICANS_any
time_CRS[is.na(time_CRS)] <- 'missing'
time_ICANS[is.na(time_ICANS)] <- 'missing'
print(fisher.test(table(time_CRS, time_ICANS), simulate.p.value = TRUE, B = 1e5))


# =============================================================================
# 4. EXPANSION × TOXICITY × RESPONSE (Fig 5)
# =============================================================================
cat("\n=== Expansion × toxicity × response ===\n")

# Subset to patients with expansion data
pm_exp <- pm[!is.na(pm$expansion_above_median), ]
cat(sprintf("  Patients with expansion data: %d\n", nrow(pm_exp)))

# Expansion × CRS
cat("\n--- Expansion × CRS (Grade 2-4 vs 0-1) ---\n")
addmargins(table(pm_exp$expansion_above_median, pm_exp$CRS_grade_2_4))
print(fisher.test(table(pm_exp$expansion_above_median, pm_exp$CRS_grade_2_4)))

# Expansion × ICANS
cat("\n--- Expansion × ICANS (Grade 1-4 vs 0) ---\n")
addmargins(table(pm_exp$expansion_above_median, pm_exp$ICANS_grade_1_4))
print(fisher.test(table(pm_exp$expansion_above_median, pm_exp$ICANS_grade_1_4)))

# Expansion × Response
cat("\n--- Expansion × Response ---\n")
addmargins(table(pm_exp$expansion_above_median, pm_exp$response_3mo))
print(fisher.test(table(pm_exp$expansion_above_median, pm_exp$response_3mo)))

# CRS × Response (within expansion cohort)
cat("\n--- CRS × Response (expansion cohort) ---\n")
addmargins(table(pm_exp$CRS_grade_2_4, pm_exp$response_3mo))
print(fisher.test(table(pm_exp$CRS_grade_2_4, pm_exp$response_3mo)))

# ICANS × Response (within expansion cohort)
cat("\n--- ICANS × Response (expansion cohort) ---\n")
addmargins(table(pm_exp$ICANS_grade_1_4, pm_exp$response_3mo))
print(fisher.test(table(pm_exp$ICANS_grade_1_4, pm_exp$response_3mo)))

# --- Combined expansion + response variable ---
cat("\n--- Expansion × Response combined groups ---\n")
pm_exp$expans_resp <- interaction(
  ifelse(pm_exp$expansion_above_median == 1, "HighExp", "LowExp"),
  pm_exp$response_3mo, sep = "_"
)
table(pm_exp$expans_resp)

# Combined × CRS
cat("\n  Combined × CRS:\n")
addmargins(table(pm_exp$expans_resp, pm_exp$CRS_grade_2_4))
print(fisher.test(table(pm_exp$expans_resp, pm_exp$CRS_grade_2_4)))

# Combined × ICANS
cat("\n  Combined × ICANS:\n")
addmargins(table(pm_exp$expans_resp, pm_exp$ICANS_grade_1_4))
print(fisher.test(table(pm_exp$expans_resp, pm_exp$ICANS_grade_1_4)))

# --- Wilcoxon: continuous expansion by response ---
cat("\n--- Wilcoxon: continuous expansion × response ---\n")
print(wilcox.test(expansion_cells_per_ul_d6_9 ~ factor(response_3mo), data = pm_exp))

# --- Wilcoxon: continuous expansion by CRS ---
cat("\n--- Wilcoxon: continuous expansion × response ---\n")
print(wilcox.test(expansion_cells_per_ul_d6_9 ~ factor(CRS_grade_2_4), data = pm_exp))

# --- Wilcoxon: continuous expansion by response ---
cat("\n--- Wilcoxon: continuous expansion × response ---\n")
print(wilcox.test(expansion_cells_per_ul_d6_9 ~ factor(ICANS_grade_1_4), data = pm_exp))


cat("\n", strrep("=", 60), "\n")
cat("s06_clinical_associations.R complete.\n")
cat(strrep("=", 60), "\n")
