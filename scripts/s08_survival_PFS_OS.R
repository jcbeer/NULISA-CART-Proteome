################################################################################
# NULISA-CART-Proteome
# s08_survival_PFS_OS.R
#
# Progression-free survival (PFS) and overall survival (OS) analyses.
#
# Kaplan-Meier survival curves stratified by CAR T-cell expansion at D6-9
# (above vs. below median), using both absolute cell counts and AUC.
#
# Original analysis by Emily M. Blaum; integrated by Joanne C. Beer.
#
# Generates:
#   - Fig 5F: PFS by absolute CAR T-cell number at D6-9
#   - Fig 5G: OS by absolute CAR T-cell number at D6-9
#   - PFS and OS by AUC (for reference; not shown in main figures)
#
# Inputs:
#   - data/private/survival_analysis.xlsx
#     This file contains protected health information (PHI) including dates
#     of infusion, progression, death, and last contact. It cannot be shared
#     publicly and is .gitignored. Contact the corresponding author
#     for access.
#
#     Expected columns:
#       ID               — anonymized patient identifier
#       leu_date         — date of leukapheresis
#       infusion_date    — date of CAR T-cell infusion
#       last_contact_date — date of last contact
#       progression_date — date of disease progression (= last_contact_date
#                          if censored / lost to follow-up)
#       PFS_status       — progression event indicator (1 = progressed, 0 = censored)
#       death_date       — date of death (= last_contact_date if censored)
#       OS_status        — death event indicator (1 = died, 0 = censored)
#       AUC              — CAR T-cell expansion by AUC through D6-9
#                          (1 = above median, 0 = below median)
#       absol            — CAR T-cell expansion by absolute count at D6-9
#                          (1 = above median, 0 = below median)
################################################################################

source(file.path("scripts", "s00_setup.R"))


# =============================================================================
# 1. LOAD AND PREPARE SURVIVAL DATA
# =============================================================================
cat("=== Loading survival data ===\n")

surv_data <- read_excel(file.path(dataDir_private, "survival_analysis.xlsx"))

# Calculate time intervals (in days) from infusion to event/censoring
surv_data <- surv_data %>%
  mutate(
    os_days  = as.duration(infusion_date %--% death_date) / ddays(1),
    pfs_days = as.duration(infusion_date %--% progression_date) / ddays(1)
  )

# Create labeled factor for expansion groups
surv_data$absol_label <- factor(
  surv_data$absol,
  levels = c(0, 1),
  labels = c("Below median", "Above median")
)
surv_data$AUC_label <- factor(
  surv_data$AUC,
  levels = c(0, 1),
  labels = c("Below median", "Above median")
)

cat(sprintf("Loaded %d patients\n", nrow(surv_data)))
cat(sprintf("PFS events: %d / %d\n", sum(surv_data$PFS_status), nrow(surv_data)))
cat(sprintf("OS events:  %d / %d\n", sum(surv_data$OS_status), nrow(surv_data)))


# =============================================================================
# 2. FIG 5F: PFS BY ABSOLUTE CAR T-CELL NUMBER AT D6-9
# =============================================================================
cat("\n=== Fig 5F: PFS by absolute CAR T-cell count ===\n")

fig5f <- survfit2(Surv(pfs_days, PFS_status) ~ absol_label, data = surv_data) %>%
  ggsurvfit() +
  labs(x = "Days", y = "Probability of progression-free survival") +
  add_confidence_interval() +
  add_risktable() +
  add_pvalue(location = "caption") +
  scale_fill_manual(values = c("blue", "pink"),
                    labels = c("Below median", "Above median")) +
  scale_color_manual(values = c("blue", "pink"),
                     labels = c("Below median", "Above median")) +
  add_legend_title("Absolute CAR T-cell number/ul") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 2000))

pdf(file.path(figureDir_main, "Fig5F_PFS_absolute_expansion.pdf"),
    width = 8, height = 6)
print(fig5f)
dev.off()
cat("Saved: Fig5F_PFS_absolute_expansion.pdf\n")


# =============================================================================
# 3. FIG 5G: OS BY ABSOLUTE CAR T-CELL NUMBER AT D6-9
# =============================================================================
cat("\n=== Fig 5G: OS by absolute CAR T-cell count ===\n")

fig5g <- survfit2(Surv(os_days, OS_status) ~ absol_label, data = surv_data) %>%
  ggsurvfit() +
  labs(x = "Days", y = "Probability of overall survival") +
  add_confidence_interval() +
  add_risktable() +
  add_pvalue(location = "caption") +
  scale_fill_manual(values = c("blue", "pink"),
                    labels = c("Below median", "Above median")) +
  scale_color_manual(values = c("blue", "pink"),
                     labels = c("Below median", "Above median")) +
  add_legend_title("Absolute CAR T-cells/ul blood") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 2000))

pdf(file.path(figureDir_main, "Fig5G_OS_absolute_expansion.pdf"),
    width = 8, height = 6)
print(fig5g)
dev.off()
cat("Saved: Fig5G_OS_absolute_expansion.pdf\n")


# =============================================================================
# 4. PFS BY AUC (reference — not in main figures)
# =============================================================================
cat("\n=== PFS by AUC ===\n")

pfs_auc <- survfit2(Surv(pfs_days, PFS_status) ~ AUC_label, data = surv_data) %>%
  ggsurvfit() +
  labs(x = "Days", y = "Probability of progression-free survival") +
  add_confidence_interval() +
  add_risktable() +
  add_pvalue(location = "caption") +
  scale_fill_manual(values = c("blue", "pink"),
                    labels = c("Below median", "Above median")) +
  scale_color_manual(values = c("blue", "pink"),
                     labels = c("Below median", "Above median")) +
  add_legend_title("AUC") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 2000))

pdf(file.path(figureDir_main, "Fig5_PFS_AUC_expansion.pdf"),
    width = 8, height = 6)
print(pfs_auc)
dev.off()
cat("Saved: Fig5_PFS_AUC_expansion.pdf\n")


# =============================================================================
# 5. OS BY AUC (reference — not in main figures)
# =============================================================================
cat("\n=== OS by AUC ===\n")

os_auc <- survfit2(Surv(os_days, OS_status) ~ AUC_label, data = surv_data) %>%
  ggsurvfit() +
  labs(x = "Days", y = "Probability of overall survival") +
  add_confidence_interval() +
  add_risktable() +
  add_pvalue(location = "caption") +
  scale_fill_manual(values = c("blue", "pink"),
                    labels = c("Below median", "Above median")) +
  scale_color_manual(values = c("blue", "pink"),
                     labels = c("Below median", "Above median")) +
  add_legend_title("AUC") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 2000))

pdf(file.path(figureDir_main, "Fig5_OS_AUC_expansion.pdf"),
    width = 8, height = 6)
print(os_auc)
dev.off()
cat("Saved: Fig5_OS_AUC_expansion.pdf\n")


cat("\n", strrep("=", 60), "\n")
cat("s08_survival_PFS_OS.R complete.\n")
cat(strrep("=", 60), "\n")
