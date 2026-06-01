################################################################################
# NULISA-CART-Proteome
# s07_cox_models.R
#
# Cox proportional hazards models for time-to-CRS and time-to-ICANS,
# using protein expression at Day of Apheresis (DA) and Day 0 (D0).
# All models adjusted for age and sex.
#
# Generates:
#   - Forest plots (Fig 4A-D)
#   - Adjusted cumulative incidence curves for significant proteins
#   - Supplementary Table 5: Cox PH model results
#
# Inputs:
#   - data/private/cleaned_data.RData
#
# Four analysis sets:
#   1. CRS  @ DA:  Surv(time_CRS_2, event_CRS_2)  ~ protein + age + sex
#   2. CRS  @ D0:  Surv(time_CRS_2, event_CRS_2)  ~ protein + age + sex
#   3. ICANS @ DA: Surv(time_ICANS_1, event_ICANS_1) ~ protein + age + sex
#   4. ICANS @ D0: Surv(time_ICANS_1, event_ICANS_1) ~ protein + age + sex
################################################################################

source(file.path("scripts", "s00_setup.R"))

library(survival)
library(survminer)
library(ggforestplot)
library(patchwork)

# Load cleaned data
load(file.path(dataDir_private, "cleaned_data.RData"))


# =============================================================================
# 1. COX MODEL FITTING FUNCTION
# =============================================================================

fit_cox_models <- function(NPQ_data, sample_metadata, patient_metadata,
                            timepoint_day_cat, time_var, event_var) {
  
  # Subset to the specified timepoint
  tp_samples <- sample_metadata[sample_metadata$day_cat == timepoint_day_cat, ]
  tp_data <- NPQ_data[, tp_samples$SampleName]
  
  # Merge time-to-event data from patient_metadata
  tp_samples <- tp_samples %>%
    left_join(patient_metadata[, c("patientID", "age", "sex", time_var, event_var)],
              by = "patientID", suffix = c("", ".pm"))
  
  # Use age/sex from patient_metadata if present (avoid .pm suffix issues)
  if ("age.pm" %in% colnames(tp_samples)) {
    tp_samples$age <- tp_samples$age.pm
    tp_samples$sex <- tp_samples$sex.pm
  }
  
  cat(sprintf("  %s: %d samples, %d events\n",
              timepoint_day_cat, nrow(tp_samples),
              sum(tp_samples[[event_var]], na.rm = TRUE)))
  
  # Fit Cox model for each protein
  results <- data.frame(
    target = rownames(tp_data),
    logHR = NA_real_, logHR_se = NA_real_,
    HR = NA_real_, lower95 = NA_real_, upper95 = NA_real_,
    pval = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(tp_data))) {
    expr <- as.numeric(tp_data[i, ])
    
    fit <- tryCatch(
      coxph(Surv(tp_samples[[time_var]], tp_samples[[event_var]]) ~ 
              expr + tp_samples$age + tp_samples$sex),
      error = function(e) NULL
    )
    
    if (!is.null(fit)) {
      s <- summary(fit)
      results$logHR[i]    <- coef(fit)["expr"]
      results$logHR_se[i] <- s$coefficients["expr", "se(coef)"]
      results$HR[i]       <- exp(coef(fit)["expr"])
      results$lower95[i]  <- exp(confint(fit)["expr", 1])
      results$upper95[i]  <- exp(confint(fit)["expr", 2])
      results$pval[i]     <- s$coefficients["expr", "Pr(>|z|)"]
    }
  }
  
  # Add FDR correction
  results$pval_FDR <- p.adjust(results$pval, method = "BH")
  
  return(list(results = results, samples = tp_samples, data = tp_data))
}


# =============================================================================
# 2. FOREST PLOT FUNCTION
# =============================================================================

create_forest_plot <- function(cox_results, title, sig_threshold = 0.05) {
  
  sig <- cox_results[cox_results$pval < sig_threshold & !is.na(cox_results$pval), ]
  sig <- sig[order(sig$logHR, decreasing = TRUE), ]
  
  if (nrow(sig) == 0) {
    cat(sprintf("  No significant proteins for: %s\n", title))
    return(NULL)
  }
  
  p <- forestplot(
    df = sig,
    name = target,
    estimate = logHR,
    se = logHR_se,
    logodds = TRUE,
    title = title,
    xlab = 'Hazard Ratio'
  )
  
  return(p)
}


# =============================================================================
# 3. ADJUSTED CUMULATIVE INCIDENCE CURVE FUNCTION
# =============================================================================

create_cuminc_plot <- function(target_name, protein_data, sample_data,
                                time_var, event_var,
                                cox_results,
                                outcome_label, timepoint_label,
                                max_time = 16, show_ci = TRUE) {
  
  # Get protein expression and create binary group for visualization
  expr <- as.numeric(protein_data[target_name, ])
  expr_median <- median(expr, na.rm = TRUE)
  target_group <- factor(ifelse(expr >= expr_median, "High", "Low"),
                          levels = c("Low", "High"))
  
  # Build data frame
  cox_data <- data.frame(
    time = sample_data[[time_var]],
    event = sample_data[[event_var]],
    target_group = target_group,
    age = sample_data$age,
    sex = sample_data$sex
  )
  cox_data <- cox_data[complete.cases(cox_data), ]
  
  # Fit binary Cox model for visualization curves
  cox_fit <- coxph(Surv(time, event) ~ target_group + age + sex, data = cox_data)
  
  # Get HR annotation from continuous model results
  target_row <- cox_results[cox_results$target == target_name, ]
  hr <- target_row$HR
  hr_lower <- target_row$lower95
  hr_upper <- target_row$upper95
  pval <- target_row$pval
  
  # Predict adjusted survival curves
  mean_age <- mean(cox_data$age, na.rm = TRUE)
  mode_sex <- names(sort(table(cox_data$sex), decreasing = TRUE))[1]
  
  newdata <- data.frame(
    target_group = factor(c("Low", "High"), levels = c("Low", "High")),
    age = mean_age, sex = mode_sex
  )
  
  surv_fits <- lapply(1:2, function(j) survfit(cox_fit, newdata = newdata[j, ]))
  
  plot_data <- do.call(rbind, lapply(1:2, function(j) {
    sf <- surv_fits[[j]]
    data.frame(
      time = sf$time, cuminc = 1 - sf$surv,
      lower = 1 - sf$upper, upper = 1 - sf$lower,
      group = c("Low", "High")[j]
    )
  }))
  
  # Add time 0 and filter
  plot_data <- rbind(
    data.frame(time = 0, cuminc = 0, lower = 0, upper = 0, group = "Low"),
    data.frame(time = 0, cuminc = 0, lower = 0, upper = 0, group = "High"),
    plot_data
  )
  plot_data <- plot_data[plot_data$time <= max_time, ]
  
  # Dynamic legend/annotation positioning
  max_cuminc <- max(plot_data$cuminc, na.rm = TRUE)
  if (max_cuminc < 0.6) {
    legend_pos <- c(0.75, 0.95); annot_y <- 0.15; annot_vjust <- 0
  } else {
    legend_pos <- c(0.75, 0.15); annot_y <- 0.95; annot_vjust <- 1
  }
  
  p <- ggplot(plot_data, aes(x = time, y = cuminc, color = group, fill = group)) +
    geom_step(linewidth = 1.2) +
    scale_color_manual(values = c("Low" = "#2E9FDF", "High" = "#E7B800"),
                        labels = c(paste0(target_name, " Low"), paste0(target_name, " High"))) +
    scale_fill_manual(values = c("Low" = "#2E9FDF", "High" = "#E7B800"),
                       labels = c(paste0(target_name, " Low"), paste0(target_name, " High"))) +
    labs(x = "Days from CAR-T Infusion",
         y = paste0("Cumulative Incidence of ", outcome_label),
         title = paste0(target_name, " (", timepoint_label, ")"),
         color = NULL, fill = NULL) +
    coord_cartesian(xlim = c(0, max_time), ylim = c(0, 1)) +
    scale_x_continuous(breaks = seq(0, max_time, by = 4)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                        labels = scales::percent_format(accuracy = 1)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", alpha = 0.3, color = "gray50") +
    annotate("text", x = max_time * 0.05, y = annot_y,
             label = sprintf("HR = %.2f (95%% CI: %.2f - %.2f)\np = %.3f",
                             hr, hr_lower, hr_upper, pval),
             hjust = 0, vjust = annot_vjust, size = 3.5, fontface = "bold") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = legend_pos,
      legend.justification = c(1, if (max_cuminc < 0.6) 1 else 0),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.title = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  if (show_ci) {
    p <- p + geom_ribbon(aes(ymin = lower, ymax = upper),
                          alpha = 0.2, color = NA, show.legend = FALSE)
  }
  
  return(p)
}


# =============================================================================
# 4. RUN ALL 4 ANALYSIS SETS
# =============================================================================
cat("=== Fitting Cox PH models ===\n")

# Define the 4 analysis configurations
analyses <- list(
  CRS_DA = list(day_cat = "DA",      time_var = "time_CRS_2",   event_var = "event_CRS_2",
                outcome = "Severe CRS",  timepoint = "Day of Apheresis"),
  CRS_D0 = list(day_cat = "DAY_0",   time_var = "time_CRS_2",   event_var = "event_CRS_2",
                outcome = "Severe CRS",  timepoint = "Day 0"),
  ICANS_DA = list(day_cat = "DA",     time_var = "time_ICANS_1", event_var = "event_ICANS_1",
                  outcome = "ICANS",     timepoint = "Day of Apheresis"),
  ICANS_D0 = list(day_cat = "DAY_0",  time_var = "time_ICANS_1", event_var = "event_ICANS_1",
                  outcome = "ICANS",     timepoint = "Day 0")
)

cox_results <- list()

for (aname in names(analyses)) {
  a <- analyses[[aname]]
  cat(sprintf("\n--- %s ---\n", aname))
  
  cox_results[[aname]] <- fit_cox_models(
    NPQ_data, sample_metadata, patient_metadata,
    timepoint_day_cat = a$day_cat,
    time_var = a$time_var,
    event_var = a$event_var
  )
  
  n_sig <- sum(cox_results[[aname]]$results$pval < 0.05, na.rm = TRUE)
  n_fdr <- sum(cox_results[[aname]]$results$pval_FDR < 0.05, na.rm = TRUE)
  cat(sprintf("  Significant: %d (unadj), %d (FDR)\n", n_sig, n_fdr))
}


# =============================================================================
# 5. FOREST PLOTS (Fig 4)
# =============================================================================
cat("\n=== Generating forest plots ===\n")

forest_plots <- lapply(names(analyses), function(aname) {
  a <- analyses[[aname]]
  create_forest_plot(
    cox_results[[aname]]$results,
    title = paste0(a$outcome, " - ", a$timepoint)
  )
})
names(forest_plots) <- names(analyses)

pdf(file.path(figureDir_main, 'Fig4_forest_plots_cox.pdf'), width = 10, height = 10)
grid.arrange(
  forest_plots$CRS_DA, forest_plots$ICANS_DA,
  forest_plots$CRS_D0, forest_plots$ICANS_D0,
  ncol = 2
)
dev.off()
cat("Saved: Fig4_forest_plots_cox.pdf\n")


# =============================================================================
# 6. ADJUSTED CUMULATIVE INCIDENCE CURVES (for significant proteins)
# =============================================================================
cat("\n=== Generating cumulative incidence curves ===\n")

for (aname in names(analyses)) {
  a <- analyses[[aname]]
  res <- cox_results[[aname]]
  
  sig_results <- res$results[res$results$pval < 0.05 & !is.na(res$results$pval), ]
  sig_results <- sig_results[order(sig_results$HR, decreasing = TRUE), ]
  sig_targets <- sig_results$target
  
  if (length(sig_targets) == 0) {
    cat(sprintf("  %s: no significant targets\n", aname))
    next
  }
  
  cat(sprintf("  %s: %d significant targets\n", aname, length(sig_targets)))
  
  cuminc_plots <- list()
  for (target in sig_targets) {
    tryCatch({
      cuminc_plots[[target]] <- create_cuminc_plot(
        target_name = target,
        protein_data = res$data,
        sample_data = res$samples,
        time_var = a$time_var,
        event_var = a$event_var,
        cox_results = res$results,
        outcome_label = a$outcome,
        timepoint_label = a$timepoint
      )
    }, error = function(e) {
      cat(sprintf("    ERROR for %s: %s\n", target, e$message))
    })
  }
  
  if (length(cuminc_plots) > 0) {
    filename <- sprintf("ExtDataFig_cuminc_%s.pdf", aname)
    pdf(file.path(figureDir_extended, filename),
        width = 32, height = ceiling(length(cuminc_plots) / 8) * 4)
    print(wrap_plots(cuminc_plots, ncol = 8))
    dev.off()
    cat(sprintf("  Saved: %s\n", filename))
  }
}


# =============================================================================
# 7. SAVE RESULTS
# =============================================================================
cat("\n=== Saving results ===\n")

cox_results_output <- lapply(cox_results, function(x) x$results)
save(cox_results_output,
     file = file.path(dataDir_public, "cox_model_results.RData"))
cat("Saved: cox_model_results.RData\n")


# =============================================================================
# 8. SUPPLEMENTARY TABLE 5: Cox PH model results
# =============================================================================
cat("\n=== Generating Supplementary Table 5 ===\n")

write_supp_table(
  list(
    "CRS_DayOfApheresis" = cox_results$CRS_DA$results,
    "CRS_Day0"           = cox_results$CRS_D0$results,
    "ICANS_DayOfApheresis" = cox_results$ICANS_DA$results,
    "ICANS_Day0"           = cox_results$ICANS_D0$results
  ),
  "SuppTable5_cox_PH_results.xlsx"
)


cat("\n", strrep("=", 60), "\n")
cat("s07_cox_models.R complete.\n")
cat(strrep("=", 60), "\n")
for (aname in names(cox_results)) {
  res <- cox_results[[aname]]$results
  n_sig <- sum(res$pval < 0.05, na.rm = TRUE)
  n_fdr <- sum(res$pval_FDR < 0.05, na.rm = TRUE)
  cat(sprintf("  %s: %d sig (unadj), %d sig (FDR) of %d proteins\n",
              aname, n_sig, n_fdr, nrow(res)))
}
