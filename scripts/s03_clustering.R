################################################################################
# NULISA-CART-Proteome
# s03_clustering.R
#
# K-means clustering of LME model coefficients to identify groups of proteins
# with similar temporal dynamics. 
#
# Generates:
#   - Trajectory plots for each model set (Figs 2B, 3A, 3D, 5B)
#   - Cluster assignment tables
#   - Supplementary Table 3: Cluster assignments
#
# Inputs:
#   - data/public/LME_coefs_for_clustering.RData (from s02)
#
# Outputs:
#   - data/public/cluster_results.RData (cluster assignments + named mappings)
#   - supplementary_tables/SuppTable3_cluster_assignments.xlsx
################################################################################

source(file.path("scripts", "s00_setup.R"))

# Load LME coefficients from s02
load(file.path(dataDir_public, "LME_coefs_for_clustering.RData"))
# Contains: LME_coefs, LME_sig_proteins, LME_LRT_stats


# =============================================================================
# 1. CLUSTER NAME CONFIGURATION
# =============================================================================
# After the first run, inspect the trajectory plots and fill in the biological
# names for each k-means cluster number. The 'kmeans_order' vector maps 
# k-means cluster numbers to the desired display order (1 = first displayed).
#
# Example: if k-means cluster 4 should be displayed first and named 
# "Myeloid Alarm", set kmeans_order[4] = 1 and names[4] = "1. Myeloid Alarm"
#
# Set APPLY_NAMES = TRUE once the mapping is verified.

APPLY_NAMES <- TRUE  # Set to FALSE for first run to inspect clusters

cluster_config <- list(
  
  # --- Overall time trend ---
  time = list(
    names = c(
      "Myeloid Alarm", 
      "Vascular Collapse",
      "Niche Collapse",  
      "The Brake", 
      "Peak Adaptive Activation",   
      "The Cytokine Storm"          
    ),
    display_order = c(1, 3, 2, 5, 6, 4)  
  ),
  
  # --- CRS × time interaction ---
  CRS = list(
    names = c(
      "Innate-to-Adaptive Bridge",  
      "Vascular & Neuro-Supportive Failure", 
      "Dysregulated Repair",  
      "Exhaustion & Chronic Inflammation",   
      "Myeloid Resolution Failure",           
      "Th1/Th2 Effector Hyper-Activation"         
    ),
    display_order = c(1, 6, 4, 2, 5, 3)
  ),
  
  # --- ICANS × time interaction ---
  ICANS = list(
    names = c(
      "Exhaustion & Chronic Inflammation",
      "Effector Cytokine & Chemokine Activation",
      "Immune Checkpoint & Myeloid Remodeling",
      "Neuroinflammation",
      "Niche Homeostasis",
      "Neuro-Recovery Failure"
    ),
    display_order = c(2, 1, 3, 4, 5, 6)
  ),
  
  # --- Expansion × time interaction ---
  expans = list(
    names = c(
      "Exhaustion Management",                    
      "The Effector Engine",                     
      "Baseline Niche Displacement",         
      "Interferons & Cytotoxicity",    
      "Metabolic & Growth Support",             
      "Early Chemotaxis & Innate Priming"       
    ),
    display_order = c(6, 2, 4, 1, 3, 5)
  )
)


# =============================================================================
# 2. CLUSTERING FUNCTION
# =============================================================================

do_clustering <- function(coef_df, lrt_stats, 
                          scale_data = FALSE,
                          p_threshold = 0.05,
                          title_label = '') {
  
  # 1. Prepare coefficient matrix
  coef_matrix <- coef_df
  rownames(coef_matrix) <- coef_df$target
  coef_matrix$target <- NULL
  
  # Remove the LRT p-value column if present
  if ("Chisq_test_pval" %in% colnames(coef_matrix)) {
    coef_matrix$Chisq_test_pval <- NULL
  }
  
  # 2. Filter to significant proteins (LRT chi-squared test)
  sig_targets <- lrt_stats$target[lrt_stats$Chisq_test_pval < p_threshold]
  coef_matrix <- coef_matrix[rownames(coef_matrix) %in% sig_targets, , drop = FALSE]
  
  # 3. Remove zero-variance rows
  row_vars <- apply(coef_matrix, 1, var, na.rm = TRUE)
  coef_matrix <- coef_matrix[row_vars > 0, , drop = FALSE]
  
  cat(sprintf("  %s: %d proteins after filtering (LRT p < %.2f)\n", 
              title_label, nrow(coef_matrix), p_threshold))
  
  # 4. Z-score if requested (time model only)
  if (scale_data) {
    coef_matrix <- t(scale(t(coef_matrix)))
  }
  
  # 5. Determine k based on number of proteins
  k <- get_n_clusters(nrow(coef_matrix))
  cat(sprintf("  Using k = %d clusters\n", k))
  
  # 6. K-means clustering
  set.seed(123)
  km <- kmeans(as.matrix(coef_matrix), centers = k, nstart = 25)
  
  # 7. Create cluster assignment data frame
  cluster_df <- data.frame(
    target = rownames(coef_matrix),
    cluster = km$cluster,
    stringsAsFactors = FALSE
  )
  
  return(list(
    df = cluster_df, 
    coef_matrix = coef_matrix,
    k = k
  ))
}


# ---------------------------------------------------------------------------
# traj_plot: trajectory plot with cluster naming and reordering
# Produces a single-row faceted plot (ncol = k) matching manuscript figures.
# ---------------------------------------------------------------------------
traj_plot <- function(cluster_result, config = NULL, 
                       center_baseline = FALSE,
                       title = "", 
                       ylab = "Coefficient",
                       xlab = NULL,
                       label_wrap_width = 25) {
  
  coef_matrix <- cluster_result$coef_matrix
  cluster_df  <- cluster_result$df
  k           <- cluster_result$k
  
  # Build plot data: join coefficients with cluster assignments (only target + cluster)
  plot_data <- as.data.frame(coef_matrix) %>%
    rownames_to_column("target") %>%
    left_join(cluster_df[, c("target", "cluster")], by = "target")
  
  # Apply naming and reordering if config is provided
  if (!is.null(config) && !all(grepl("^Cluster \\d+$", config$names))) {
    # Create mapping: k-means cluster → display number and name
    mapping <- data.frame(
      cluster = config$display_order,
      display_number = seq_along(config$display_order),
      stringsAsFactors = FALSE
    )
    # Build display labels from names in display order
    clean_names <- gsub("^\\d+\\.\\s*", "", config$names)
    mapping$cluster_label <- paste0(mapping$display_number, ". ", 
                                     clean_names[config$display_order])
    
    plot_data <- plot_data %>% left_join(mapping, by = "cluster")
  } else {
    # Default: use k-means numbers as-is
    plot_data$display_number <- plot_data$cluster
    plot_data$cluster_label <- paste("Cluster", plot_data$cluster)
  }
  
  # Pivot to long format
  plot_data <- plot_data %>%
    pivot_longer(
      cols = -c(target, cluster, display_number, cluster_label),
      names_to = "time",
      values_to = "coefficient"
    )
  
  # Optional: center at baseline (subtract cluster-mean DA value)
  if (center_baseline) {
    da_col <- colnames(coef_matrix)[1]  # first column = DA/baseline
    baseline_means <- plot_data %>%
      filter(time == da_col) %>%
      group_by(display_number) %>%
      summarize(cluster_baseline = mean(coefficient, na.rm = TRUE), .groups = "drop")
    
    plot_data <- plot_data %>%
      left_join(baseline_means, by = "display_number") %>%
      mutate(coefficient = coefficient - cluster_baseline) %>%
      dplyr::select(-cluster_baseline)
  }
  
  # Set factor levels for proper ordering
  plot_data$time <- factor(plot_data$time, levels = colnames(coef_matrix))
  plot_data$cluster_label <- factor(plot_data$cluster_label,
                                     levels = unique(plot_data$cluster_label[
                                       order(plot_data$display_number)]))
  
  # Color palette (matches s18)
  cluster_colors <- c("1" = "#F8766D", "2" = "#B79F00", "3" = "#00BA38",
                       "4" = "#00BFC4", "5" = "#619CFF", "6" = "#F564E3")
  
  # Create plot
  p <- ggplot(plot_data, aes(x = time, y = coefficient, group = target, 
                              color = factor(display_number))) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    geom_line(alpha = 0.5) +
    facet_wrap(~cluster_label, scales = "fixed", ncol = k,
               labeller = labeller(cluster_label = label_wrap_gen(width = label_wrap_width))) +
    scale_color_manual(values = cluster_colors[1:k]) +
    theme_bw() +
    scale_x_discrete(labels = TIMEPOINT_LABELS_SHORT) +
    labs(title = title, y = ylab, x = xlab) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      strip.text = element_text(size = 9)
    )
  
  return(list(plot_data = plot_data, plot = p))
}


# =============================================================================
# 3. RUN CLUSTERING FOR ALL 4 MODEL SETS
# =============================================================================
cat("=== Running k-means clustering ===\n")

cluster_results <- list(
  time = do_clustering(
    LME_coefs$time_coefs, LME_LRT_stats$time,
    scale_data = TRUE,  # z-score for time model only
    title_label = "Overall Time Trend"
  ),
  CRS = do_clustering(
    LME_coefs$CRS_coefs, LME_LRT_stats$CRS,
    scale_data = FALSE,
    title_label = "Severe CRS × Time Interaction"
  ),
  ICANS = do_clustering(
    LME_coefs$ICANS_coefs, LME_LRT_stats$ICANS,
    scale_data = FALSE,
    title_label = "Severe ICANS × Time Interaction"
  ),
  expans = do_clustering(
    LME_coefs$expans_coefs, LME_LRT_stats$expans,
    scale_data = FALSE,
    title_label = "CAR T-cell Expansion × Time Interaction"
  )
)


# =============================================================================
# 4. PRINT CLUSTER SUMMARIES (for verifying name assignments)
# =============================================================================
cat("\n=== Cluster summaries ===\n")

for (analysis in names(cluster_results)) {
  cat(sprintf("\n--- %s (k=%d) ---\n", analysis, cluster_results[[analysis]]$k))
  
  summary_df <- cluster_results[[analysis]]$df %>%
    group_by(cluster) %>%
    summarise(
      n_proteins = n(),
      proteins = paste(sort(target), collapse = ", "),
      .groups = "drop"
    )
  
  for (i in 1:nrow(summary_df)) {
    cat(sprintf("  Cluster %d (%d proteins): %s\n", 
                summary_df$cluster[i], summary_df$n_proteins[i],
                substr(summary_df$proteins[i], 1, 80)))
    if (nchar(summary_df$proteins[i]) > 80) cat("    ...\n")
  }
}


# =============================================================================
# 5. APPLY BIOLOGICAL NAMES AND REORDER (if APPLY_NAMES = TRUE)
# =============================================================================

apply_cluster_names <- function(cluster_df, config) {
  # Map k-means cluster numbers to biological names
  cluster_df$cluster_name <- config$names[cluster_df$cluster]
  
  # Map to display order
  # display_order[i] = which k-means cluster is displayed at position i
  # So k-means cluster display_order[1] gets display number 1, etc.
  km_to_display <- setNames(seq_along(config$display_order), config$display_order)
  cluster_df$display_number <- km_to_display[as.character(cluster_df$cluster)]
  
  # Create final display label: "1. Biological Name"
  # Extract the name without any existing number prefix
  clean_names <- gsub("^\\d+\\.\\s*", "", config$names)
  named_by_display <- setNames(clean_names, seq_along(config$display_order))
  # Map: for each protein, get its display number, then look up the name
  # The name at display position i comes from k-means cluster display_order[i]
  display_name_map <- setNames(
    paste0(seq_along(config$display_order), ". ", clean_names[config$display_order]),
    config$display_order
  )
  cluster_df$cluster_label <- display_name_map[as.character(cluster_df$cluster)]
  
  return(cluster_df)
}

if (APPLY_NAMES) {
  cat("\n=== Applying biological cluster names ===\n")
  
  for (analysis in names(cluster_results)) {
    config <- cluster_config[[analysis]]
    has_real_names <- !all(grepl("^Cluster \\d+$", config$names))
    
    if (has_real_names) {
      cluster_results[[analysis]]$df <- apply_cluster_names(
        cluster_results[[analysis]]$df, 
        config
      )
      cat(sprintf("  %s: %d clusters named\n", analysis,
                  length(unique(cluster_results[[analysis]]$df$cluster_label))))
    } else {
      # Default placeholders for analyses not yet named (e.g., ICANS)
      cluster_results[[analysis]]$df$cluster_name <- paste("Cluster", cluster_results[[analysis]]$df$cluster)
      cluster_results[[analysis]]$df$display_number <- cluster_results[[analysis]]$df$cluster
      cluster_results[[analysis]]$df$cluster_label <- paste("Cluster", cluster_results[[analysis]]$df$cluster)
      cat(sprintf("  %s: using default cluster numbers (names not yet assigned)\n", analysis))
    }
  }
} else {
  cat("\n=== APPLY_NAMES = FALSE: using default cluster numbers ===\n")
  cat("  Inspect trajectory plots, then update cluster_config and set APPLY_NAMES = TRUE\n")
  
  # Add placeholder labels
  for (analysis in names(cluster_results)) {
    cluster_results[[analysis]]$df$cluster_name <- paste("Cluster", cluster_results[[analysis]]$df$cluster)
    cluster_results[[analysis]]$df$display_number <- cluster_results[[analysis]]$df$cluster
    cluster_results[[analysis]]$df$cluster_label <- paste("Cluster", cluster_results[[analysis]]$df$cluster)
  }
}


# =============================================================================
# 6. PROTEIN-TO-GENE TRANSLATION (for pathway analysis in s04)
# =============================================================================
cat("\n=== Translating protein names to gene symbols ===\n")

# Load target_metadata for the translation map
load(file.path(dataDir_private, "cleaned_data.RData"))  # for target_metadata

translate_proteins <- function(protein_vector, target_meta) {
  gene_list <- c()
  for (p in protein_vector) {
    idx <- which(target_meta$target == p)
    if (length(idx) > 0) {
      # Split multi-gene symbols (heterodimers: "INHBA; INHBB")
      genes <- trimws(unlist(strsplit(target_meta$gene_symbol[idx], ";")))
      gene_list <- c(gene_list, genes)
    } else {
      gene_list <- c(gene_list, p)  # fallback: use target name
    }
  }
  return(sort(unique(gene_list)))
}

# All panel genes (background for enrichment)
all_panel_genes <- translate_proteins(target_metadata$target, target_metadata)
cat(sprintf("  Panel background: %d unique genes\n", length(all_panel_genes)))

# Per-cluster gene lists
cluster_gene_lists <- lapply(names(cluster_results), function(analysis) {
  df <- cluster_results[[analysis]]$df
  clusters <- split(df$target, df$cluster)
  lapply(clusters, translate_proteins, target_meta = target_metadata)
})
names(cluster_gene_lists) <- names(cluster_results)


# =============================================================================
# 7. GENERATE AND SAVE TRAJECTORY PLOTS
# =============================================================================
cat("\n=== Generating trajectory plots ===\n")

# Use config if names are applied, otherwise NULL for defaults
get_config <- function(analysis) {
  if (APPLY_NAMES) cluster_config[[analysis]] else NULL
}

time_traj <- traj_plot(
  cluster_results$time, config = get_config("time"),
  center_baseline = TRUE,
  title = "Overall Time Trend",
  ylab = "Standardized coefficient"
)

CRS_traj <- traj_plot(
  cluster_results$CRS, config = get_config("CRS"),
  center_baseline = FALSE,
  title = "Severe (Grade 2-4) vs. Non-Severe (Grade 0-1) CRS",
  ylab = expression(Log[2]~Fold~Change)
)

ICANS_traj <- traj_plot(
  cluster_results$ICANS, config = get_config("ICANS"),
  center_baseline = FALSE,
  title = "Severe (Grade 1-4) vs. Non-Severe (Grade 0) ICANS",
  ylab = expression(Log[2]~Fold~Change)
)

expans_traj <- traj_plot(
  cluster_results$expans, config = get_config("expans"),
  center_baseline = FALSE,
  title = "High vs. Low CAR T-cell Expansion at Day 6-9",
  ylab = expression(Log[2]~Fold~Change)
)

# Save plots (single-row layout: width=14, height=3)
traj_plots <- list(
  time   = list(traj = time_traj,   file = "Fig2B_clusters_time.pdf"),
  CRS    = list(traj = CRS_traj,    file = "Fig3A_clusters_CRS.pdf"),
  ICANS  = list(traj = ICANS_traj,  file = "Fig3D_clusters_ICANS.pdf"),
  expans = list(traj = expans_traj, file = "Fig5B_clusters_expansion.pdf")
)

for (analysis in names(traj_plots)) {
  spec <- traj_plots[[analysis]]
  pdf(file.path(figureDir_main, spec$file), width = 14, height = 3)
  print(spec$traj$plot)
  dev.off()
  cat(sprintf("  Saved: %s\n", spec$file))
}


# =============================================================================
# 8. SAVE RESULTS
# =============================================================================
cat("\n=== Saving results ===\n")

save(cluster_results, cluster_gene_lists, all_panel_genes, cluster_config,
     file = file.path(dataDir_public, "cluster_results.RData"))
cat("Saved: cluster_results.RData\n")


# =============================================================================
# 9. SUPPLEMENTARY TABLE 3: Cluster assignments
# =============================================================================
cat("\n=== Generating Supplementary Table 3 ===\n")

supp3_sheets <- lapply(names(cluster_results), function(analysis) {
  cluster_results[[analysis]]$df %>%
    arrange(display_number, target) %>%
    dplyr::select(target, cluster, cluster_label)
})
names(supp3_sheets) <- c("Overall_Time", "CRS_Interaction", 
                          "ICANS_Interaction", "Expansion_Interaction")

write_supp_table(supp3_sheets, "SuppTable3_cluster_assignments.xlsx")


# =============================================================================
# 10. SUMMARY
# =============================================================================
cat("\n", strrep("=", 60), "\n")
cat("s03_clustering.R complete.\n")
cat(strrep("=", 60), "\n")
cat(sprintf("  APPLY_NAMES = %s\n", APPLY_NAMES))
for (analysis in names(cluster_results)) {
  n_proteins <- nrow(cluster_results[[analysis]]$df)
  k <- cluster_results[[analysis]]$k
  has_names <- !all(grepl("^Cluster \\d+$", cluster_config[[analysis]]$names))
  status <- if (has_names) "named" else "DEFAULT"
  cat(sprintf("  %s: %d proteins → %d clusters [%s]\n", analysis, n_proteins, k, status))
}

# Check for any analyses still needing names
unnamed <- names(cluster_config)[sapply(cluster_config, function(cfg) 
  all(grepl("^Cluster \\d+$", cfg$names)))]
if (length(unnamed) > 0) {
  cat(sprintf("\n  TODO: Assign biological names for: %s\n", paste(unnamed, collapse = ", ")))
}
