################################################################################
# NULISA-CART-Proteome
# s05_heatmaps.R
#
# Annotated coefficient heatmaps with cluster assignments and pathway
# enrichment annotations for each model set.
#
# Generates:
#   - Fig 2C: Time effect coefficient heatmap
#   - Fig 3B/E: CRS/ICANS interaction heatmaps 
#   - Fig 5C: Expansion interaction heatmap
#
# Inputs:
#   - data/public/cluster_results.RData (from s03)
#   - data/public/pathway_enrichment_results.RData (from s04)
#   - data/private/cleaned_data.RData (for target_metadata)
################################################################################

source(file.path("scripts", "s00_setup.R"))

library(rcartocolor)

# Load data
load(file.path(dataDir_public, "cluster_results.RData"))
load(file.path(dataDir_public, "pathway_enrichment_results.RData"))
load(file.path(dataDir_private, "cleaned_data.RData"))  # for target_metadata


# =============================================================================
# 1. PREPARE PATHWAY ANNOTATIONS
# =============================================================================
cat("=== Preparing pathway annotations ===\n")

# --- Gene-to-protein back-translation ---
# Pathway geneID fields contain gene symbols. We need to map these back to
# NULISAseq protein names for the heatmap annotations.

# Build reverse translation: gene symbol → protein target name(s)
build_reverse_map <- function(target_meta) {
  rev_map <- list()
  for (i in seq_len(nrow(target_meta))) {
    genes <- trimws(unlist(strsplit(target_meta$gene_symbol[i], ";")))
    target <- target_meta$target[i]
    for (g in genes) {
      rev_map[[g]] <- c(rev_map[[g]], target)
    }
  }
  rev_map
}

reverse_map <- build_reverse_map(target_metadata)

# Convert geneID string to protein names, handling heterodimers
genes_to_proteins <- function(gene_string, rev_map, valid_targets) {
  if (is.na(gene_string) || gene_string == "") return(character(0))
  
  genes <- unlist(strsplit(as.character(gene_string), "/"))
  proteins <- unique(unlist(lapply(genes, function(g) rev_map[[g]])))
  
  # Keep only proteins that are in the cluster
  proteins <- proteins[proteins %in% valid_targets]
  return(proteins)
}


# --- Select top 5 pathways per cluster and expand to protein-level ---
prepare_pathway_matrix <- function(master_results, analysis_name,
                                    cluster_df, top_n = 5,
                                    pval_cutoff = 0.1) {
  
  valid_targets <- cluster_df$target
  
  # Select top pathways
  top_pathways <- master_results %>%
    filter(Analysis == analysis_name, pvalue < pval_cutoff) %>%
    group_by(Cluster) %>%
    arrange(pvalue) %>%
    slice_head(n = top_n) %>%
    ungroup()
  
  if (nrow(top_pathways) == 0) {
    # Return empty matrix if no significant pathways
    empty_mat <- matrix(FALSE, nrow = length(valid_targets), ncol = 0)
    rownames(empty_mat) <- valid_targets
    return(as.data.frame(empty_mat))
  }
  
  # Create pathway display name
  top_pathways$pathway_name <- paste0(top_pathways$Database, ": ", top_pathways$Description)
  
  # Expand to protein-level: one row per protein × pathway
  # IMPORTANT: keep Cluster column so we only mark proteins in the cluster
  # where the pathway was enriched
  pathway_proteins <- top_pathways %>%
    rowwise() %>%
    mutate(target = list(genes_to_proteins(geneID, reverse_map, valid_targets))) %>%
    tidyr::unnest(target) %>%
    ungroup() %>%
    dplyr::select(target, pathway_name, Cluster) %>%
    distinct()
  
  # Filter: only keep protein-pathway pairs where the protein belongs to
  # the cluster in which the pathway was enriched
  pathway_proteins <- pathway_proteins %>%
    left_join(cluster_df[, c("target", "cluster_label")], by = "target") %>%
    filter(Cluster == cluster_label) %>%
    dplyr::select(target, pathway_name) %>%
    distinct()
  
  # Pivot to wide: proteins × pathways (TRUE/FALSE)
  if (nrow(pathway_proteins) == 0) {
    empty_mat <- matrix(FALSE, nrow = length(valid_targets), ncol = 0)
    rownames(empty_mat) <- valid_targets
    return(as.data.frame(empty_mat))
  }
  
  pathway_matrix <- pathway_proteins %>%
    mutate(is_present = TRUE) %>%
    pivot_wider(names_from = pathway_name, values_from = is_present, values_fill = FALSE) %>%
    as.data.frame()
  
  rownames(pathway_matrix) <- pathway_matrix$target
  pathway_matrix$target <- NULL
  
  # Add missing targets (proteins in cluster but not in any pathway)
  missing <- setdiff(valid_targets, rownames(pathway_matrix))
  if (length(missing) > 0) {
    missing_rows <- matrix(FALSE, nrow = length(missing), ncol = ncol(pathway_matrix),
                            dimnames = list(missing, colnames(pathway_matrix)))
    pathway_matrix <- rbind(pathway_matrix, missing_rows)
  }
  
  return(pathway_matrix)
}

# Analysis name mapping
analysis_names <- c(time = "TimeEffect", CRS = "SevereCRS",
                     ICANS = "SevereICANS", expans = "CARTExpansion")

# Build pathway matrices for all analyses
pathway_matrices <- lapply(names(analysis_names), function(key) {
  aname <- analysis_names[key]
  cat(sprintf("  %s: ", aname))
  mat <- prepare_pathway_matrix(master_results, aname,
                                 cluster_results[[key]]$df,
                                 top_n = 5,
                                 pval_cutoff = PATHWAY_P_THRESHOLD)
  cat(sprintf("%d pathways × %d proteins\n", ncol(mat), nrow(mat)))
  mat
})
names(pathway_matrices) <- names(analysis_names)


# =============================================================================
# 2. HEATMAP FUNCTION
# =============================================================================

create_coef_heatmap <- function(cluster_result, pathway_matrix, config,
                                 center_DA = FALSE,
                                 cluster_within = TRUE,
                                 col_breaks = c(-4, 0, 4),
                                 legend_title = "Standardized\ncoefficient",
                                 protein_fontsize = 4,
                                 protein_name_height = 1.25) {
  
  coef_matrix <- cluster_result$coef_matrix
  cluster_df  <- cluster_result$df
  
  # --- Cluster colors ---
  cluster_colors <- c("1" = "#F8766D", "2" = "#B79F00", "3" = "#00BA38",
                       "4" = "#00BFC4", "5" = "#619CFF", "6" = "#F564E3")
  
  # --- Pathway annotation colors ---
  pathway_colors_pool <- c(
    carto_pal(n = 12, "Safe"),
    palette.colors(palette = "Okabe-Ito")[-1],
    "#8c564b", "#e377c2", "#bcbd22", "#17becf", "#000000",
    "#332288", "#117733", "#999933", "#882255"
  )
  pathway_colors_pool <- unique(pathway_colors_pool)
  
  # --- Build merged data: coefs + cluster assignment ---
  coefs_df <- as.data.frame(coef_matrix) %>%
    rownames_to_column("target") %>%
    left_join(cluster_df[, c("target", "cluster", "display_number")], by = "target")
  
  # Sort by display number (reversed for landscape orientation)
  coefs_df <- coefs_df[order(coefs_df$display_number, decreasing = TRUE), ]
  rownames(coefs_df) <- coefs_df$target
  
  # Extract coefficient matrix
  time_cols <- setdiff(colnames(coef_matrix), c("target"))
  coefs <- coefs_df[, time_cols]
  colnames(coefs) <- TIMEPOINT_LABELS
  
  # --- Optional: center by cluster-mean DA ---
  if (center_DA) {
    for (d in unique(coefs_df$display_number)) {
      idx <- which(coefs_df$display_number == d)
      baseline_mean <- mean(coefs[idx, 1], na.rm = TRUE)
      coefs[idx, ] <- coefs[idx, ] - baseline_mean
    }
  }
  
  # --- Pathway annotation ---
  # Align pathway_matrix rows to coefs_df order
  pathway_annot <- pathway_matrix[rownames(coefs_df), , drop = FALSE]
  pathway_annot[is.na(pathway_annot)] <- FALSE
  pathway_annot <- apply(pathway_annot, 2, function(x) as.character(as.numeric(x)))
  
  # Sort pathways by primary cluster
  pathway_cluster <- sapply(colnames(pathway_annot), function(pname) {
    present <- as.logical(as.numeric(pathway_annot[, pname]))
    clusters <- coefs_df$display_number[present]
    if (length(clusters) > 0) {
      as.numeric(names(sort(table(clusters), decreasing = TRUE)[1]))
    } else NA
  })
  pathway_annot <- pathway_annot[, order(pathway_cluster, decreasing = FALSE, na.last = TRUE), drop = FALSE]
  
  # --- Annotation objects ---
  cluster_annot <- data.frame(cluster = as.character(coefs_df$display_number))
  pathway_annot_df <- as.data.frame(pathway_annot)
  
  # Cluster color list
  cluster_col_list <- list("cluster" = cluster_colors)
  
  # Pathway color list
  pathway_col_list <- list()
  for (i in seq_along(colnames(pathway_annot))) {
    ci <- ((i - 1) %% length(pathway_colors_pool)) + 1
    pathway_col_list[[colnames(pathway_annot)[i]]] <- c("0" = "grey90", "1" = pathway_colors_pool[ci])
  }
  
  # Top annotation (cluster bar)
  top_annot <- HeatmapAnnotation(
    df = cluster_annot,
    col = cluster_col_list,
    show_legend = FALSE,
    annotation_name_gp = gpar(fontsize = 10),
    annotation_name_side = 'left',
    simple_anno_size = unit(4, "mm"),
    which = "column"
  )
  
  # Bottom annotation (protein names + pathway bars)
  bottom_annot <- HeatmapAnnotation(
    proteins = anno_text(
      rownames(coefs_df),
      gp = gpar(fontsize = protein_fontsize),
      rot = 90, just = "right", location = 1
    ),
    df = pathway_annot_df,
    col = pathway_col_list,
    show_legend = FALSE,
    annotation_name_gp = gpar(fontsize = 10),
    annotation_name_side = "left",
    simple_anno_size = unit(4, "mm"),
    which = "column",
    annotation_height = unit(
      c(protein_name_height, rep(4, ncol(pathway_annot_df))),
      c("cm", rep("mm", ncol(pathway_annot_df)))
    )
  )
  
  # --- Heatmap ---
  col_fun <- colorRamp2(col_breaks, c("blue", "gray90", "red"))
  
  ht <- Heatmap(
    t(as.matrix(coefs)),
    col = col_fun,
    column_title = NULL,
    show_row_names = TRUE,
    show_column_names = FALSE,
    cluster_columns = cluster_within,
    show_column_dend = FALSE,
    cluster_rows = FALSE,
    column_split = cluster_annot$cluster,
    cluster_column_slices = FALSE,
    clustering_method_columns = 'ward.D2',
    row_names_gp = gpar(fontsize = 10),
    top_annotation = top_annot,
    bottom_annotation = bottom_annot,
    show_heatmap_legend = TRUE,
    heatmap_legend_param = list(title = legend_title),
    border = TRUE,
    column_gap = unit(3, "mm"),
    row_gap = unit(0, "mm"),
    row_names_side = 'left'
  )
  
  # --- Dynamic height ---
  n_pathways <- ncol(pathway_annot)
  total_height <- 3 + (protein_name_height / 2.54) + (n_pathways * 0.15)
  
  return(list(heatmap = ht, height = total_height))
}


# =============================================================================
# 3. GENERATE HEATMAPS
# =============================================================================
cat("\n=== Generating heatmaps ===\n")

heatmap_specs <- list(
  time = list(
    file = "Fig2C_heatmap_time.pdf",
    center_DA = TRUE, cluster_within = TRUE,
    col_breaks = c(-4, 0, 4),
    legend_title = "Standardized\ncoefficient"
  ),
  CRS = list(
    file = "Fig3B_heatmap_CRS.pdf",
    center_DA = FALSE, cluster_within = TRUE,
    col_breaks = c(-2, 0, 2),
    legend_title = "Log\u2082 Fold\nChange"
  ),
  ICANS = list(
    file = "Fig3E_heatmap_ICANS.pdf",
    center_DA = FALSE, cluster_within = TRUE,
    col_breaks = c(-2, 0, 2),
    legend_title = "Log\u2082 Fold\nChange"
  ),
  expans = list(
    file = "Fig5C_heatmap_expansion.pdf",
    center_DA = FALSE, cluster_within = TRUE,
    col_breaks = c(-2, 0, 2),
    legend_title = "Log\u2082 Fold\nChange"
  )
)

for (key in names(heatmap_specs)) {
  spec <- heatmap_specs[[key]]
  cat(sprintf("  %s...", key))
  
  result <- create_coef_heatmap(
    cluster_result = cluster_results[[key]],
    pathway_matrix = pathway_matrices[[key]],
    config = cluster_config[[key]],
    center_DA = spec$center_DA,
    cluster_within = spec$cluster_within,
    col_breaks = spec$col_breaks,
    legend_title = spec$legend_title
  )
  
  cairo_pdf(file.path(figureDir_main, spec$file), width = 18, height = result$height)
  draw(result$heatmap, padding = unit(c(2, 10, 2, 4), "mm"), merge_legend = TRUE)
  dev.off()
  
  cat(sprintf(" saved (%s, h=%.1f in)\n", spec$file, result$height))
}


cat("\n", strrep("=", 60), "\n")
cat("s05_heatmaps.R complete.\n")
cat(strrep("=", 60), "\n")
