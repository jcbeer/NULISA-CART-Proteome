################################################################################
# NULISA-CART-Proteome
# s04_pathway_enrichment.R
#
# Over-representation analysis (ORA) for each protein cluster across five
# pathway databases: GO:BP, GO:MF, KEGG, Reactome, MSigDB Hallmark.
#
# Generates:
#   - Enrichment dot plots per model set (Extended Data Figs)
#   - Master enrichment results table
#   - Supplementary Table 4: Pathway enrichment results
#
# Inputs:
#   - data/public/cluster_results.RData (from s03)
#
# Outputs:
#   - data/public/pathway_enrichment_results.RData
#   - supplementary_tables/SuppTable4_pathway_enrichment.xlsx
################################################################################

source(file.path("scripts", "s00_setup.R"))

library(org.Hs.eg.db)
library(ReactomePA)
library(msigdbr)

# Load cluster results from s03
load(file.path(dataDir_public, "cluster_results.RData"))
# Contains: cluster_results, cluster_gene_lists, all_panel_genes, cluster_config


# =============================================================================
# 1. SETUP: BACKGROUND AND HALLMARK GENE SETS
# =============================================================================
cat("=== Setting up enrichment analysis ===\n")

# Convert background genes to Entrez IDs (for KEGG and Reactome)
all_entrez_universe <- bitr(all_panel_genes, fromType = "SYMBOL",
                             toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

cat(sprintf("  Background: %d gene symbols → %d Entrez IDs\n",
            length(all_panel_genes), length(all_entrez_universe)))

# Build Hallmark term-to-gene mapping (uses gene symbols directly)
hallmark_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  mutate(gs_name = gsub("^HALLMARK_", "", gs_name),
         gs_name = gsub("_", " ", gs_name),
         gs_name = tools::toTitleCase(tolower(gs_name)))

cat(sprintf("  Hallmark gene sets: %d\n", n_distinct(hallmark_t2g$gs_name)))


# =============================================================================
# 2. ENRICHMENT FUNCTION (single cluster)
# =============================================================================

run_enrichment <- function(gene_symbols, cluster_name = "",
                            pval_cutoff = 0.1, min_gs = 10, max_gs = 200) {
  
  cat(sprintf("  %s (%d genes, minGS=%d)...", cluster_name, length(gene_symbols), min_gs))
  
  results <- list()
  
  # Convert to Entrez
  entrez_genes <- tryCatch(
    bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
    error = function(e) data.frame(SYMBOL = character(), ENTREZID = character())
  )
  entrez_ids <- entrez_genes$ENTREZID
  
  # --- GO: Biological Process ---
  tryCatch({
    ego_bp <- enrichGO(gene = gene_symbols, universe = all_panel_genes,
                        OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                        pAdjustMethod = "BH", pvalueCutoff = pval_cutoff,
                        qvalueCutoff = 1, minGSSize = min_gs, maxGSSize = max_gs)
    if (!is.null(ego_bp) && nrow(ego_bp@result) > 1) {
      ego_bp <- clusterProfiler::simplify(ego_bp, cutoff = 0.7, by = "p.adjust", select_fun = min)
    }
    if (!is.null(ego_bp) && nrow(ego_bp@result) > 0) {
      results$GO_BP <- ego_bp@result %>% filter(pvalue < pval_cutoff) %>%
        mutate(Database = "GO_BP", Cluster = cluster_name)
    }
  }, error = function(e) NULL)
  
  # --- GO: Molecular Function ---
  tryCatch({
    ego_mf <- enrichGO(gene = gene_symbols, universe = all_panel_genes,
                        OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "MF",
                        pAdjustMethod = "BH", pvalueCutoff = pval_cutoff,
                        qvalueCutoff = 1, minGSSize = min_gs, maxGSSize = max_gs)
    if (!is.null(ego_mf) && nrow(ego_mf@result) > 1) {
      ego_mf <- clusterProfiler::simplify(ego_mf, cutoff = 0.7, by = "p.adjust", select_fun = min)
    }
    if (!is.null(ego_mf) && nrow(ego_mf@result) > 0) {
      results$GO_MF <- ego_mf@result %>% filter(pvalue < pval_cutoff) %>%
        mutate(Database = "GO_MF", Cluster = cluster_name)
    }
  }, error = function(e) NULL)
  
  # --- KEGG ---
  tryCatch({
    kegg_res <- enrichKEGG(gene = entrez_ids, universe = all_entrez_universe,
                            organism = "hsa", pvalueCutoff = pval_cutoff,
                            qvalueCutoff = 1, minGSSize = min_gs, maxGSSize = max_gs)
    if (!is.null(kegg_res) && nrow(kegg_res@result) > 0) {
      kegg_res <- setReadable(kegg_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      results$KEGG <- kegg_res@result %>% filter(pvalue < pval_cutoff) %>%
        mutate(Database = "KEGG", Cluster = cluster_name)
    }
  }, error = function(e) NULL)
  
  # --- Reactome ---
  tryCatch({
    react_res <- enrichPathway(gene = entrez_ids, universe = all_entrez_universe,
                                pvalueCutoff = pval_cutoff, qvalueCutoff = 1,
                                minGSSize = min_gs, maxGSSize = max_gs, readable = TRUE)
    if (!is.null(react_res) && nrow(react_res@result) > 0) {
      results$Reactome <- react_res@result %>% filter(pvalue < pval_cutoff) %>%
        mutate(Database = "Reactome", Cluster = cluster_name)
    }
  }, error = function(e) NULL)
  
  # --- Hallmark (MSigDB) ---
  tryCatch({
    hallmark_res <- enricher(gene = gene_symbols, universe = all_panel_genes,
                              TERM2GENE = hallmark_t2g, pvalueCutoff = pval_cutoff,
                              qvalueCutoff = 1, minGSSize = min_gs, maxGSSize = max_gs)
    if (!is.null(hallmark_res) && nrow(hallmark_res@result) > 0) {
      results$Hallmark <- hallmark_res@result %>% filter(pvalue < pval_cutoff) %>%
        mutate(Database = "Hallmark", Cluster = cluster_name)
    }
  }, error = function(e) NULL)
  
  n_total <- sum(sapply(results, nrow))
  cat(sprintf(" %d terms\n", n_total))
  
  return(results)
}


# =============================================================================
# 3. WRAPPER: RUN ENRICHMENT FOR ALL CLUSTERS IN ONE MODEL SET
# =============================================================================

run_analysis <- function(cluster_list, analysis_name, pval_cutoff = 0.1) {
  
  cat(sprintf("\n%s\n%s\n", strrep("=", 60), analysis_name))
  
  all_results <- list()
  
  for (cname in names(cluster_list)) {
    genes <- cluster_list[[cname]]
    
    # Dynamic minGSSize based on cluster size
    min_gs <- ifelse(length(genes) < 15, 3, ifelse(length(genes) <= 30, 5, 10))
    
    enrichment <- run_enrichment(genes, cluster_name = cname,
                                  pval_cutoff = pval_cutoff, min_gs = min_gs)
    
    for (db in names(enrichment)) {
      if (!is.null(enrichment[[db]]) && nrow(enrichment[[db]]) > 0) {
        enrichment[[db]]$Analysis <- analysis_name
        all_results[[paste(cname, db, sep = "_")]] <- enrichment[[db]]
      }
    }
  }
  
  combined <- bind_rows(all_results)
  
  # Print summary
  if (nrow(combined) > 0) {
    cat(sprintf("\n  Total terms: %d\n", nrow(combined)))
    summary_stats <- combined %>%
      group_by(Database) %>%
      summarise(n = n(), .groups = "drop")
    for (i in 1:nrow(summary_stats)) {
      cat(sprintf("    %s: %d\n", summary_stats$Database[i], summary_stats$n[i]))
    }
  }
  
  return(combined)
}


# =============================================================================
# 4. APPLY CLUSTER LABELS FROM cluster_config
# =============================================================================
# Rename cluster_gene_lists from "1", "2", ... to biological names 
# using cluster_config (same config as s03)

analysis_names <- c(time = "TimeEffect", CRS = "SevereCRS",
                     ICANS = "SevereICANS", expans = "CARTExpansion")

labeled_gene_lists <- list()

for (analysis in names(cluster_gene_lists)) {
  config <- cluster_config[[analysis]]
  gene_list <- cluster_gene_lists[[analysis]]
  
  has_real_names <- !all(grepl("^Cluster \\d+$", config$names))
  
  if (has_real_names) {
    # Build display labels: "1. Name" in display order
    clean_names <- gsub("^\\d+\\.\\s*", "", config$names)
    display_labels <- paste0(seq_along(config$display_order), ". ",
                              clean_names[config$display_order])
    
    # Reorder gene list: display_order[i] = k-means cluster at position i
    reordered <- gene_list[as.character(config$display_order)]
    names(reordered) <- display_labels
  } else {
    # Default: keep "Cluster 1", etc.
    names(gene_list) <- paste("Cluster", names(gene_list))
    reordered <- gene_list
  }
  
  labeled_gene_lists[[analysis_names[analysis]]] <- reordered
}

# Verify labeling
for (aname in names(labeled_gene_lists)) {
  cat(sprintf("\n%s clusters:\n", aname))
  for (cname in names(labeled_gene_lists[[aname]])) {
    cat(sprintf("  %s (%d genes)\n", cname, length(labeled_gene_lists[[aname]][[cname]])))
  }
}


# =============================================================================
# 5. RUN ENRICHMENT FOR ALL 4 ANALYSES
# =============================================================================
cat("\n=== Running pathway enrichment ===\n")

all_enrichment_results <- list()

for (aname in names(labeled_gene_lists)) {
  all_enrichment_results[[aname]] <- run_analysis(
    labeled_gene_lists[[aname]], 
    analysis_name = aname,
    pval_cutoff = PATHWAY_P_THRESHOLD  # 0.1, from s00
  )
}


# =============================================================================
# 6. COMBINE INTO MASTER RESULTS
# =============================================================================
cat("\n=== Combining results ===\n")

master_results <- bind_rows(all_enrichment_results)

# Standardize columns (Hallmark lacks some optional columns)
cols_keep <- c("Analysis", "Cluster", "Database", "ID", "Description",
               "geneID", "Count", "GeneRatio", "BgRatio",
               "pvalue", "p.adjust", "qvalue")
cols_present <- intersect(cols_keep, colnames(master_results))
master_results <- master_results[, cols_present]

cat(sprintf("  Total enrichment terms: %d\n", nrow(master_results)))
cat(sprintf("  Terms with p < 0.05: %d\n", sum(master_results$pvalue < 0.05)))

# Summary table
cat("\n  Terms per analysis × database:\n")
print(table(master_results$Analysis, master_results$Database))


# =============================================================================
# 7. VISUALIZATION: ENRICHMENT DOT PLOTS
# =============================================================================
cat("\n=== Generating enrichment plots ===\n")

create_enrichment_plot <- function(results, analysis_name,
                                    pvalue_cutoff = 0.1,
                                    top_n_per_db = 3,
                                    plot_width = 16, plot_height = 10,
                                    label_wrap_width = 25) {
  
  plot_data <- results %>%
    filter(Analysis == analysis_name, pvalue < pvalue_cutoff) %>%
    group_by(Cluster, Database) %>%
    arrange(pvalue) %>%
    slice_head(n = top_n_per_db) %>%
    ungroup()
  
  if (nrow(plot_data) == 0) {
    cat(sprintf("  %s: no terms to plot\n", analysis_name))
    return(NULL)
  }
  
  # Parse GeneRatio to numeric
  plot_data <- plot_data %>%
    mutate(GeneRatio_num = sapply(GeneRatio, function(x) {
      parts <- as.numeric(unlist(strsplit(as.character(x), "/")))
      parts[1] / parts[2]
    })) %>%
    mutate(Description = ifelse(nchar(Description) > 45,
                                 paste0(substr(Description, 1, 42), "..."),
                                 Description))
  
  # Sort by gene ratio within each cluster
  plot_data <- plot_data %>%
    mutate(Description = reorder(Description, GeneRatio_num))
  
  p <- ggplot(plot_data,
              aes(x = GeneRatio_num, y = Description, size = Count, color = pvalue)) +
    geom_point(alpha = 0.8) +
    scale_color_gradient(low = "#B2182B", high = "#D1E5F0", name = "P-value",
                          trans = "log10",
                          breaks = 10^(-6:-1),
                          labels = function(x) parse(text = sprintf("10^%d", round(log10(x))))) +
    scale_size_continuous(name = "Gene\nCount", range = c(3, 8)) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    facet_grid(Database ~ Cluster, scales = "free", space = "free_y",
               labeller = labeller(Cluster = label_wrap_gen(width = label_wrap_width))) +
    theme_bw() +
    labs(title = paste(analysis_name, ": Pathway Enrichment"),
         subtitle = sprintf("Top %d terms per database (p < %.1f)", top_n_per_db, pvalue_cutoff),
         x = "Gene Ratio", y = NULL) +
    theme(
      strip.text.x = element_text(size = 8, face = "bold"),
      strip.text.y = element_text(size = 9, angle = 0, face = "italic"),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 8),
      legend.position = "bottom",
      panel.spacing = unit(0.5, "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  return(p)
}

# Generate and save plots
plot_files <- c(TimeEffect = "ExtDataFig2B_enrichment_time.pdf",
                SevereCRS = "ExtDataFig3B_enrichment_CRS.pdf",
                SevereICANS = "ExtDataFig5B_enrichment_ICANS.pdf",
                CARTExpansion = "ExtDataFig6B_enrichment_expansion.pdf")

for (aname in names(plot_files)) {
  p <- create_enrichment_plot(master_results, aname,
                               pvalue_cutoff = PATHWAY_P_THRESHOLD,
                               top_n_per_db = 3)
  if (!is.null(p)) {
    pdf(file.path(figureDir_extended, plot_files[aname]), width = 16, height = 10)
    print(p)
    dev.off()
    cat(sprintf("  Saved: %s\n", plot_files[aname]))
  }
}


# =============================================================================
# 8. SAVE RESULTS
# =============================================================================
cat("\n=== Saving results ===\n")

save(master_results, labeled_gene_lists, all_panel_genes,
     file = file.path(dataDir_public, "pathway_enrichment_results.RData"))
cat("Saved: pathway_enrichment_results.RData\n")


# =============================================================================
# 9. SUPPLEMENTARY TABLE 4: Pathway enrichment results
# =============================================================================
cat("\n=== Generating Supplementary Table 4 ===\n")

# Sheet per analysis: top 10 per cluster, all databases, p < 0.05
supp4_sheets <- list()

for (aname in names(all_enrichment_results)) {
  sheet_data <- all_enrichment_results[[aname]] %>%
    filter(pvalue < 0.05) %>%
    group_by(Cluster) %>%
    arrange(pvalue) %>%
    slice_head(n = 10) %>%
    ungroup() %>%
    dplyr::select(Cluster, Database, Description, GeneRatio,
                  pvalue, p.adjust, qvalue, geneID, Count) %>%
    arrange(Cluster, pvalue)
  
  supp4_sheets[[aname]] <- sheet_data
}

# Add a master sheet with all results
supp4_sheets[["All_Results"]] <- master_results %>%
  arrange(Analysis, Cluster, pvalue)

write_supp_table(supp4_sheets, "SuppTable4_pathway_enrichment.xlsx")


# =============================================================================
# 10. CLUSTER-PATHWAY SUMMARY TABLES (for supplementary PDF)
# =============================================================================
cat("\n=== Generating cluster-pathway summary tables ===\n")

# For each analysis, create a table with:
#   Cluster | Protein Names | Top Pathways (pooled across databases)

create_summary_table <- function(cluster_df, enrichment_results, analysis_name,
                                  top_n = 5, pval_cutoff = 0.05) {
  
  # Get cluster labels in display order
  cluster_labels <- cluster_df %>%
    arrange(display_number) %>%
    pull(cluster_label) %>%
    unique()
  
  rows <- lapply(cluster_labels, function(clabel) {
    # Get protein (target) names for this cluster
    proteins <- cluster_df %>%
      filter(cluster_label == clabel) %>%
      pull(target) %>%
      sort()
    
    # Get top pathways pooled across databases
    top_paths <- enrichment_results %>%
      filter(Analysis == analysis_name, Cluster == clabel, pvalue < pval_cutoff) %>%
      arrange(pvalue) %>%
      slice_head(n = top_n) %>%
      mutate(entry = paste0(Database, ": ", Description))
    
    pathway_text <- if (nrow(top_paths) > 0) {
      paste(top_paths$entry, collapse = "\n")
    } else {
      ""
    }
    
    data.frame(
      Cluster = clabel,
      N_Proteins = length(proteins),
      Protein_Names = paste(proteins, collapse = ", "),
      Top_Pathways = pathway_text,
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows(rows)
}

summary_sheets <- list()

# Map analysis names back to cluster_results keys
results_keys <- c(TimeEffect = "time", SevereCRS = "CRS",
                   SevereICANS = "ICANS", CARTExpansion = "expans")

for (aname in names(labeled_gene_lists)) {
  rkey <- results_keys[aname]
  summary_sheets[[aname]] <- create_summary_table(
    cluster_results[[rkey]]$df,
    master_results,
    analysis_name = aname,
    top_n = 5, pval_cutoff = PATHWAY_P_THRESHOLD  # 0.1, matching enrichment threshold
  )
}

write_supp_table(summary_sheets,
                  "cluster_pathway_summary.xlsx",
                  dir = suppDir_working)


cat("\n", strrep("=", 60), "\n")
cat("s04_pathway_enrichment.R complete.\n")
cat(strrep("=", 60), "\n")
cat(sprintf("  Total enrichment terms: %d\n", nrow(master_results)))
cat(sprintf("  Databases: GO_BP, GO_MF, KEGG, Reactome, Hallmark\n"))
cat(sprintf("  P-value threshold: %.1f\n", PATHWAY_P_THRESHOLD))
cat(sprintf("  GO simplify cutoff: 0.7\n"))
