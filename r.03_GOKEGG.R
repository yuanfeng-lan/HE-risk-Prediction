rm(list = ls()); gc()
ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "03_GOKEGG")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(ORIGINAL_DIR)

library(data.table)
library(ggplot2)
library(stringr)
library(scales)
library(cowplot)
library(clusterProfiler)
library(enrichplot)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyverse)

# Setup output directories
script.dir <- file.path(output)
files <- list.files(script.dir, pattern = "DEG\\.csv$", full.names = TRUE, ignore.case = TRUE)
output_dir <- file.path(output, "enrich_results")
if (!dir.exists(output_dir)) dir.create(output_dir)

# Helper function to calculate gene ratio
calculate_gene_ratio <- function(ratio_str) {
  if (is.character(ratio_str) && grepl("/", ratio_str)) {
    parts <- strsplit(ratio_str, "/")[[1]]
    num <- as.numeric(parts[1])
    den <- as.numeric(parts[2])
    return(ifelse(den != 0, num / den, 0))
  } else {
    return(as.numeric(ratio_str))
  }
}

# Helper function for enrichment analysis
run_enrichment <- function(genes, sample_name, output_dir) {
  message("Processing: ", sample_name)
  
  # Convert gene symbols to Entrez IDs
  entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # GO enrichment
  go <- enrichGO(gene = entrez$ENTREZID,
                 OrgDb = org.Hs.eg.db,
                 keyType = "ENTREZID",
                 ont = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 1,
                 qvalueCutoff = 1,
                 readable = TRUE)
  
  # KEGG enrichment
  kegg <- enrichKEGG(gene = entrez$ENTREZID,
                     keyType = "kegg",
                     organism = "hsa",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 1)
  
  # Combine results
  go_results <- as.data.frame(go)
  kegg_results <- as.data.frame(kegg)
  kegg_results$ONTOLOGY <- "KEGG"
  
  cols <- c("ONTOLOGY", "ID", "Description", "GeneRatio", "BgRatio", 
            "pvalue", "p.adjust", "qvalue", "geneID", "Count")
  go_results <- go_results[, cols]
  kegg_results <- kegg_results[, cols]
  go_results$ONTOLOGY <- paste0("GO", go_results$ONTOLOGY)
  
  all_results <- rbind(go_results, kegg_results)
  
  # Save all results
  out_path <- file.path(output_dir, sample_name)
  if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)
  write.csv(all_results, file.path(out_path, paste0(sample_name, "_all_results.csv")))
  
  # Filter significant results
  sig_results <- subset(all_results, pvalue <= 0.05)
  
  if (nrow(sig_results) > 0) {
    # Create bubble plot
    create_bubble_plot(sig_results, sample_name, out_path)
  } else {
    message("No significant enrichment found for: ", sample_name)
  }
  
  return(all_results)
}

# Helper function to create bubble plot
create_bubble_plot <- function(results_df, sample_name, out_path) {
  results_df <- results_df %>%
    mutate(GeneRatio_numeric = sapply(GeneRatio, calculate_gene_ratio))
  
  # Select top 5 per ontology
  filtered_df <- results_df %>%
    group_by(ONTOLOGY) %>%
    arrange(desc(GeneRatio_numeric)) %>%
    slice_head(n = 5) %>%
    ungroup()
  
  filtered_df <- filtered_df %>%
    mutate(log10_p_value = -log10(pvalue)) %>%
    mutate(Description_wrapped = str_wrap(Description, width = 40))
  
  # Order factors
  filtered_df <- filtered_df[order(match(filtered_df$ONTOLOGY, unique(filtered_df$ONTOLOGY)), 
                                   filtered_df$GeneRatio_numeric), ]
  filtered_df$Description_wrapped <- factor(filtered_df$Description_wrapped, 
                                           levels = filtered_df$Description_wrapped)
  filtered_df$ONTOLOGY <- gsub("GO", "GO ", filtered_df$ONTOLOGY)
  
  # Create plot
  p <- ggplot(filtered_df, aes(x = GeneRatio_numeric, y = Description_wrapped, 
                               size = Count, color = log10_p_value)) +
    geom_point(alpha = 0.7) +
    scale_color_gradientn(
      name = expression(-log[10](italic(p))),
      colors = c("#3C5488", "#F1C40F", "#DC0000"),
      limits = c(min(filtered_df$log10_p_value), max(filtered_df$log10_p_value)),
      na.value = "grey50"
    ) +
    labs(x = "GeneRatio", y = "Terms", title = "") +
    facet_wrap(~ ONTOLOGY, ncol = 1, scales = "free_y") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 16),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      plot.title = element_text(hjust = 0.5),
      strip.text = element_text(size = 15, face = "bold")
    )
  
  ggsave(file.path(out_path, paste0(sample_name, "_enrichment_bubble_plot.png")), 
         p, width = 8, height = 11)
  ggsave(file.path(out_path, paste0(sample_name, "_enrichment_bubble_plot.pdf")), 
         p, width = 8, height = 11)
}

#### Process DEG files from 01_DEG directory ####
deg_files <- list.files("./01_DEG", pattern = "DEG.*\\.csv$", full.names = TRUE)

for (fl in deg_files) {
  sample <- tools::file_path_sans_ext(basename(fl))
  
  # Read DEG file
  deg_data <- read.csv(fl)
  
  # Extract gene column
  if ("X" %in% names(deg_data)) {
    genes <- deg_data$X
  } else if ("x" %in% names(deg_data)) {
    genes <- deg_data$x
  } else if ("Gene" %in% names(deg_data)) {
    genes <- deg_data$Gene
  } else {
    genes <- rownames(deg_data)
  }
  
  # Run enrichment if genes exist
  if (length(genes) > 0) {
    run_enrichment(genes, sample, output_dir)
  }
}

#### Process GSE15654 specific DEG subsets ####
gse15654_deg <- read.csv(file.path("01_DEG", "01.DEG_GSE15654_Poor_prognosis-Good_prognosis_all.csv"))
gse15654_up <- subset(gse15654_deg, P.Value <= 0.05 & logFC >= 0.5)
gse15654_down <- subset(gse15654_deg, P.Value <= 0.05 & logFC <= -0.5)
gse15654_all <- rbind(gse15654_up, gse15654_down)

# Process up-regulated genes
if (nrow(gse15654_up) > 0) {
  run_enrichment(gse15654_up$X, "GSE15654_DEG_UP", output_dir)
}

# Process down-regulated genes
if (nrow(gse15654_down) > 0) {
  run_enrichment(gse15654_down$X, "GSE15654_DEG_DOWN", output_dir)
}

# Process all genes
if (nrow(gse15654_all) > 0) {
  run_enrichment(gse15654_all$X, "GSE15654_DEG_ALL", output_dir)
}

#### Process common HE genes from previous analysis ####
he_genes_files <- list.files("./02_Mfuzz", pattern = "HE_.*_gene\\.csv$", full.names = TRUE)

for (fl in he_genes_files) {
  sample <- tools::file_path_sans_ext(basename(fl))
  
  # Read gene list
  gene_data <- read.csv(fl)
  if ("genes" %in% names(gene_data)) {
    genes <- gene_data$genes
  } else if ("x" %in% names(gene_data)) {
    genes <- gene_data$x
  } else {
    genes <- gene_data[, 1]
  }
  
  # Run enrichment if genes exist
  if (length(genes) > 0) {
    # Skip if sample already processed
    if (!sample %in% c("GSE41919_GSE57193_UP_gene", "GSE41919_GSE57193_DOWN_gene", 
                       "GSE41919_GSE57193_ALL_gene")) {
      run_enrichment(genes, sample, output_dir)
    }
  }
}

message("Enrichment analysis completed successfully!")
