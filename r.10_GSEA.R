rm(list = ls()); gc()
ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "10_GSEA")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(ORIGINAL_DIR)

library(tidyverse)
library(data.table)
library(msigdbr)
library(fgsea)
library(ggplot2)
library(ggridges)
library(tidyr)

# Load gene sets
c5bp <- tryCatch(
  msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP"),
  error = function(e) NULL
)
save(c5bp, file = file.path(output, "c5bp.Rdata"))

# Prepare pathways
prepare_pathways <- function(c5bp) {
  pathways_list <- list()
  if (!is.null(c5bp)) {
    c5_sets <- split(c5bp$gene_symbol, c5bp$gs_name)
    names(c5_sets) <- paste0("C5_BP|", names(c5_sets))
    pathways_list <- c(pathways_list, c5_sets)
  }
  return(pathways_list)
}

pathways_list <- prepare_pathways(c5bp)

# Load data
exp_df <- read.csv(file.path("09_SEM", "09_SEM__exp_corrected.csv"), row.names = 1)
group_merged <- read.csv(file.path("09_SEM", "09_SEM__group_merged.csv"), row.names = 1)

# Load models
final_model_1 <- readRDS(file.path("07_xgboost", "07_xgboost_GSE139602_final_xgboost_model.Rdata"))
final_model_2 <- readRDS(file.path("07_xgboost", "07_xgboost_GSE15654_final_xgboost_model_cox.Rdata"))

# Load model feature genes
gene1 <- read.csv(file.path("06_RF", "06_RF_GSE139602_final_cogenes.csv"))$x
gene2 <- read.csv(file.path("06_RF", "06_RF_GSE15654_final_cogenes.csv"))$x

# Prepare expression matrix
exp_t <- as.data.frame(t(exp_df))
gene1_avail <- gene1[gene1 %in% colnames(exp_t)]
gene2_avail <- gene2[gene2 %in% colnames(exp_t)]

cat("Model 1 features available:", length(gene1_avail), "/", length(gene1), "\n")
cat("Model 2 features available:", length(gene2_avail), "/", length(gene2), "\n")

# Calculate model prediction scores
exp_t_scaled <- scale(exp_t)
group_merged$risk_score_1 <- predict(final_model_1, as.matrix(exp_t_scaled[, gene1_avail]))
group_merged$risk_score_2 <- predict(final_model_2, as.matrix(exp_t_scaled[, gene2_avail]))

cat("Model 1 score range:", range(group_merged$risk_score_1), "\n")
cat("Model 2 score range:", range(group_merged$risk_score_2), "\n")

# Core function: GSEA based on model scores with gradient permutation testing
run_score_based_gsea <- function(score_vector,
                                 score_name,
                                 expr_matrix,
                                 pathways_list,
                                 output_dir,
                                 labels = "",
                                 nperm_gradients = c(1000, 10000, 100000),
                                 gradient_labels = c("1k", "10k", "100k"),
                                 pval_threshold = 0.01,
                                 top_n = 5) {
  
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("Processing GSEA for:", score_name, "\n")
  cat(rep("=", 60), "\n", sep = "")
  
  # Calculate correlations and rank genes
  expr <- as.matrix(expr_matrix)
  mode(expr) <- "numeric"
  
  cors <- apply(expr, 1, function(x) {
    suppressWarnings(cor(x, score_vector, method = "spearman", use = "pairwise.complete.obs"))
  })
  
  cors <- cors[!is.na(cors)]
  ranked <- sort(cors, decreasing = TRUE)
  
  cat("  Number of genes ranked:", length(ranked), "\n")
  
  # Gradient permutation testing
  gsea_results_list <- list()
  
  for (i in seq_along(nperm_gradients)) {
    cat("  Running GSEA with nperm =", nperm_gradients[i], "(", gradient_labels[i], ")\n")
    
    set.seed(123)
    fgsea_res <- fgsea(
      pathways = pathways_list,
      stats = ranked,
      minSize = 15,
      maxSize = 500,
      nperm = nperm_gradients[i]
    )
    
    fgsea_res <- fgsea_res %>%
      as.data.frame() %>%
      arrange(padj, pval, desc(NES)) %>%
      mutate(
        nperm_setting = nperm_gradients[i],
        gradient_label = gradient_labels[i]
      )
    
    gsea_results_list[[i]] <- fgsea_res
  }
  
  # Combine all gradient results
  all_gradient_results <- bind_rows(gsea_results_list)
  
  all_gradient_results_clean <- all_gradient_results
  list_cols <- names(all_gradient_results_clean)[sapply(all_gradient_results_clean, is.list)]
  for (col in list_cols) {
    all_gradient_results_clean[[col]] <- sapply(all_gradient_results_clean[[col]], function(x) {
      if (is.null(x) || length(x) == 0) return(NA_character_)
      paste(as.character(x), collapse = ";")
    })
  }
  
  write.csv(all_gradient_results_clean, 
            file.path(output_dir, paste0(labels, score_name, "_GSEA_gradient_comparison.csv")), 
            row.names = FALSE)
  
  # Identify significant pathways per gradient
  significant_pathways <- all_gradient_results %>%
    group_by(gradient_label) %>%
    filter(pval < pval_threshold) %>%
    summarise(
      pathways = list(pathway),
      count = n()
    ) %>%
    ungroup()
  
  cat("\n  Significant pathways (pval <", pval_threshold, "):\n")
  print(significant_pathways %>% select(gradient_label, count))
  
  # Calculate overlap rates
  if (nrow(significant_pathways) >= 2) {
    sets_list <- list()
    for (j in 1:nrow(significant_pathways)) {
      label <- significant_pathways$gradient_label[j]
      sets_list[[label]] <- unlist(significant_pathways$pathways[j])
    }
    
    cat("\n  Overlap rates between gradients:\n")
    labels_vec <- significant_pathways$gradient_label
    for (i in 1:(length(labels_vec) - 1)) {
      for (j in (i + 1):length(labels_vec)) {
        set1 <- sets_list[[labels_vec[i]]]
        set2 <- sets_list[[labels_vec[j]]]
        inter <- length(intersect(set1, set2))
        union <- length(union(set1, set2))
        rate <- ifelse(union == 0, NA, inter / union * 100)
        cat("    ", labels_vec[i], "vs", labels_vec[j], ":", round(rate, 2), "%\n")
      }
    }
  }
  
  # Select final results (highest precision)
  fgsea_res <- gsea_results_list[[length(gsea_results_list)]]
  
  fgsea_res_clean <- fgsea_res
  list_cols_final <- names(fgsea_res_clean)[sapply(fgsea_res_clean, is.list)]
  for (col in list_cols_final) {
    fgsea_res_clean[[col]] <- sapply(fgsea_res_clean[[col]], function(x) {
      if (is.null(x) || length(x) == 0) return(NA_character_)
      paste(as.character(x), collapse = ";")
    })
  }
  
  write.csv(fgsea_res_clean, 
            file.path(output_dir, paste0(labels, score_name, "_GSEA_results_100k.csv")), 
            row.names = FALSE)
  
  # Visualization
  sig <- fgsea_res %>% filter(pval < pval_threshold)
  if (nrow(sig) == 0) {
    warning("No pathways with pval < ", pval_threshold, " for ", score_name)
    sig <- fgsea_res %>% head(10)
  }
  
  top_pos <- sig %>% arrange(desc(NES)) %>% head(top_n)
  top_neg <- sig %>% arrange(NES) %>% head(top_n)
  selected <- bind_rows(top_pos, top_neg) %>% distinct(pathway, .keep_all = TRUE)
  selected_paths <- selected$pathway
  
  display_names <- sub("^C5_BP\\|GOBP_", "", selected_paths)
  display_names <- sub("^C5_BP\\|", "", display_names)
  display_names <- str_replace_all(str_to_lower(display_names), "_", " ")
  names(display_names) <- selected_paths
  
  calc_running_es <- function(stats, geneset) {
    N <- length(stats)
    selected_idx <- which(names(stats) %in% geneset)
    Nh <- length(selected_idx)
    if (Nh == 0) return(rep(0, N))
    weights <- abs(stats)^1
    hit <- integer(N)
    hit[selected_idx] <- 1
    Phit <- cumsum(hit * weights) / sum(weights[hit == 1])
    Pmiss <- cumsum((1 - hit) / (N - Nh))
    Phit - Pmiss
  }
  
  df_list <- lapply(selected_paths, function(pw) {
    es <- calc_running_es(ranked, pathways_list[[pw]])
    data.frame(
      rank = seq_along(ranked),
      ES = es,
      ESpos = es - min(es),
      pathway = pw,
      pathway_display = display_names[pw],
      NES = selected$NES[match(pw, selected$pathway)],
      neglog10pval = -log10(selected$pval[match(pw, selected$pathway)] + 1e-300),
      stringsAsFactors = FALSE
    )
  })
  df <- bind_rows(df_list)
  
  pos_order <- selected %>% arrange(desc(NES)) %>% pull(pathway)
  df$pathway_display <- factor(df$pathway_display, levels = display_names[pos_order])
  
  max_chars <- 40
  df$pathway_display_wrapped <- str_wrap(df$pathway_display, width = max_chars)
  
  p <- ggplot(df, aes(x = rank, y = pathway_display_wrapped, height = ESpos, fill = neglog10pval)) +
    geom_ridgeline(stat = "identity", scale = 1, colour = "black", size = 0.2) +
    scale_fill_gradient(
      low = "#4DBBD5", 
      high = "#E64B35", 
      name = expression(-log[10](italic(p) - value))
    ) +
    labs(
      x = "Gene rank (by correlation with model score)", 
      y = "", 
      title = bquote(.(score_name) ~ " GSEA (100k permutations)")
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10, face = "bold"),
      plot.title = element_text(hjust = 0.5),
      strip.text = element_text(size = 15, face = "bold")
    )
  
  ggsave(file.path(output_dir, paste0(labels, score_name, "_GSEA_results_100k.png")), 
         p, width = 6.5, height = 5)
  ggsave(file.path(output_dir, paste0(labels, score_name, "_GSEA_results_100k.pdf")), 
         p, width = 6.5, height = 5)
  
  cat("\n  Results saved for:", score_name, "\n")
  
  return(list(
    results = fgsea_res,
    gradient_results = all_gradient_results,
    plot = p,
    selected_pathways = selected,
    ranked_genes = ranked
  ))
}

# Run GSEA for both models
result_model1 <- run_score_based_gsea(
  score_vector = group_merged$risk_score_1,
  score_name = "Model1_Progression",
  expr_matrix = exp_df,
  pathways_list = pathways_list,
  output_dir = output,
  labels = "10_GSEA_",
  nperm_gradients = c(1000, 10000, 100000),
  gradient_labels = c("1k", "10k", "100k"),
  pval_threshold = 0.01,
  top_n = 5
)

result_model2 <- run_score_based_gsea(
  score_vector = group_merged$risk_score_2,
  score_name = "Model2_Prognostic",
  expr_matrix = exp_df,
  pathways_list = pathways_list,
  output_dir = output,
  labels = "10_GSEA_",
  nperm_gradients = c(1000, 10000, 100000),
  gradient_labels = c("1k", "10k", "100k"),
  pval_threshold = 0.01,
  top_n = 5
)

# Generate gradient stability summary
target_genes <- c("PRB2", "TUBA1C", "NPC2", "LRRC32", "TLN1", "SOX9", "SERPINA3", "RNASE4")

gradient_summary <- data.frame()

# Single gene results
for (gene in target_genes) {
  grad_file <- file.path(output, paste0("10_GSEA_", gene, "_GSEA_gradient_comparison.csv"))
  
  if (file.exists(grad_file)) {
    grad_res <- read.csv(grad_file)
    
    if (!"gradient_label" %in% colnames(grad_res)) {
      cat("Warning: gradient_label column not found in", gene, "\n")
      next
    }
    
    sig_by_gradient <- grad_res %>%
      filter(pval < 0.01) %>%
      group_by(gradient_label) %>%
      summarise(
        pathways = list(pathway),
        count = n()
      ) %>%
      ungroup()
    
    if (nrow(sig_by_gradient) >= 2) {
      sets_list <- list()
      for (j in 1:nrow(sig_by_gradient)) {
        label <- sig_by_gradient$gradient_label[j]
        sets_list[[label]] <- unlist(sig_by_gradient$pathways[j])
      }
      
      labels_vec <- sig_by_gradient$gradient_label
      
      overlap_1k_vs_10k <- NA
      overlap_1k_vs_100k <- NA
      overlap_10k_vs_100k <- NA
      
      for (i in 1:(length(labels_vec) - 1)) {
        for (j in (i + 1):length(labels_vec)) {
          set1 <- sets_list[[labels_vec[i]]]
          set2 <- sets_list[[labels_vec[j]]]
          inter <- length(intersect(set1, set2))
          union <- length(union(set1, set2))
          rate <- ifelse(union == 0, NA, inter / union * 100)
          
          pair <- paste(labels_vec[i], "vs", labels_vec[j])
          if (pair == "1k vs 10k" || pair == "10k vs 1k") {
            overlap_1k_vs_10k <- rate
          } else if (pair == "1k vs 100k" || pair == "100k vs 1k") {
            overlap_1k_vs_100k <- rate
          } else if (pair == "10k vs 100k" || pair == "100k vs 10k") {
            overlap_10k_vs_100k <- rate
          }
        }
      }
      
      gradient_summary <- bind_rows(gradient_summary, data.frame(
        Gene = gene,
        Type = "Single Gene",
        n_1k = ifelse("1k" %in% sig_by_gradient$gradient_label, 
                      sig_by_gradient$count[sig_by_gradient$gradient_label == "1k"], 0),
        n_10k = ifelse("10k" %in% sig_by_gradient$gradient_label, 
                       sig_by_gradient$count[sig_by_gradient$gradient_label == "10k"], 0),
        n_100k = ifelse("100k" %in% sig_by_gradient$gradient_label, 
                        sig_by_gradient$count[sig_by_gradient$gradient_label == "100k"], 0),
        overlap_1k_vs_10k = round(overlap_1k_vs_10k, 2),
        overlap_1k_vs_100k = round(overlap_1k_vs_100k, 2),
        overlap_10k_vs_100k = round(overlap_10k_vs_100k, 2)
      ))
    }
  }
}

# Model 1 results
grad_file <- file.path(output, "10_GSEA_Model1_Progression_GSEA_gradient_comparison.csv")
if (file.exists(grad_file)) {
  grad_res <- read.csv(grad_file)
  
  sig_by_gradient <- grad_res %>%
    filter(pval < 0.01) %>%
    group_by(gradient_label) %>%
    summarise(pathways = list(pathway), count = n()) %>%
    ungroup()
  
  if (nrow(sig_by_gradient) >= 2) {
    sets_list <- list()
    for (j in 1:nrow(sig_by_gradient)) {
      label <- sig_by_gradient$gradient_label[j]
      sets_list[[label]] <- unlist(sig_by_gradient$pathways[j])
    }
    labels_vec <- sig_by_gradient$gradient_label
    
    overlap_1k_vs_10k <- NA
    overlap_1k_vs_100k <- NA
    overlap_10k_vs_100k <- NA
    
    for (i in 1:(length(labels_vec) - 1)) {
      for (j in (i + 1):length(labels_vec)) {
        set1 <- sets_list[[labels_vec[i]]]
        set2 <- sets_list[[labels_vec[j]]]
        inter <- length(intersect(set1, set2))
        union <- length(union(set1, set2))
        rate <- ifelse(union == 0, NA, inter / union * 100)
        pair <- paste(labels_vec[i], "vs", labels_vec[j])
        if (pair == "1k vs 10k" || pair == "10k vs 1k") {
          overlap_1k_vs_10k <- rate
        } else if (pair == "1k vs 100k" || pair == "100k vs 1k") {
          overlap_1k_vs_100k <- rate
        } else if (pair == "10k vs 100k" || pair == "100k vs 10k") {
          overlap_10k_vs_100k <- rate
        }
      }
    }
    
    gradient_summary <- bind_rows(gradient_summary, data.frame(
      Gene = "Model1_Progression",
      Type = "Model Score",
      n_1k = ifelse("1k" %in% sig_by_gradient$gradient_label, 
                    sig_by_gradient$count[sig_by_gradient$gradient_label == "1k"], 0),
      n_10k = ifelse("10k" %in% sig_by_gradient$gradient_label, 
                     sig_by_gradient$count[sig_by_gradient$gradient_label == "10k"], 0),
      n_100k = ifelse("100k" %in% sig_by_gradient$gradient_label, 
                      sig_by_gradient$count[sig_by_gradient$gradient_label == "100k"], 0),
      overlap_1k_vs_10k = round(overlap_1k_vs_10k, 2),
      overlap_1k_vs_100k = round(overlap_1k_vs_100k, 2),
      overlap_10k_vs_100k = round(overlap_10k_vs_100k, 2)
    ))
  }
}

# Model 2 results
grad_file <- file.path(output, "10_GSEA_Model2_Prognostic_GSEA_gradient_comparison.csv")
if (file.exists(grad_file)) {
  grad_res <- read.csv(grad_file)
  
  sig_by_gradient <- grad_res %>%
    filter(pval < 0.01) %>%
    group_by(gradient_label) %>%
    summarise(pathways = list(pathway), count = n()) %>%
    ungroup()
  
  if (nrow(sig_by_gradient) >= 2) {
    sets_list <- list()
    for (j in 1:nrow(sig_by_gradient)) {
      label <- sig_by_gradient$gradient_label[j]
      sets_list[[label]] <- unlist(sig_by_gradient$pathways[j])
    }
    labels_vec <- sig_by_gradient$gradient_label
    
    overlap_1k_vs_10k <- NA
    overlap_1k_vs_100k <- NA
    overlap_10k_vs_100k <- NA
    
    for (i in 1:(length(labels_vec) - 1)) {
      for (j in (i + 1):length(labels_vec)) {
        set1 <- sets_list[[labels_vec[i]]]
        set2 <- sets_list[[labels_vec[j]]]
        inter <- length(intersect(set1, set2))
        union <- length(union(set1, set2))
        rate <- ifelse(union == 0, NA, inter / union * 100)
        pair <- paste(labels_vec[i], "vs", labels_vec[j])
        if (pair == "1k vs 10k" || pair == "10k vs 1k") {
          overlap_1k_vs_10k <- rate
        } else if (pair == "1k vs 100k" || pair == "100k vs 1k") {
          overlap_1k_vs_100k <- rate
        } else if (pair == "10k vs 100k" || pair == "100k vs 10k") {
          overlap_10k_vs_100k <- rate
        }
      }
    }
    
    gradient_summary <- bind_rows(gradient_summary, data.frame(
      Gene = "Model2_Prognostic",
      Type = "Model Score",
      n_1k = ifelse("1k" %in% sig_by_gradient$gradient_label, 
                    sig_by_gradient$count[sig_by_gradient$gradient_label == "1k"], 0),
      n_10k = ifelse("10k" %in% sig_by_gradient$gradient_label, 
                     sig_by_gradient$count[sig_by_gradient$gradient_label == "10k"], 0),
      n_100k = ifelse("100k" %in% sig_by_gradient$gradient_label, 
                      sig_by_gradient$count[sig_by_gradient$gradient_label == "100k"], 0),
      overlap_1k_vs_10k = round(overlap_1k_vs_10k, 2),
      overlap_1k_vs_100k = round(overlap_1k_vs_100k, 2),
      overlap_10k_vs_100k = round(overlap_10k_vs_100k, 2)
    ))
  }
}

# Save summary table
gene_order <- c("PRB2", "TUBA1C", "NPC2", "LRRC32", "TLN1", "SOX9", "SERPINA3", "RNASE4", 
                "Model1_Progression", "Model2_Prognostic")
gradient_summary <- gradient_summary %>%
  mutate(Gene = factor(Gene, levels = gene_order)) %>%
  arrange(Gene, Type)

write.csv(gradient_summary, 
          file.path(output, "10_GSEA_GSEA_gradient_stability_summary_with_models.csv"), 
          row.names = FALSE)

cat("\n=== GSEA Gradient Stability Summary (with Model Scores) ===\n")
print(gradient_summary)

# Calculate average overlap rates
avg_overall <- gradient_summary %>%
  summarise(
    Avg_1k_vs_10k = round(mean(overlap_1k_vs_10k, na.rm = TRUE), 2),
    Avg_1k_vs_100k = round(mean(overlap_1k_vs_100k, na.rm = TRUE), 2),
    Avg_10k_vs_100k = round(mean(overlap_10k_vs_100k, na.rm = TRUE), 2)
  )

cat("\n=== Average Overlap Rates Across All Genes/Models ===\n")
print(avg_overall)

# Create gradient stability plot
plot_data <- gradient_summary %>%
  select(Gene, Type, overlap_1k_vs_10k, overlap_1k_vs_100k, overlap_10k_vs_100k) %>%
  pivot_longer(
    cols = starts_with("overlap"),
    names_to = "Comparison",
    values_to = "Overlap_Rate"
  ) %>%
  mutate(
    Comparison = case_when(
      Comparison == "overlap_1k_vs_10k" ~ "1k vs 10k",
      Comparison == "overlap_1k_vs_100k" ~ "1k vs 100k",
      Comparison == "overlap_10k_vs_100k" ~ "10k vs 100k"
    ),
    Comparison = factor(Comparison, levels = c("1k vs 10k", "1k vs 100k", "10k vs 100k"))
  )

plot_data$Gene_Type <- ifelse(plot_data$Gene %in% c("Model1_Progression", "Model2_Prognostic"), 
                              "Model Score", "Single Gene")

# Determine which genes should be bold (model scores)
is_bold <- plot_data$Gene %in% c("Model1_Progression", "Model2_Prognostic")

p_stability <- ggplot(plot_data, aes(x = Gene, y = Overlap_Rate, fill = Comparison)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "red", alpha = 0.6) +
  geom_hline(yintercept = 50, linetype = "dotted", color = "orange", alpha = 0.6) +
  scale_fill_manual(
    values = c("1k vs 10k" = "#4DBBD5", "1k vs 100k" = "#F39B7F", "10k vs 100k" = "#E64B35"),
    name = "Permutation Comparison"
  ) +
  labs(
    x = "Gene / Model Score", 
    y = "Overlap Rate of Significant Pathways (%)",
    title = "GSEA Gradient Stability Across Permutation Counts",
    subtitle = "Single genes vs. model prediction scores"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      face = ifelse(unique(plot_data$Gene) %in% c("Model1_Progression", "Model2_Prognostic"), 
                    "bold", "italic"), 
      size = 9, angle = 45, hjust = 1
    ),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "bottom",
    legend.text = element_text(size = 9)
  )

ggsave(file.path(output, "10_GSEA_GSEA_gradient_stability_with_models.png"), 
       p_stability, width = 12, height = 6)
ggsave(file.path(output, "10_GSEA_GSEA_gradient_stability_with_models.pdf"), 
       p_stability, width = 12, height = 6)

cat("\n=== All Analyses Complete ===\n")
cat("Results saved to:", output, "\n")
