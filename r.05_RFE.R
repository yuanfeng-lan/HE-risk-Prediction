rm(list = ls()); gc()
ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "05_RFE")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(ORIGINAL_DIR)

library(tidyverse)
library(caret)
library(ggplot2)
library(ggplotify)
library(gridExtra)
library(grid)
library(ggtext)

# Helper function for RFE analysis
run_rfe_analysis <- function(exp_data, gene_list, response_data, 
                             response_col, output_dir, prefix,
                             seed = 123) {
  
  # Filter expression data
  exp_filtered <- exp_data[gene_list, ]
  exp_data_t <- as.data.frame(t(exp_filtered))
  
  # Prepare data
  response_var <- response_data[, response_col, drop = FALSE]
  data <- cbind(exp_data_t, response_var)
  dat <- data
  num_features <- ncol(dat) - 1
  
  # Run RFE
  set.seed(seed)
  control <- rfeControl(functions = lmFuncs, method = "repeatedcv", 
                        number = 5, repeats = 3)
  results <- rfe(dat[, 1:num_features],
                 dat[, response_col],
                 sizes = c(1:num_features),
                 rfeControl = control,
                 method = "RMSE")
  
  # Save selected predictors
  predictors_selected <- predictors(results)
  write.csv(data.frame(genes = predictors_selected), 
            file.path(output_dir, paste0(prefix, "results.csv")),
            row.names = FALSE)
  
  # Create RFE performance plot
  p <- create_rfe_plot(results, prefix)
  
  # Save plot
  ggsave(file.path(output_dir, paste0(prefix, "RFE_results.png")), 
         p, width = 6, height = 6)
  ggsave(file.path(output_dir, paste0(prefix, "RFE_results.pdf")), 
         p, width = 6, height = 6)
  
  return(list(results = results, predictors = predictors_selected))
}

# Helper function for RFE with gene label plot
run_rfe_analysis_with_labels <- function(exp_data, gene_list, response_data, 
                                         response_col, output_dir, prefix,
                                         seed = 123) {
  
  # Filter expression data
  exp_filtered <- exp_data[gene_list, ]
  exp_data_t <- as.data.frame(t(exp_filtered))
  
  # Prepare data
  response_var <- response_data[, response_col, drop = FALSE]
  data <- cbind(exp_data_t, response_var)
  dat <- data
  num_features <- ncol(dat) - 1
  
  # Run RFE
  set.seed(seed)
  control <- rfeControl(functions = lmFuncs, method = "repeatedcv", 
                        number = 5, repeats = 3)
  results <- rfe(dat[, 1:num_features],
                 dat[, response_col],
                 sizes = c(1:num_features),
                 rfeControl = control,
                 method = "RMSE")
  
  # Save selected predictors
  predictors_selected <- predictors(results)
  write.csv(data.frame(genes = predictors_selected), 
            file.path(output_dir, paste0(prefix, "results.csv")),
            row.names = FALSE)
  
  # Create RFE performance plot
  p1 <- create_rfe_plot(results, prefix)
  
  # Create gene label plot
  p2 <- create_gene_label_plot(gene_list, predictors_selected)
  
  # Combine plots
  combined_plot <- grid.arrange(
    p1, p2,
    ncol = 2,
    widths = c(4, 1),
    layout_matrix = rbind(c(1, 2))
  )
  
  # Save combined plot
  ggsave(file.path(output_dir, paste0(prefix, "RFE_results_with_labels.png")), 
         combined_plot, width = 6, height = 6)
  ggsave(file.path(output_dir, paste0(prefix, "RFE_results_with_labels.pdf")), 
         combined_plot, width = 6, height = 6)
  
  # Save individual plot
  ggsave(file.path(output_dir, paste0(prefix, "RFE_results.png")), 
         p1, width = 6, height = 6)
  ggsave(file.path(output_dir, paste0(prefix, "RFE_results.pdf")), 
         p1, width = 6, height = 6)
  
  return(list(results = results, predictors = predictors_selected))
}

# Helper function to create RFE plot
create_rfe_plot <- function(results, prefix) {
  df <- results$results
  df$Variables <- as.numeric(as.character(df$Variables))
  best_n <- as.numeric(results[["bestSubset"]])
  best_row <- which(df$Variables == best_n)[1]
  best_rmse <- df$RMSE[best_row]
  
  p <- ggplot(df, aes(x = Variables, y = RMSE)) +
    geom_line(color = "#2c7fb8", size = 1) +
    geom_point(color = "#2c7fb8", size = 3) +
    { if ("RMSESD" %in% names(df)) 
      geom_errorbar(aes(ymin = RMSE - RMSESD, ymax = RMSE + RMSESD), 
                    width = 0.2, color = "#2c7fb8", alpha = 0.6)
      else NULL } +
    geom_point(data = df[best_row, , drop = FALSE], 
               aes(x = Variables, y = RMSE), color = "#e31a1c", size = 4) +
    geom_vline(xintercept = best_n, linetype = "dashed", colour = "grey40") +
    annotate("text", x = best_n, y = best_rmse, 
             label = paste0("Best: ", best_n, "\nRMSE=", round(best_rmse, 3)),
             vjust = -1.2, hjust = 0.5, color = "#e31a1c", size = 6) +
    scale_x_continuous(breaks = df$Variables) +
    labs(x = "Number of Features", y = "RMSE (Repeated Cross-Validation)",
         title = "RFE Performance", subtitle = "", caption = "") +
    theme_minimal(base_size = 15) +
    theme(
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_line(color = "grey95"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
      axis.title = element_text(face = "bold", size = 15),
      axis.text = element_text(color = "black", size = 15),
      panel.border = element_rect(colour = "black", fill = NA, size = 1)
    )
  
  return(p)
}

# Helper function to create gene label plot
create_gene_label_plot <- function(gene_list, selected_genes) {
  label_df <- data.frame(
    x = rep(1, length(gene_list)),
    y = seq_along(gene_list),
    label = sapply(gene_list, function(gene) {
      if (gene %in% selected_genes) {
        paste0("<span style='color:red;'><i>", gene, "</i></span>")
      } else {
        paste0("<i>", gene, "</i>")
      }
    })
  )
  
  p <- ggplot(label_df, aes(x = x, y = y, label = label)) +
    geom_richtext(aes(label = label), hjust = 0, size = 4, 
                  fill = NA, label.color = NA) +
    scale_y_reverse(breaks = seq_along(gene_list), labels = NULL) +
    theme_void() +
    theme(
      plot.margin = margin(0, 0, 0, -2, "cm"),
      panel.grid = element_blank()
    )
  
  return(p)
}

#### GSE139602 RFE Analysis ####
# Load data
gse139602_up <- read.csv(file.path("02_Mfuzz", "HE_GSE139602_UP_gene.csv"), row.names = 1)
gse139602_down <- read.csv(file.path("02_Mfuzz", "HE_GSE139602_DOWN_gene.csv"), row.names = 1)
top_genes_gse139602 <- c(gse139602_up$x, gse139602_down$x)

exp_gse139602 <- read.csv(file.path("00_rawdata", "00.rawdata_GSE139602_exp.csv"), row.names = 1)
group_gse139602 <- read.csv(file.path("00_rawdata", "00.rawdata_GSE139602_group.csv"), row.names = 1)

# Match samples
group_gse139602 <- group_gse139602[colnames(exp_gse139602), , drop = FALSE]

# Convert group labels to numeric
group_gse139602$group <- ifelse(group_gse139602$characteristics_ch1 == "disease state: Healthy", 0,
                         ifelse(group_gse139602$characteristics_ch1 == "disease state: eCLD", 1,
                         ifelse(group_gse139602$characteristics_ch1 == "disease state: Compensated Cirrhosis", 2,
                         ifelse(group_gse139602$characteristics_ch1 == "disease state: Decompesated Cirrhosis", 3,
                         ifelse(group_gse139602$characteristics_ch1 == "disease state: Acute-on-chronic liver failure", 4, NA)))))

# Run RFE with labels
result_gse139602 <- run_rfe_analysis_with_labels(
  exp_data = exp_gse139602,
  gene_list = top_genes_gse139602,
  response_data = group_gse139602,
  response_col = "group",
  output_dir = output,
  prefix = "GSE139602_",
  seed = 47
)

#### GSE15654 RFE Analysis ####
# Load data
gse15654_genes <- read.csv(file.path("02_Mfuzz", "HE_GSE15654_ALL_gene_DEG.csv"))
top_genes_gse15654 <- gse15654_genes$x

exp_gse15654 <- read.csv(file.path("00_rawdata", "00.rawdata_GSE15654_exp.csv"), row.names = 1)
surv_gse15654 <- read.csv(file.path("00_rawdata", "00.rawdata_GSE15654_sur.csv"), row.names = 1)

# Match samples
surv_gse15654 <- surv_gse15654[colnames(exp_gse15654), , drop = FALSE]

# Run RFE
result_gse15654 <- run_rfe_analysis(
  exp_data = exp_gse15654,
  gene_list = top_genes_gse15654,
  response_data = surv_gse15654,
  response_col = "death",
  output_dir = output,
  prefix = "GSE15654_",
  seed = 18
)

message("RFE analysis completed successfully!")

# Print summary
cat("\n=== Summary ===\n")
cat("GSE139602 selected genes (", length(result_gse139602$predictors), "):\n")
print(result_gse139602$predictors)
cat("\nGSE15654 selected genes (", length(result_gse15654$predictors), "):\n")
print(result_gse15654$predictors)

# Save combined results
combined_genes <- list(
  GSE139602 = result_gse139602$predictors,
  GSE15654 = result_gse15654$predictors
)

write.csv(data.frame(
  Dataset = rep(names(combined_genes), sapply(combined_genes, length)),
  Gene = unlist(combined_genes)
), file.path(output, "RFE_selected_genes_summary.csv"), row.names = FALSE)

# Find overlapping genes
overlap_genes <- intersect(result_gse139602$predictors, result_gse15654$predictors)
if (length(overlap_genes) > 0) {
  write.csv(data.frame(genes = overlap_genes), 
            file.path(output, "RFE_overlap_genes.csv"), 
            row.names = FALSE)
  cat("\nOverlapping genes between datasets (", length(overlap_genes), "):\n")
  print(overlap_genes)
} else {
  cat("\nNo overlapping genes found between datasets.\n")
}
