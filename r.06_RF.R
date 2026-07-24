rm(list = ls()); gc()
ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "06_RF")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(ORIGINAL_DIR)

library(randomForest)
library(randomForestSRC)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggvenn)

# Helper function for Random Forest analysis (regression)
run_rf_regression <- function(exp_data, gene_list, response_data, 
                              response_col, output_dir, prefix,
                              seed = 123, top_n = 10) {
  
  # Filter expression data
  exp_filtered <- exp_data[gene_list, ]
  exp_data_t <- as.data.frame(t(exp_filtered))
  
  # Prepare data
  response_var <- response_data[, response_col, drop = FALSE]
  data <- cbind(exp_data_t, response_var)
  
  # Run Random Forest
  set.seed(seed)
  forest_result <- randomForest(as.formula(paste(response_col, "~ .")), 
                                data = data, importance = TRUE)
  
  # Extract importance
  imp <- as.data.frame(importance(forest_result))
  imp$feature <- rownames(imp)
  
  # Get column names
  mse_col <- "%IncMSE"
  node_col <- "IncNodePurity"
  
  # Calculate standardized importance
  imp2 <- imp %>%
    mutate(mse_z = as.numeric(scale(.data[[mse_col]])),
           node_z = as.numeric(scale(.data[[node_col]])),
           mean_z = (mse_z + node_z) / 2) %>%
    arrange(mean_z)
  
  imp2$feature <- factor(imp2$feature, levels = imp2$feature)
  
  # Select top features
  top_feats <- imp2 %>% slice_max(mean_z, n = top_n) %>% pull(feature)
  imp2 <- imp2 %>% mutate(top10 = feature %in% top_feats)
  
  # Create importance plot
  p <- create_importance_plot(imp2, top_n, prefix)
  
  # Save results
  ggsave(file.path(output_dir, paste0(prefix, "RF_results.png")), 
         p, width = 6, height = 6)
  ggsave(file.path(output_dir, paste0(prefix, "RF_results.pdf")), 
         p, width = 6, height = 6)
  
  write.csv(imp2, file.path(output_dir, paste0(prefix, "RF_results_table.csv")), 
            row.names = FALSE)
  
  # Save selected genes
  selected_genes <- imp2 %>% filter(top10 == TRUE) %>% pull(feature)
  write.csv(data.frame(genes = selected_genes), 
            file.path(output_dir, paste0(prefix, "RF_results.csv")), 
            row.names = FALSE)
  
  return(list(model = forest_result, importance = imp2, selected_genes = selected_genes))
}

# Helper function for Random Forest survival analysis
run_rf_survival <- function(exp_data, gene_list, surv_data, 
                            time_col, event_col, output_dir, prefix,
                            seed = 123, top_n = 10) {
  
  # Filter expression data
  exp_filtered <- exp_data[gene_list, ]
  exp_data_t <- as.data.frame(t(exp_filtered))
  
  # Prepare data
  surv_subset <- surv_data[, c(time_col, event_col), drop = FALSE]
  data <- cbind(exp_data_t, surv_subset)
  
  # Create formula
  form <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ ."))
  
  # Determine mtry
  mf <- model.frame(form, data = data)
  n_xvar <- ncol(mf) - 1
  mtry_val <- floor(sqrt(n_xvar))
  if (mtry_val < 1) mtry_val <- 1
  
  # Run Random Forest survival
  set.seed(seed)
  rf_surv <- randomForestSRC::rfsrc(
    formula = form,
    data = data,
    ntree = 1000,
    mtry = mtry_val,
    importance = TRUE,
    block.size = 1
  )
  
  # Extract importance
  var_importance <- rf_surv$importance
  imp <- as.data.frame(var_importance)
  
  if (ncol(imp) == 1) {
    imp$importance <- as.numeric(imp[[1]])
  } else {
    imp$importance <- rowMeans(imp, na.rm = TRUE)
  }
  imp$feature <- rownames(imp)
  
  # Calculate standardized importance
  imp2 <- imp %>%
    mutate(importance_z = as.numeric(scale(importance))) %>%
    arrange(importance_z)
  
  imp2$feature <- factor(imp2$feature, levels = imp2$feature)
  
  # Select top features
  top_feats <- imp2 %>% slice_max(importance_z, n = top_n) %>% pull(feature)
  imp2 <- imp2 %>% mutate(top10 = feature %in% top_feats)
  
  # Create importance plot
  p <- create_importance_plot_survival(imp2, top_n, prefix)
  
  # Save results
  ggsave(file.path(output_dir, paste0(prefix, "RF_results.png")), 
         p, width = 4, height = 6)
  ggsave(file.path(output_dir, paste0(prefix, "RF_results.pdf")), 
         p, width = 4, height = 6)
  
  write.csv(imp2, file.path(output_dir, paste0(prefix, "RF_results_table.csv")), 
            row.names = FALSE)
  
  # Save selected genes
  selected_genes <- imp2 %>% filter(top10 == TRUE) %>% pull(feature)
  write.csv(data.frame(genes = selected_genes), 
            file.path(output_dir, paste0(prefix, "RF_results.csv")), 
            row.names = FALSE)
  
  return(list(model = rf_surv, importance = imp2, selected_genes = selected_genes))
}

# Helper function to create importance plot
create_importance_plot <- function(imp2, top_n, prefix) {
  label_x <- -3
  rng <- range(imp2$mean_z, na.rm = TRUE)
  
  imp_long <- imp2[, c("feature", "mse_z", "node_z", "mean_z", "top10")] %>%
    pivot_longer(cols = c(mse_z, node_z), names_to = "metric", values_to = "zvalue") %>%
    mutate(
      metric = ifelse(metric == "mse_z", "%IncMSE (z)", "IncNodePurity (z)"),
      feature = factor(feature, levels = levels(imp2$feature))
    )
  
  p <- ggplot() +
    geom_col(data = imp2, aes(x = mean_z, y = feature, fill = top10), 
             width = 0.6, alpha = 0.6) +
    geom_line(data = imp_long, aes(x = zvalue, y = feature, group = feature),
              color = "#7F8C8D", linetype = "dashed", size = 0.3) +
    geom_point(data = imp_long, aes(x = zvalue, y = feature, color = metric, shape = metric),
               size = 3) +
    geom_label(data = imp2,
               aes(x = label_x, y = feature, 
                   label = paste0("italic('", feature, "')")),
               fill = "white",
               color = ifelse(imp2$top10, "#E64B35", "#7F8C8D"),
               hjust = 0,
               size = 4.5,
               parse = TRUE,
               linewidth = 0.15,
               label.r = grid::unit(0.12, "lines")) +
    scale_fill_manual(values = c("FALSE" = "#7F8C8D", "TRUE" = "#E64B35"), guide = "none") +
    scale_color_manual(values = c("%IncMSE (z)" = "#3C5488", "IncNodePurity (z)" = "#00A087")) +
    labs(x = "Standardized importance (z)", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title = element_blank(),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 12)
    ) +
    coord_cartesian(xlim = c(label_x, rng[2] + 0.02 * (rng[2] - rng[1])))
  
  return(p)
}

# Helper function to create survival importance plot
create_importance_plot_survival <- function(imp2, top_n, prefix) {
  label_x <- -1.5
  rng <- range(imp2$importance_z, na.rm = TRUE)
  
  p <- ggplot() +
    geom_col(data = imp2, aes(x = importance_z, y = feature, fill = top10), 
             width = 0.6, alpha = 0.6) +
    geom_point(data = imp2, aes(x = importance_z, y = feature, color = top10), 
               size = 2) +
    geom_label(data = imp2,
               aes(x = label_x, y = feature, 
                   label = paste0("italic('", feature, "')")),
               fill = "white",
               color = ifelse(imp2$top10, "#E64B35", "#7F8C8D"),
               hjust = 0,
               size = 4.5,
               parse = TRUE,
               linewidth = 0.15,
               label.r = grid::unit(0.12, "lines")) +
    scale_fill_manual(values = c("FALSE" = "#7F8C8D", "TRUE" = "#E64B35"), guide = "none") +
    scale_color_manual(values = c("FALSE" = "#7F8C8D", "TRUE" = "#E64B35"), guide = "none") +
    labs(x = "Standardized importance (z)", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title = element_blank(),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 12)
    ) +
    coord_cartesian(xlim = c(label_x, rng[2] + 0.02 * (rng[2] - rng[1])))
  
  return(p)
}

# Helper function for Venn diagram
create_venn_diagram <- function(venn_list, output_dir, prefix, 
                                colors = c("#00A087", "#3C5488", "#F39B7F")) {
  p <- ggvenn(
    venn_list,
    names(venn_list),
    text_size = 8,
    fill_color = colors,
    fill_alpha = 0.7,
    stroke_color = "black",
    stroke_size = 1,
    show_percentage = FALSE
  )
  
  ggsave(file.path(output_dir, paste0(prefix, "venn_final_cogenes.pdf")), 
         p, width = 4, height = 4)
  ggsave(file.path(output_dir, paste0(prefix, "venn_final_cogenes.png")), 
         p, width = 4, height = 4)
  
  return(p)
}

#### GSE139602 Random Forest Analysis ####
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

# Run Random Forest
result_gse139602 <- run_rf_regression(
  exp_data = exp_gse139602,
  gene_list = top_genes_gse139602,
  response_data = group_gse139602,
  response_col = "group",
  output_dir = output,
  prefix = "GSE139602_",
  seed = 5,
  top_n = 10
)

# Load results from other methods
lasso_gse139602 <- read.csv(file.path("04_LASSO", "GSE139602_result.csv"))
rfe_gse139602 <- read.csv(file.path("05_RFE", "GSE139602_results.csv"))
rf_gse139602 <- read.csv(file.path("06_RF", "GSE139602_RF_results.csv"))

# Create Venn diagram
venn_gse139602 <- list(
  LASSO = lasso_gse139602$x,
  RFE = rfe_gse139602$x,
  RandomForest = rf_gse139602$x
)

create_venn_diagram(venn_gse139602, output, "GSE139602_")

# Get common genes
common_genes_gse139602 <- Reduce(intersect, venn_gse139602)
write.csv(data.frame(genes = common_genes_gse139602), 
          file.path(output, "GSE139602_final_cogenes.csv"), 
          row.names = FALSE)

#### GSE15654 Random Forest Survival Analysis ####
# Load data
gse15654_genes <- read.csv(file.path("02_Mfuzz", "HE_GSE15654_ALL_gene_DEG.csv"))
top_genes_gse15654 <- gse15654_genes$x

exp_gse15654 <- read.csv(file.path("00_rawdata", "00.rawdata_GSE15654_exp.csv"), row.names = 1)
surv_gse15654 <- read.csv(file.path("00_rawdata", "00.rawdata_GSE15654_sur.csv"), row.names = 1)

# Match samples
surv_gse15654 <- surv_gse15654[colnames(exp_gse15654), , drop = FALSE]

# Run Random Forest survival
result_gse15654 <- run_rf_survival(
  exp_data = exp_gse15654,
  gene_list = top_genes_gse15654,
  surv_data = surv_gse15654,
  time_col = "days_to_death",
  event_col = "death",
  output_dir = output,
  prefix = "GSE15654_",
  seed = 2,
  top_n = 10
)

# Load results from other methods
lasso_gse15654 <- read.csv(file.path("04_LASSO", "GSE15654_result.csv"))
rfe_gse15654 <- read.csv(file.path("05_RFE", "GSE15654_results.csv"))
rf_gse15654 <- read.csv(file.path("06_RF", "GSE15654_RF_results.csv"))

# Create Venn diagram
venn_gse15654 <- list(
  LASSO = lasso_gse15654$x,
  RFE = rfe_gse15654$x,
  RandomForest = rf_gse15654$x
)

create_venn_diagram(venn_gse15654, output, "GSE15654_")

# Get common genes
common_genes_gse15654 <- Reduce(intersect, venn_gse15654)
write.csv(data.frame(genes = common_genes_gse15654), 
          file.path(output, "GSE15654_final_cogenes.csv"), 
          row.names = FALSE)

message("Random Forest analysis completed successfully!")

# Print summary
cat("\n=== Summary ===\n")
cat("GSE139602 selected genes (", length(result_gse139602$selected_genes), "):\n")
print(result_gse139602$selected_genes)
cat("\nGSE139602 common genes across all methods (", length(common_genes_gse139602), "):\n")
print(common_genes_gse139602)

cat("\nGSE15654 selected genes (", length(result_gse15654$selected_genes), "):\n")
print(result_gse15654$selected_genes)
cat("\nGSE15654 common genes across all methods (", length(common_genes_gse15654), "):\n")
print(common_genes_gse15654)

# Save combined summary
combined_summary <- data.frame(
  Dataset = c(rep("GSE139602", length(common_genes_gse139602)),
              rep("GSE15654", length(common_genes_gse15654))),
  Gene = c(common_genes_gse139602, common_genes_gse15654)
)
write.csv(combined_summary, 
          file.path(output, "final_common_genes_summary.csv"), 
          row.names = FALSE)
