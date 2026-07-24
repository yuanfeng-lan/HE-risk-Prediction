rm(list = ls()); gc()
ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "11_PCA")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(ORIGINAL_DIR)

library(scatterplot3d)
library(ggplot2)
library(dendextend)
library(tidyr)
library(dplyr)
library(ggdendro)

# Helper function for PCA and clustering analysis
perform_pca_clustering <- function(exp_file, group_file, group_subset = NULL,
                                   GSE_name, output_dir, labels = "",
                                   group_colors = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F")) {
  
  # Create output directory
  outdir <- file.path(output_dir, GSE_name)
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  # Load data
  exp <- read.csv(exp_file, header = TRUE, row.names = 1)
  group <- read.csv(group_file, header = TRUE, row.names = 1)
  
  # Subset if specified
  if (!is.null(group_subset)) {
    group <- subset(group, group != group_subset)
  }
  
  # Match samples
  exp <- exp[, rownames(group)]
  exp <- na.omit(exp)
  
  # Scale expression
  exp <- scale(exp)
  exp_data <- as.data.frame(t(exp))
  
  # Define group type based on available columns
  if ("group" %in% colnames(group)) {
    if (GSE_name == "GSE139602") {
      group$type <- ifelse(group$characteristics_ch1 == "disease state: Healthy", "Healthy",
                    ifelse(group$characteristics_ch1 == "disease state: eCLD", "eCLD",
                    ifelse(group$characteristics_ch1 == "disease state: Compensated Cirrhosis", "CC",
                    ifelse(group$characteristics_ch1 == "disease state: Decompesated Cirrhosis", "DC",
                    ifelse(group$characteristics_ch1 == "disease state: Acute-on-chronic liver failure", "ACLF", NA)))))
      group$type <- factor(group$type, levels = c("Healthy", "eCLD", "CC", "DC", "ACLF"))
    } else if (GSE_name == "GSE15654") {
      group$type <- factor(group$group, levels = c("Good_prognosis", "Poor_prognosis"))
    } else if (GSE_name == "GSE41919" || GSE_name == "GSE57193") {
      group$type <- factor(group$group, levels = c("cirrhosis with HE", "cirrhosis without HE"))
    } else {
      group$type <- factor(group$group)
    }
  } else {
    group$type <- factor(group[, 1])
  }
  
  # Perform PCA
  pca <- prcomp(exp_data, center = FALSE, scale. = FALSE)
  scores <- as.data.frame(pca$x[, 1:3])
  colnames(scores) <- c("PC1", "PC2", "PC3")
  explained <- round(100 * (pca$sdev^2) / sum(pca$sdev^2), 1)
  
  # Prepare colors
  groups <- group$type
  cols <- group_colors[as.integer(groups)]
  cols <- as.character(cols)
  
  # Create 3D PCA plot
  png(file.path(outdir, paste0(labels, GSE_name, "_pca3d_plot.png")), 
      width = 5, height = 5, units = "in", res = 300)
  s3d <- scatterplot3d(scores$PC1, scores$PC2, scores$PC3,
                       color = cols, pch = 19, cex.symbols = 1.2,
                       xlab = paste0("PC1 (", explained[1], "%)"),
                       ylab = paste0("PC2 (", explained[2], "%)"),
                       zlab = paste0("PC3 (", explained[3], "%)"),
                       main = paste0(GSE_name, " PCA 3D"))
  op <- par(xpd = NA)
  
  # Adjust legend position based on GSE
  if (GSE_name == "GSE139602") {
    legend("topright", legend = levels(groups), col = group_colors[1:length(levels(groups))],
           pch = 19, bty = "n", inset = c(-0.12, 0))
  } else if (GSE_name == "GSE15654") {
    legend("topright", legend = levels(groups), col = group_colors[1:length(levels(groups))],
           pch = 19, bty = "n", inset = c(0, -0.12))
  } else {
    legend("topright", legend = levels(groups), col = group_colors[1:length(levels(groups))],
           pch = 19, bty = "n", inset = c(-0.12, -0.14))
  }
  par(op)
  dev.off()
  
  pdf(file.path(outdir, paste0(labels, GSE_name, "_pca3d_plot.pdf")), width = 5, height = 5)
  s3d <- scatterplot3d(scores$PC1, scores$PC2, scores$PC3,
                       color = cols, pch = 19, cex.symbols = 1.2,
                       xlab = paste0("PC1 (", explained[1], "%)"),
                       ylab = paste0("PC2 (", explained[2], "%)"),
                       zlab = paste0("PC3 (", explained[3], "%)"),
                       main = paste0(GSE_name, " PCA 3D"))
  op <- par(xpd = NA)
  if (GSE_name == "GSE139602") {
    legend("topright", legend = levels(groups), col = group_colors[1:length(levels(groups))],
           pch = 19, bty = "n", inset = c(-0.12, 0))
  } else if (GSE_name == "GSE15654") {
    legend("topright", legend = levels(groups), col = group_colors[1:length(levels(groups))],
           pch = 19, bty = "n", inset = c(0, -0.12))
  } else {
    legend("topright", legend = levels(groups), col = group_colors[1:length(levels(groups))],
           pch = 19, bty = "n", inset = c(-0.12, -0.14))
  }
  par(op)
  dev.off()
  
  # Create hierarchical clustering dendrogram
  hc <- hclust(dist(exp_data))
  known_groups <- group$type
  names(known_groups) <- rownames(group)
  dend <- as.dendrogram(hc)
  
  sample_colors <- group_colors[known_groups[order.dendrogram(dend)]]
  dend_colored <- dend %>%
    set("labels_col", sample_colors) %>%
    set("labels_cex", ifelse(GSE_name == "GSE15654", 0.2, 0.8)) %>%
    set("branches_lwd", ifelse(GSE_name == "GSE15654", 1, 2))
  
  png(file.path(outdir, paste0(labels, GSE_name, "_hclust_plot.png")), 
      width = 6, height = 5, units = "in", res = 300)
  plot(dend_colored, main = paste0(GSE_name, " Cluster Dendrogram"))
  unique_groups <- unique(known_groups)
  
  # Adjust legend position
  if (GSE_name == "GSE139602" || GSE_name == "GSE15654") {
    legend("topright", legend = unique_groups, fill = group_colors[1:length(unique_groups)], 
           title = "Groups", cex = 0.8)
  } else if (GSE_name == "GSE41919") {
    legend("topleft", legend = unique_groups, fill = group_colors[1:length(unique_groups)], 
           title = "Groups", cex = 0.8)
  } else {
    legend(x = 5.8, y = 42, legend = unique_groups, fill = group_colors[1:length(unique_groups)], 
           title = "Groups", cex = 0.8)
  }
  dev.off()
  
  pdf(file.path(outdir, paste0(labels, GSE_name, "_hclust_plot.pdf")), width = 6, height = 5)
  plot(dend_colored, main = paste0(GSE_name, " Cluster Dendrogram"))
  unique_groups <- unique(known_groups)
  if (GSE_name == "GSE139602" || GSE_name == "GSE15654") {
    legend("topright", legend = unique_groups, fill = group_colors[1:length(unique_groups)], 
           title = "Groups", cex = 0.8)
  } else if (GSE_name == "GSE41919") {
    legend("topleft", legend = unique_groups, fill = group_colors[1:length(unique_groups)], 
           title = "Groups", cex = 0.8)
  } else {
    legend(x = 5.8, y = 42, legend = unique_groups, fill = group_colors[1:length(unique_groups)], 
           title = "Groups", cex = 0.8)
  }
  dev.off()
  
  return(list(pca = pca, scores = scores, dendrogram = dend_colored, group = group))
}

#### Run PCA and clustering for all datasets ####

# GSE139602
result_gse139602 <- perform_pca_clustering(
  exp_file = "./00_rawdata/00.rawdata_GSE139602_exp.csv",
  group_file = "./00_rawdata/00.rawdata_GSE139602_group.csv",
  group_subset = NULL,
  GSE_name = "GSE139602",
  output_dir = output,
  labels = "11_PCA_"
)

# GSE15654
result_gse15654 <- perform_pca_clustering(
  exp_file = "./00_rawdata/00.rawdata_GSE15654_exp.csv",
  group_file = "./00_rawdata/00.rawdata_GSE15654_group.csv",
  group_subset = NULL,
  GSE_name = "GSE15654",
  output_dir = output,
  labels = "11_PCA_"
)

# GSE41919
result_gse41919 <- perform_pca_clustering(
  exp_file = "./00_rawdata/00.rawdata_GSE41919_exp.csv",
  group_file = "./00_rawdata/00.rawdata_GSE41919_group.csv",
  group_subset = "non-cirrhotic control",
  GSE_name = "GSE41919",
  output_dir = output,
  labels = "11_PCA_"
)

# GSE57193
result_gse57193 <- perform_pca_clustering(
  exp_file = "./00_rawdata/00.rawdata_GSE57193_exp.csv",
  group_file = "./00_rawdata/00.rawdata_GSE57193_group.csv",
  group_subset = "healthy",
  GSE_name = "GSE57193",
  output_dir = output,
  labels = "11_PCA_"
)

message("PCA and clustering analysis completed successfully!")

# Print summary
cat("\n=== Analysis Summary ===\n")
cat("GSE139602: PCA completed,", nrow(result_gse139602$scores), "samples\n")
cat("GSE15654: PCA completed,", nrow(result_gse15654$scores), "samples\n")
cat("GSE41919: PCA completed,", nrow(result_gse41919$scores), "samples\n")
cat("GSE57193: PCA completed,", nrow(result_gse57193$scores), "samples\n")
cat("\nResults saved to:", output, "\n")
