rm(list = ls()); gc()
ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "02_Mfuzz")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(ORIGINAL_DIR)

library(tidyverse)
library(stringr)
library(Mfuzz)
library(ggvenn)

# Helper function for Mfuzz clustering
run_mfuzz_clustering <- function(exp_file, group_file, cluster_num = 12, 
                                 output_dir, color = "#E64B35") {
  
  exp_data <- read.csv(exp_file, row.names = 1)
  group <- read.csv(group_file, row.names = 1)
  group$group <- sub("^[^:]*:\\s?", "", group$characteristics_ch1)
  
  # Prepare expression matrix
  exp_df <- exp_data[, rownames(group)]
  exp_df <- t(exp_df) %>% as.data.frame()
  exp_df$group <- group$group
  exp_df$group <- factor(exp_df$group, 
                         levels = c("Healthy", "eCLD", "Compensated Cirrhosis",
                                    "Decompesated Cirrhosis", "Acute-on-chronic liver failure"))
  
  # Aggregate by group means
  sample_agg <- aggregate(exp_df[, 1:nrow(exp_data)], 
                          by = list(exp_df$group), mean, na.rm = TRUE)
  rownames(sample_agg) <- sample_agg[, 1]
  sample_agg <- data.frame(t(sample_agg[, -1]))
  colnames(sample_agg) <- c("Healthy", "eCLD", "CC", "DC", "ACLF")
  
  # Create ExpressionSet and preprocess
  eset <- ExpressionSet(assayData = as.matrix(sample_agg))
  eset <- filter.NA(eset, thres = 0.25)
  eset <- fill.NA(eset, mode = 'mean')
  eset <- filter.std(eset, min.std = 0)
  eset <- standardise(eset)
  
  # Perform Mfuzz clustering
  set.seed(123)
  mfuzz_cluster <- mfuzz(eset, c = cluster_num, m = mestimate(eset))
  
  # Save cluster plots
  plot_clusters(eset, mfuzz_cluster, cluster_num, output_dir, color)
  
  # Save cluster memberships
  save_cluster_memberships(mfuzz_cluster, cluster_num, output_dir)
  
  return(list(eset = eset, cluster = mfuzz_cluster))
}

plot_clusters <- function(eset, mfuzz_cluster, cluster_num, output_dir, color) {
  png(file.path(output_dir, "02_Mfuzz_results_plot.png"), 
      width = 12, height = 8, res = 300, units = "in")
  mfuzz.plot2(eset, cl = mfuzz_cluster, mfrow = c(3, 4),
              time.labels = colnames(eset), centre = TRUE, x11 = FALSE,
              col = color, xlab = "", ylab = "Expression",
              cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.0)
  dev.off()
  
  pdf(file.path(output_dir, "02_Mfuzz_results_plot.pdf"), width = 12, height = 8)
  mfuzz.plot2(eset, cl = mfuzz_cluster, mfrow = c(3, 4),
              time.labels = colnames(eset), centre = TRUE, x11 = FALSE,
              col = color, xlab = "", ylab = "Expression",
              cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.0)
  dev.off()
}

save_cluster_memberships <- function(mfuzz_cluster, cluster_num, output_dir) {
  mfu_output <- file.path(output_dir, "mfuzz")
  if (!dir.exists(mfu_output)) {
    dir.create(mfu_output, recursive = TRUE)
  }
  
  for(i in 1:cluster_num) {
    cluster_genes <- names(mfuzz_cluster$cluster[unname(mfuzz_cluster$cluster) == i])
    write.csv(mfuzz_cluster[[4]][cluster_genes, i], 
              file.path(mfu_output, paste0("mfuzz_", i, ".csv")))
  }
  
  # Extract clusters 2 and 8 with coefficient threshold
  extract_cluster_genes <- function(cluster_id, coef_threshold = 0.4) {
    df <- read.csv(file.path(mfu_output, paste0("mfuzz_", cluster_id, ".csv")))
    df <- df[order(df$x, decreasing = TRUE), ]
    colnames(df) <- c("x", "Coefficient")
    df_filtered <- subset(df, Coefficient >= coef_threshold)
    write.csv(df_filtered, file.path(output_dir, paste0("mfuzz_", cluster_id, "_DEG.csv")))
    return(df_filtered)
  }
  
  cluster2_up <- extract_cluster_genes(2)
  cluster8_down <- extract_cluster_genes(8)
  
  # Combine clusters
  all_genes <- rbind(cluster2_up, cluster8_down)
  write.csv(all_genes, file.path(output_dir, "GSE139602_gene_all_DEG.csv"))
  write.csv(cluster2_up, file.path(output_dir, "GSE139602_cluster2_UP_DEG.csv"))
  write.csv(cluster8_down, file.path(output_dir, "GSE139602_cluster8_DOWN_DEG.csv"))
  
  return(list(cluster2_up = cluster2_up, cluster8_down = cluster8_down, all = all_genes))
}

# Helper function for Venn diagrams
create_venn <- function(df1, df2, plot_name, output_dir, 
                        label1 = NULL, label2 = NULL, 
                        fill_colors = c("#00A087", "#3C5488")) {
  
  # Extract gene names from dataframes
  get_genes <- function(df) {
    if ("X" %in% names(df)) {
      return(df$X)
    } else if ("x" %in% names(df)) {
      return(df$x)
    } else if (is.vector(df)) {
      return(df)
    } else {
      return(rownames(df))
    }
  }
  
  gene1 <- get_genes(df1)
  gene2 <- get_genes(df2)
  
  # Set labels if not provided
  if (is.null(label1)) label1 <- deparse(substitute(df1))
  if (is.null(label2)) label2 <- deparse(substitute(df2))
  
  venn_list <- setNames(list(gene1, gene2), c(label1, label2))
  
  p <- ggvenn(venn_list,
              text_size = 6.5,
              fill_color = fill_colors,
              fill_alpha = 0.7,
              stroke_color = "black",
              stroke_size = 1)
  
  ggsave(file.path(output_dir, paste0(plot_name, "_venn_diagram.pdf")), 
         p, width = 6, height = 5)
  ggsave(file.path(output_dir, paste0(plot_name, "_venn_diagram.png")), 
         p, width = 6, height = 5)
  
  # Save overlapping genes
  overlap_genes <- intersect(gene1, gene2)
  write.csv(data.frame(genes = overlap_genes), 
            file.path(output_dir, paste0(plot_name, "_gene.csv")), 
            row.names = FALSE)
  
  return(overlap_genes)
}

#### Run Mfuzz clustering on GSE139602 ####
mfuzz_results <- run_mfuzz_clustering(
  exp_file = "./00_rawdata/00.rawdata_GSE139602_exp.csv",
  group_file = "./00_rawdata/00.rawdata_GSE139602_group.csv",
  cluster_num = 12,
  output_dir = output,
  color = "#E64B35"
)

# Extract cluster results
cluster2_up <- read.csv(file.path(output, "GSE139602_cluster2_UP_DEG.csv"), row.names = 1)
cluster8_down <- read.csv(file.path(output, "GSE139602_cluster8_DOWN_DEG.csv"), row.names = 1)
gse139602_all <- read.csv(file.path(output, "GSE139602_gene_all_DEG.csv"), row.names = 1)

#### Load DEG data from other datasets ####
load_deg_data <- function(file_path, direction = "all") {
  deg <- read.csv(file_path, row.names = 1)
  
  if (direction == "up") {
    return(subset(deg, P.Value <= 0.05 & logFC >= 0.5))
  } else if (direction == "down") {
    return(subset(deg, P.Value <= 0.05 & logFC <= -0.5))
  } else {
    up <- subset(deg, P.Value <= 0.05 & logFC >= 0.5)
    down <- subset(deg, P.Value <= 0.05 & logFC <= -0.5)
    return(rbind(up, down))
  }
}

gse41919_up <- load_deg_data("./01_DEG/01.DEG_GSE41919_with_HE-without_HE_all.csv", "up")
gse41919_down <- load_deg_data("./01_DEG/01.DEG_GSE41919_with_HE-without_HE_all.csv", "down")
gse41919_all <- load_deg_data("./01_DEG/01.DEG_GSE41919_with_HE-without_HE_all.csv", "all")

gse57193_up <- load_deg_data("./01_DEG/01.DEG_GSE57193_cirrhosis_with_HE-cirrhosis_all.csv", "up")
gse57193_down <- load_deg_data("./01_DEG/01.DEG_GSE57193_cirrhosis_with_HE-cirrhosis_all.csv", "down")
gse57193_all <- load_deg_data("./01_DEG/01.DEG_GSE57193_cirrhosis_with_HE-cirrhosis_all.csv", "all")

gse15654 <- read.csv("./01_DEG/01.DEG_GSE15654_Poor_prognosis-Good_prognosis_all.csv", row.names = 1)
gse15654_up <- subset(gse15654, P.Value <= 0.05 & logFC >= 0.5)
gse15654_down <- subset(gse15654, P.Value <= 0.05 & logFC <= -0.5)
gse15654_all <- rbind(gse15654_up, gse15654_down)

#### Venn diagrams between GSE41919 and GSE57193 ####
venn_41919_57193_up <- create_venn(gse41919_up, gse57193_up, 
                                   "GSE41919_GSE57193_UP", output,
                                   "GSE41919_UP", "GSE57193_UP")

venn_41919_57193_down <- create_venn(gse41919_down, gse57193_down,
                                     "GSE41919_GSE57193_DOWN", output,
                                     "GSE41919_DOWN", "GSE57193_DOWN")

venn_41919_57193_all <- create_venn(gse41919_all, gse57193_all,
                                    "GSE41919_GSE57193_ALL", output,
                                    "GSE41919_ALL", "GSE57193_ALL")

# Read common genes
common_up <- read.csv(file.path(output, "GSE41919_GSE57193_UP_gene.csv"), row.names = 1)
common_down <- read.csv(file.path(output, "GSE41919_GSE57193_DOWN_gene.csv"), row.names = 1)
common_all <- read.csv(file.path(output, "GSE41919_GSE57193_ALL_gene.csv"), row.names = 1)

#### Venn diagrams with GSE139602 clusters ####
venn_he_up <- create_venn(common_up, cluster2_up, 
                          "HE_GSE139602_UP", output,
                          "HE_Common_UP", "GSE139602_Cluster2_UP")

venn_he_down <- create_venn(common_down, cluster8_down,
                            "HE_GSE139602_DOWN", output,
                            "HE_Common_DOWN", "GSE139602_Cluster8_DOWN")

venn_he_all <- create_venn(common_all, gse139602_all,
                           "HE_GSE139602_ALL", output,
                           "HE_Common_ALL", "GSE139602_ALL")

#### Venn diagrams with GSE15654 ####
venn_15654_up <- create_venn(common_up, gse15654_up,
                             "HE_GSE15654_UP", output,
                             "HE_Common_UP", "GSE15654_UP")

venn_15654_down <- create_venn(common_down, gse15654_down,
                               "HE_GSE15654_DOWN", output,
                               "HE_Common_DOWN", "GSE15654_DOWN")

venn_15654_all <- create_venn(common_all, gse15654_all,
                              "HE_GSE15654_ALL", output,
                              "HE_Common_ALL", "GSE15654_ALL")
