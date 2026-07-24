rm(list = ls()); gc()
ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "09_SEM")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(ORIGINAL_DIR)

library(tidyverse)
library(mediation)
library(limma)
library(sva)
library(corrplot)
library(pheatmap)
library(lavaan)
library(semPlot)
library(dplyr)
library(ggplot2)
library(scales)
library(grid)
library(ggrepel)

# Helper function for normality testing
test_normality <- function(exp_data, by_group = FALSE, group_info = NULL) {
  if (is.data.frame(exp_data)) exp_data <- as.matrix(exp_data)
  
  if (!by_group) {
    pvals <- apply(exp_data, 1, function(x) {
      x <- as.numeric(x)
      if (all(is.na(x)) || sum(!is.na(x)) < 3) return(NA_real_)
      res <- tryCatch(shapiro.test(x), error = function(e) NULL)
      if (is.null(res)) NA_real_ else as.numeric(res$p.value)
    })
    
    df <- data.frame(gene = rownames(exp_data), p.value = as.numeric(pvals), stringsAsFactors = FALSE)
    df$adj.p <- p.adjust(df$p.value, method = "BH")
    return(df[order(df$adj.p), ])
  } else {
    if (is.null(group_info)) stop("group_info required for by_group = TRUE")
    groups <- levels(factor(group_info))
    res_list <- list()
    
    for (g in groups) {
      cols <- which(group_info == g)
      if (length(cols) < 3) {
        tmp <- data.frame(gene = rownames(exp_data), group = g, p.value = NA_real_, stringsAsFactors = FALSE)
      } else {
        pvals <- apply(exp_data[, cols, drop = FALSE], 1, function(x) {
          x <- as.numeric(x)
          if (all(is.na(x)) || sum(!is.na(x)) < 3) return(NA_real_)
          res <- tryCatch(shapiro.test(x), error = function(e) NULL)
          if (is.null(res)) NA_real_ else as.numeric(res$p.value)
        })
        tmp <- data.frame(gene = rownames(exp_data), group = g, p.value = as.numeric(pvals), stringsAsFactors = FALSE)
        tmp$adj.p <- p.adjust(tmp$p.value, method = "BH")
      }
      res_list[[g]] <- tmp
    }
    
    df_all <- do.call(rbind, lapply(res_list, function(x) {
      if (!"adj.p" %in% colnames(x)) x$adj.p <- NA_real_
      x
    }))
    rownames(df_all) <- NULL
    return(df_all[order(df_all$group, df_all$adj.p), ])
  }
}

# Helper function for batch correction
perform_batch_correction <- function(exp_data, batch_info, group_info = NULL) {
  sample_order <- colnames(exp_data)
  batch_factor <- batch_info[match(sample_order, names(batch_info))]
  unique_batches <- unique(batch_factor)
  
  cat("Unique batches:", unique_batches, "\n")
  cat("Number of unique batches:", length(unique_batches), "\n")
  
  if (length(unique_batches) > 1 && !any(is.na(unique_batches))) {
    if (!is.null(group_info)) {
      group_factor <- group_info[match(sample_order, names(group_info))]
      mod <- model.matrix(~ group_factor)
      exp_corrected <- ComBat(dat = exp_data, batch = batch_factor, mod = mod)
    } else {
      exp_corrected <- ComBat(dat = exp_data, batch = batch_factor)
    }
  } else {
    warning("Batch correction not performed - insufficient batches")
    exp_corrected <- exp_data
  }
  
  return(exp_corrected)
}

# Helper function for PCA plot
create_pca_plot <- function(exp_data, group_info, project_info = NULL, title = "PCA Plot") {
  pca <- prcomp(t(exp_data), scale. = TRUE)
  pca_data <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    Sample = colnames(exp_data),
    Group = group_info
  )
  
  if (!is.null(project_info)) {
    pca_data$Project <- project_info
  }
  
  var_explained <- pca$sdev^2 / sum(pca$sdev^2)
  pc1_var <- round(var_explained[1] * 100, 2)
  pc2_var <- round(var_explained[2] * 100, 2)
  
  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = factor(Group))) +
    geom_point(size = 2, alpha = 0.7) +
    labs(
      x = paste0("PC1 (", pc1_var, "%)"),
      y = paste0("PC2 (", pc2_var, "%)"),
      title = title,
      color = "Group"
    ) +
    theme_minimal()
  
  return(list(plot = p, pca = pca, data = pca_data))
}

# Helper function for correlation heatmap
create_correlation_heatmap <- function(exp_data, output_dir, filename_prefix) {
  cor_matrix <- cor(t(exp_data), use = "pairwise.complete.obs")
  p_matrix <- cor.mtest(cor_matrix)$p
  
  genes <- colnames(cor_matrix)
  
  cor_df <- as.data.frame(as.table(cor_matrix), stringsAsFactors = FALSE) %>%
    rename(gene_y = Var1, gene_x = Var2, r = Freq)
  
  p_df <- as.data.frame(as.table(p_matrix), stringsAsFactors = FALSE) %>%
    rename(gene_y = Var1, gene_x = Var2, p = Freq)
  
  plot_df <- cor_df %>%
    left_join(p_df, by = c("gene_y", "gene_x")) %>%
    mutate(
      sig = case_when(
        p < 0.001 ~ "***",
        p < 0.01  ~ "**",
        p < 0.05  ~ "*",
        TRUE      ~ ""
      ),
      label = ifelse(gene_x == gene_y, sprintf("%.2f", r), sprintf("%.2f%s", r, sig))
    )
  
  plot_df$gene_x <- factor(plot_df$gene_x, levels = genes)
  plot_df$gene_y <- factor(plot_df$gene_y, levels = rev(genes))
  
  p <- ggplot(plot_df, aes(x = gene_x, y = gene_y, fill = r)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = label), size = 2.8, color = "black") +
    scale_fill_gradient2(
      low = "#8491B4", mid = "white", high = "#F39B7F",
      midpoint = 0, limits = c(-1, 1), name = "Correlation"
    ) +
    scale_x_discrete(position = "top") +
    coord_fixed() +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text.x.top = element_text(face = "italic", color = "black", 
                                     angle = 45, hjust = 0, vjust = 0,
                                     margin = margin(b = 0)),
      axis.text.y.left = element_text(face = "italic", color = "black",
                                      margin = margin(r = 0))
    )
  
  ggsave(file.path(output_dir, paste0(filename_prefix, "_corrplot_full.pdf")), 
         p, width = 4.5, height = 4.5, dpi = 300)
  ggsave(file.path(output_dir, paste0(filename_prefix, "_corrplot_full.png")), 
         p, width = 4.5, height = 4.5, dpi = 300)
  
  return(list(plot = p, cor_matrix = cor_matrix, p_matrix = p_matrix))
}

#### Load and prepare data ####
# Load expression data
exp01 <- read.csv(file.path("00_rawdata", "00.rawdata_GSE139602_exp.csv"), row.names = 1)
group01 <- read.csv(file.path("00_rawdata", "00.rawdata_GSE139602_group.csv"), row.names = 1)
group01 <- group01[colnames(exp01), , drop = FALSE]

# Convert group labels
group01$group <- ifelse(group01$characteristics_ch1 == "disease state: Healthy", 0,
                 ifelse(group01$characteristics_ch1 == "disease state: eCLD", 1,
                 ifelse(group01$characteristics_ch1 == "disease state: Compensated Cirrhosis", 2,
                 ifelse(group01$characteristics_ch1 == "disease state: Decompesated Cirrhosis", 3,
                 ifelse(group01$characteristics_ch1 == "disease state: Acute-on-chronic liver failure", 4, NA)))))

exp02 <- read.csv(file.path("00_rawdata", "00.rawdata_GSE15654_exp.csv"), row.names = 1)
group02 <- read.csv(file.path("00_rawdata", "00.rawdata_GSE15654_group.csv"), row.names = 1)

# Load models and gene lists
model1 <- readRDS(file.path("07_xgboost", "07_xgboost_GSE139602_final_xgboost_model.Rdata"))
model2 <- readRDS(file.path("07_xgboost", "07_xgboost_GSE15654_final_xgboost_model_cox.Rdata"))

model1_gene <- read.csv(file.path("07_xgboost", "07_xgboost_GSE139602_shap_mat.csv"), row.names = 1)
model2_gene <- read.csv(file.path("07_xgboost", "07_xgboost_GSE15654_shap_mat.csv"), row.names = 1)

all_features <- c(colnames(model1_gene), colnames(model2_gene))

# Find common genes
common_genes <- intersect(rownames(exp01), rownames(exp02))
exp01 <- exp01[common_genes, ]
exp02 <- exp02[common_genes, ]
exp_merged <- cbind(exp01, exp02)

# Prepare group information
group01$project <- "GSE139602"
group02$project <- "GSE15654"
group_merged <- data.frame(
  row.names = c(rownames(group01), rownames(group02)),
  group = c(group01$group, group02$group),
  project = c(group01$project, group02$project)
)
group_merged <- group_merged[colnames(exp_merged), ]

#### Batch correction ####
# Pre-correction visualization
png(file.path(output, "09_SEM_pre_merged.png"), width = 5, height = 4, res = 300, units = "in")
boxplot(exp_merged, xaxt = "n", col = "lightblue", 
        main = "Gene Expression Distribution (Pre-batch Correction)", 
        ylab = "Expression Value")
dev.off()

pdf(file.path(output, "09_SEM_pre_merged.pdf"), width = 5, height = 4)
boxplot(exp_merged, xaxt = "n", col = "lightblue", 
        main = "Gene Expression Distribution (Pre-batch Correction)", 
        ylab = "Expression Value")
dev.off()

# Perform batch correction
group_merged$sample <- rownames(group_merged)
sample_order <- colnames(exp_merged)
batch_info <- group_merged$project[match(sample_order, group_merged$sample)]
group_info <- group_merged$group[match(sample_order, group_merged$sample)]
names(batch_info) <- sample_order
names(group_info) <- sample_order

exp_corrected <- perform_batch_correction(exp_merged, batch_info, group_info)

# Post-correction visualization
png(file.path(output, "09_SEM_post_merged.png"), width = 5, height = 4, res = 300, units = "in")
boxplot(exp_corrected, xaxt = "n", col = "lightblue", 
        main = "Gene Expression Distribution (Post-batch Correction)", 
        ylab = "Expression Value")
dev.off()

pdf(file.path(output, "09_SEM_post_merged.pdf"), width = 5, height = 4)
boxplot(exp_corrected, xaxt = "n", col = "lightblue", 
        main = "Gene Expression Distribution (Post-batch Correction)", 
        ylab = "Expression Value")
dev.off()

# Save corrected data
write.csv(group_merged, file.path(output, "09_SEM_group_merged.csv"))
write.csv(exp_merged, file.path(output, "09_SEM_exp_merged.csv"))
write.csv(exp_corrected, file.path(output, "09_SEM_exp_corrected.csv"))

# PCA plots
pca_pre <- create_pca_plot(exp_merged, group_merged$group, group_merged$project, 
                           "PCA of Pre-corrected Expression Data")
ggsave(file.path(output, "09_SEM_pca_pre.png"), pca_pre$plot, width = 6, height = 5)
ggsave(file.path(output, "09_SEM_pca_pre.pdf"), pca_pre$plot, width = 6, height = 5)

pca_post <- create_pca_plot(exp_corrected, group_merged$group, group_merged$project,
                            "PCA of Post-corrected Expression Data")
ggsave(file.path(output, "09_SEM_pca_post.png"), pca_post$plot, width = 6, height = 5)
ggsave(file.path(output, "09_SEM_pca_post.pdf"), pca_post$plot, width = 6, height = 5)

#### Normality testing ####
exp_data <- exp_corrected[all_features, ]
norm_res <- test_normality(exp_data, by_group = FALSE)
write.csv(norm_res, file.path(output, "09_SEM_normality_shapiro_overall.csv"))

#### Correlation analysis ####
cor_results <- create_correlation_heatmap(exp_data, output, "09_SEM")

#### Mediation analysis ####
# Prepare data for mediation
X <- as.data.frame(t(exp_corrected[all_features, ]))
X <- scale(X)
X <- as.data.frame(X)

# Predict progression and survival
progression <- predict(model1, scale(t(exp_corrected[colnames(model1_gene), ])))
survival <- predict(model2, scale(t(exp_corrected[colnames(model2_gene), ])))

X$progression <- progression
X$survival <- survival

# Define mediator and outcome models
mediator_vars <- c("NPC2", "TLN1", "TUBA1C", "LRRC32", "PRB2")
outcome_vars <- c("SOX9", "RNASE4", "SERPINA3")

# Linear regression models for mediation
model_a <- lm(progression ~ NPC2 + TLN1 + TUBA1C + LRRC32 + PRB2, data = X)
model_bc <- lm(survival ~ NPC2 + TLN1 + TUBA1C + LRRC32 + PRB2 + SOX9 + RNASE4 + SERPINA3 + progression, data = X)

# Perform mediation analysis
mediation_results <- list()
for (treat in mediator_vars) {
  med_result <- mediate(model.m = model_a,
                        model.y = model_bc,
                        treat = treat,
                        mediator = "progression",
                        boot = TRUE,
                        sims = 1000)
  mediation_results[[treat]] <- summary(med_result)
}

#### Structural Equation Modeling (SEM) ####
# Define SEM model
sem_model <- '
  # Direct effects on progression
  progression ~ a_NPC2*NPC2 + a_TLN1*TLN1 + a_TUBA1C*TUBA1C + a_LRRC32*LRRC32 + a_PRB2*PRB2
  
  # Effects on outcome mediators
  SOX9 ~ b_SOX9*progression
  RNASE4 ~ b_RNASE4*progression
  SERPINA3 ~ b_SERPINA3*progression
  
  # Effects on survival
  survival ~ c_SOX9*SOX9 + c_RNASE4*RNASE4 + c_SERPINA3*SERPINA3 + c_progression*progression
  
  # Indirect effects
  indirect_NPC2_SOX9 := a_NPC2 * b_SOX9 * c_SOX9
  indirect_NPC2_RNASE4 := a_NPC2 * b_RNASE4 * c_RNASE4
  indirect_NPC2_SERPINA3 := a_NPC2 * b_SERPINA3 * c_SERPINA3
  
  indirect_TLN1_SOX9 := a_TLN1 * b_SOX9 * c_SOX9
  indirect_TLN1_RNASE4 := a_TLN1 * b_RNASE4 * c_RNASE4
  indirect_TLN1_SERPINA3 := a_TLN1 * b_SERPINA3 * c_SERPINA3
  
  indirect_TUBA1C_SOX9 := a_TUBA1C * b_SOX9 * c_SOX9
  indirect_TUBA1C_RNASE4 := a_TUBA1C * b_RNASE4 * c_RNASE4
  indirect_TUBA1C_SERPINA3 := a_TUBA1C * b_SERPINA3 * c_SERPINA3
  
  indirect_LRRC32_SOX9 := a_LRRC32 * b_SOX9 * c_SOX9
  indirect_LRRC32_RNASE4 := a_LRRC32 * b_RNASE4 * c_RNASE4
  indirect_LRRC32_SERPINA3 := a_LRRC32 * b_SERPINA3 * c_SERPINA3
  
  indirect_PRB2_SOX9 := a_PRB2 * b_SOX9 * c_SOX9
  indirect_PRB2_RNASE4 := a_PRB2 * b_RNASE4 * c_RNASE4
  indirect_PRB2_SERPINA3 := a_PRB2 * b_SERPINA3 * c_SERPINA3
'

# Fit SEM model
set.seed(123)
fit <- sem(sem_model, data = X, se = "bootstrap", bootstrap = 2000, missing = "FIML")

# Extract results
fit_summary <- summary(fit, standardized = TRUE, rsquare = TRUE)
fit_measures <- fitmeasures(fit, c("cfi", "tli", "rmsea", "srmr"))
parameter_estimates <- parameterEstimates(fit, standardized = TRUE)

# Save results
write.csv(parameter_estimates, file.path(output, "09_SEM_sem_results.csv"))
write.csv(fit_measures, file.path(output, "09_SEM_fit_measures.csv"))

# Extract indirect effects
indirect_effects <- parameter_estimates[grepl("^indirect_", parameter_estimates$label), ]
indirect_effects <- indirect_effects %>%
  mutate(p_adj = p.adjust(pvalue, method = "fdr"))

write.csv(indirect_effects, file.path(output, "09_SEM_indirect_effects.csv"))
write.csv(indirect_effects[indirect_effects$pvalue < 0.05, ], 
          file.path(output, "09_SEM_significant_indirect.csv"))

#### Create SEM path diagram ####
# Prepare node positions
nodes <- bind_rows(
  tibble(node = c("NPC2", "TLN1", "TUBA1C", "LRRC32", "PRB2"),
         x = c(-4, -2, 0, 2, 4), y = 4),
  tibble(node = "progression", x = 0, y = 3),
  tibble(node = c("SOX9", "RNASE4", "SERPINA3"),
         x = c(-2, 0, 2), y = 2),
  tibble(node = "survival", x = 0, y = 1)
)

# Prepare edges
pe <- parameter_estimates %>%
  filter(op == "~") %>%
  transmute(
    from = rhs,
    to = lhs,
    est = std.all,
    pval = pvalue
  )

edges <- pe %>%
  inner_join(nodes %>% rename(from = node, x = x, y = y), by = "from") %>%
  inner_join(nodes %>% rename(to = node, xend = x, yend = y), by = "to") %>%
  mutate(
    dx = xend - x,
    curve_type = case_when(
      dx > 0 ~ "right",
      dx < 0 ~ "left",
      TRUE ~ "vertical"
    ),
    xm = (x + xend) / 2 + ifelse(dx == 0, 0.18, 0),
    ym = (y + yend) / 2 + ifelse(dx == 0, 0, 0.12 * sign(dx))
  )

# Prepare node labels with italics for genes
gene_nodes <- c("NPC2", "TLN1", "TUBA1C", "LRRC32", "PRB2", "SOX9", "RNASE4", "SERPINA3")
nodes <- nodes %>%
  mutate(
    label_expr = ifelse(
      node %in% gene_nodes,
      paste0("italic('", node, "')"),
      paste0("'", node, "'")
    )
  )

# Create SEM path diagram
p_sem <- ggplot() +
  geom_curve(
    data = edges %>% filter(curve_type == "right"),
    aes(x = x, y = y, xend = xend, yend = yend, color = est, linewidth = abs(est)),
    curvature = 0.25, alpha = 0.9, lineend = "round",
    arrow = arrow(length = unit(0.18, "cm"), type = "closed")
  ) +
  geom_curve(
    data = edges %>% filter(curve_type == "left"),
    aes(x = x, y = y, xend = xend, yend = yend, color = est, linewidth = abs(est)),
    curvature = -0.25, alpha = 0.9, lineend = "round",
    arrow = arrow(length = unit(0.18, "cm"), type = "closed")
  ) +
  geom_curve(
    data = edges %>% filter(curve_type == "vertical"),
    aes(x = x, y = y, xend = xend, yend = yend, color = est, linewidth = abs(est)),
    curvature = 0.20, alpha = 0.9, lineend = "round",
    arrow = arrow(length = unit(0.18, "cm"), type = "closed")
  ) +
  geom_point(
    data = nodes,
    aes(x = x, y = y),
    shape = 21, size = 14, stroke = 0.9,
    fill = alpha("#F2F2F2", 0.45), color = "#4D4D4D"
  ) +
  geom_label_repel(
    data = edges,
    aes(x = xm, y = ym, label = sprintf("%.2f", est)),
    size = 3.0, color = "black", fill = "#F7F7F7",
    alpha = 0.95, label.size = 0.15, label.r = unit(0.08, "lines"),
    box.padding = 0.12, point.padding = 0.05,
    min.segment.length = 0, segment.color = NA,
    segment.size = 0.25, force = 1.2,
    max.overlaps = Inf, show.legend = FALSE
  ) +
  geom_text(
    data = nodes,
    aes(x = x, y = y, label = label_expr),
    parse = TRUE, size = 4.3, fontface = "bold", color = "#1F1F1F"
  ) +
  scale_color_gradient2(
    low = "#2C7BB6", mid = "#BDBDBD", high = "#D7191C", midpoint = 0,
    name = "Standardized coefficient"
  ) +
  scale_linewidth(range = c(0.7, 2.8), name = "|Standardized coefficient|") +
  coord_cartesian(xlim = c(-5, 5), ylim = c(0.6, 4.4), clip = "off") +
  theme_void(base_size = 13) +
  theme(
    legend.position = "right",
    plot.margin = margin(10, 20, 10, 20)
  )

ggsave(file.path(output, "09_SEM_semplot.png"), p_sem, width = 7.5, height = 6, dpi = 300)
ggsave(file.path(output, "09_SEM_semplot.pdf"), p_sem, width = 7.5, height = 6, dpi = 300)

message("SEM analysis completed successfully!")

# Print summary
cat("\n=== SEM Fit Measures ===\n")
print(fit_measures)

cat("\n=== Significant Indirect Effects ===\n")
print(indirect_effects[indirect_effects$pvalue < 0.05, c("label", "est", "pvalue", "p_adj")])
