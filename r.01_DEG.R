rm(list = ls()); gc()
ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "01_DEG")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(ORIGINAL_DIR)

library(limma)
library(ggplot2)
library(ggrepel)
library(ggvenn)

# Helper function for DEG analysis
run_deg_analysis <- function(exp_file, group_file, group_subset, 
                             group_names, contrast_formula, 
                             output_prefix, logFC_cutoff = 0.5, p_cutoff = 0.05) {
  
  exp <- read.csv(exp_file, header = TRUE, row.names = 1)
  group <- read.csv(group_file, header = TRUE, row.names = 1)
  
  # Subset groups if specified
  if (!is.null(group_subset)) {
    group <- subset(group, group != group_subset)
  }
  exp <- exp[, rownames(group)]
  
  # Design matrix
  design <- model.matrix(~0 + factor(group$group))
  rownames(design) <- rownames(group)
  colnames(design) <- group_names
  
  # Contrast matrix
  contrast.matrix <- makeContrasts(contrast_formula, levels = design)
  
  # Fit linear model
  set.seed(123)
  fit <- lmFit(exp, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  # Extract results
  tempOutput <- topTable(fit2, coef = 1, n = Inf)
  nrDEG <- na.omit(tempOutput)
  
  # Save all DEG results
  write.csv(nrDEG, file = file.path(output, paste0(output_prefix, "_all.csv")))
  
  # Filter significant DEGs
  nrDEG_sig <- subset(nrDEG, abs(logFC) > logFC_cutoff & P.Value < p_cutoff)
  gene_list <- rownames(nrDEG_sig)
  
  # Add change column
  nrDEG$change <- ifelse(nrDEG$P.Value < p_cutoff & abs(nrDEG$logFC) >= logFC_cutoff,
                         ifelse(nrDEG$logFC > 0, 'Up', 'Down'),
                         'NS')
  
  write.csv(nrDEG, file = file.path(output, paste0(output_prefix, "_all_2.csv")))
  
  # Volcano plot
  nrDEG$Label <- ""
  nrDEG <- nrDEG[order(nrDEG$P.Value), ]
  nrDEG$Gene <- rownames(nrDEG)
  
  p <- ggplot(nrDEG, aes(x = logFC, y = -log10(P.Value), colour = change)) +
    geom_point(alpha = 0.4, size = 2) +
    scale_color_manual(values = c("#4DBBD5", "#d2dae2", "#E64B35")) +
    geom_vline(xintercept = 0, lty = 4, col = "black", lwd = 0.8) +
    geom_vline(xintercept = logFC_cutoff, lty = 2, col = "grey50", lwd = 0.4) +
    geom_vline(xintercept = -logFC_cutoff, lty = 2, col = "grey50", lwd = 0.4) +
    geom_hline(yintercept = -log10(p_cutoff), lty = 4, col = "black", lwd = 0.8) +
    labs(title = output_prefix,
         x = expression(log[2]~FC),
         y = expression(-log[10]~italic(p))) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text = element_text(size = 12),
      legend.position = "right",
      legend.text = element_text(size = 12),
      legend.title = element_blank()
    ) +
    geom_text_repel(data = nrDEG, aes(label = Label),
                    size = 3.5,
                    box.padding = unit(0.5, "lines"),
                    point.padding = unit(0.8, "lines"),
                    segment.color = "black",
                    show.legend = FALSE,
                    max.overlaps = 10000)
  
  ggsave(file.path(output, paste0(output_prefix, "_volcano.png")), 
         plot = p, width = 5, height = 5, units = "in")
  ggsave(file.path(output, paste0(output_prefix, "_volcano.pdf")), 
         plot = p, width = 5, height = 5, units = "in")
  
  return(list(deg_all = nrDEG, deg_sig = nrDEG_sig, gene_list = gene_list))
}

#### GSE41919 ####
result1 <- run_deg_analysis(
  exp_file = './00_rawdata/00.rawdata_GSE41919_exp.csv',
  group_file = './00_rawdata/00.rawdata_GSE41919_group.csv',
  group_subset = "non-cirrhotic control",
  group_names = c('with_HE', 'without_HE'),
  contrast_formula = "with_HE-without_HE",
  output_prefix = "01.DEG_GSE41919_with_HE-without_HE"
)

#### GSE57193 ####
result2 <- run_deg_analysis(
  exp_file = './00_rawdata/00.rawdata_GSE57193_exp.csv',
  group_file = './00_rawdata/00.rawdata_GSE57193_group.csv',
  group_subset = "healthy",
  group_names = c('cirrhosis', 'cirrhosis_with_HE'),
  contrast_formula = "cirrhosis_with_HE-cirrhosis",
  output_prefix = "01.DEG_GSE57193_cirrhosis_with_HE-cirrhosis"
)

#### GSE15654 ####
exp <- read.csv('./00_rawdata/00.rawdata_GSE15654_exp.csv', header = TRUE, row.names = 1)
group <- read.csv('./00_rawdata/00.rawdata_GSE15654_group.csv', header = TRUE, row.names = 1)
exp <- exp[, rownames(group)]
exp <- scale(exp)

design <- model.matrix(~0 + factor(group$group))
rownames(design) <- rownames(group)
colnames(design) <- c('Good_prognosis', 'Poor_prognosis')
contrast.matrix <- makeContrasts("Poor_prognosis-Good_prognosis", levels = design)

set.seed(123)
fit <- lmFit(exp, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

tempOutput <- topTable(fit2, coef = 1, n = Inf)
nrDEG <- na.omit(tempOutput)
write.csv(nrDEG, file = "./01_DEG/01.DEG_GSE15654_Poor_prognosis-Good_prognosis_all.csv")

nrDEG_sig <- subset(nrDEG, abs(logFC) > 0.5 & P.Value < 0.05)
gene_list <- rownames(nrDEG_sig)

nrDEG$change <- ifelse(nrDEG$P.Value < 0.05 & abs(nrDEG$logFC) >= 0.5,
                       ifelse(nrDEG$logFC > 0, 'Up', 'Down'),
                       'NS')
write.csv(nrDEG, file = "./01_DEG/01.DEG_GSE15654_Poor_prognosis-Good_prognosis_all_2.csv")

# Volcano plot for GSE15654
nrDEG$Label <- ""
nrDEG <- nrDEG[order(nrDEG$P.Value), ]
nrDEG$Gene <- rownames(nrDEG)

p <- ggplot(nrDEG, aes(x = logFC, y = -log10(P.Value), colour = change)) +
  geom_point(alpha = 0.4, size = 2) +
  scale_color_manual(values = c("#4DBBD5", "#d2dae2", "#E64B35")) +
  geom_vline(xintercept = 0, lty = 4, col = "black", lwd = 0.8) +
  geom_vline(xintercept = 0.5, lty = 2, col = "grey50", lwd = 0.4) +
  geom_vline(xintercept = -0.5, lty = 2, col = "grey50", lwd = 0.4) +
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.8) +
  labs(title = "GSE15654",
       x = expression(log[2]~FC),
       y = expression(-log[10]~italic(p))) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.title = element_blank()
  ) +
  geom_text_repel(data = nrDEG, aes(label = Label),
                  size = 3.5,
                  box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"),
                  segment.color = "black",
                  show.legend = FALSE,
                  max.overlaps = 10000)

ggsave("./01_DEG/01.DEG_GSE15654_DEG_volcano.png", plot = p, width = 5, height = 5, units = "in")
ggsave("./01_DEG/01.DEG_GSE15654_DEG_volcano.pdf", plot = p, width = 5, height = 5, units = "in")

#### Venn Diagram ####
# Extract gene lists from all three datasets
gene_lists <- list(
  GSE41919 = result1$gene_list,
  GSE57193 = result2$gene_list,
  GSE15654 = gene_list
)

# Create Venn diagram
venn_plot <- ggvenn(
  gene_lists,
  fill_color = c("#E64B35", "#4DBBD5", "#00A087"),
  stroke_size = 0.5,
  set_name_size = 4
)

ggsave(file.path(output, "01.DEG_venn_plot.png"), 
       plot = venn_plot, width = 6, height = 6, units = "in")
ggsave(file.path(output, "01.DEG_venn_plot.pdf"), 
       plot = venn_plot, width = 6, height = 6, units = "in")

# Save overlapping genes
overlap_genes <- Reduce(intersect, gene_lists)
write.csv(data.frame(genes = overlap_genes), 
          file = file.path(output, "01.DEG_overlap_genes.csv"), 
          row.names = FALSE)
