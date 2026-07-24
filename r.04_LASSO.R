rm(list = ls()); gc()
ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "04_LASSO")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(ORIGINAL_DIR)

library(glmnet)
library(ggplot2)
library(reshape2)
library(ggsci)
library(survival)
library(survminer)
library(ggtext)

# Helper function for LASSO analysis
run_lasso_analysis <- function(exp_data, gene_list, group_data = NULL, 
                               surv_data = NULL, output_dir, prefix,
                               family_type = "gaussian", lambda_label = NULL) {
  
  # Filter expression data
  exp_filtered <- exp_data[gene_list, ]
  x <- as.matrix(t(exp_filtered))
  
  # Set response variable
  if (family_type == "cox") {
    y <- Surv(surv_data$days_to_death, surv_data$death)
  } else {
    y <- group_data$group
  }
  
  # Cross-validation
  set.seed(123)
  cvfit <- cv.glmnet(x, y, family = family_type, nlambda = 100, alpha = 1)
  fit <- glmnet(x, y, family = family_type, nlambda = 100, alpha = 1)
  
  # Get optimal lambda
  lambda_min <- round(cvfit$lambda.min, digits = 4)
  message("Optimal lambda for ", prefix, ": ", lambda_min)
  
  # Extract coefficients
  coef_matrix <- coef(fit)
  tmp <- as.data.frame(as.matrix(coef_matrix))
  tmp$coef <- rownames(tmp)
  tmp <- reshape::melt(tmp, id = "coef")
  tmp$variable <- as.numeric(gsub("s", "", tmp$variable))
  tmp$coef <- gsub('_', '-', tmp$coef)
  tmp$lambda <- fit$lambda[tmp$variable + 1]
  tmp$norm <- apply(abs(coef_matrix[-1, ]), 2, sum)[tmp$variable + 1]
  
  # Fit with optimal lambda
  fit2 <- glmnet(x = x, y = y, alpha = 1, family = family_type, lambda = cvfit$lambda.min)
  selected_genes <- rownames(fit2$beta)[as.numeric(fit2$beta) != 0]
  
  # Create coefficient plot
  p1 <- create_coefficient_plot(tmp, cvfit, selected_genes, lambda_min, prefix)
  
  # Create CV plot
  p2 <- create_cv_plot(cvfit, lambda_min, prefix)
  
  # Save plots
  ggsave(file.path(output_dir, paste0(prefix, "result01.pdf")), p1, w = 6, h = 6)
  ggsave(file.path(output_dir, paste0(prefix, "result01.png")), p1, w = 6, h = 6)
  ggsave(file.path(output_dir, paste0(prefix, "result02.pdf")), p2, w = 6, h = 6)
  ggsave(file.path(output_dir, paste0(prefix, "result02.png")), p2, w = 6, h = 6)
  
  # Save selected genes
  write.csv(data.frame(genes = selected_genes), 
            file.path(output_dir, paste0(prefix, "result.csv")), 
            row.names = FALSE)
  
  return(list(selected_genes = selected_genes, lambda_min = lambda_min, cvfit = cvfit))
}

# Helper function for coefficient plot
create_coefficient_plot <- function(tmp, cvfit, selected_genes, lambda_min, prefix) {
  p <- ggplot(tmp, aes(log(lambda), value, color = coef)) +
    geom_vline(xintercept = log(cvfit$lambda.min),
               size = 0.8, color = 'grey60',
               alpha = 0.8, linetype = 2) +
    geom_line(size = 1) +
    xlab(expression(paste(lambda, " (log scale)"))) +
    ylab('Coefficients') +
    theme_bw(base_rect_size = 2) +
    scale_color_manual(
      name = "Coefficient",
      values = c(pal_npg()(9), pal_d3()(7)),
      breaks = unique(tmp$coef),
      labels = function(x) {
        sapply(x, function(gene) {
          if (gene %in% selected_genes) {
            paste0("<span style='color:red;'><i>", gene, "</i></span>")
          } else {
            paste0("<i>", gene, "</i>")
          }
        })
      }
    ) +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = c(0.01, 0.01)) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 15, color = 'black'),
      axis.text = element_text(size = 12, color = 'black'),
      legend.title = element_blank(),
      legend.text = element_markdown(size = 14),
      legend.position = "right",
      legend.direction = "vertical",
      legend.key.height = grid::unit(6, "mm"),
      plot.margin = ggplot2::margin(t = 2, b = 2, l = 2, r = 2, unit = "pt")
    ) +
    annotate('text', x = -3.0, y = -0.8,
             label = paste0("Optimal lambda = ", lambda_min),
             color = 'black', size = 6) +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 2)))
  
  return(p)
}

# Helper function for cross-validation plot
create_cv_plot <- function(cvfit, lambda_min, prefix) {
  xx <- data.frame(
    lambda = cvfit[["lambda"]],
    cvm = cvfit[["cvm"]],
    cvsd = cvfit[["cvsd"]],
    cvup = cvfit[["cvup"]],
    cvlo = cvfit[["cvlo"]],
    nozezo = cvfit[["nzero"]]
  )
  xx$ll <- log(xx$lambda)
  xx$NZERO <- paste0(xx$nozezo, ' vars')
  
  p <- ggplot(xx, aes(ll, cvm, color = NZERO)) +
    geom_errorbar(aes(x = ll, ymin = cvlo, ymax = cvup),
                  width = 0.05, size = 1) +
    geom_vline(xintercept = xx$ll[which.min(xx$cvm)],
               size = 0.8, color = 'grey60', alpha = 0.8,
               linetype = 2) +
    geom_point(size = 2) +
    xlab(expression(paste(lambda, " (log scale)"))) +
    ylab('Partial Likelihood Deviance') +
    theme_bw(base_rect_size = 1.5) +
    scale_color_manual(values = c(pal_npg()(10), pal_d3()(10),
                                  pal_lancet()(10), pal_aaas()(10),
                                  pal_jco()(10))) +
    scale_x_continuous(expand = c(0.02, 0.02)) +
    scale_y_continuous(expand = c(0.02, 0.02)) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 15, color = 'black'),
      axis.text = element_text(size = 12, color = 'black'),
      legend.title = element_blank(),
      legend.text = element_text(size = 12, color = 'black'),
      legend.position = 'none'
    ) +
    annotate('text', x = -4.0, y = 1.8,
             label = paste0("Optimal lambda = ", lambda_min),
             color = 'black', size = 6) +
    guides(col = guide_legend(ncol = 6))
  
  return(p)
}

#### GSE139602 LASSO Analysis ####
# Load data
gse139602_up <- read.csv(file.path("02_Mfuzz", "HE_GSE139602_UP_gene.csv"), row.names = 1)
gse139602_down <- read.csv(file.path("02_Mfuzz", "HE_GSE139602_DOWN_gene.csv"), row.names = 1)
top_genes <- c(gse139602_up$x, gse139602_down$x)

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

# Run LASSO
result_gse139602 <- run_lasso_analysis(
  exp_data = exp_gse139602,
  gene_list = top_genes,
  group_data = group_gse139602,
  output_dir = output,
  prefix = "GSE139602_",
  family_type = "gaussian",
  lambda_label = 0.0386
)

#### GSE15654 LASSO Analysis (Cox) ####
# Load data
gse15654_genes <- read.csv(file.path("02_Mfuzz", "HE_GSE15654_ALL_gene_DEG.csv"))
top_genes_cox <- gse15654_genes$x

exp_gse15654 <- read.csv(file.path("00_rawdata", "00.rawdata_GSE15654_exp.csv"), row.names = 1)
group_gse15654 <- read.csv(file.path("00_rawdata", "00.rawdata_GSE15654_group.csv"), row.names = 1)
surv_gse15654 <- read.csv(file.path("00_rawdata", "00.rawdata_GSE15654_sur.csv"), row.names = 1)
surv_gse15654 <- surv_gse15654[colnames(exp_gse15654), , drop = FALSE]

# Run Cox LASSO
result_gse15654 <- run_lasso_analysis(
  exp_data = exp_gse15654,
  gene_list = top_genes_cox,
  surv_data = surv_gse15654,
  output_dir = output,
  prefix = "GSE15654_",
  family_type = "cox",
  lambda_label = 0.0159
)

message("LASSO analysis completed successfully!")

# Print summary
cat("\n=== Summary ===\n")
cat("GSE139602 selected genes (", length(result_gse139602$selected_genes), "):\n")
print(result_gse139602$selected_genes)
cat("\nGSE15654 selected genes (", length(result_gse15654$selected_genes), "):\n")
print(result_gse15654$selected_genes)

# Save combined results
combined_genes <- list(
  GSE139602 = result_gse139602$selected_genes,
  GSE15654 = result_gse15654$selected_genes
)
write.csv(data.frame(
  Dataset = rep(names(combined_genes), sapply(combined_genes, length)),
  Gene = unlist(combined_genes)
), file.path(output, "LASSO_selected_genes_summary.csv"), row.names = FALSE)
