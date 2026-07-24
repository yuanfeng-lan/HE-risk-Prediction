rm(list = ls()); gc()
ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "08_model_validation")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(ORIGINAL_DIR)

library(shapviz)
library(xgboost)
library(ggplot2)
library(tidyverse)
library(caret)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(viridis)
library(cowplot)
library(scales)
library(ggbeeswarm)
library(patchwork)
library(ggplotify)
library(grid)
library(gridExtra)
library(ggpubr)
library(survival)
library(survminer)

# Helper function for boxplot with p-value
create_boxplot <- function(data, x_var, y_var, fill_var, color_var, 
                           x_labels = NULL, title = "", y_label = "", 
                           p_value = NULL, output_dir, filename,
                           colors = NULL, width = 4, height = 4) {
  
  p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], 
                        fill = .data[[fill_var]], colour = .data[[color_var]])) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.5, 
                 colour = "black", size = 0.4) +
    geom_jitter(position = position_jitter(width = 0.15, height = 0), 
                size = 1.8, alpha = 0.9) +
    theme_minimal(base_size = 15) +
    labs(x = "", y = y_label, title = title) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 14, face = "bold", angle = 45, 
                                 vjust = 1, hjust = 1, margin = margin(t = 2)),
      axis.text.y = element_text(size = 14),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    ) +
    coord_cartesian(clip = "off")
  
  # Add custom x labels if provided
  if (!is.null(x_labels)) {
    p <- p + scale_x_discrete(labels = x_labels)
  }
  
  # Add custom colors if provided
  if (!is.null(colors)) {
    p <- p + scale_fill_manual(values = colors) +
      scale_colour_manual(values = colors)
  }
  
  # Add p-value annotation if provided
  if (!is.null(p_value)) {
    y_max <- max(data[[y_var]], na.rm = TRUE)
    p <- p + annotation_custom(
      grob = textGrob(
        label = bquote(italic(p) ~ "=" ~ .(p_value)),
        gp = gpar(fontsize = 15),
        hjust = 0.5
      ),
      xmin = 0.2,
      xmax = 1.5,
      ymin = y_max * 0.85,
      ymax = y_max * 1.05
    )
  }
  
  ggsave(file.path(output_dir, paste0(filename, ".png")), p, width = width, height = height, dpi = 300)
  ggsave(file.path(output_dir, paste0(filename, ".pdf")), p, width = width, height = height)
  
  return(p)
}

# Helper function for gene expression boxplot
create_gene_boxplot <- function(data, gene, group_col, output_dir, 
                                colors = NULL, title = NULL, 
                                width = 4, height = 4) {
  
  if (is.null(colors)) {
    colors <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F")
    names(colors) <- unique(data[[group_col]])
  }
  
  p_value <- compare_means(as.formula(paste(gene, "~", group_col)), 
                           data = data, method = "kruskal.test")$p.format
  
  p <- ggplot(data, aes(x = .data[[group_col]], y = .data[[gene]], 
                        fill = .data[[group_col]], colour = .data[[group_col]])) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.5, 
                 colour = "black", size = 0.4) +
    geom_jitter(position = position_jitter(width = 0.15, height = 0), 
                size = 1.8, alpha = 0.9) +
    scale_fill_manual(values = colors) +
    scale_colour_manual(values = colors) +
    theme_minimal(base_size = 15) +
    labs(x = "", y = "Relative gene expression", title = ifelse(is.null(title), gene, title)) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 14, face = "bold", angle = 45, 
                                 vjust = 1, hjust = 1, margin = margin(t = 2)),
      axis.text.y = element_text(size = 14),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 16, face = "italic")
    ) +
    coord_cartesian(clip = "off") +
    annotation_custom(
      grob = textGrob(
        label = bquote(italic(p) ~ "=" ~ .(p_value)),
        gp = gpar(fontsize = 15),
        hjust = 0.5
      ),
      xmin = 1.1,
      xmax = 1.5,
      ymin = max(data[[gene]], na.rm = TRUE) * 0.85,
      ymax = max(data[[gene]], na.rm = TRUE) * 1.05
    )
  
  ggsave(file.path(output_dir, paste0(gene, "_plot.png")), p, width = width, height = height, dpi = 300)
  ggsave(file.path(output_dir, paste0(gene, "_plot.pdf")), p, width = width, height = height)
  
  return(p)
}

# Helper function for KM survival plot
create_km_plot <- function(data, time_col, status_col, group_col, 
                           output_dir, filename, title = NULL,
                           xlab = "Time (years)", ylab = "Survival probability",
                           colors = c("#3C5488", "#DC0000"),
                           width = 5, height = 5) {
  
  fit <- survfit(as.formula(paste0("Surv(", time_col, ", ", status_col, ") ~ ", group_col)), 
                 data = data)
  
  cox_uni <- coxph(as.formula(paste0("Surv(", time_col, ", ", status_col, ") ~ ", group_col)), 
                   data = data)
  
  hr <- summary(cox_uni)$coefficients[, "exp(coef)"]
  ci_low <- summary(cox_uni)$conf.int[, "lower .95"]
  ci_high <- summary(cox_uni)$conf.int[, "upper .95"]
  pval_cox <- summary(cox_uni)$coefficients[, "Pr(>|z|)"]
  p_txt <- if (is.na(pval_cox)) "NA" else format(signif(pval_cox, 3), scientific = TRUE)
  
  p_note_expr <- sprintf(
    "atop('HR=%.2f (95%% CI %.2f-%.2f)', italic(p)==%s)",
    hr, ci_low, ci_high, p_txt
  )
  
  p_km <- ggsurvplot(
    fit,
    data = data,
    risk.table = TRUE,
    pval = FALSE,
    conf.int = FALSE,
    surv.median.line = "hv",
    palette = colors,
    xlab = xlab,
    ylab = ylab,
    legend.title = ifelse(is.null(title), group_col, title),
    legend.labs = c("Low", "High"),
    risk.table.height = 0.25,
    ggtheme = theme_minimal(base_size = 13)
  )
  
  p_km$plot <- p_km$plot +
    annotate(
      "text",
      x = 10, y = 0.05,
      label = p_note_expr,
      parse = TRUE,
      hjust = 1.02, vjust = 0,
      size = 4
    ) +
    labs(color = ifelse(is.null(title), group_col, title), 
         fill = ifelse(is.null(title), group_col, title)) +
    theme(legend.title = element_text(face = "italic"))
  
  ggsave(file.path(output_dir, paste0(filename, ".png")), 
         print(p_km), width = width, height = height, dpi = 300)
  ggsave(file.path(output_dir, paste0(filename, ".pdf")), 
         print(p_km), width = width, height = height)
  
  return(list(plot = p_km, cox = cox_uni))
}

#### GSE139602 Validation ####
GSE <- "GSE139602_"

# Load data
gene <- read.csv(file.path("06_RF", paste0("06_RF_", GSE, "final_cogenes.csv")))$x

exp <- read.csv(file.path("00_rawdata", "00.rawdata_GSE139602_exp.csv"), row.names = 1)
group <- read.csv(file.path("00_rawdata", "00.rawdata_GSE139602_group.csv"), row.names = 1)
group <- group[colnames(exp), , drop = FALSE]

# Convert group labels
group$group <- ifelse(group$characteristics_ch1 == "disease state: Healthy", 0,
               ifelse(group$characteristics_ch1 == "disease state: eCLD", 1,
               ifelse(group$characteristics_ch1 == "disease state: Compensated Cirrhosis", 2,
               ifelse(group$characteristics_ch1 == "disease state: Decompesated Cirrhosis", 3,
               ifelse(group$characteristics_ch1 == "disease state: Acute-on-chronic liver failure", 4, NA)))))

# Prepare expression data
exp_filtered <- exp[gene, ]
exp_data <- as.data.frame(t(exp_filtered))
exp_data <- as.data.frame(scale(exp_data))
group_data <- group[, 2, drop = FALSE]
data <- cbind(exp_data, group_data)
X_mat <- as.matrix(data[, 1:(ncol(data) - 1)])

# Load model and make predictions
final_model <- readRDS(file.path("07_xgboost", paste0("07_xgboost_", GSE, "final_xgboost_model.Rdata")))
pred_train <- predict(final_model, X_mat)
data$pred <- pred_train

# Add group labels
data$group_label <- ifelse(group$characteristics_ch1 == "disease state: Healthy", "Healthy",
                    ifelse(group$characteristics_ch1 == "disease state: eCLD", "eCLD",
                    ifelse(group$characteristics_ch1 == "disease state: Compensated Cirrhosis", "CC",
                    ifelse(group$characteristics_ch1 == "disease state: Decompesated Cirrhosis", "DC",
                    ifelse(group$characteristics_ch1 == "disease state: Acute-on-chronic liver failure", "ACLF", NA)))))

data$group_label <- factor(data$group_label, levels = c("Healthy", "eCLD", "CC", "DC", "ACLF"))

# Create output directory
outdir <- file.path(output, "plot_01_GSE139602")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Plot individual genes
for (gene_name in gene) {
  create_gene_boxplot(data, gene_name, "group_label", outdir)
}

# Plot risk score across groups
colors <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F")
names(colors) <- levels(data$group_label)

p_value <- compare_means(pred ~ group_label, data = data, method = "kruskal.test")$p.format

p <- create_boxplot(
  data = data,
  x_var = "group_label",
  y_var = "pred",
  fill_var = "group_label",
  color_var = "group_label",
  title = "",
  y_label = "Risk Score",
  p_value = p_value,
  output_dir = outdir,
  filename = paste0("08_model_validation_", GSE, "risk_score_plot"),
  colors = colors,
  width = 4,
  height = 4
)

#### Validate on GSE41919 ####
# Load data
gse41919_exp <- read.csv('./00_rawdata/00.rawdata_GSE41919_exp.csv', header = TRUE, row.names = 1)
gse41919_group <- read.csv('./00_rawdata/00.rawdata_GSE41919_group.csv', header = TRUE, row.names = 1)
gse41919_group <- subset(gse41919_group, group != "non-cirrhotic control")
gse41919_exp <- gse41919_exp[, rownames(gse41919_group)]
gse41919_exp <- as.data.frame(t(gse41919_exp))

# Predict risk scores
gse41919_df <- gse41919_exp[, gene]
gse41919_df <- scale(gse41919_df)
gse41919_risk_score <- predict(final_model, gse41919_df)
gse41919_group$risk_score <- gse41919_risk_score
gse41919_group$group <- factor(gse41919_group$group, levels = c("cirrhosis without HE", "cirrhosis with HE"))

# Create risk score plot
p_value <- compare_means(risk_score ~ group, data = gse41919_group, method = "wilcox.test")$p.format
p <- create_boxplot(
  data = gse41919_group,
  x_var = "group",
  y_var = "risk_score",
  fill_var = "group",
  color_var = "group",
  x_labels = c("Cirrhosis\nwithout HE", "Cirrhosis\nwith HE"),
  title = "GSE41919",
  y_label = "Risk Score from Model",
  p_value = p_value,
  output_dir = outdir,
  filename = paste0("08_model_validation_", GSE, "Risk_score_GSE41919_plot"),
  colors = c("#4DBBD5", "#E64B35"),
  width = 4,
  height = 4
)

write.csv(gse41919_group, file.path(outdir, paste0("08_model_validation_", GSE, "Risk_score_GSE41919_data.csv")))

#### Validate on GSE57193 ####
gse57193_exp <- read.csv('./00_rawdata/00.rawdata_GSE57193_exp.csv', header = TRUE, row.names = 1)
gse57193_group <- read.csv('./00_rawdata/00.rawdata_GSE57193_group.csv', header = TRUE, row.names = 1)
gse57193_group <- subset(gse57193_group, group != "healthy")
gse57193_exp <- gse57193_exp[, rownames(gse57193_group)]
gse57193_exp <- as.data.frame(t(gse57193_exp))

# Predict risk scores
gse57193_df <- gse57193_exp[, gene]
gse57193_df <- scale(gse57193_df)
gse57193_risk_score <- predict(final_model, gse57193_df)
gse57193_group$risk_score <- gse57193_risk_score
gse57193_group$group <- ifelse(gse57193_group$group == "cirrhosis", "cirrhosis without HE", "cirrhosis with HE")
gse57193_group$group <- factor(gse57193_group$group, levels = c("cirrhosis without HE", "cirrhosis with HE"))

# Create risk score plot
p_value <- compare_means(risk_score ~ group, data = gse57193_group, method = "wilcox.test")$p.format
p <- create_boxplot(
  data = gse57193_group,
  x_var = "group",
  y_var = "risk_score",
  fill_var = "group",
  color_var = "group",
  x_labels = c("Cirrhosis\nwithout HE", "Cirrhosis\nwith HE"),
  title = "GSE57193",
  y_label = "Risk Score from Model",
  p_value = p_value,
  output_dir = outdir,
  filename = paste0("08_model_validation_", GSE, "Risk_score_GSE57193_plot"),
  colors = c("#4DBBD5", "#E64B35"),
  width = 4,
  height = 4
)

write.csv(gse57193_group, file.path(outdir, paste0("08_model_validation_", GSE, "Risk_score_GSE57193_data.csv")))

#### GSE15654 Validation ####
GSE <- "GSE15654_"

# Load data
gene <- read.csv(file.path("06_RF", "06_RF_GSE15654_final_cogenes.csv"))$x

exp <- read.csv(file.path("00_rawdata", "00.rawdata_GSE15654_exp.csv"), row.names = 1)
surv <- read.csv(file.path("00_rawdata", "00.rawdata_GSE15654_sur.csv"), row.names = 1)
surv <- surv[colnames(exp), , drop = FALSE]

# Prepare expression data
exp_filtered <- exp[gene, ]
exp_data <- as.data.frame(t(exp_filtered))
exp_data <- scale(exp_data)
Xmat <- as.matrix(exp_data)

# Load model and make predictions
final_model <- readRDS(file.path("07_xgboost", "07_xgboost_GSE15654_final_xgboost_model_cox.Rdata"))
pred_train <- predict(final_model, Xmat)

# Convert time to years
time_train <- surv$days_to_death / 365
status_train <- surv$death
exp_data <- as.data.frame(exp_data)

# Create KM plots for each gene
for (gene_name in gene) {
  data_km <- data.frame(
    time = time_train,
    status = status_train,
    group = ifelse(exp_data[[gene_name]] > median(exp_data[[gene_name]], na.rm = TRUE), "High", "Low")
  )
  data_km$group <- factor(data_km$group, levels = c("Low", "High"))
  
  create_km_plot(
    data = data_km,
    time_col = "time",
    status_col = "status",
    group_col = "group",
    output_dir = output,
    filename = paste0("08_model_validation_", GSE, gene_name, "_KM_plot"),
    title = gene_name,
    ylab = "Survival probability (OS)"
  )
}

# Decompensation-free survival
time_train <- surv$days_to_decomp / 365
status_train <- surv$decomp

for (gene_name in gene) {
  data_km <- data.frame(
    time = time_train,
    status = status_train,
    group = ifelse(exp_data[[gene_name]] > median(exp_data[[gene_name]], na.rm = TRUE), "High", "Low")
  )
  data_km$group <- factor(data_km$group, levels = c("Low", "High"))
  
  create_km_plot(
    data = data_km,
    time_col = "time",
    status_col = "status",
    group_col = "group",
    output_dir = output,
    filename = paste0("08_model_validation_", GSE, gene_name, "_Decompensation_free_survival"),
    title = gene_name,
    ylab = "Cumulative survival free of decompensation"
  )
}

# Child endpoint-free survival
time_train <- surv$days_to_child / 365
status_train <- surv$child

for (gene_name in gene) {
  data_km <- data.frame(
    time = time_train,
    status = status_train,
    group = ifelse(exp_data[[gene_name]] > median(exp_data[[gene_name]], na.rm = TRUE), "High", "Low")
  )
  data_km$group <- factor(data_km$group, levels = c("Low", "High"))
  
  create_km_plot(
    data = data_km,
    time_col = "time",
    status_col = "status",
    group_col = "group",
    output_dir = output,
    filename = paste0("08_model_validation_", GSE, gene_name, "_Child_free_survival"),
    title = gene_name,
    ylab = "Cumulative survival free of Child endpoint"
  )
}

# HCC progression-free survival
time_train <- surv$days_to_hcc / 365
status_train <- surv$hcc

for (gene_name in gene) {
  data_km <- data.frame(
    time = time_train,
    status = status_train,
    group = ifelse(exp_data[[gene_name]] > median(exp_data[[gene_name]], na.rm = TRUE), "High", "Low")
  )
  data_km$group <- factor(data_km$group, levels = c("Low", "High"))
  
  create_km_plot(
    data = data_km,
    time_col = "time",
    status_col = "status",
    group_col = "group",
    output_dir = output,
    filename = paste0("08_model_validation_", GSE, gene_name, "_HCC_free_survival"),
    title = gene_name,
    ylab = "Cumulative survival free of HCC progression"
  )
}

# Risk score KM plots
risk_score <- pred_train
outdir <- file.path(output, GSE)

for (event_type in c("decomp", "child", "hcc")) {
  time_col <- switch(event_type,
                     "decomp" = "days_to_decomp",
                     "child" = "days_to_child",
                     "hcc" = "days_to_hcc")
  
  time_train <- surv[[time_col]] / 365
  status_train <- surv[[event_type]]
  
  data_km <- data.frame(
    time = time_train,
    status = status_train,
    group = ifelse(risk_score > median(risk_score, na.rm = TRUE), "High", "Low")
  )
  data_km$group <- factor(data_km$group, levels = c("Low", "High"))
  
  ylab_text <- switch(event_type,
                      "decomp" = "Cumulative survival free of decompensation",
                      "child" = "Cumulative survival free of Child endpoint",
                      "hcc" = "Cumulative survival free of HCC progression")
  
  create_km_plot(
    data = data_km,
    time_col = "time",
    status_col = "status",
    group_col = "group",
    output_dir = output,
    filename = paste0("08_model_validation_", GSE, event_type, "_free_survival"),
    title = "Risk Score",
    ylab = ylab_text
  )
}

message("Model validation completed successfully!")
