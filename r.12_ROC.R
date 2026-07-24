rm(list = ls()); gc()
ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "12_ROC")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(ORIGINAL_DIR)

library(ggplot2)
library(pROC)
library(xgboost)
library(timeROC)
library(survival)

# Helper function for multi-class ROC
perform_multiclass_roc <- function(pred_data, group_col, score_col, output_dir, 
                                   GSE_name, labels = "", 
                                   colors = c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000","#7E6148","#B09C85")) {
  
  # Get unique groups
  groups <- unique(pred_data[[group_col]])
  
  # Generate all pairwise combinations
  group_combinations <- combn(groups, 2, simplify = FALSE)
  
  roc_results <- list()
  auc_values <- data.frame(Group1 = character(), Group2 = character(), AUC = numeric())
  
  for (comb in group_combinations) {
    group1 <- comb[1]
    group2 <- comb[2]
    
    data_subset <- pred_data %>% filter(.data[[group_col]] %in% c(group1, group2))
    data_subset$BinaryLabel <- ifelse(data_subset[[group_col]] == group1, 0, 1)
    
    roc_obj <- roc(data_subset$BinaryLabel, data_subset[[score_col]], 
                   levels = c(0, 1), direction = "<", quiet = TRUE)
    
    roc_results[[paste(group1, group2, sep = " vs ")]] <- roc_obj
    auc_values <- rbind(auc_values, data.frame(Group1 = group1, Group2 = group2, AUC = auc(roc_obj)))
  }
  
  # Plot ROC curves
  plot_colors <- adjustcolor(colors, alpha.f = 0.5)
  
  pdf(file.path(output_dir, paste0(labels, GSE_name, "_ROC_results.pdf")), width = 6, height = 6)
  plot(roc_results[[1]], col = plot_colors[1], main = "", lwd = 2, 
       cex.lab = 1.5, xlab = "1 - Specificity", ylab = "Sensitivity",
       xlim = c(0, 1), ylim = c(0, 1))
  
  for (i in 2:length(roc_results)) {
    plot(roc_results[[i]], col = plot_colors[i], add = TRUE, lwd = 2)
  }
  
  legend_labels <- sapply(names(roc_results), function(name) {
    auc_value <- auc(roc_results[[name]])
    paste(name, sprintf("(AUC = %.2f)", auc_value))
  })
  legend("bottomright", legend = legend_labels, col = plot_colors, lwd = 2)
  dev.off()
  
  png(file.path(output_dir, paste0(labels, GSE_name, "_ROC_results.png")), 
      width = 6, height = 6, units = "in", res = 300)
  plot(roc_results[[1]], col = plot_colors[1], main = "", lwd = 2, 
       cex.lab = 1.5, xlab = "1 - Specificity", ylab = "Sensitivity",
       xlim = c(0, 1), ylim = c(0, 1))
  
  for (i in 2:length(roc_results)) {
    plot(roc_results[[i]], col = plot_colors[i], add = TRUE, lwd = 2)
  }
  
  legend_labels <- sapply(names(roc_results), function(name) {
    auc_value <- auc(roc_results[[name]])
    paste(name, sprintf("(AUC = %.2f)", auc_value))
  })
  legend("bottomright", legend = legend_labels, col = plot_colors, lwd = 2)
  dev.off()
  
  write.csv(auc_values, file.path(output_dir, paste0(labels, GSE_name, "_ROC_results.csv")), row.names = FALSE)
  
  return(list(roc_results = roc_results, auc_values = auc_values))
}

# Helper function for time-dependent ROC
perform_timeroc <- function(pred_data, time_col, event_col, score_col, 
                            times = c(3, 5, 10), output_dir, GSE_name, labels = "",
                            colors = c("#E64B35","#4DBBD5","#00A087")) {
  
  roc_results <- list()
  auc_values <- c()
  
  for (t in times) {
    roc_obj <- timeROC(T = pred_data[[time_col]], 
                       delta = pred_data[[event_col]], 
                       marker = pred_data[[score_col]], 
                       cause = 1, times = t, iid = TRUE)
    roc_results[[paste0(t, "yr")]] <- roc_obj
    auc_values <- c(auc_values, roc_obj$AUC[2])
    cat(sprintf("%d-year AUC: %.2f\n", t, roc_obj$AUC[2]))
  }
  
  # Plot ROC curves
  plot_colors <- adjustcolor(colors, alpha.f = 0.85)
  
  pdf(file.path(output_dir, paste0(labels, GSE_name, "_ROC_results.pdf")), width = 6, height = 6)
  
  plot(roc_results[[1]]$FP[, 1], roc_results[[1]]$TP[, 2], 
       type = "l", col = plot_colors[1], lwd = 2,
       xlab = "1 - Specificity", ylab = "Sensitivity", cex.lab = 1.5,
       main = "", xlim = c(0, 1), ylim = c(0, 1))
  
  abline(a = 0, b = 1, lty = 2, col = "gray")
  
  for (i in 2:length(roc_results)) {
    lines(roc_results[[i]]$FP[, 1], roc_results[[i]]$TP[, 2], 
          col = plot_colors[i], lwd = 2)
  }
  
  legend_labels <- sapply(seq_along(times), function(i) {
    sprintf("%d years (AUC = %.2f)", times[i], auc_values[i])
  })
  legend("bottomright", legend = legend_labels, col = plot_colors, lwd = 2, cex = 1.2)
  dev.off()
  
  png(file.path(output_dir, paste0(labels, GSE_name, "_ROC_results.png")), 
      width = 6, height = 6, units = "in", res = 300)
  
  plot(roc_results[[1]]$FP[, 1], roc_results[[1]]$TP[, 2], 
       type = "l", col = plot_colors[1], lwd = 2,
       xlab = "1 - Specificity", ylab = "Sensitivity", cex.lab = 1.5,
       main = "", xlim = c(0, 1), ylim = c(0, 1))
  
  abline(a = 0, b = 1, lty = 2, col = "gray")
  
  for (i in 2:length(roc_results)) {
    lines(roc_results[[i]]$FP[, 1], roc_results[[i]]$TP[, 2], 
          col = plot_colors[i], lwd = 2)
  }
  
  legend_labels <- sapply(seq_along(times), function(i) {
    sprintf("%d years (AUC = %.2f)", times[i], auc_values[i])
  })
  legend("bottomright", legend = legend_labels, col = plot_colors, lwd = 2, cex = 1.2)
  dev.off()
  
  return(list(roc_results = roc_results, auc_values = auc_values))
}

# Helper function for model fusion
fuse_models <- function(exp_data, group_data, gene1, gene2, model1, model2,
                        output_dir, GSE_name, labels = "",
                        colors = c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F")) {
  
  # Prepare expression data
  exp_t <- as.data.frame(t(exp_data))
  
  # Predict with model 1
  df1 <- exp_t[, gene1[gene1 %in% colnames(exp_t)]]
  df1 <- scale(df1)
  group_data$risk_score_1 <- predict(model1, df1)
  
  # Predict with model 2
  df2 <- exp_t[, gene2[gene2 %in% colnames(exp_t)]]
  df2 <- scale(df2)
  group_data$risk_score_2 <- predict(model2, df2)
  
  # Calculate AUC for each model
  roc_1 <- roc(group_data$BinaryLabel, group_data$risk_score_1, 
               levels = c(0, 1), direction = "<", quiet = TRUE)
  roc_2 <- roc(group_data$BinaryLabel, group_data$risk_score_2, 
               levels = c(0, 1), direction = "<", quiet = TRUE)
  
  auc_1 <- as.numeric(auc(roc_1))
  auc_2 <- as.numeric(auc(roc_2))
  
  cat("\n=== Model Performance in", GSE_name, "===\n")
  cat("Model 1 (Progression) AUC:", round(auc_1, 4), "\n")
  cat("Model 2 (Prognostic) AUC:", round(auc_2, 4), "\n")
  
  # Calculate fusion weights based on AUC
  w1 <- auc_1 / (auc_1 + auc_2)
  w2 <- auc_2 / (auc_1 + auc_2)
  
  cat("\n=== Fusion Weights ===\n")
  cat("AUC-proportional: w1 =", round(w1, 4), ", w2 =", round(w2, 4), "\n")
  
  # Scale risk scores and compute combined score
  group_data$risk_score_1_scaled <- scale(group_data$risk_score_1)
  group_data$risk_score_2_scaled <- scale(group_data$risk_score_2)
  group_data$risk_score <- w1 * group_data$risk_score_1_scaled + 
                           w2 * group_data$risk_score_2_scaled
  
  # Evaluate combined score
  roc_combined <- roc(group_data$BinaryLabel, as.numeric(group_data$risk_score), 
                      levels = c(0, 1), direction = "<", quiet = TRUE)
  auc_combined <- as.numeric(auc(roc_combined))
  
  cat("\n=== Combined Score Performance ===\n")
  cat("Combined AUC:", round(auc_combined, 4), "\n")
  
  return(list(
    group_data = group_data,
    roc_1 = roc_1,
    roc_2 = roc_2,
    roc_combined = roc_combined,
    auc_1 = auc_1,
    auc_2 = auc_2,
    auc_combined = auc_combined,
    w1 = w1,
    w2 = w2
  ))
}

#### GSE139602 Multi-class ROC ####
GSE <- "GSE139602"

exp <- read.csv(file.path("00_rawdata", "00.rawdata_GSE139602_exp.csv"), row.names = 1)
group <- read.csv(file.path("00_rawdata", "00.rawdata_GSE139602_group.csv"), row.names = 1)
group <- group[colnames(exp), , drop = FALSE]

group$group <- ifelse(group$characteristics_ch1 == "disease state: Healthy", "Healthy",
               ifelse(group$characteristics_ch1 == "disease state: eCLD", "eCLD",
               ifelse(group$characteristics_ch1 == "disease state: Compensated Cirrhosis", "CC",
               ifelse(group$characteristics_ch1 == "disease state: Decompesated Cirrhosis", "DC",
               ifelse(group$characteristics_ch1 == "disease state: Acute-on-chronic liver failure", "ACLF", NA)))))

pred_data <- read.csv(file.path("07_xgboost", paste0("07_xgboost_", GSE, "_pred_data.csv")), row.names = 1)
group <- group[rownames(pred_data), ]
pred_data$group <- group$group

roc_results <- perform_multiclass_roc(
  pred_data = pred_data,
  group_col = "group",
  score_col = "pred",
  output_dir = output,
  GSE_name = GSE,
  labels = "12_ROC_"
)

#### GSE15654 Time-dependent ROC ####
GSE <- "GSE15654"

gene <- read.csv(file.path("06_RF", paste0("06_RF_", GSE, "_final_cogenes.csv")))$x
exp <- read.csv(file.path("00_rawdata", "00.rawdata_GSE15654_exp.csv"), row.names = 1)
group <- read.csv(file.path("00_rawdata", "00.rawdata_GSE15654_group.csv"), row.names = 1)
surv <- read.csv(file.path("00_rawdata", "00.rawdata_GSE15654_sur.csv"), row.names = 1)
surv <- surv[colnames(exp), , drop = FALSE]

final_model <- readRDS(file.path("07_xgboost", paste0("07_xgboost_", GSE, "_final_xgboost_model_cox.Rdata")))

exp <- exp[gene, ]
exp_data <- as.data.frame(t(exp))
exp_data <- scale(exp_data)
Xmat <- as.matrix(exp_data)

pred_train <- predict(final_model, Xmat)
group$pred <- pred_train
surv <- surv[rownames(group), ]
group$death <- surv$death
group$time <- surv$days_to_death / 365

timeroc_results <- perform_timeroc(
  pred_data = group,
  time_col = "time",
  event_col = "death",
  score_col = "pred",
  times = c(3, 5, 10),
  output_dir = output,
  GSE_name = GSE,
  labels = "12_ROC_"
)

#### GSE41919 Model Fusion ####
GSE <- "GSE41919"

# Load genes and models
gene1 <- read.csv(file.path("06_RF", "06_RF_GSE139602_final_cogenes.csv"))$x
gene2 <- read.csv(file.path("06_RF", "06_RF_GSE15654_final_cogenes.csv"))$x

model1 <- readRDS(file.path("07_xgboost", "07_xgboost_GSE139602_final_xgboost_model.Rdata"))
model2 <- readRDS(file.path("07_xgboost", "07_xgboost_GSE15654_final_xgboost_model_cox.Rdata"))

# Load GSE41919 data
exp_data <- read.csv('./00_rawdata/00.rawdata_GSE41919_exp.csv', header = TRUE, row.names = 1)
group_data <- read.csv('./00_rawdata/00.rawdata_GSE41919_group.csv', header = TRUE, row.names = 1)
group_data <- subset(group_data, group != "non-cirrhotic control")
exp_data <- exp_data[, rownames(group_data)]
exp_data <- as.data.frame(t(exp_data))

group_data$BinaryLabel <- ifelse(group_data$group == "cirrhosis without HE", 0, 1)

fusion_gse41919 <- fuse_models(
  exp_data = exp_data,
  group_data = group_data,
  gene1 = gene1,
  gene2 = gene2,
  model1 = model1,
  model2 = model2,
  output_dir = output,
  GSE_name = "GSE41919",
  labels = "12_ROC_"
)

#### GSE57193 Model Fusion ####
GSE <- "GSE57193"

# Load GSE57193 data
exp_data <- read.csv('./00_rawdata/00.rawdata_GSE57193_exp.csv', header = TRUE, row.names = 1)
group_data <- read.csv('./00_rawdata/00.rawdata_GSE57193_group.csv', header = TRUE, row.names = 1)
group_data <- subset(group_data, group != "healthy")
exp_data <- exp_data[, rownames(group_data)]
exp_data <- as.data.frame(t(exp_data))

group_data$BinaryLabel <- ifelse(group_data$group == "cirrhosis with HE", 1, 0)

fusion_gse57193 <- fuse_models(
  exp_data = exp_data,
  group_data = group_data,
  gene1 = gene1,
  gene2 = gene2,
  model1 = model1,
  model2 = model2,
  output_dir = output,
  GSE_name = "GSE57193",
  labels = "12_ROC_"
)

#### Combined ROC plot for GSE41919 and GSE57193 ####
plot_colors <- adjustcolor(c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", 
                              "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"), alpha.f = 0.5)

pdf(file.path(output, "12_ROC_GSE41919_GSE57193_combined_ROC.pdf"), width = 6, height = 6)
plot(fusion_gse41919$roc_combined, col = plot_colors[1], 
     main = "Cirrhosis without HE vs with HE", lwd = 2, 
     cex.lab = 1.5, xlab = "1 - Specificity", ylab = "Sensitivity",
     xlim = c(1, 0), ylim = c(0, 1))
plot(fusion_gse57193$roc_combined, col = plot_colors[2], add = TRUE, lwd = 2,
     xlim = c(1, 0), ylim = c(0, 1))
legend("bottomright", legend = c(
  sprintf("GSE41919 (AUC = %.2f)", fusion_gse41919$auc_combined),
  sprintf("GSE57193 (AUC = %.2f)", fusion_gse57193$auc_combined)
), col = plot_colors, lwd = 2, cex = 1.2)
dev.off()

png(file.path(output, "12_ROC_GSE41919_GSE57193_combined_ROC.png"), 
    width = 6, height = 6, units = "in", res = 300)
plot(fusion_gse41919$roc_combined, col = plot_colors[1], 
     main = "Cirrhosis without HE vs with HE", lwd = 2, 
     cex.lab = 1.5, xlab = "1 - Specificity", ylab = "Sensitivity",
     xlim = c(1, 0), ylim = c(0, 1))
plot(fusion_gse57193$roc_combined, col = plot_colors[2], add = TRUE, lwd = 2,
     xlim = c(1, 0), ylim = c(0, 1))
legend("bottomright", legend = c(
  sprintf("GSE41919 (AUC = %.2f)", fusion_gse41919$auc_combined),
  sprintf("GSE57193 (AUC = %.2f)", fusion_gse57193$auc_combined)
), col = plot_colors, lwd = 2, cex = 1.2)
dev.off()

message("ROC analysis completed successfully!")

# Print summary
cat("\n=== ROC Analysis Summary ===\n")
cat("GSE139602: Multi-class ROC completed\n")
cat("GSE15654: Time-dependent ROC completed\n")
cat("GSE41919: Model fusion AUC =", round(fusion_gse41919$auc_combined, 4), "\n")
cat("GSE57193: Model fusion AUC =", round(fusion_gse57193$auc_combined, 4), "\n")
cat("\nResults saved to:", output, "\n")
