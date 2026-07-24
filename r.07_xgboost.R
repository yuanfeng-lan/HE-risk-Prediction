rm(list = ls()); gc()
ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "07_xgboost")

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
library(survival)
library(survminer)

# Helper function for XGBoost regression
run_xgboost_regression <- function(exp_data, gene_list, response_data, 
                                   response_col, output_dir, prefix,
                                   seed = 123) {
  
  # Prepare data
  exp_filtered <- exp_data[gene_list, ]
  exp_data_t <- as.data.frame(t(exp_filtered))
  response_var <- response_data[, response_col, drop = FALSE]
  data <- cbind(exp_data_t, response_var)
  
  # Scale features
  X_mat <- as.matrix(data[, 1:(ncol(data) - 1)])
  X_mat <- scale(X_mat)
  y_vec <- as.numeric(data[, response_col])
  
  # Determine nfold
  n <- nrow(X_mat)
  nfold <- if (n <= 30) n else min(5, floor(n / 5))
  
  # Create DMatrix
  dtrain <- xgb.DMatrix(data = X_mat, label = y_vec)
  
  # Set parameters
  params <- list(
    objective = "reg:squarederror",
    eval_metric = "rmse",
    learning_rate = 0.05,
    max_depth = 3,
    min_child_weight = 5,
    subsample = 0.7,
    colsample_bytree = 0.7,
    lambda = 5,
    alpha = 0.1
  )
  
  # Cross-validation
  set.seed(seed)
  cv <- xgb.cv(
    params = params,
    data = dtrain,
    nrounds = 2000,
    nfold = nfold,
    early_stopping_rounds = 50,
    verbose = 1,
    showsd = TRUE,
    stratified = FALSE
  )
  
  # Get best nrounds
  elog <- cv$evaluation_log
  test_col <- grep("^test_.*_mean$", names(elog), value = TRUE)[1]
  best_nrounds <- if (!is.null(cv$best_iteration)) {
    cv$best_iteration
  } else if (!is.na(test_col)) {
    which.min(elog[[test_col]])
  } else {
    nrow(elog)
  }
  message("best_nrounds = ", best_nrounds)
  
  # Train final model
  final_model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = best_nrounds,
    watchlist = list(train = dtrain),
    verbose = 1
  )
  
  # Make predictions
  pred_train <- predict(final_model, X_mat)
  pred_train <- as.numeric(pred_train)
  y_true <- as.numeric(y_vec)
  
  # Calculate metrics
  train_rmse <- sqrt(mean((pred_train - y_true)^2, na.rm = TRUE))
  train_mae <- mean(abs(pred_train - y_true), na.rm = TRUE)
  train_r2 <- 1 - sum((pred_train - y_true)^2, na.rm = TRUE) / 
    sum((y_true - mean(y_true, na.rm = TRUE))^2, na.rm = TRUE)
  
  cat("Train RMSE:", train_rmse, "\n")
  cat("Train MAE :", train_mae, "\n")
  cat("Train R2  :", train_r2, "\n")
  
  # Save predictions
  response_data$pred <- pred_train
  write.csv(response_data, file.path(output_dir, paste0(prefix, "pred_data.csv")))
  
  # Create calibration plot
  create_calibration_plot(pred_train, y_true, train_r2, output_dir, prefix)
  
  # Save model
  saveRDS(final_model, file.path(output_dir, paste0(prefix, "final_xgboost_model.Rdata")))
  
  # SHAP analysis
  shap_results <- perform_shap_analysis(final_model, X_mat, gene_list, output_dir, prefix)
  
  return(list(
    model = final_model,
    pred = pred_train,
    y_true = y_true,
    rmse = train_rmse,
    mae = train_mae,
    r2 = train_r2,
    shap = shap_results
  ))
}

# Helper function for XGBoost Cox regression
run_xgboost_cox <- function(exp_data, gene_list, surv_data, 
                            time_col, event_col, output_dir, prefix,
                            seed = 123) {
  
  # Prepare data
  exp_filtered <- exp_data[gene_list, ]
  exp_data_t <- as.data.frame(t(exp_filtered))
  exp_scaled <- scale(exp_data_t)
  X_mat <- as.matrix(exp_scaled)
  
  time_train <- surv_data[, time_col]
  status_train <- surv_data[, event_col]
  
  # Create DMatrix
  dtrain <- xgb.DMatrix(data = X_mat, label = as.numeric(time_train))
  
  # Set parameters for Cox regression
  params <- list(
    objective = "survival:cox",
    eval_metric = "cox-nloglik",
    eta = 0.05,
    max_depth = 6,
    subsample = 0.9,
    colsample_bytree = 1,
    colsample_bylevel = 1,
    colsample_bynode = 1,
    min_child_weight = 1,
    gamma = 0,
    lambda = 0.5,
    alpha = 0
  )
  
  # Cross-validation
  set.seed(seed)
  cv <- xgb.cv(
    params = params,
    data = dtrain,
    nrounds = 2000,
    nfold = 5,
    early_stopping_rounds = 50,
    verbose = 1,
    showsd = TRUE
  )
  
  # Get best nrounds
  elog <- cv$evaluation_log
  test_col <- grep("^test_.*_mean$", names(elog), value = TRUE)[1]
  best_nrounds <- if (!is.null(cv$best_iteration)) {
    cv$best_iteration
  } else if (!is.na(test_col)) {
    which.min(elog[[test_col]])
  } else {
    nrow(elog)
  }
  message("best_nrounds = ", best_nrounds)
  
  # Train final model
  final_model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = best_nrounds,
    watchlist = list(train = dtrain),
    verbose = 1
  )
  
  # Make predictions
  pred_train <- predict(final_model, X_mat)
  
  # Save model and predictions
  saveRDS(final_model, file.path(output_dir, paste0(prefix, "final_xgboost_model_cox.Rdata")))
  
  # Create survival analysis
  surv_results <- create_survival_analysis(pred_train, time_train, status_train, 
                                           output_dir, prefix)
  
  # SHAP analysis
  shap_results <- perform_shap_analysis(final_model, X_mat, gene_list, output_dir, prefix)
  
  return(list(
    model = final_model,
    pred = pred_train,
    survival = surv_results,
    shap = shap_results
  ))
}

# Helper function for calibration plot
create_calibration_plot <- function(pred_train, y_true, train_r2, output_dir, prefix) {
  df_cal <- data.frame(pred = as.numeric(pred_train), obs = as.numeric(y_true))
  bins <- min(10, nrow(df_cal))
  df_cal <- df_cal %>% mutate(bin = ntile(pred, bins))
  
  calib_summary <- df_cal %>%
    group_by(bin) %>%
    summarise(
      mean_pred = mean(pred, na.rm = TRUE),
      mean_obs = mean(obs, na.rm = TRUE),
      n = n(),
      sd_obs = sd(obs, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      se = ifelse(n > 1, sd_obs / sqrt(n), 0),
      ci_lo = mean_obs - 1.96 * se,
      ci_hi = mean_obs + 1.96 * se
    )
  
  p_cal <- ggplot(df_cal, aes(x = pred, y = obs)) +
    geom_point(alpha = 0.4, size = 3, color = "#4DBBD5") +
    geom_point(data = calib_summary,
               mapping = aes(x = mean_pred, y = mean_obs),
               inherit.aes = FALSE,
               color = "#E64B35", size = 4.5) +
    geom_errorbar(data = calib_summary,
                  mapping = aes(x = mean_pred, ymin = ci_lo, ymax = ci_hi),
                  inherit.aes = FALSE,
                  width = 0.02 * diff(range(df_cal$pred, na.rm = TRUE)),
                  color = "#E64B35", alpha = 0.8) +
    geom_smooth(method = "loess", se = TRUE, color = "#F39B7F", 
                fill = alpha("#F39B7F", 0.2)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#3C5488") +
    labs(x = "Predicted value", y = "Observed value", title = "") +
    annotate("text", x = Inf, y = Inf, 
             label = paste0("R² = ", formatC(train_r2, format = "f", digits = 3)),
             hjust = 1.1, vjust = 1.5, size = 8) +
    theme_minimal(base_size = 12) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.background = element_rect(fill = "white"),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text = element_text(size = 12)
    ) +
    coord_cartesian(expand = FALSE)
  
  ggsave(file.path(output_dir, paste0(prefix, "calibration.png")), 
         p_cal, width = 6, height = 6, dpi = 300)
  ggsave(file.path(output_dir, paste0(prefix, "calibration.pdf")), 
         p_cal, width = 6, height = 6)
}

# Helper function for SHAP analysis
perform_shap_analysis <- function(model, X_mat, gene_list, output_dir, prefix, top_n = 12) {
  # Calculate SHAP values
  shap_contrib <- predict(model, X_mat, predcontrib = TRUE)
  shap_values <- shap_contrib[, -ncol(shap_contrib), drop = FALSE]
  baseline_val <- mean(shap_contrib[, ncol(shap_contrib)])
  
  # Create shapviz object
  sv <- shapviz(object = model, x = X_mat, X_pred = X_mat)
  
  # Save SHAP matrix
  shap_mat <- sv$S
  X_df <- sv$X
  colnames(shap_mat) <- colnames(X_df)
  rownames(shap_mat) <- rownames(X_df)
  write.csv(shap_mat, file.path(output_dir, paste0(prefix, "shap_mat.csv")))
  
  # Create importance bar plot
  mean_abs <- colMeans(abs(shap_mat), na.rm = TRUE)
  ord <- order(mean_abs, decreasing = TRUE)
  top_features <- colnames(shap_mat)[ord][seq_len(min(top_n, length(colnames(shap_mat))))]
  mean_abs_top <- mean_abs[top_features]
  
  df_bar <- data.frame(
    feature = factor(top_features, levels = rev(top_features)),
    mean_abs = mean_abs_top
  )
  
  p_bar <- ggplot(df_bar, aes(x = mean_abs, y = feature, fill = mean_abs)) +
    geom_col(width = 0.7) +
    scale_fill_gradientn(colors = c("#4DBBD5", "#F1C40F", "#E64B35"), 
                         guide = guide_colorbar(title = "mean |SHAP|")) +
    labs(x = "Mean |SHAP value|", y = NULL, title = "Features - SHAP importance") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text = element_text(size = 12),
      axis.text.y = element_text(face = "italic"),
      legend.text = element_text(size = 12),
      legend.title = element_blank()
    ) +
    geom_text(aes(label = sprintf("%.3f", mean_abs)), hjust = -0.1, size = 5) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2)))
  
  ggsave(file.path(output_dir, paste0(prefix, "shap_bar.png")), 
         p_bar, width = 6, height = 4, dpi = 300)
  ggsave(file.path(output_dir, paste0(prefix, "shap_bar.pdf")), 
         p_bar, width = 6, height = 4)
  
  # Create beeswarm plot
  shap_df <- as.data.frame(shap_mat)
  shap_df$sample <- rownames(shap_df)
  long_shap <- shap_df %>%
    pivot_longer(cols = -sample, names_to = "feature", values_to = "shap")
  
  X_long <- X_df %>%
    mutate(sample = rownames(.)) %>%
    pivot_longer(cols = -sample, names_to = "feature", values_to = "feat_value")
  
  plot_df <- left_join(long_shap, X_long, by = c("sample", "feature")) %>%
    filter(feature %in% top_features) %>%
    mutate(feature = factor(feature, levels = rev(top_features)))
  
  p_beeswarm <- ggplot(plot_df, aes(x = shap, y = feature, color = feat_value)) +
    ggbeeswarm::geom_quasirandom(groupOnX = FALSE, size = 2, alpha = 0.8, width = 0.3) +
    scale_color_gradient(low = "#91D1C2", high = "#DC0000", name = "Feature value") +
    geom_vline(xintercept = 0, color = "grey40", linetype = "dashed") +
    labs(x = "SHAP value", y = NULL, title = "SHAP beeswarm") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text = element_text(size = 12),
      axis.text.y = element_text(face = "italic"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14)
    )
  
  ggsave(file.path(output_dir, paste0(prefix, "shap_beeswarm.png")), 
         p_beeswarm, width = 6, height = 4, dpi = 300)
  ggsave(file.path(output_dir, paste0(prefix, "shap_beeswarm.pdf")), 
         p_beeswarm, width = 6, height = 4)
  
  # Create waterfall plots for specific samples
  create_waterfall_plots(sv, output_dir, prefix)
  
  return(list(
    shap_mat = shap_mat,
    sv = sv,
    top_features = top_features,
    importance = mean_abs
  ))
}

# Helper function for waterfall plots
create_waterfall_plots <- function(sv, output_dir, prefix, samples = c(2, 38, 96)) {
  for (sample_id in samples) {
    p <- sv_waterfall(sv, row_id = sample_id, fill_colors = c("#8491B4", "#F39B7F")) +
      theme(
        axis.text = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16)
      )
    
    p_force <- sv_force(sv, row_id = sample_id, fill_colors = c("#8491B4", "#F39B7F")) +
      theme(
        axis.text = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)
      )
    
    ggsave(file.path(output_dir, paste0(prefix, "waterfall_sample", sample_id, ".png")), 
           p, width = 5.4, height = 4, dpi = 300)
    ggsave(file.path(output_dir, paste0(prefix, "waterfall_sample", sample_id, ".pdf")), 
           p, width = 5.4, height = 4)
    ggsave(file.path(output_dir, paste0(prefix, "waterfall_sample", sample_id, "_force.png")), 
           p_force, width = 4, height = 2, dpi = 300)
    ggsave(file.path(output_dir, paste0(prefix, "waterfall_sample", sample_id, "_force.pdf")), 
           p_force, width = 4, height = 2)
  }
}

# Helper function for survival analysis
create_survival_analysis <- function(pred_train, time_train, status_train, output_dir, prefix) {
  # Calculate concordance index
  cf <- survival::concordancefit(Surv(time_train, status_train), pred_train)
  cidx_train <- cf$concordance
  print(paste("Concordance index:", cidx_train))
  
  # Split into high and low risk groups
  group_train <- ifelse(pred_train > median(pred_train, na.rm = TRUE), "High", "Low")
  df_plot <- data.frame(
    time = time_train / 365,
    status = status_train,
    group = factor(group_train, levels = c("Low", "High"))
  )
  
  # Fit survival curve
  fit <- survfit(Surv(time, status) ~ group, data = df_plot)
  
  # Cox regression
  cox_uni <- coxph(Surv(time, status) ~ group, data = df_plot)
  summary(cox_uni)
  
  hr <- summary(cox_uni)$coefficients[, "exp(coef)"]
  ci_low <- summary(cox_uni)$conf.int[, "lower .95"]
  ci_high <- summary(cox_uni)$conf.int[, "upper .95"]
  pval_cox <- summary(cox_uni)$coefficients[, "Pr(>|z|)"]
  p_txt <- if (is.na(pval_cox)) "NA" else format(signif(pval_cox, 3), scientific = TRUE)
  
  # Create annotation
  p_note_expr <- sprintf(
    "atop('HR=%.2f (95%% CI %.2f-%.2f)', italic(p)==%s)",
    hr, ci_low, ci_high, p_txt
  )
  
  # Create KM plot
  p_km <- ggsurvplot(
    fit,
    data = df_plot,
    risk.table = TRUE,
    pval = FALSE,
    conf.int = FALSE,
    surv.median.line = "hv",
    palette = c("#3C5488", "#DC0000"),
    xlab = "Time (years)",
    legend.title = "Risk Score",
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
    )
  
  ggsave(file.path(output_dir, paste0(prefix, "KM_plot.png")), 
         print(p_km), width = 5, height = 5, dpi = 300)
  ggsave(file.path(output_dir, paste0(prefix, "KM_plot.pdf")), 
         print(p_km), width = 5, height = 5)
  
  return(list(
    concordance = cidx_train,
    cox = cox_uni,
    km_plot = p_km
  ))
}

#### GSE139602 XGBoost Regression ####
# Load data
gse139602_genes <- read.csv(file.path("06_RF", "GSE139602_final_cogenes.csv"))
gene_list_gse139602 <- gse139602_genes$x

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

# Run XGBoost regression
result_gse139602 <- run_xgboost_regression(
  exp_data = exp_gse139602,
  gene_list = gene_list_gse139602,
  response_data = group_gse139602,
  response_col = "group",
  output_dir = output,
  prefix = "GSE139602_",
  seed = 123
)

#### GSE15654 XGBoost Cox Regression ####
# Load data
gse15654_genes <- read.csv(file.path("06_RF", "GSE15654_final_cogenes.csv"))
gene_list_gse15654 <- gse15654_genes$x

exp_gse15654 <- read.csv(file.path("00_rawdata", "00.rawdata_GSE15654_exp.csv"), row.names = 1)
surv_gse15654 <- read.csv(file.path("00_rawdata", "00.rawdata_GSE15654_sur.csv"), row.names = 1)

# Match samples
surv_gse15654 <- surv_gse15654[colnames(exp_gse15654), , drop = FALSE]

# Run XGBoost Cox regression
result_gse15654 <- run_xgboost_cox(
  exp_data = exp_gse15654,
  gene_list = gene_list_gse15654,
  surv_data = surv_gse15654,
  time_col = "days_to_death",
  event_col = "death",
  output_dir = output,
  prefix = "GSE15654_",
  seed = 123
)

message("XGBoost analysis completed successfully!")

# Print summary
cat("\n=== Summary ===\n")
cat("GSE139602:\n")
cat("  RMSE:", result_gse139602$rmse, "\n")
cat("  MAE:", result_gse139602$mae, "\n")
cat("  R²:", result_gse139602$r2, "\n")

cat("\nGSE15654:\n")
cat("  Concordance index:", result_gse15654$survival$concordance, "\n")
