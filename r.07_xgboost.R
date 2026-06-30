rm(list = ls()); gc()
ORIGINAL_DIR <- "~/data_hdd/R/Yclub/2026.3.4_HE_Cirrhosis"
output <- file.path(ORIGINAL_DIR, "07_xgboost")
labels <- "07_xgboost_"

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
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(cowplot)     
library(scales)
library(ggbeeswarm)
library(patchwork)
library(ggplotify)
library(grid)
library(gridExtra)
library(survival); library(survminer)

####GSE139602####
GSE <- "GSE139602_"
gene = read.csv(file.path(ORIGINAL_DIR,"06_RF",paste0("06_RF_",GSE,"final_cogenes.csv"))) 
gene <- gene$x
HE_GSE139602_UP_gene <- read.csv(file.path(ORIGINAL_DIR,"02_Mfuzz","HE_GSE139602_UP_gene.csv"),row.names = 1)
HE_GSE139602_DOWN_gene <- read.csv(file.path(ORIGINAL_DIR,"02_Mfuzz","HE_GSE139602_DOWN_gene.csv"),row.names = 1)
exp <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE139602_exp.csv"),row.names = 1)
group <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE139602_group.csv"),row.names = 1)
group <- group[colnames(exp),,drop=F]
table(group$characteristics_ch1)
group$group <- ifelse(group$characteristics_ch1=="disease state: Healthy",0,
                      ifelse(group$characteristics_ch1=="disease state: eCLD",1,
                             ifelse(group$characteristics_ch1=="disease state: Compensated Cirrhosis",2,
                                    ifelse(group$characteristics_ch1=="disease state: Decompesated Cirrhosis",3,
                                           ifelse(group$characteristics_ch1=="disease state: Acute-on-chronic liver failure",4,NA)))))
exp <- exp[gene,]
exp_data <- as.data.frame(t(exp))
group_data <- group[,2,drop=F]
data <- cbind(exp_data,group_data)
X_mat <- as.matrix(data[,c(1:(ncol(data)-1))])
X_mat <- scale(X_mat)
y_vec <- as.numeric(data$group)
n <- nrow(X_mat)
if (n <= 30) {
  nfold <- n   
} else {
  nfold <- min(5, floor(n/5)) 
}
dtrain <- xgb.DMatrix(data = X_mat, label = y_vec)

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

set.seed(123)
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

elog <- cv$evaluation_log
test_col <- grep("^test_.*_mean$", names(elog), value = TRUE)[1]
if (!is.null(cv$best_iteration)) {
  best_nrounds <- cv$best_iteration
} else if (!is.na(test_col)) {
  best_nrounds <- which.min(elog[[test_col]])
} else {
  best_nrounds <- nrow(elog)
}
message("best_nrounds = ", best_nrounds)

final_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = best_nrounds,
  watchlist = list(train = dtrain),
  verbose = 1
)

pred_train <- predict(final_model, X_mat)
pred_train <- as.numeric(pred_train)
y_true <- as.numeric(y_vec)
group$pred <- pred_train
write.csv(group,file.path(output,paste0(labels,GSE,"pred_data.csv")))

train_rmse <- sqrt(mean((pred_train - y_true)^2, na.rm = TRUE))
train_mae  <- mean(abs(pred_train - y_true), na.rm = TRUE)
train_r2   <- 1 - sum((pred_train - y_true)^2, na.rm = TRUE) / sum((y_true - mean(y_true, na.rm = TRUE))^2, na.rm = TRUE)

cat("Train RMSE:", train_rmse, "\n")
cat("Train MAE :", train_mae, "\n")
cat("Train R2  :", train_r2, "\n")


df_cal <- data.frame(
  pred = as.numeric(pred_train),
  obs  = as.numeric(y_true)
)
train_rmse <- sqrt(mean((df_cal$pred - df_cal$obs)^2, na.rm = TRUE))
bins <- min(10, nrow(df_cal))
df_cal <- df_cal %>% mutate(bin = ntile(pred, bins))

if (!"bin" %in% colnames(df_cal)) df_cal <- df_cal %>% mutate(bin = ntile(pred, 10))

calib_summary <- df_cal %>%
  dplyr::group_by(bin) %>%
  dplyr::summarise(
    mean_pred = mean(pred, na.rm = TRUE),
    mean_obs  = mean(obs,  na.rm = TRUE),
    n         = n(),
    sd_obs    = sd(obs, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
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
  geom_smooth(method = "loess", se = TRUE, color = "#F39B7F", fill = alpha("#F39B7F", 0.2)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#3C5488") +
  labs(x = "Predicted value", y = "Observed value",
       title = "") +
  theme_minimal(base_size = 12) +
  annotate("text", x = Inf, y = Inf, label = paste0("R² = ", formatC(train_r2, format = "f", digits = 3)),
           hjust = 1.1, vjust = 1.5, size = 8)
p_cal <- p_cal +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_rect(fill = "white"),
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16),  
    axis.text = element_text(size = 12)

  ) +
  coord_cartesian(expand = FALSE)  
print(p_cal)
ggsave(file.path(output, paste0(labels,GSE, "calibration.png")), p_cal, width = 6, height = 6, dpi = 300)
ggsave(file.path(output, paste0(labels,GSE, "calibration.pdf")), p_cal, width = 6, height = 6)
saveRDS(final_model,file.path(output,paste0(labels,GSE,"final_xgboost_model.Rdata")))

pred_train <- predict(final_model, X_mat)
shap_contrib <- predict(final_model, as.matrix(X_mat), predcontrib = TRUE)
shap_values <- shap_contrib[, -ncol(shap_contrib), drop = FALSE]
baseline_val <- mean(shap_contrib[, ncol(shap_contrib)])
sv <- shapviz(object = final_model, x = as.matrix(X_mat), X_pred = X_mat)
xvars <- gene
p <- sv_dependence(sv, v = xvars, share_y = TRUE)
ggsave(file.path(output,paste0(labels,GSE, "dependence_all.png")), p, width = 12, height = 8, dpi = 300)
ggsave(file.path(output,paste0(labels,GSE, "dependence_all.pdf")), p, width = 12, height = 8)
top_n <- 12
out_prefix <- file.path(output, paste0(labels,GSE))
pal_name <- "Spectral"

shap_mat <- NULL; X_df <- NULL; baseline_val <- NULL
shap_mat <- sv$S
X_df <- sv$X
baseline_val <- sv$baseline
colnames(shap_mat) <- colnames(X_df)
rownames(shap_mat) <- rownames(X_df)
features <- colnames(shap_mat)

write.csv(shap_mat,file.path(output,paste0(labels,GSE,"shap_mat.csv")))

mean_abs <- colMeans(abs(shap_mat), na.rm = TRUE)
ord <- order(mean_abs, decreasing = TRUE)
top_features <- features[ord][seq_len(min(top_n, length(features)))]
mean_abs_top <- mean_abs[top_features]

df_bar <- data.frame(feature = factor(top_features, levels = rev(top_features)),
                     mean_abs = mean_abs_top)
colors_bar <- c("#4DBBD5","#F1C40F","#E64B35")
p_bar <- ggplot(df_bar, aes(x = mean_abs, y = feature, fill = mean_abs)) +
  geom_col(width = 0.7) +
  scale_fill_gradientn(colors = colors_bar, guide = guide_colorbar(title = "mean |SHAP|")) +
  labs(x = "Mean |SHAP value|", y = NULL, title = paste0(" Features - SHAP importance")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5,size = 20),
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16),  
        axis.text = element_text(size = 12),
        axis.text.y = element_text(face = "italic"),
        legend.text = element_text(size = 12),
        legend.title = element_blank()
        ) +
  geom_text(aes(label = sprintf("%.3f", mean_abs)), 
            hjust = -0.1, size = 5) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2)))
ggsave(paste0(out_prefix, "shap_bar.png"), p_bar, width = 6, height = 4, dpi = 300)
ggsave(paste0(out_prefix, "shap_bar.pdf"), p_bar, width = 6, height = 4)

shap_df <- as.data.frame(shap_mat)
shap_df$sample <- rownames(shap_df)
long_shap <- shap_df %>%
  pivot_longer(cols = -sample, names_to = "feature", values_to = "shap")

X_long <- X_df %>% mutate(sample = rownames(.)) %>%
  pivot_longer(cols = -sample, names_to = "feature", values_to = "feat_value")

plot_df <- left_join(long_shap, X_long, by = c("sample", "feature")) %>%
  filter(feature %in% top_features) %>%
  mutate(feature = factor(feature, levels = rev(top_features)))

p_beeswarm <- ggplot(plot_df, aes(x = shap, y = feature, color = feat_value)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE, size = 2, alpha = 0.8, width = 0.3) +
  scale_color_gradient(low = "#91D1C2", high = "#DC0000", name = "Feature value")+
  geom_vline(xintercept = 0, color = "grey40", linetype = "dashed") +
  labs(x = "SHAP value", y = NULL, title = paste0("SHAP beeswarm")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5,size = 20),
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16),  
        axis.text = element_text(size = 12),
        axis.text.y = element_text(face = "italic"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
ggsave(paste0(out_prefix, "shap_beeswarm.png"), p_beeswarm, width = 6, height = 4, dpi = 300)
ggsave(paste0(out_prefix, "shap_beeswarm.pdf"), p_beeswarm, width = 6, height = 4)


p1 <- sv_waterfall(sv, row_id = 2,fill_colors = c("#8491B4","#F39B7F")) +
  theme(axis.text = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),

        )
p1_2 <- sv_force(sv, row_id = 2,fill_colors = c("#8491B4","#F39B7F"))+
  theme(axis.text = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),

  )
ggsave(paste0(out_prefix, "waterfall_sample2.png"), p1, width = 5.4, height = 4, dpi = 300)
ggsave(paste0(out_prefix, "waterfall_sample2.pdf"), p1, width = 5.4, height = 4)
ggsave(paste0(out_prefix, "waterfall_sample2_2.png"), p1_2, width = 4, height = 2, dpi = 300)
ggsave(paste0(out_prefix, "waterfall_sample2_2.pdf"), p1_2, width = 4, height = 2)

p2 <- sv_waterfall(sv, row_id = 38,fill_colors = c("#8491B4","#F39B7F")) +
  theme(axis.text = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
  )
p2_2 <- sv_force(sv, row_id = 38,fill_colors = c("#8491B4","#F39B7F"))+
  theme(axis.text = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
  )

ggsave(paste0(out_prefix, "waterfall_sample38.png"), p2, width = 5.4, height = 4, dpi = 300)
ggsave(paste0(out_prefix, "waterfall_sample38.pdf"), p2, width = 5.4, height = 4)
ggsave(paste0(out_prefix, "waterfall_sample2_38.png"), p2_2, width = 4, height = 2, dpi = 300)
ggsave(paste0(out_prefix, "waterfall_sample2_38.pdf"), p2_2, width = 4, height = 2)


####GSE15654####
GSE <- "GSE15654_"
gene = read.csv(file.path(ORIGINAL_DIR,"06_RF",paste0("06_RF_",GSE,"final_cogenes.csv"))) 
gene <- gene$x

exp <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE15654_exp.csv"),row.names = 1)
group <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE15654_group.csv"),row.names = 1)
surv <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE15654_sur.csv"),row.names = 1)
surv <- surv[colnames(exp),,drop=F]
exp <- exp[gene,]
exp_data <- as.data.frame(t(exp))
exp_data <- scale(exp_data)
Xmat <- as.matrix(exp_data)
time_train <- surv$days_to_death
status_train <- surv$death
dtrain <- xgb.DMatrix(data = Xmat, label = as.numeric(time_train))
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
set.seed(123)
cv <- xgb.cv(
  params = params,
  data = dtrain,
  nrounds = 2000,
  nfold = 5,
  early_stopping_rounds = 50,
  verbose = 1,
  showsd = TRUE
)
elog <- cv$evaluation_log
test_col <- grep("^test_.*_mean$", names(elog), value = TRUE)[1]
if (!is.null(cv$best_iteration)) {
  best_nrounds <- cv$best_iteration
} else if (!is.na(test_col)) {
  best_nrounds <- which.min(elog[[test_col]])
} else {
  best_nrounds <- nrow(elog)
}
message("best_nrounds = ", best_nrounds)
final_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = best_nrounds,
  watchlist = list(train = dtrain),
  verbose = 1
)
xgb.importance(model = final_model) %>% head(50)
saveRDS(final_model,file.path(output,paste0(labels,GSE,"final_xgboost_model_cox.Rdata")))
pred_train <- predict(final_model, Xmat)
group$pred <- pred_train
write.csv(group,file.path(output,paste0(labels,GSE,"pred_data.csv")))

table(pred_train)
pred_train
cf <- survival::concordancefit(Surv(time_train, status_train), pred_train)
str(cf)
cidx_train <- cf$concordance
print(cidx_train)
group_train <- ifelse(pred_train > median(pred_train, na.rm = TRUE), "High", "Low")
df_plot <- data.frame(time = time_train, status = status_train, group = factor(group_train, levels = c("Low","High")))
df_plot$time <- df_plot$time/365
fit <- survfit(Surv(time, status) ~ group, data = df_plot)
cox_uni <- coxph(Surv(time, status) ~ group, data = df_plot)
summary(cox_uni)
hr <- summary(cox_uni)$coefficients[,"exp(coef)"]
ci_low <- summary(cox_uni)$conf.int[,"lower .95"]
ci_high <- summary(cox_uni)$conf.int[,"upper .95"]
pval_cox <- summary(cox_uni)$coefficients[,"Pr(>|z|)"]
p_txt <- if (is.na(pval_cox)) "NA" else format(signif(pval_cox, 3), scientific = TRUE)

# plotmath 字符串：第一行 HR，第二行 italic(p)
p_note_expr <- sprintf(
  "atop('HR=%.2f (95%% CI %.2f-%.2f)', italic(p)==%s)",
  hr, ci_low, ci_high, p_txt
)
p_km <- ggsurvplot(
  fit,
  data = df_plot,
  risk.table = TRUE,    
  pval = F,           
  conf.int = FALSE,      
  surv.median.line = "hv",
  palette = c("#3C5488","#DC0000"),
  xlab = "Time (years)",
  legend.title = "Risk Score",
  legend.labs = c("Low","High") ,
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
png(file = file.path(output,paste0(labels,GSE,"KM_plot.png")),width =5,height = 5,res= 300,units = "in")
print(p_km)
dev.off()
pdf(file = file.path(output,paste0(labels,GSE,"KM_plot.pdf")),width =5,height = 5)
print(p_km)
dev.off()
X_mat <- Xmat
shap_mat <- predict(final_model, Xmat, predcontrib = TRUE)
pred_train <- predict(final_model, X_mat)
shap_contrib <- predict(final_model, as.matrix(X_mat), predcontrib = TRUE)
shap_values <- shap_contrib[, -ncol(shap_contrib), drop = FALSE]
baseline_val <- mean(shap_contrib[, ncol(shap_contrib)])
sv <- shapviz(object = final_model, x = as.matrix(X_mat), X_pred = X_mat)
xvars <- gene
p <- sv_dependence(sv, v = xvars, share_y = TRUE)
ggsave(file.path(output,paste0(labels,GSE, "dependence_all.png")), p, width = 12, height = 8, dpi = 300)
ggsave(file.path(output,paste0(labels,GSE, "dependence_all.pdf")), p, width = 12, height = 8)


top_n <- 12
out_prefix <- file.path(output, paste0(labels,GSE))
pal_name <- "Spectral"
shap_mat <- NULL; X_df <- NULL; baseline_val <- NULL
shap_mat <- sv$S
X_df <- sv$X
baseline_val <- sv$baseline
colnames(shap_mat) <- colnames(X_df)
rownames(shap_mat) <- rownames(X_df)
features <- colnames(shap_mat)

write.csv(shap_mat,file.path(output,paste0(labels,GSE,"shap_mat.csv")))
mean_abs <- colMeans(abs(shap_mat), na.rm = TRUE)
ord <- order(mean_abs, decreasing = TRUE)
top_features <- features[ord][seq_len(min(top_n, length(features)))]
mean_abs_top <- mean_abs[top_features]
df_bar <- data.frame(feature = factor(top_features, levels = rev(top_features)),
                     mean_abs = mean_abs_top)
colors_bar <- c("#4DBBD5","#F1C40F","#E64B35")
p_bar <- ggplot(df_bar, aes(x = mean_abs, y = feature, fill = mean_abs)) +
  geom_col(width = 0.7) +
  scale_fill_gradientn(colors = colors_bar, guide = guide_colorbar(title = "mean |SHAP|")) +
  labs(x = "Mean |SHAP value|", y = NULL, title = paste0(" Features - SHAP importance")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5,size = 20),
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16),  
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.text.y = element_text(face = "italic"),
        legend.title = element_blank()
  ) +
  geom_text(aes(label = sprintf("%.3f", mean_abs)), 
            hjust = -0.1, size = 5) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2)))
ggsave(paste0(out_prefix, "shap_bar.png"), p_bar, width = 6, height = 4, dpi = 300)
ggsave(paste0(out_prefix, "shap_bar.pdf"), p_bar, width = 6, height = 4)
shap_df <- as.data.frame(shap_mat)
shap_df$sample <- rownames(shap_df)
long_shap <- shap_df %>%
  pivot_longer(cols = -sample, names_to = "feature", values_to = "shap")

X_long <- X_df %>% mutate(sample = rownames(.)) %>%
  pivot_longer(cols = -sample, names_to = "feature", values_to = "feat_value")

plot_df <- left_join(long_shap, X_long, by = c("sample", "feature")) %>%
  filter(feature %in% top_features) %>%
  mutate(feature = factor(feature, levels = rev(top_features)))
p_beeswarm <- ggplot(plot_df, aes(x = shap, y = feature, color = feat_value)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE, size = 2, alpha = 0.8, width = 0.3) +
  scale_color_gradient(low = "#91D1C2", high = "#DC0000", name = "Feature value")+
  geom_vline(xintercept = 0, color = "grey40", linetype = "dashed") +
  labs(x = "SHAP value", y = NULL, title = paste0("SHAP beeswarm")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5,size = 20),
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16),  
        axis.text = element_text(size = 12),
        axis.text.y = element_text(face = "italic"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
ggsave(paste0(out_prefix, "shap_beeswarm.png"), p_beeswarm, width = 6, height = 4, dpi = 300)
ggsave(paste0(out_prefix, "shap_beeswarm.pdf"), p_beeswarm, width = 6, height = 4)

p1 <- sv_waterfall(sv, row_id = 2,fill_colors = c("#8491B4","#F39B7F")) +
  theme(axis.text = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
  )
p1_2 <- sv_force(sv, row_id = 2,fill_colors = c("#8491B4","#F39B7F"))+
  theme(axis.text = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
  )
ggsave(paste0(out_prefix, "waterfall_sample2.png"), p1, width = 5.4, height = 4, dpi = 300)
ggsave(paste0(out_prefix, "waterfall_sample2.pdf"), p1, width = 5.4, height = 4)
ggsave(paste0(out_prefix, "waterfall_sample2_2.png"), p1_2, width = 4, height = 2, dpi = 300)
ggsave(paste0(out_prefix, "waterfall_sample2_2.pdf"), p1_2, width = 4, height = 2)

p2 <- sv_waterfall(sv, row_id = 96,fill_colors = c("#8491B4","#F39B7F")) +
  theme(axis.text = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
  )
p2_2 <- sv_force(sv, row_id = 96,fill_colors = c("#8491B4","#F39B7F"))+
  theme(axis.text = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
  )

ggsave(paste0(out_prefix, "waterfall_sample96.png"), p2, width = 5.4, height = 4, dpi = 300)
ggsave(paste0(out_prefix, "waterfall_sample96.pdf"), p2, width = 5.4, height = 4)
ggsave(paste0(out_prefix, "waterfall_sample2_96.png"), p2_2, width = 4, height = 2, dpi = 300)
ggsave(paste0(out_prefix, "waterfall_sample2_96.pdf"), p2_2, width = 4, height = 2)

