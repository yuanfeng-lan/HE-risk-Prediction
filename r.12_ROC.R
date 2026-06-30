rm(list = ls()); gc()
ORIGINAL_DIR <- "~/data_hdd/R/Yclub/2026.3.4_HE_Cirrhosis"
output <- file.path(ORIGINAL_DIR, "12_ROC")
labels <- "12_ROC_"
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(ORIGINAL_DIR)
redame_tex <- "# 此步骤目的是对模型的风险评分通过ROC曲线来评其预测的性能"

library(ggplot2)
library(pROC)
library(xgboost)
library(timeROC)
library(survival)
####GSE139602####
GSE <- "GSE139602"
exp <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE139602_exp.csv"),row.names = 1)
group <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE139602_group.csv"),row.names = 1)
group <- group[colnames(exp),,drop=F]
table(group$characteristics_ch1)
group$group <- ifelse(group$characteristics_ch1=="disease state: Healthy","Healtehy",
                      ifelse(group$characteristics_ch1=="disease state: eCLD","eCLD",
                             ifelse(group$characteristics_ch1=="disease state: Compensated Cirrhosis","CC",
                                    ifelse(group$characteristics_ch1=="disease state: Decompesated Cirrhosis","DC",
                                           ifelse(group$characteristics_ch1=="disease state: Acute-on-chronic liver failure","ACLF",NA)))))
#读取风险评分表
pred_data <- read.csv(file.path(ORIGINAL_DIR,"07_xgboost",paste0("07_xgboost_",GSE,"_pred_data.csv")),row.names = 1)
group <- group[rownames(pred_data),]
pred_data$group <- group$group

group_combinations <- combn(unique(pred_data$group), 2, simplify = FALSE)

# 初始化存储结果的列表
roc_results <- list()
auc_values <- data.frame(Group1 = character(), Group2 = character(), AUC = numeric())

# 进行两两分组的ROC分析
for (comb in group_combinations) {
  group1 <- comb[1]
  group2 <- comb[2]
  
  # 筛选当前分组的数据
  data_subset <- pred_data %>% filter(group %in% c(group1, group2))
  
  # 创建二分类标签（1表示group1，0表示group2）
  data_subset$BinaryLabel <- ifelse(data_subset$group == group1, 0, 1)
  
  # 计算ROC曲线
  roc_obj <- roc(data_subset$BinaryLabel, data_subset$pred, levels = c(0, 1), direction = "<")
  
  # 保存ROC对象
  roc_results[[paste(group1, group2, sep = " vs ")]] <- roc_obj
  
  # 保存AUC值
  auc_values <- rbind(auc_values, data.frame(Group1 = group1, Group2 = group2, AUC = auc(roc_obj)))
}
auc_values
# 打印AUC值
plot_colors <- adjustcolor(c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000","#7E6148","#B09C85"), alpha.f =0.5) # 设置透明度

pdf(file.path(output,paste0(labels,GSE,"_ROC_results.pdf")),width = 6,height = 6)
plot(roc_results[[1]], col = plot_colors[1], main = "", lwd = 2, cex.lab = 1.5,xlab = "1 - Specificity", ylab = "Sensitivity",
     xlim = c(0, 1), ylim = c(0, 1))
for (i in 2:length(roc_results)) {
  plot(roc_results[[i]], col = plot_colors[i], add = TRUE, lwd = 2)
}

# 添加图例，并标注AUC值
legend_labels <- sapply(names(roc_results), function(name) {
  auc_value <- auc(roc_results[[name]])
  paste(name, sprintf("(AUC = %.2f)", auc_value))
})
legend("bottomright", legend = legend_labels, col = plot_colors, lwd = 2)

dev.off()

png(file.path(output,paste0(labels,GSE,"_ROC_results.png")),width = 6,height = 6,units = "in",res = 300)
plot(roc_results[[1]], col = plot_colors[1], main = "", lwd = 2, cex.lab = 1.5,xlab = "1 - Specificity", ylab = "Sensitivity",
     xlim = c(0, 1), ylim = c(0, 1))
for (i in 2:length(roc_results)) {
  plot(roc_results[[i]], col = plot_colors[i], add = TRUE, lwd = 2)
}

# 添加图例，并标注AUC值
legend_labels <- sapply(names(roc_results), function(name) {
  auc_value <- auc(roc_results[[name]])
  paste(name, sprintf("(AUC = %.2f)", auc_value))
})
legend("bottomright", legend = legend_labels, col = plot_colors, lwd = 2)

dev.off()
write.csv(auc_values,file.path(output,paste0(labels,GSE,"ROC_results.csv")))

####GSE15654####
GSE <- 'GSE15654'
gene = read.csv(file.path(ORIGINAL_DIR,"06_RF",paste0("06_RF_",GSE,"_final_cogenes.csv"))) 
gene <- gene$x
exp <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE15654_exp.csv"),row.names = 1)
group <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE15654_group.csv"),row.names = 1)
surv <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE15654_sur.csv"),row.names = 1)
surv <- surv[colnames(exp),,drop=F]
final_model <- readRDS(file.path("07_xgboost",paste0("07_xgboost_",GSE,"_final_xgboost_model_cox.Rdata")))
exp <- exp[gene,]
exp_data <- as.data.frame(t(exp))
exp_data <- scale(exp_data)
Xmat <- as.matrix(exp_data)

pred_train <- predict(final_model, Xmat)
group$pred <- pred_train
surv <- surv[rownames(group),]
group$group <- surv$death
group$time <- surv$days_to_death
group$time <- group$time/365

set.seed(123)
pred_data <- group
# 计算3年、5年和10年的生存ROC曲线
roc_3yr <- timeROC(T = pred_data$time, delta = pred_data$group, marker = pred_data$pred, 
                   cause = 1, times = 3, iid = TRUE)
roc_5yr <- timeROC(T = pred_data$time, delta = pred_data$group, marker = pred_data$pred, 
                   cause = 1, times = 5, iid = TRUE)
roc_10yr <- timeROC(T = pred_data$time, delta = pred_data$group, marker = pred_data$pred, 
                    cause = 1, times = 10, iid = TRUE)

# 提取AUC值
auc_3yr <- roc_3yr$AUC[2]
auc_5yr <- roc_5yr$AUC[2]
auc_10yr <- roc_10yr$AUC[2]

# 打印AUC值
cat(sprintf("3-year AUC: %.2f\n", auc_3yr))
cat(sprintf("5-year AUC: %.2f\n", auc_5yr))
cat(sprintf("10-year AUC: %.2f\n", auc_10yr))

# Prepare ROC data for plotting
roc_data <- data.frame(
  FPR = c(roc_3yr$FP[, 1], roc_5yr$FP[, 1], roc_10yr$FP[, 1]),
  TPR = c(roc_3yr$TP[, 1], roc_5yr$TP[, 1], roc_10yr$TP[, 1]),
  Time = factor(rep(c("3 years", "5 years", "10 years"), each = length(roc_3yr$FP[, 1])))
)

# Custom colors
plot_colors <- adjustcolor(c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000","#7E6148","#B09C85"), alpha.f = 0.85)

pdf(file.path(output,paste0(labels,GSE,"_ROC_results.pdf")),width = 6,height = 6)

# 绘制3年生存ROC曲线
plot(roc_3yr$FP[, 1], roc_3yr$TP[, 2], type = "l", col = plot_colors[1], lwd = 2,
     xlab = "1 - Specificity", ylab = "Sensitivity", cex.lab = 1.5,
     main = "",
     xlim = c(0, 1), ylim = c(0, 1))

# 添加对角线（y = x）
abline(a = 0, b = 1, lty = 2, col = "gray") # lty = 2 表示虚线，col = "gray" 表示灰色

# 添加5年生存ROC曲线
lines(roc_5yr$FP[, 1], roc_5yr$TP[, 2], col = plot_colors[2], lwd = 2)

# 添加10年生存ROC曲线
lines(roc_10yr$FP[, 1], roc_10yr$TP[, 2], col = plot_colors[3], lwd = 2)

# 添加图例，并标注AUC值
legend("bottomright", legend = c(
  sprintf("3 years (AUC = %.2f)", auc_3yr),
  sprintf("5 years (AUC = %.2f)", auc_5yr),
  sprintf("10 years (AUC = %.2f)", auc_10yr)
), col = plot_colors, lwd = 2, cex = 1.2)
dev.off()

png(file.path(output,paste0(labels,GSE,"_ROC_results.png")),width = 6,height = 6,units = "in",res = 300)
plot(roc_3yr$FP[, 1], roc_3yr$TP[, 2], type = "l", col = plot_colors[1], lwd = 2,
     xlab = "1 - Specificity", ylab = "Sensitivity", cex.lab = 1.5,
     main = "",
     xlim = c(0, 1), ylim = c(0, 1))

# 添加对角线（y = x）
abline(a = 0, b = 1, lty = 2, col = "gray") # lty = 2 表示虚线，col = "gray" 表示灰色

# 添加5年生存ROC曲线
lines(roc_5yr$FP[, 1], roc_5yr$TP[, 2], col = plot_colors[2], lwd = 2)

# 添加10年生存ROC曲线
lines(roc_10yr$FP[, 1], roc_10yr$TP[, 2], col = plot_colors[3], lwd = 2)

# 添加图例，并标注AUC值
legend("bottomright", legend = c(
  sprintf("3 years (AUC = %.2f)", auc_3yr),
  sprintf("5 years (AUC = %.2f)", auc_5yr),
  sprintf("10 years (AUC = %.2f)", auc_10yr)
), col = plot_colors, lwd = 2, cex = 1.2)
dev.off()


####GSE41919####
GSE <- "GSE41919"
gene = read.csv(file.path(ORIGINAL_DIR,"06_RF",paste0("06_RF_","GSE139602_","final_cogenes.csv"))) 
gene <- gene$x

GSE41919_exp = read.csv(file = './00_rawdata/00.rawdata_GSE41919_exp.csv', header = TRUE,row.names = 1) 
GSE41919_group = read.csv(file = './00_rawdata/00.rawdata_GSE41919_group.csv',  header = TRUE,row.names = 1) 
GSE41919_group <- subset(GSE41919_group,group!="non-cirrhotic control")
GSE41919_exp <- GSE41919_exp[,rownames(GSE41919_group)]
GSE41919_exp <- as.data.frame(t(GSE41919_exp))
GSE41919_df <- GSE41919_exp[,gene]
GSE41919_df <- scale(GSE41919_df)
final_model <- readRDS(file.path(ORIGINAL_DIR, "07_xgboost",paste0("07_xgboost_","GSE139602","_final_xgboost_model.Rdata")))

GSE41919_risk_score <- predict(final_model,GSE41919_df)
GSE41919_group$risk_score_1 <- GSE41919_risk_score

gene = read.csv(file.path(ORIGINAL_DIR,"06_RF","06_RF_GSE15654_final_cogenes.csv")) 
gene <- gene$x

GSE41919_df_2 <- GSE41919_exp[,gene]
GSE41919_df_2 <- scale(GSE41919_df_2)
final_model_2 <- readRDS(file.path(ORIGINAL_DIR, "07_xgboost",paste0("07_xgboost_","GSE15654","_final_xgboost_model_cox.Rdata")))
GSE41919_risk_score <- predict(final_model_2,GSE41919_df_2)
GSE41919_group$risk_score_2 <- GSE41919_risk_score
GSE41919_group$risk_score <- (GSE41919_group$risk_score_1 + GSE41919_group$risk_score_2)/2

GSE41919_group$BinaryLabel <- ifelse(GSE41919_group$group == "cirrhosis without HE", 0, 1)

# 计算ROC曲线
roc_obj_01 <- roc(GSE41919_group$BinaryLabel, GSE41919_group$risk_score, levels = c(0, 1), direction = "<")

####GSE57193
GSE57193_exp = read.csv(file = './00_rawdata/00.rawdata_GSE57193_exp.csv', header = TRUE,row.names = 1) 
GSE57193_group = read.csv(file = './00_rawdata/00.rawdata_GSE57193_group.csv',  header = TRUE,row.names = 1) 
GSE57193_group <- subset(GSE57193_group,group!="healthy")
GSE57193_exp <- GSE57193_exp[,rownames(GSE57193_group)]
GSE57193_exp <- as.data.frame(t(GSE57193_exp))
gene = read.csv(file.path(ORIGINAL_DIR,"06_RF",paste0("06_RF_","GSE139602_","final_cogenes.csv"))) 
gene <- gene$x

GSE57193_df <- GSE57193_exp[,gene]
GSE57193_df <- scale(GSE57193_df)
GSE57193_risk_score <- predict(final_model,GSE57193_df)
GSE57193_group$risk_score_1 <- GSE57193_risk_score
GSE57193_group$group <- ifelse(GSE57193_group$group=="cirrhosis","cirrhosis without HE","cirrhosis with HE")
GSE57193_group$group <- factor(GSE57193_group$group,levels = c("cirrhosis without HE","cirrhosis with HE"))
gene = read.csv(file.path(ORIGINAL_DIR,"06_RF","06_RF_GSE15654_final_cogenes.csv")) 
gene <- gene$x
GSE57193_df_2 <- GSE57193_exp[,gene]
GSE57193_df_2 <- scale(GSE57193_df_2)

final_model_2 <- readRDS(file.path(ORIGINAL_DIR, "07_xgboost",paste0("07_xgboost_","GSE15654","_final_xgboost_model_cox.Rdata")))
GSE57193_risk_score <- predict(final_model_2,GSE57193_df_2)
GSE57193_group$risk_score_2 <- GSE57193_risk_score
GSE57193_group$risk_score <- (GSE57193_group$risk_score_1 + GSE57193_group$risk_score_2)/2

GSE57193_group$BinaryLabel <- ifelse(GSE57193_group$group == "cirrhosis without HE", 1, 0)

# 计算ROC曲线
roc_obj_02 <- roc(GSE57193_group$BinaryLabel, GSE57193_group$risk_score, levels = c(0,1), direction = ">")
pdf(file.path(output,paste0(labels,GSE,"_ROC_results.pdf")),width = 6,height = 6)
plot(roc_obj_01, col = plot_colors[1], main = "Cirrhosis without HE vs with HE", lwd = 2, cex.lab = 1.5,xlab = "1 - Specificity", ylab = "Sensitivity",
     xlim = c(1, 0), ylim = c(0, 1))
plot(roc_obj_02, col = plot_colors[2], add = TRUE, lwd = 2,
     xlim = c(1, 0), ylim = c(0, 1))
legend("bottomright", legend = c(
  sprintf("GSE41919 (AUC = %.2f)", auc(roc_obj_01)),
  sprintf("GSE57193 (AUC = %.2f)", auc(roc_obj_02))
), col = plot_colors, lwd = 2, cex = 1.2)
dev.off()

png(file.path(output,paste0(labels,GSE,"_ROC_results.png")),width = 6,height = 6,units = "in",res = 300)
plot(roc_obj_01, col = plot_colors[1], main = "Cirrhosis without HE vs with HE", lwd = 2, cex.lab = 1.5,xlab = "1 - Specificity", ylab = "Sensitivity",
     xlim = c(1, 0), ylim = c(0, 1))
plot(roc_obj_02, col = plot_colors[2], add = TRUE, lwd = 2,
     xlim = c(1, 0), ylim = c(0, 1))
legend("bottomright", legend = c(
  sprintf("GSE41919 (AUC = %.2f)", auc(roc_obj_01)),
  sprintf("GSE57193 (AUC = %.2f)", auc(roc_obj_02))
), col = plot_colors, lwd = 2, cex = 1.2)
dev.off()

