rm(list = ls()); gc()
ORIGINAL_DIR <- "~/data_hdd/R/Yclub/2026.3.4_HE_Cirrhosis"
output <- file.path(ORIGINAL_DIR, "08_model_validation")
labels <- "08_model_validation_"

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
exp_data <- as.data.frame(scale(exp_data))
group_data <- group[,2,drop=F]
data <- cbind(exp_data,group_data)
X_mat <- as.matrix(data[,c(1:(ncol(data)-1))])
y_vec <- as.numeric(data$group)
final_model <- readRDS(file.path(ORIGINAL_DIR, "07_xgboost",paste0("07_xgboost_",GSE,"final_xgboost_model.Rdata")))
pred_train <- predict(final_model, X_mat)
data$pred <- pred_train
data_df <- cbind(exp_data,pred_train)
data_df$group <- ifelse(group$characteristics_ch1=="disease state: Healthy","Healthy",
                        ifelse(group$characteristics_ch1=="disease state: eCLD","eCLD",
                               ifelse(group$characteristics_ch1=="disease state: Compensated Cirrhosis","CC",
                                      ifelse(group$characteristics_ch1=="disease state: Decompesated Cirrhosis","DC",
                                             ifelse(group$characteristics_ch1=="disease state: Acute-on-chronic liver failure","ACLF",NA)))))
data_df$group <- factor(data_df$group,levels = c("Healthy","eCLD","CC","DC","ACLF"))
outdir <- file.path(output,"plot_01_GSE139602")
cols <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F")
names(cols) <- unique(data_df$group)
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}
test_plot1 <- function(gene, gene_label, palette = "Set2"){
  library(ggplot2); library(ggpubr)
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")
  library(RColorBrewer)
  stopifnot(exists("data_df"))
  stopifnot(dir.exists(outdir))
  data_df$group <- factor(data_df$group)
  groups <- levels(data_df$group)
  cols <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F")
  names(cols) <- groups
  test_formula <- as.formula(paste(gene, "~ group"))
  p_value <- compare_means(test_formula, data = data_df, method = "kruskal.test")$p.format
  p <- ggplot(data_df, aes(x = group, y = .data[[gene]], fill = group, colour = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.5, colour = "black", size = 0.4) +
    geom_jitter(position = position_jitter(width = 0.15, height = 0), size = 1.8, alpha = 0.9) +
    scale_fill_manual(values = cols) +
    scale_colour_manual(values = cols) +
    theme_minimal(base_size = 15) +
    labs(x = "", y = paste0("Relative gene expression"),title = gene_label) +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 14),
          axis.text.x  = element_text(size = 14, face = "bold",angle = 45,vjust = 1, hjust = 1, margin = margin(t = 2)),
          axis.text.y  = element_text(size = 14),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 16, face = "italic"))+
    coord_cartesian(clip = "off") +
    annotation_custom(
      grob = textGrob(
        label = bquote(italic(p) ~ "=" ~ .(p_value)),
        gp = gpar(fontsize = 15),
        hjust = 0.5
      ),
      xmin = 1.1,  
      xmax = 1.5,
      ymin = max(data_df[[gene]], na.rm = TRUE) * 0.85,  
      ymax = max(data_df[[gene]], na.rm = TRUE) * 1.05
    )
  
  
  ggsave(file.path(outdir, paste0(labels, gene_label, "_plot.png")), p, width = 4, height = 4, dpi = 300)
  ggsave(file.path(outdir, paste0(labels, gene_label, "_plot.pdf")), p, width = 4, height = 4)
}
colnames(exp_data)
test_plot1("NPC2","NPC2")
test_plot1("TLN1","TLN1")
test_plot1("TUBA1C","TUBA1C")
test_plot1("LRRC32","LRRC32")
test_plot1("PRB2","PRB2")
#GSE41919
GSE41919_exp = read.csv(file = './00_rawdata/00.rawdata_GSE41919_exp.csv', header = TRUE,row.names = 1) 
GSE41919_group = read.csv(file = './00_rawdata/00.rawdata_GSE41919_group.csv',  header = TRUE,row.names = 1) 
GSE41919_group <- subset(GSE41919_group,group!="non-cirrhotic control")
GSE41919_exp <- GSE41919_exp[,rownames(GSE41919_group)]
GSE41919_exp <- as.data.frame(t(GSE41919_exp))
GSE41919_df <- GSE41919_exp[,gene]
GSE41919_df <- scale(GSE41919_df)
GSE41919_risk_score <- predict(final_model,GSE41919_df)
GSE41919_group$risk_score <- GSE41919_risk_score
GSE41919_group$group <- factor(GSE41919_group$group,levels = c("cirrhosis without HE","cirrhosis with HE"))
p_value <- compare_means(risk_score ~ group, data = GSE41919_group, method = "wilcox.test")$p.format
p <- ggplot(GSE41919_group, aes(x = group, y = risk_score, fill = group, colour = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.5, colour = "black", size = 0.4) +
  geom_jitter(position = position_jitter(width = 0.15, height = 0), size = 1.8, alpha = 0.9) +
  scale_fill_manual(values = c("#4DBBD5","#E64B35")) +
  scale_colour_manual(values = c("#4DBBD5","#E64B35")) +
  theme_minimal(base_size = 15) +
  scale_x_discrete(labels = c(paste0("Cirrhosis","\n","without HE"),paste0("Cirrhosis","\n","with HE")))+
  labs(x = "", y = "Risk Score from Model",title = "GSE41919") +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 14),
        axis.text.x  = element_text(size = 14, face = "bold",angle = 45,vjust = 1, hjust = 1, margin = margin(t = 2)),
        axis.text.y  = element_text(size = 14),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))+
  coord_cartesian(clip = "off") +
  annotation_custom(
    grob = textGrob(
      label = bquote(italic(p) ~ "=" ~ .(p_value)), 
      gp = gpar(fontsize = 15),  
      hjust = 0.5
    ),
    xmin = 0.2, 
    xmax = 1.5,
    ymin = max(GSE41919_group$risk_score, na.rm = TRUE) * 0.85,  
    ymax = max(GSE41919_group$risk_score, na.rm = TRUE) * 1.05
  ) 

ggsave(file.path(outdir, paste0(labels,GSE,"Risk_score_GSE41919_plot.png")), p, width = 4, height = 4)
ggsave(file.path(outdir, paste0(labels,GSE,"Risk_score_GSE41919_plot.pdf")), p, width = 4, height = 4)
write.csv(GSE41919_group,file.path(outdir,paste0(labels,GSE,"Risk_score_GSE41919_data.csv")))
#GSE57193
GSE57193_exp = read.csv(file = './00_rawdata/00.rawdata_GSE57193_exp.csv', header = TRUE,row.names = 1) 
GSE57193_group = read.csv(file = './00_rawdata/00.rawdata_GSE57193_group.csv',  header = TRUE,row.names = 1) 
GSE57193_group <- subset(GSE57193_group,group!="healthy")
GSE57193_exp <- GSE57193_exp[,rownames(GSE57193_group)]
GSE57193_exp <- as.data.frame(t(GSE57193_exp))
GSE57193_df <- GSE57193_exp[,gene]
GSE57193_df <- scale(GSE57193_df)
GSE57193_risk_score <- predict(final_model,GSE57193_df)
GSE57193_group$risk_score <- GSE57193_risk_score
GSE57193_group$group <- ifelse(GSE57193_group$group=="cirrhosis","cirrhosis without HE","cirrhosis with HE")
GSE57193_group$group <- factor(GSE57193_group$group,levels = c("cirrhosis without HE","cirrhosis with HE"))
p_value <- compare_means(risk_score ~ group, data = GSE57193_group, method = "wilcox.test")$p.format
p <- ggplot(GSE57193_group, aes(x = group, y = risk_score, fill = group, colour = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.5, colour = "black", size = 0.4) +
  geom_jitter(position = position_jitter(width = 0.15, height = 0), size = 1.8, alpha = 0.9) +
  scale_fill_manual(values = c("#4DBBD5","#E64B35")) +
  scale_colour_manual(values = c("#4DBBD5","#E64B35")) +
  theme_minimal(base_size = 15) +
  scale_x_discrete(labels = c(paste0("Cirrhosis","\n","without HE"),paste0("Cirrhosis","\n","with HE")))+
  labs(x = "", y = "Risk Score from Model",title = "GSE57193") +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 14),
        axis.text.x  = element_text(size = 14, face = "bold",angle = 45,vjust = 1, hjust = 1, margin = margin(t = 2)),
        axis.text.y  = element_text(size = 14),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))+
  coord_cartesian(clip = "off") +
  annotation_custom(
    grob = textGrob(
      label = bquote(italic(p) ~ "=" ~ .(p_value)), 
      gp = gpar(fontsize = 15),  
      hjust = 0.5
    ),
    xmin = 0.2,  
    xmax = 1.5,
    ymin = max(GSE41919_group$risk_score, na.rm = TRUE) * 0.85,  
    ymax = max(GSE41919_group$risk_score, na.rm = TRUE) * 1.05
  ) 

ggsave(file.path(outdir, paste0(labels,GSE,"Risk_score_GSE57193_plot.png")), p, width = 4, height = 4)
ggsave(file.path(outdir, paste0(labels,GSE,"Risk_score_GSE57193_plot.pdf")), p, width = 4, height = 4)
write.csv(GSE57193_group,file.path(outdir,paste0(labels,GSE,"Risk_score_GSE57193_data.csv")))

groups <- levels(data_df$group)
cols <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F")
names(cols) <- groups

p_value <- compare_means(pred_train ~ group, data = data_df, method = "kruskal.test")$p.format
p <- ggplot(data_df, aes(x = group, y = pred_train, fill = group, colour = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.5, colour = "black", size = 0.4) +
  geom_jitter(position = position_jitter(width = 0.15, height = 0), size = 1.8, alpha = 0.9) +
  scale_fill_manual(values = cols) +
  scale_colour_manual(values = cols) +
  theme_minimal(base_size = 15) +
  labs(x = "", y = paste0("Risk Score"),title = "") +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 14),
        axis.text.x  = element_text(size = 14, face = "bold",angle = 45,vjust = 1, hjust = 1, margin = margin(t = 2)),
        axis.text.y  = element_text(size = 14),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))+
  coord_cartesian(clip = "off") +
  annotation_custom(
    grob = textGrob(
      label = bquote(italic(p) ~ "=" ~ .(p_value)), 
      gp = gpar(fontsize = 15),  
      hjust = 0.5
    ),
    xmin = 1.1, 
    xmax = 1.5,
    ymin = max(data_df$pred_train, na.rm = TRUE) * 0.85,  
    ymax = max(data_df$pred_train, na.rm = TRUE) * 1.05
  )

ggsave(file.path(outdir, paste0(labels,GSE,"risk_score_plot.png")), p, width = 4, height = 4, dpi = 300)
ggsave(file.path(outdir, paste0(labels,GSE,"risk_score_plot.pdf")), p, width = 4, height = 4)
#





####GSE15654####
GSE <- "GSE15654_"
gene = read.csv(file.path(ORIGINAL_DIR,"06_RF","06_RF_GSE15654_final_cogenes.csv")) 
gene <- gene$x
gene
exp <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE15654_exp.csv"),row.names = 1)
group <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE15654_group.csv"),row.names = 1)
surv <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE15654_sur.csv"),row.names = 1)
surv <- surv[colnames(exp),,drop=F]
exp <- exp[gene,]
exp_data <- as.data.frame(t(exp))
exp_data <- scale(exp_data)
Xmat <- as.matrix(exp_data)
time_train <- surv$days_to_death
time_train <- time_train/365
status_train <- surv$death
exp_data <- as.data.frame(exp_data)
final_model <- readRDS(file.path(ORIGINAL_DIR,"07_xgboost","07_xgboost_GSE15654_final_xgboost_model_cox.Rdata"))
pred_train <- predict(final_model, Xmat)

for (i in gene) {
  group_train <- ifelse(exp_data[[i]] > median(exp_data[[i]], na.rm = TRUE), "High", "Low")
  df_plot <- data.frame(time = time_train, status = status_train, group = factor(group_train, levels = c("Low","High")))
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
    pval = FALSE,                  # 关闭默认 p 值文本
    conf.int = FALSE,
    surv.median.line = "hv",
    palette = c("#3C5488", "#DC0000"),
    xlab = "Time (years)",
    ylab = "Survival probability (OS)",
    legend.title = i,
    legend.labs = c("Low", "High"),
    risk.table.height = 0.25,
    ggtheme = theme_minimal(base_size = 13)
  )
  
  # 手动添加：p 斜体
  p_km$plot <- p_km$plot +
    annotate(
      "text",
      x = 10, y = 0.05,
      label = p_note_expr,
      parse = TRUE,
      hjust = 1.02, vjust = 0,
      size = 4
    )
  p_km$plot <- p_km$plot +
    labs(color = i, fill = i) +
    theme(legend.title = element_text(face = "italic"))
  
  png(file = file.path(output,paste0(labels,GSE,i,"_KM_plot.png")),width =5,height = 5,res= 300,units = "in")
  print(p_km)
  dev.off()
  pdf(file = file.path(output,paste0(labels,GSE,i,"_KM_plot.pdf")),width =5,height = 5)
  print(p_km)
  dev.off()
}

time_train <- surv$days_to_decomp
time_train <- time_train/365
status_train <- surv$decomp

for (i in gene) {
  group_train <- ifelse(exp_data[[i]] > median(exp_data[[i]], na.rm = TRUE), "High", "Low")
  df_plot <- data.frame(time = time_train, status = status_train, group = factor(group_train, levels = c("Low","High")))
  fit <- survfit(Surv(time, status) ~ group, data = df_plot)
  cox_uni <- coxph(Surv(time, status) ~ group, data = df_plot)
  summary(cox_uni)
  hr <- summary(cox_uni)$coefficients[,"exp(coef)"]
  ci_low <- summary(cox_uni)$conf.int[,"lower .95"]
  ci_high <- summary(cox_uni)$conf.int[,"upper .95"]
  pval_cox <- summary(cox_uni)$coefficients[,"Pr(>|z|)"]
  p_txt <- if (is.na(pval_cox)) "NA" else format(signif(pval_cox, 3), scientific = TRUE)
  
  p_note_expr <- sprintf(
    "atop('HR=%.2f (95%% CI %.2f-%.2f)', italic(p)==%s)",
    hr, ci_low, ci_high, p_txt
  )
  p_km <- ggsurvplot(
    fit,
    data = df_plot,
    risk.table = TRUE,    
    pval = FALSE,           
    conf.int = FALSE,      
    surv.median.line = "hv",
    palette = c("#3C5488","#DC0000"),
    xlab = "Time (years)",
    ylab = "Cumulative survival free of decompensation",
    legend.title = i,
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
  p_km$plot <- p_km$plot +
    labs(color = i, fill = i) +
    theme(legend.title = element_text(face = "italic"))
  png(file = file.path(output,paste0(labels,GSE,i,"_Decompensation−free survival_plot.png")),width =5,height = 5,res= 300,units = "in")
  print(p_km)
  dev.off()
  pdf(file = file.path(output,paste0(labels,GSE,i,"_Decompensation−free survival_plot.pdf")),width =5,height = 5)
  print(p_km)
  dev.off()
}



time_train <- surv$days_to_child
time_train <- time_train/365
status_train <- surv$child

for (i in gene) {
  group_train <- ifelse(exp_data[[i]] > median(exp_data[[i]], na.rm = TRUE), "High", "Low")
  df_plot <- data.frame(time = time_train, status = status_train, group = factor(group_train, levels = c("Low","High")))
  fit <- survfit(Surv(time, status) ~ group, data = df_plot)
  cox_uni <- coxph(Surv(time, status) ~ group, data = df_plot)
  summary(cox_uni)
  hr <- summary(cox_uni)$coefficients[,"exp(coef)"]
  ci_low <- summary(cox_uni)$conf.int[,"lower .95"]
  ci_high <- summary(cox_uni)$conf.int[,"upper .95"]
  pval_cox <- summary(cox_uni)$coefficients[,"Pr(>|z|)"]
  p_txt <- if (is.na(pval_cox)) "NA" else format(signif(pval_cox, 3), scientific = TRUE)
  
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
    ylab = "Cumulative survival free of Child endpoint",
    legend.title = i,
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
  p_km$plot <- p_km$plot +
    labs(color = i, fill = i) +
    theme(legend.title = element_text(face = "italic"))
  png(file = file.path(output,paste0(labels,GSE,i,"_Cumulative survival free of Child endpoint_plot.png")),width =5,height = 5,res= 300,units = "in")
  print(p_km)
  dev.off()
  pdf(file = file.path(output,paste0(labels,GSE,i,"_Cumulative survival free of Child endpoint_plot.pdf")),width =5,height = 5)
  print(p_km)
  dev.off()
}

time_train <- surv$days_to_hcc
time_train <- time_train/365
status_train <- surv$hcc

for (i in gene) {
  group_train <- ifelse(exp_data[[i]] > median(exp_data[[i]], na.rm = TRUE), "High", "Low")
  df_plot <- data.frame(time = time_train, status = status_train, group = factor(group_train, levels = c("Low","High")))
  fit <- survfit(Surv(time, status) ~ group, data = df_plot)
  cox_uni <- coxph(Surv(time, status) ~ group, data = df_plot)
  summary(cox_uni)
  hr <- summary(cox_uni)$coefficients[,"exp(coef)"]
  ci_low <- summary(cox_uni)$conf.int[,"lower .95"]
  ci_high <- summary(cox_uni)$conf.int[,"upper .95"]
  pval_cox <- summary(cox_uni)$coefficients[,"Pr(>|z|)"]
  p_txt <- if (is.na(pval_cox)) "NA" else format(signif(pval_cox, 3), scientific = TRUE)
  
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
    ylab = "Cumulative survival free of HCC progression",
    legend.title = i,
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
  p_km$plot <- p_km$plot +
    labs(color = i, fill = i) +
    theme(legend.title = element_text(face = "italic"))
  png(file = file.path(output,paste0(labels,GSE,i,"_Cumulative survival free of HCC progression_plot.png")),width =5,height = 5,res= 300,units = "in")
  print(p_km)
  dev.off()
  pdf(file = file.path(output,paste0(labels,GSE,i,"_Cumulative survival free of HCC progression_plot.pdf")),width =5,height = 5)
  print(p_km)
  dev.off()
}
#
group_train <- ifelse(pred_train > median(pred_train, na.rm = TRUE), "High", "Low")
time_train <- surv$days_to_decomp
time_train <- time_train/365
status_train <- surv$decomp
df_plot <- data.frame(time = time_train, status = status_train, group = factor(group_train, levels = c("Low","High")))
fit <- survfit(Surv(time, status) ~ group, data = df_plot)
cox_uni <- coxph(Surv(time, status) ~ group, data = df_plot)
summary(cox_uni)
hr <- summary(cox_uni)$coefficients[,"exp(coef)"]
ci_low <- summary(cox_uni)$conf.int[,"lower .95"]
ci_high <- summary(cox_uni)$conf.int[,"upper .95"]
pval_cox <- summary(cox_uni)$coefficients[,"Pr(>|z|)"]
p_txt <- if (is.na(pval_cox)) "NA" else format(signif(pval_cox, 3), scientific = TRUE)

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
  ylab = "Cumulative survival free of decompensation",
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
png(file = file.path(output,paste0(labels,GSE,"_Decompensation−free survival_plot.png")),width =5,height = 5,res= 300,units = "in")
print(p_km)
dev.off()
pdf(file = file.path(output,paste0(labels,GSE,"_Decompensation−free survival_plot.pdf")),width =5,height = 5)
print(p_km)
dev.off()
#
group_train <- ifelse(pred_train > median(pred_train, na.rm = TRUE), "High", "Low")
time_train <- surv$days_to_child
time_train <- time_train/365
status_train <- surv$child
df_plot <- data.frame(time = time_train, status = status_train, group = factor(group_train, levels = c("Low","High")))
fit <- survfit(Surv(time, status) ~ group, data = df_plot)
cox_uni <- coxph(Surv(time, status) ~ group, data = df_plot)
summary(cox_uni)
hr <- summary(cox_uni)$coefficients[,"exp(coef)"]
ci_low <- summary(cox_uni)$conf.int[,"lower .95"]
ci_high <- summary(cox_uni)$conf.int[,"upper .95"]
pval_cox <- summary(cox_uni)$coefficients[,"Pr(>|z|)"]
p_txt <- if (is.na(pval_cox)) "NA" else format(signif(pval_cox, 3), scientific = TRUE)

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
  ylab = "Cumulative survival free of Child endpoint",
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
png(file = file.path(output,paste0(labels,GSE,"_Cumulative survival free of Child endpoint_plot.png")),width =5,height = 5,res= 300,units = "in")
print(p_km)
dev.off()
pdf(file = file.path(output,paste0(labels,GSE,"_Cumulative survival free of Child endpoint_plot.pdf")),width =5,height = 5)
print(p_km)
dev.off()

#
group_train <- ifelse(pred_train > median(pred_train, na.rm = TRUE), "High", "Low")
time_train <- surv$days_to_hcc
time_train <- time_train/365
status_train <- surv$hcc
df_plot <- data.frame(time = time_train, status = status_train, group = factor(group_train, levels = c("Low","High")))
fit <- survfit(Surv(time, status) ~ group, data = df_plot)
cox_uni <- coxph(Surv(time, status) ~ group, data = df_plot)
summary(cox_uni)
hr <- summary(cox_uni)$coefficients[,"exp(coef)"]
ci_low <- summary(cox_uni)$conf.int[,"lower .95"]
ci_high <- summary(cox_uni)$conf.int[,"upper .95"]
pval_cox <- summary(cox_uni)$coefficients[,"Pr(>|z|)"]
p_txt <- if (is.na(pval_cox)) "NA" else format(signif(pval_cox, 3), scientific = TRUE)

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
  ylab = "Cumulative survival free of HCC progression",
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
png(file = file.path(output,paste0(labels,GSE,"_Cumulative survival free of HCC progression_plot.png")),width =5,height = 5,res= 300,units = "in")
print(p_km)
dev.off()
pdf(file = file.path(output,paste0(labels,GSE,"_Cumulative survival free of HCC progression_plot.pdf")),width =5,height = 5)
print(p_km)
dev.off()


#GSE41919
outdir <- file.path(output,paste0(GSE,"GSE41919"))
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}
GSE41919_exp = read.csv(file = './00_rawdata/00.rawdata_GSE41919_exp.csv', header = TRUE,row.names = 1) 
GSE41919_group = read.csv(file = './00_rawdata/00.rawdata_GSE41919_group.csv',  header = TRUE,row.names = 1) 
GSE41919_group <- subset(GSE41919_group,group!="non-cirrhotic control")
GSE41919_exp <- GSE41919_exp[,rownames(GSE41919_group)]
GSE41919_exp <- as.data.frame(t(GSE41919_exp))
GSE41919_df <- GSE41919_exp[,gene]
GSE41919_df <- scale(GSE41919_df)
GSE41919_risk_score <- predict(final_model,GSE41919_df)
GSE41919_group$risk_score <- GSE41919_risk_score
GSE41919_group$group <- factor(GSE41919_group$group,levels = c("cirrhosis without HE","cirrhosis with HE"))
p_value <- compare_means(risk_score ~ group, data = GSE41919_group, method = "wilcox.test")$p.format
p <- ggplot(GSE41919_group, aes(x = group, y = risk_score, fill = group, colour = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.5, colour = "black", size = 0.4) +
  geom_jitter(position = position_jitter(width = 0.15, height = 0), size = 1.8, alpha = 0.9) +
  scale_fill_manual(values = c("#4DBBD5","#E64B35")) +
  scale_colour_manual(values = c("#4DBBD5","#E64B35")) +
  theme_minimal(base_size = 15) +
  scale_x_discrete(labels = c(paste0("Cirrhosis","\n","without HE"),paste0("Cirrhosis","\n","with HE")))+
  labs(x = "", y = "Risk Score from Model",title = "GSE41919") +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 14),
        axis.text.x  = element_text(size = 14, face = "bold",angle = 45,vjust = 1, hjust = 1, margin = margin(t = 2)),
        axis.text.y  = element_text(size = 14),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))+
  coord_cartesian(clip = "off") +
  annotation_custom(
    grob = textGrob(
      label = bquote(italic(p) ~ "=" ~ .(p_value)), 
      gp = gpar(fontsize = 15),  
      hjust = 0.5
    ),
    xmin = 0.2,  
    xmax = 1.5,
    ymin = max(GSE41919_group$risk_score, na.rm = TRUE) * 0.85,  
    ymax = max(GSE41919_group$risk_score, na.rm = TRUE) * 1.05
  ) 
ggsave(file.path(outdir, paste0(labels,GSE,"Risk_score_GSE41919_plot.png")), p, width = 4, height = 4)
ggsave(file.path(outdir, paste0(labels,GSE,"Risk_score_GSE41919_plot.pdf")), p, width = 4, height = 4)
write.csv(GSE41919_group,file.path(outdir,paste0(labels,GSE,"Risk_score_GSE41919_data.csv")))

data_df <- GSE41919_exp[,gene]
for (i in gene) {
  plot_df <- data_df[,i,drop = F]
  plot_df$group <- ifelse(plot_df[,i]>=median(plot_df[,i]),"High","Low")
  pre <- as.formula(paste0(i," ~ ","group"))
  p_value <- compare_means(pre, data = plot_df, method = "wilcox.test")$p.format
  p <- ggplot(plot_df, aes(x = group, y = plot_df[,i], fill = group, colour = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.5, colour = "black", size = 0.4) +
    geom_jitter(position = position_jitter(width = 0.15, height = 0), size = 1.8, alpha = 0.9) +
    scale_fill_manual(values = c("#4DBBD5","#E64B35")) +
    scale_colour_manual(values = c("#4DBBD5","#E64B35")) +
    theme_minimal(base_size = 15) +
    scale_x_discrete(labels = c(paste0("Cirrhosis","\n","without HE"),paste0("Cirrhosis","\n","with HE")))+
    labs(x = "", y = "Relative gene expression",title = "GSE41919",subtitle = i) +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 14),
          axis.text.x  = element_text(size = 14, face = "bold",angle = 45,vjust = 1, hjust = 1, margin = margin(t = 2)),
          axis.text.y  = element_text(size = 14),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 16,face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 14),face = "italic")+
    coord_cartesian(clip = "off") +
    annotation_custom(
      grob = textGrob(
        label = bquote(italic(p) ~ "=" ~ .(p_value)), 
        gp = gpar(fontsize = 15), 
        hjust = 0.5
      ),
      xmin = 2,  
      xmax = 1.5,
      ymin = max(plot_df[,i], na.rm = TRUE) * 0.85,  
      ymax = max(plot_df[,i], na.rm = TRUE) * 1.05
    ) 
  ggsave(file.path(outdir, paste0(labels,GSE,i,"_GSE41919_plot.png")), p, width = 4, height = 4)
  ggsave(file.path(outdir, paste0(labels,GSE,i,"_GSE41919_plot.pdf")), p, width = 4, height = 4)
}


#GSE57193
outdir <- file.path(output,paste0(GSE,"GSE57193"))
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}
GSE57193_exp = read.csv(file = './00_rawdata/00.rawdata_GSE57193_exp.csv', header = TRUE,row.names = 1) 
GSE57193_group = read.csv(file = './00_rawdata/00.rawdata_GSE57193_group.csv',  header = TRUE,row.names = 1) 
GSE57193_group <- subset(GSE57193_group,group!="healthy")
GSE57193_exp <- GSE57193_exp[,rownames(GSE57193_group)]
GSE57193_exp <- as.data.frame(t(GSE57193_exp))
GSE57193_df <- GSE57193_exp[,gene]
GSE57193_df <- scale(GSE57193_df)
GSE57193_risk_score <- predict(final_model,GSE57193_df)
GSE57193_group$risk_score <- GSE57193_risk_score
GSE57193_group$group <- ifelse(GSE57193_group$group=="cirrhosis","cirrhosis without HE","cirrhosis with HE")
GSE57193_group$group <- factor(GSE57193_group$group,levels = c("cirrhosis without HE","cirrhosis with HE"))
p_value <- compare_means(risk_score ~ group, data = GSE57193_group, method = "wilcox.test")$p.format
p <- ggplot(GSE57193_group, aes(x = group, y = risk_score, fill = group, colour = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.5, colour = "black", size = 0.4) +
  geom_jitter(position = position_jitter(width = 0.15, height = 0), size = 1.8, alpha = 0.9) +
  scale_fill_manual(values = c("#4DBBD5","#E64B35")) +
  scale_colour_manual(values = c("#4DBBD5","#E64B35")) +
  theme_minimal(base_size = 15) +
  scale_x_discrete(labels = c(paste0("Cirrhosis","\n","without HE"),paste0("Cirrhosis","\n","with HE")))+
  labs(x = "", y = "Risk Score from Model",title = "GSE57193") +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 14),
        axis.text.x  = element_text(size = 14, face = "bold",angle = 45,vjust = 1, hjust = 1, margin = margin(t = 2)),
        axis.text.y  = element_text(size = 14),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))+
  coord_cartesian(clip = "off") +
  annotation_custom(
    grob = textGrob(
      label = bquote(italic(p) ~ "=" ~ .(p_value)), 
      gp = gpar(fontsize = 15),  
      hjust = 0.5
    ),
    xmin = 2,  
    xmax = 1.5,
    ymin = max(GSE41919_group$risk_score, na.rm = TRUE) * 0.85,  
    ymax = max(GSE41919_group$risk_score, na.rm = TRUE) * 1.05
  ) 

ggsave(file.path(outdir, paste0(labels,GSE,"Risk_score_GSE57193_plot.png")), p, width = 4, height = 4)
ggsave(file.path(outdir, paste0(labels,GSE,"Risk_score_GSE57193_plot.pdf")), p, width = 4, height = 4)
write.csv(GSE57193_group,file.path(outdir,paste0(labels,GSE,"Risk_score_GSE57193_data.csv")))

data_df <- GSE57193_exp[,gene]
for (i in gene) {
  plot_df <- data_df[,i,drop = F]
  plot_df$group <- ifelse(plot_df[,i]>=median(plot_df[,i]),"High","Low")
  pre <- as.formula(paste0(i," ~ ","group"))
  p_value <- compare_means(pre, data = plot_df, method = "wilcox.test")$p.format
  p <- ggplot(plot_df, aes(x = group, y = plot_df[,i], fill = group, colour = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.5, colour = "black", size = 0.4) +
    geom_jitter(position = position_jitter(width = 0.15, height = 0), size = 1.8, alpha = 0.9) +
    scale_fill_manual(values = c("#4DBBD5","#E64B35")) +
    scale_colour_manual(values = c("#4DBBD5","#E64B35")) +
    theme_minimal(base_size = 15) +
    scale_x_discrete(labels = c(paste0("Cirrhosis","\n","without HE"),paste0("Cirrhosis","\n","with HE")))+
    labs(x = "", y = "Relative gene expression",title = "GSE57193",subtitle = i) +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 14),
          axis.text.x  = element_text(size = 14, face = "bold",angle = 45,vjust = 1, hjust = 1, margin = margin(t = 2)),
          axis.text.y  = element_text(size = 14),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 16,face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 14,face = "italic"))+
    coord_cartesian(clip = "off") +
    annotation_custom(
      grob = textGrob(
        label = bquote(italic(p) ~ "=" ~ .(p_value)), 
        gp = gpar(fontsize = 15), 
        hjust = 0.5
      ),
      xmin = 2,  
      xmax = 1.5,
      ymin = max(plot_df[,i], na.rm = TRUE) * 0.85,  
      ymax = max(plot_df[,i], na.rm = TRUE) * 1.05
    ) 
  ggsave(file.path(outdir, paste0(labels,GSE,i,"_GSE41919_plot.png")), p, width = 4, height = 4)
  ggsave(file.path(outdir, paste0(labels,GSE,i,"_GSE41919_plot.pdf")), p, width = 4, height = 4)
}



