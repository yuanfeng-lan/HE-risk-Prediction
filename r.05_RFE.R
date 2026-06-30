rm(list = ls()); gc()
ORIGINAL_DIR <- "~/data_hdd/R/Yclub/2026.3.4_HE_Cirrhosis"
output <- file.path(ORIGINAL_DIR, "05_RFE")
labels <- "05_RFE_"

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(ORIGINAL_DIR)
library(tidyverse)
library(caret)
library(ggplot2)
library(ggplotify)
library(gridExtra) # 用于拼图
library(grid)   
library(ggtext)
####GSE139602####
GSE <- "GSE139602_"
HE_GSE139602_UP_gene <- read.csv(file.path(ORIGINAL_DIR,"02_Mfuzz","HE_GSE139602_UP_gene.csv"),row.names = 1)
HE_GSE139602_DOWN_gene <- read.csv(file.path(ORIGINAL_DIR,"02_Mfuzz","HE_GSE139602_DOWN_gene.csv"),row.names = 1)
topgene <- c(HE_GSE139602_UP_gene$x,HE_GSE139602_DOWN_gene$x)
topgene
exp <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE139602_exp.csv"),row.names = 1)
group <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE139602_group.csv"),row.names = 1)
group <- group[colnames(exp),,drop=F]
table(group$characteristics_ch1)
group$group <- ifelse(group$characteristics_ch1=="disease state: Healthy",0,
                      ifelse(group$characteristics_ch1=="disease state: eCLD",1,
                             ifelse(group$characteristics_ch1=="disease state: Compensated Cirrhosis",2,
                                    ifelse(group$characteristics_ch1=="disease state: Decompesated Cirrhosis",3,
                                           ifelse(group$characteristics_ch1=="disease state: Acute-on-chronic liver failure",4,NA)))))
exp <- exp[topgene,]
exp_data <- as.data.frame(t(exp))
group_data <- group[,2,drop=F]
data <- cbind(exp_data,group_data)
dat <- data
num <- ncol(dat)-1
set.seed(047)
control <- rfeControl(functions = lmFuncs, method = "repeatedcv", number = 5, repeats = 3) 
results <- rfe(dat[,1:num],            
               dat$group,  
               sizes = c(1:num),
               rfeControl = control,
               method = "RMSE")
predictors_results <- predictors(results)
write.csv(predictors(results), file.path(output,paste0(labels,GSE,"results.csv") ))

df <- results$results
df$Variables <- as.numeric(as.character(df$Variables))
best_n <- as.numeric(results[["bestSubset"]])
best_row <- which(df$Variables == best_n)[1]
best_rmse <- df$RMSE[best_row]
p1 <- ggplot(df, aes(x = Variables, y = RMSE)) +
  geom_line(color = "#2c7fb8", size = 1) +
  geom_point(color = "#2c7fb8", size = 3) +
  { if ("RMSESD" %in% names(df)) 
    geom_errorbar(aes(ymin = RMSE - RMSESD, ymax = RMSE + RMSESD), width = 0.2, color = "#2c7fb8", alpha = 0.6)
    else NULL } +
  geom_point(data = df[best_row, , drop = FALSE], aes(x = Variables, y = RMSE), color = "#e31a1c", size = 4) +
  geom_vline(xintercept = best_n, linetype = "dashed", colour = "grey40") +
  annotate("text", x = best_n, y = best_rmse, label = paste0("Best: ", best_n, "\nRMSE=", round(best_rmse, 3)),
           vjust = -1.2, hjust = 0.5, color = "#e31a1c", size = 6) +
  scale_x_continuous(breaks = df$Variables) +
  labs(x = "Number of Features", y = "RMSE (Repeated Cross-Validation)",
       title = "RFE Performance", subtitle = "",
       caption = "") +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(color = "black", size = 15),
    panel.border = element_rect(colour = "black", fill = NA, size = 1) 
  )

# 创建列标签
topgene_labels <- sapply(topgene, function(gene) {
  if (gene %in% predictors_results) {
    paste0("<span style='color:red;'><i>", gene, "</i></span>") # 红色高亮并斜体
  } else {
    paste0("<i>", gene, "</i>") # 仅斜体
  }
})

# 转换为数据框以便绘图
label_df <- data.frame(
  x = rep(1, length(topgene)),
  y = seq_along(topgene),
  label = sapply(topgene, function(gene) {
    if (gene %in% predictors_results) {
      paste0("<span style='color:red;'><i>", gene, "</i></span>") # 红色高亮并斜体
    } else {
      paste0("<i>", gene, "</i>") # 仅斜体
    }
  })
)

# 绘制列标签
p2 <- ggplot(label_df, aes(x = x, y = y, label = label)) +
  geom_richtext(aes(label = label), hjust = 0, size = 4, fill = NA, label.color = NA) +
  scale_y_reverse(breaks = seq_along(topgene), labels = NULL) +
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, -2, "cm"),
    panel.grid = element_blank()
  )


# 拼图
pdf(file.path(output, paste0(labels, GSE, "RFE_results.pdf")), width = 6, height = 6) # 增加宽度
grid.arrange(
  p1, p2,
  ncol = 2,
  widths = c(4, 1), # 增加 p2 的宽度比例
  layout_matrix = rbind(c(1, 2)) # 确保两张图紧密排列
)
dev.off()

png(file.path(output,paste0(labels,GSE,"RFE_results.png")),width = 6,height = 6,units = "in",res = 300)
grid.arrange(
  p1, p2,
  ncol = 2,
  widths = c(4, 1), # 增加 p2 的宽度比例
  layout_matrix = rbind(c(1, 2)) # 确保两张图紧密排列
)
dev.off()
####GSE15654####
GSE <- "GSE15654_"
topgene <- read.csv(file.path(ORIGINAL_DIR,"02_Mfuzz","HE_GSE15654_ALL_gene_DEG.csv"))
topgene <- topgene$x
exp <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE15654_exp.csv"),row.names = 1)
group <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE15654_group.csv"),row.names = 1)
surv <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE15654_sur.csv"),row.names = 1)
surv <- surv[colnames(exp),,drop=F]
exp <- exp[topgene,]
exp_data <- as.data.frame(t(exp))
group_data <- surv[,2,drop=F]
data <- cbind(exp_data,group_data)
dat <- data
num <- ncol(dat)-1
set.seed(018)
control <- rfeControl(functions = lmFuncs, method = "repeatedcv", number = 5, repeats = 3) 
results <- rfe(dat[,1:num],            
               dat$death,  
               sizes = c(1:num),
               rfeControl = control,
               method = "RMSE")

predictors_results <- predictors(results)
write.csv(predictors(results), file.path(output,paste0(labels,GSE,"results.csv") ))
df <- results$results
df$Variables <- as.numeric(as.character(df$Variables))
best_n <- as.numeric(results[["bestSubset"]])
best_row <- which(df$Variables == best_n)[1]
best_rmse <- df$RMSE[best_row]

p <- ggplot(df, aes(x = Variables, y = RMSE)) +
  geom_line(color = "#2c7fb8", size = 1) +
  geom_point(color = "#2c7fb8", size = 3) +
  { if ("RMSESD" %in% names(df)) 
    geom_errorbar(aes(ymin = RMSE - RMSESD, ymax = RMSE + RMSESD), width = 0.2, color = "#2c7fb8", alpha = 0.6)
    else NULL } +
  geom_point(data = df[best_row, , drop = FALSE], aes(x = Variables, y = RMSE), color = "#e31a1c", size = 4) +
  geom_vline(xintercept = best_n, linetype = "dashed", colour = "grey40") +
  annotate("text", x = best_n, y = best_rmse, label = paste0("Best: ", best_n, "\nRMSE=", round(best_rmse, 3)),
           vjust = -1.2, hjust = 0.5, color = "#e31a1c", size = 6) +
  scale_x_continuous(breaks = df$Variables) +
  labs(x = "Number of Features", y = "RMSE (Repeated Cross-Validation)",
       title = "RFE Performance", subtitle = "",
       caption = "") +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95"),
    plot.title = element_text(face = "bold", hjust = 0.5,size = 20),
    axis.title = element_text(face = "bold",size = 15),
    axis.text = element_text(color = "black",size = 15),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)  
  )
ggsave(file.path(output,paste0(labels,GSE,"RFE_results.png")),p,width = 6,height = 6)
ggsave(file.path(output,paste0(labels,GSE,"RFE_results.pdf")),p,width = 6,height = 6)

