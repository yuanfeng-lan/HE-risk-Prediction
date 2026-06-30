rm(list = ls()); gc()
ORIGINAL_DIR <- "~/data_hdd/R/Yclub/2026.3.4_HE_Cirrhosis"
output <- file.path(ORIGINAL_DIR, "06_RF")
labels <- "06_RF_"

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(ORIGINAL_DIR)
library(randomForest)
library(randomForestSRC)
library(dplyr); library(tidyr); library(ggplot2)
library("ggvenn")

####GSE139602####
GSE <- "GSE139602_"
HE_GSE139602_UP_gene <- read.csv(file.path(ORIGINAL_DIR,"02_Mfuzz","HE_GSE139602_UP_gene.csv"),row.names = 1)
HE_GSE139602_DOWN_gene <- read.csv(file.path(ORIGINAL_DIR,"02_Mfuzz","HE_GSE139602_DOWN_gene.csv"),row.names = 1)
topgene <- c(HE_GSE139602_UP_gene$x,HE_GSE139602_DOWN_gene$x)
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
set.seed(005)
forest_result <- randomForest(group ~ ., data = data, importance = TRUE)
forest_result
imp <- as.data.frame(importance(forest_result))
imp$feature <- rownames(imp)
mse_col <- "%IncMSE"
node_col <- "IncNodePurity"

imp2 <- imp %>%
  mutate(mse_z = as.numeric(scale(.data[[mse_col]])),
         node_z = as.numeric(scale(.data[[node_col]])),
         mean_z = (mse_z + node_z) / 2) %>%
  arrange(mean_z)  

imp2$feature <- factor(imp2$feature, levels = imp2$feature)

topN <- 10
top_feats <- imp2 %>% slice_max(mean_z, n = topN) %>% pull(feature)
imp2 <- imp2 %>% mutate(top10 = as.logical(feature %in% top_feats))

imp_long <- imp2[, c("feature", "mse_z", "node_z", "mean_z", "top10")] %>%
  pivot_longer(cols = c(mse_z, node_z), names_to = "metric", values_to = "zvalue") %>%
  mutate(
    metric = ifelse(metric == "mse_z", "%IncMSE (z)", "IncNodePurity (z)"),
    feature = factor(feature, levels = levels(imp2$feature)) 
  )
rng <- range(imp2$mean_z, na.rm = TRUE)
label_x <- -3

p <- ggplot() +
  geom_col(data = imp2, aes(x = mean_z, y = feature, fill = top10), width = 0.6, alpha = 0.6) +
  geom_line(data = imp_long, aes(x = zvalue, y = feature, group = feature),
            color = "#7F8C8D", linetype = "dashed", size = 0.3) +
  geom_point(data = imp_long, aes(x = zvalue, y = feature, color = metric, shape = metric),
             size = 3) + 
  geom_label(data = imp2,
             aes(x = label_x, y = feature, 
                 label = paste0("italic('", feature, "')")), # 用引号包裹基因名
             fill = "white",                       
             color = ifelse(imp2$top10, "#E64B35", "#7F8C8D"),  
             hjust = 0,                            
             size = 4.5,
             parse = TRUE, # 启用解析以支持表达式
             linewidth = 0.15,                    
             label.r = grid::unit(0.12, "lines")) +
  scale_fill_manual(values = c("FALSE" = "#7F8C8D", "TRUE" = "#E64B35"), guide = "none") +
  scale_color_manual(values = c("%IncMSE (z)" = "#3C5488", "IncNodePurity (z)" = "#00A087")) +
  labs(x = "Standardized importance (z)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_blank(),
    axis.title = element_text(size = 15), 
    axis.text = element_text(size = 12), 
    legend.text = element_text(size = 12)
  ) +
  coord_cartesian(xlim = c(label_x, rng[2] + 0.02 * (rng[2] - rng[1])))

ggsave(file.path(output,paste0(labels,GSE,"RF_results.png")),p,width = 6,height = 6)
ggsave(file.path(output,paste0(labels,GSE,"RF_results.pdf")),p,width = 6,height = 6)
write.csv(imp2,file.path(output,paste0(labels,GSE,"RF_results_table.csv")))
results <- subset(imp2,top10=="TRUE")
choose_gene <- rownames(results)
write.csv(choose_gene,file.path(output,paste0(labels,GSE,"RF_results.csv")))

LASSO_results <- read.csv(file.path(ORIGINAL_DIR,"04_LASSO",paste0("04_LASSO_",GSE,"result.csv")))
RFE_results <- read.csv(file.path(ORIGINAL_DIR,"05_RFE",paste0("05_RFE_",GSE,"results.csv")))
RF_results <- read.csv(file.path(ORIGINAL_DIR,"06_RF",paste0("06_RF_",GSE,"RF_results.csv")))

venn <- list(LASSO=LASSO_results$x,
             RFE = RFE_results$x,
             RandomForest = RF_results$x)
p <- ggvenn(
  venn,  
  c("LASSO","RFE","RandomForest"), 
  text_size = 8,  
  fill_color = c("#00A087", "#3C5488","#F39B7F"),
  fill_alpha = 0.7,  
  stroke_color = "black",  
  stroke_size = 1 ,
  show_percentage = F
)
ggsave(file.path(output, paste0(labels,GSE,"venn_final_cogenes.pdf")),p,w=4,h=4)
ggsave(file.path(output, paste0(labels,GSE,"venn_final_cogenes.png")),p,w=4,h=4)
cogene <- intersect(LASSO_results$x,RFE_results$x)
cogene <- intersect(cogene,RF_results$x)
write.csv(cogene,file.path(output,paste0(labels,GSE,"final_cogenes.csv")))


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
group_data <- surv[,c(1,2),drop=F]
data <- cbind(exp_data,group_data)
form <- as.formula(Surv(days_to_death, death) ~ .)

mf <- model.frame(form, data = data)
n_xvar <- ncol(mf) - 1   

mtry_val <- floor(sqrt(n_xvar))
if (mtry_val < 1) mtry_val <- 1

set.seed(002)
rf_surv <- randomForestSRC::rfsrc(
  formula = form,
  data = data,
  ntree = 1000,
  mtry = mtry_val,
  importance = TRUE,
  block.size = 1
)
print(rf_surv)
var_importance <- rf_surv$importance
head(sort(var_importance, decreasing = TRUE), 20)
imp <- as.data.frame(var_importance)
if (ncol(imp) == 1) {
  imp$importance <- as.numeric(imp[[1]])
} else {
  imp$importance <- rowMeans(imp, na.rm = TRUE)
}
imp$feature <- rownames(imp)

imp2 <- imp %>%
  mutate(importance_z = as.numeric(scale(importance))) %>%
  arrange(importance_z)   

imp2$feature <- factor(imp2$feature, levels = imp2$feature)

topN <- 10
top_feats <- imp2 %>% slice_max(importance_z, n = topN) %>% pull(feature)
imp2 <- imp2 %>% mutate(top10 = feature %in% top_feats)

label_x <- -1.5
rng <- range(imp2$importance_z, na.rm = TRUE)

p <- ggplot() +
  geom_col(data = imp2, aes(x = importance_z, y = feature, fill = top10), width = 0.6, alpha = 0.6) +
  geom_point(data = imp2, aes(x = importance_z, y = feature, color = top10), size = 2) +
  geom_label(data = imp2,
             aes(x = label_x, y = feature, 
                 label = paste0("italic('", feature, "')")), # 用引号包裹基因名
             fill = "white",                       
             color = ifelse(imp2$top10, "#E64B35", "#7F8C8D"),  
             hjust = 0,                            
             size = 4.5,
             parse = TRUE, # 启用解析以支持表达式
             linewidth = 0.15,                    
             label.r = grid::unit(0.12, "lines")) +
  scale_fill_manual(values = c("FALSE" = "#7F8C8D", "TRUE" = "#E64B35"), guide = "none") +
  scale_color_manual(values = c("FALSE" = "#7F8C8D", "TRUE" = "#E64B35"), guide = "none") +
  labs(x = "Standardized importance (z)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_blank(),
    axis.title = element_text(size=15), 
    axis.text = element_text(size=12), 
    legend.text = element_text(size=12), 
  ) +
  coord_cartesian(xlim = c(label_x, rng[2] + 0.02 * (rng[2] - rng[1])))

ggsave(file.path(output,paste0(labels,GSE,"RF_results.png")),p,width = 4,height = 6)
ggsave(file.path(output,paste0(labels,GSE,"RF_results.pdf")),p,width = 4,height = 6)
write.csv(imp2,file.path(output,paste0(labels,GSE,"RF_results_table.csv")))

results <- subset(imp2,top10=="TRUE")
choose_gene <- rownames(results)
write.csv(choose_gene,file.path(output,paste0(labels,GSE,"RF_results.csv")))


LASSO_results <- read.csv(file.path(ORIGINAL_DIR,"04_LASSO",paste0("04_LASSO_",GSE,"result.csv")))
RFE_results <- read.csv(file.path(ORIGINAL_DIR,"05_RFE",paste0("05_RFE_",GSE,"results.csv")))
RF_results <- read.csv(file.path(ORIGINAL_DIR,"06_RF",paste0("06_RF_",GSE,"RF_results.csv")))

venn <- list(LASSO=LASSO_results$x,
             RFE = RFE_results$x,
             RandomForest = RF_results$x)
p <- ggvenn(
  venn,  
  c("LASSO","RFE","RandomForest"), 
  text_size = 8,  
  fill_color = c("#00A087", "#3C5488","#F39B7F"),
  fill_alpha = 0.7,  
  stroke_color = "black",  
  stroke_size = 1 ,
  show_percentage = F
)
ggsave(file.path(output, paste0(labels,GSE,"venn_final_cogenes.pdf")),p,w=4,h=4)
ggsave(file.path(output, paste0(labels,GSE,"venn_final_cogenes.png")),p,w=4,h=4)
cogene <- intersect(LASSO_results$x,RFE_results$x)
cogene <- intersect(cogene,RF_results$x)
write.csv(cogene,file.path(output,paste0(labels,GSE,"final_cogenes.csv")))
