rm(list = ls()); gc()
ORIGINAL_DIR <- "~/data_hdd/R/Yclub/2026.3.4_HE_Cirrhosis"
output <- file.path(ORIGINAL_DIR, "09_SEM")
labels <- "09_SEM_"

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(ORIGINAL_DIR)
library(tidyverse)
library(mediation)
library(limma)
library(sva)
library(corrplot)  
library(pheatmap)  
library(lavaan)
library(semPlot)
library(dplyr)
library(ggplot2)
library(scales)
library(grid)
library(ggrepel)
exp01 <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE139602_exp.csv"),row.names = 1)
group01 <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE139602_group.csv"),row.names = 1)
group01 <- group01[colnames(exp01),,drop=F]
table(group01$characteristics_ch1)
group01$group <- ifelse(group01$characteristics_ch1=="disease state: Healthy",0,
                        ifelse(group01$characteristics_ch1=="disease state: eCLD",1,
                               ifelse(group01$characteristics_ch1=="disease state: Compensated Cirrhosis",2,
                                      ifelse(group01$characteristics_ch1=="disease state: Decompesated Cirrhosis",3,
                                             ifelse(group01$characteristics_ch1=="disease state: Acute-on-chronic liver failure",4,NA)))))
exp02 <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE15654_exp.csv"),row.names = 1)
group02 <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE15654_group.csv"),row.names = 1)

model2 <- readRDS(file.path(ORIGINAL_DIR,"07_xgboost","07_xgboost_GSE15654_final_xgboost_model_cox.Rdata"))
model1 <- readRDS(file.path(ORIGINAL_DIR,"07_xgboost","07_xgboost_GSE139602_final_xgboost_model.Rdata"))

model1_gene <- read.csv(file.path(ORIGINAL_DIR,"07_xgboost","07_xgboost_GSE139602_shap_mat.csv"),row.names = 1)
model2_gene <- read.csv(file.path(ORIGINAL_DIR,"07_xgboost","07_xgboost_GSE15654_shap_mat.csv"),row.names = 1)

all_feature <- c(colnames(model1_gene),colnames(model2_gene))

intersect(all_feature,rownames(exp01))
intersect(all_feature,rownames(exp02))

samegene <- intersect(rownames(exp01),rownames(exp02))
exp01 <- exp01[samegene,]
exp02 <- exp02[samegene,]
exp_merged <- cbind(exp01,exp02)

group01$project <- "GSE139602"
group02$project <- "GSE15654"
group_merged <- data.frame(row.names = c(rownames(group01),rownames(group02)),
                           group = c(group01$group,group02$group),
                           project = c(group01$project,group02$project))
group_merged <- group_merged[colnames(exp_merged),]

pca <- prcomp(t(exp_merged), scale. = TRUE)
pca_data <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], 
                       Sample = colnames(exp_merged),
                       Group = group_merged$group)  
pca_data$project <- group_merged$project
var_explained <- pca$sdev^2 / sum(pca$sdev^2)
pc1_var <- round(var_explained[1] * 100, 2)
pc2_var <- round(var_explained[2] * 100, 2)
png(file.path(output,paste0(labels,"pre_merged.png")),width = 5,height = 4,res = 300,units = "in")
boxplot(exp_merged, 
        xaxt = "n",  
        col = "lightblue",  
        main = "Gene Expression Distribution (Pre-batch Correction)",  
        ylab = "Expression Value"  
)
dev.off()
pdf(file.path(output,paste0(labels,"pre_merged.pdf")),width = 5,height = 4)
boxplot(exp_merged, 
        xaxt = "n",  
        col = "lightblue",  
        main = "Gene Expression Distribution (Pre-batch Correction)",  
        ylab = "Expression Value" 
)
dev.off()


group_merged$sample <- rownames(group_merged)
if ("sample" %in% colnames(group_merged)) {
  sample_order <- colnames(exp_merged)
  batch_info <- group_merged$project[match(sample_order, group_merged$sample)]
  unique_batches <- unique(batch_info)
  cat("Unique batches:", unique_batches, "\n")
  cat("Number of unique batches:", length(unique_batches), "\n")
  
  if (length(unique_batches) > 1 && !any(is.na(unique_batches))) {
    if ("group" %in% colnames(group_merged)) {
      group_info <- group_merged$group[match(sample_order, group_merged$sample)]
      mod <- model.matrix(~ group_info)
      exp_corrected <- ComBat(dat = exp_merged, batch = batch_info)
    } else {
      warning("")
      exp_corrected <- ComBat(dat = exp_merged, batch = batch_info)
    }
  } else {
    warning("")
    exp_corrected <- exp_merged
  }
} else {
  stop("")
}
pca <- prcomp(t(exp_corrected), scale. = TRUE)
pca_data <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], 
                       Sample = colnames(exp_corrected),
                       Group = group_merged$group)  
pca_data$project <- group_merged$project
var_explained <- pca$sdev^2 / sum(pca$sdev^2)
pc1_var <- round(var_explained[1] * 100, 2)
pc2_var <- round(var_explained[2] * 100, 2)

ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 1,alpha = 0.5) +
  labs(x = paste0("PC1 (", pc1_var, "%)"), y = paste0("PC2 (", pc2_var, "%)"), title = "PCA of Corrected Expression Data") +
  theme_minimal()
write.csv(group_merged,file.path(output,paste0(labels,"_group_merged.csv")))
write.csv(exp_merged,file.path(output,paste0(labels,"_exp_merged.csv")))
write.csv(exp_corrected,file.path(output,paste0(labels,"_exp_corrected.csv")))

exp_corrected <- read.csv(file.path(output,paste0(labels,"_exp_corrected.csv")),row.names = 1)
png(file.path(output,paste0(labels,"post_merged.png")),width = 5,height = 4,res = 300,units = "in")
boxplot(exp_corrected, 
        xaxt = "n", 
        col = "lightblue", 
        main = "Gene Expression Distribution (Post-batch Correction)", 
        ylab = "Expression Value"  
)
dev.off()
pdf(file.path(output,paste0(labels,"post_merged.pdf")),width = 5,height = 4)
boxplot(exp_corrected, 
        xaxt = "n", 
        col = "lightblue", 
        main = "Gene Expression Distribution (Post-batch Correction)", 
        ylab = "Expression Value"
)
dev.off()



exp_data <- exp_corrected[all_feature,]


exp_normality_test_genes <- function(exp_data, by_group = FALSE){
  if (is.data.frame(exp_data)) exp_data <- as.matrix(exp_data)
  infer_group <- function(){
    if (exists("group0", envir = .GlobalEnv)) {
      gdf <- get("group0", envir = .GlobalEnv)
      if (is.data.frame(gdf) && !is.null(rownames(gdf))) {
        coln <- if ("group" %in% colnames(gdf)) "group" else if ("Group" %in% colnames(gdf)) "Group" else {
          chs <- sapply(gdf, function(x) is.character(x) || is.factor(x))
          if (any(chs)) names(chs)[which(chs)[1]] else NULL
        }
        if (!is.null(coln) && all(colnames(exp_data) %in% rownames(gdf))) {
          return(factor(as.character(gdf[colnames(exp_data), coln])))
        }
      }
    }
    cn <- colnames(exp_data)
    prefixes <- sub("^([^_\\.-]+)[_\\.-].*$", "\\1", cn)
    if (all(prefixes == cn)) prefixes <- sub("([A-Za-z]+)\\d*$", "\\1", cn)
    if (length(unique(prefixes)) > 1) return(factor(prefixes))
    return(NULL)
  }
  if (!by_group){
    pvals <- apply(exp_data, 1, function(x){
      x <- as.numeric(x)
      if (all(is.na(x)) || sum(!is.na(x)) < 3) return(NA_real_)
      res <- tryCatch(stats::shapiro.test(x), error = function(e) NULL)
      if (is.null(res)) NA_real_ else as.numeric(res$p.value)
    })
    df <- data.frame(gene = rownames(exp_data), p.value = as.numeric(pvals), stringsAsFactors = FALSE)
    df$adj.p <- p.adjust(df$p.value, method = "BH")
    df <- df[order(df$adj.p), ]
    return(df)
  } else {
    group <- infer_group()
    if (is.null(group)) stop("")
    res_list <- list()
    groups <- levels(group)
    for (g in groups){
      cols <- which(group == g)
      if (length(cols) < 3) {
        tmp <- data.frame(gene = rownames(exp_data), group = g, p.value = NA_real_, stringsAsFactors = FALSE)
      } else {
        pvals <- apply(exp_data[, cols, drop = FALSE], 1, function(x){
          x <- as.numeric(x)
          if (all(is.na(x)) || sum(!is.na(x)) < 3) return(NA_real_)
          res <- tryCatch(stats::shapiro.test(x), error = function(e) NULL)
          if (is.null(res)) NA_real_ else as.numeric(res$p.value)
        })
        tmp <- data.frame(gene = rownames(exp_data), group = g, p.value = as.numeric(pvals), stringsAsFactors = FALSE)
        tmp$adj.p <- p.adjust(tmp$p.value, method = "BH")
      }
      res_list[[g]] <- tmp
    }
    df_all <- do.call(rbind, lapply(res_list, function(x){
      if (!"adj.p" %in% colnames(x)) x$adj.p <- NA_real_
      x
    }))
    rownames(df_all) <- NULL
    return(df_all[order(df_all$group, df_all$adj.p), ])
  }
}

norm_res_overall <- exp_normality_test_genes(exp_data, by_group = FALSE)
write.csv(norm_res_overall, file.path(output, paste0(labels, "normality_shapiro_overall.csv")))


cor_matrix <- cor(t(exp_data), use = "pairwise.complete.obs")
p_matrix <- cor.mtest(cor_matrix)$p

genes <- colnames(cor_matrix)

# 2) 转成长表
cor_df <- as.data.frame(as.table(cor_matrix), stringsAsFactors = FALSE) %>%
  rename(gene_y = Var1, gene_x = Var2, r = Freq)

p_df <- as.data.frame(as.table(p_matrix), stringsAsFactors = FALSE) %>%
  rename(gene_y = Var1, gene_x = Var2, p = Freq)

plot_df <- cor_df %>%
  left_join(p_df, by = c("gene_y", "gene_x")) %>%
  mutate(
    sig = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ ""
    ),
    # 可选：对角线不显示文字
    label = ifelse(gene_x == gene_y, sprintf("%.2f", r), sprintf("%.2f%s", r, sig))
  )

# 固定顺序
plot_df$gene_x <- factor(plot_df$gene_x, levels = genes)
plot_df$gene_y <- factor(plot_df$gene_y, levels = rev(genes))

# 3) 完整热图（不再过滤 upper）
p <- ggplot(plot_df, aes(x = gene_x, y = gene_y, fill = r)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = label), size = 2.8, color = "black") +
  scale_fill_gradient2(
    low = "#8491B4", mid = "white", high = "#F39B7F",
    midpoint = 0, limits = c(-1, 1), name = "Correlation"
  ) +
  scale_x_discrete(position = "top") +
  coord_fixed() +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    # 基因名斜体
    axis.text.x.top = element_text(
      face = "italic", color = "black", angle = 45, hjust = 0, vjust = 0,
      margin = margin(b = 0)   # 上标签下移
    ),
    axis.text.y.left = element_text(
      face = "italic", color = "black",
      margin = margin(r = 0)   # 左标签右移
    )
  )

ggsave(
  filename = file.path(output, paste0(labels, "_corrplot_ggplot2_full.pdf")),
  plot = p, width = 4.5, height = 4.5, dpi = 300
)
ggsave(
  filename = file.path(output, paste0(labels, "_corrplot_ggplot2_full.png")),
  plot = p, width = 4.5, height = 4.5, dpi = 300
)



progression <- predict(model1,scale(t(exp_corrected[colnames(model1_gene),])))
survival <- predict(model2,scale(t(exp_corrected[colnames(model2_gene),])))
all_feature
X <- as.data.frame(t(exp_corrected[all_feature,]))
X <- scale(X)

X <- as.data.frame(X)
X$progression <- progression
X$survival <- survival
model_a <- lm(progression ~ NPC2 + TLN1 + TUBA1C + LRRC32 + PRB2, data = X)

model_bc <- lm(survival ~ NPC2 + TLN1 + TUBA1C + LRRC32 + PRB2 +  SOX9 + RNASE4 + SERPINA3 + progression, data = X)
treat_vars <- c("NPC2", "TLN1", "TUBA1C", "LRRC32", "PRB2")

results <- list()

for (treat in treat_vars) {
  med_result <- mediate(model.m = model_a,  
                        model.y = model_bc,
                        treat = treat,
                        mediator = "progression",
                        boot = TRUE,
                        sims = 1000)
  results[[treat]] <- summary(med_result)
}

results



#lavaan
model <- '
  # 
  progression ~ label("a_NPC2")*NPC2 + label("a_TLN1")*TLN1 + label("a_TUBA1C")*TUBA1C + label("a_LRRC32")*LRRC32 + label("a_PRB2")*PRB2
  
  # 
  SOX9 ~ label("b_SOX9")*progression
  RNASE4 ~ label("b_RNASE4")*progression
  SERPINA3 ~ label("b_SERPINA3")*progression
  
  # 
  survival ~ label("c_SOX9")*SOX9 + label("c_RNASE4")*RNASE4 + label("c_SERPINA3")*SERPINA3 + label("c_progression")*progression
  
  
  # 
  # NPC2 
  indirect_NPC2_SOX9 := a_NPC2 * b_SOX9 * c_SOX9
  indirect_NPC2_RNASE4 := a_NPC2 * b_RNASE4 * c_RNASE4
  indirect_NPC2_SERPINA3 := a_NPC2 * b_SERPINA3 * c_SERPINA3
  
  # TLN1 
  indirect_TLN1_SOX9 := a_TLN1 * b_SOX9 * c_SOX9
  indirect_TLN1_RNASE4 := a_TLN1 * b_RNASE4 * c_RNASE4
  indirect_TLN1_SERPINA3 := a_TLN1 * b_SERPINA3 * c_SERPINA3
  
  # TUBA1C 
  indirect_TUBA1C_SOX9 := a_TUBA1C * b_SOX9 * c_SOX9
  indirect_TUBA1C_RNASE4 := a_TUBA1C * b_RNASE4 * c_RNASE4
  indirect_TUBA1C_SERPINA3 := a_TUBA1C * b_SERPINA3 * c_SERPINA3
  
  # LRRC32 
  indirect_LRRC32_SOX9 := a_LRRC32 * b_SOX9 * c_SOX9
  indirect_LRRC32_RNASE4 := a_LRRC32 * b_RNASE4 * c_RNASE4
  indirect_LRRC32_SERPINA3 := a_LRRC32 * b_SERPINA3 * c_SERPINA3
  
  # PRB2 
  indirect_PRB2_SOX9 := a_PRB2 * b_SOX9 * c_SOX9
  indirect_PRB2_RNASE4 := a_PRB2 * b_RNASE4 * c_RNASE4
  indirect_PRB2_SERPINA3 := a_PRB2 * b_SERPINA3 * c_SERPINA3
'
set.seed(123)
fit <- sem(model, data = X, se = "bootstrap", bootstrap = 2000, missing = "FIML")
summary(fit, standardized = TRUE, rsquare = TRUE)
fit_results <- summary(fit, standardized = TRUE, rsquare = TRUE)
results <- fit_results$pe
write.csv(results,file.path(output,paste0(labels,"sem_results.csv")))
fitmeasures(fit, c("cfi","tli","rmsea","srmr"))
param <- parameterEstimates(fit, boot.ci.type = "perc", level = 0.95)
indirects <- param %>% filter(grepl("^indirect_", label))
indirects <- indirects %>% mutate(p_adj = p.adjust(pvalue, method = "fdr"))
print(indirects)
pe <- parameterEstimates(fit, standardized = TRUE) %>%
  filter(op == "~") %>%
  transmute(
    from = rhs,         # 自变量
    to   = lhs,         # 因变量
    est  = std.all,     # 标准化系数
    pval = pvalue
  )
gene_nodes <- c("NPC2", "TLN1", "TUBA1C", "LRRC32", "PRB2", "SOX9", "RNASE4", "SERPINA3")

# 为每个节点生成可解析标签：基因用 italic()，非基因保持普通
nodes <- nodes %>%
  mutate(
    label_expr = ifelse(
      node %in% gene_nodes,
      paste0("italic('", node, "')"),
      paste0("'", node, "'")
    )
  )

# 节点文字
geom_text(
  data = nodes,
  aes(x = x, y = y, label = label_expr),
  parse = TRUE,
  size = 4.3,
  fontface = "bold",
  color = "#1F1F1F"
)
# 2) 定义节点分层与坐标（每行居中）
nodes <- bind_rows(
  tibble(node = c("NPC2", "TLN1", "TUBA1C", "LRRC32", "PRB2"),
         x = c(-4, -2, 0, 2, 4), y = 4),
  tibble(node = "progression", x = 0, y = 3),
  tibble(node = c("SOX9", "RNASE4", "SERPINA3"),
         x = c(-2, 0, 2), y = 2),
  tibble(node = "survival", x = 0, y = 1)
)

# 3) 连接边并计算绘图参数
edges <- pe %>%
  inner_join(nodes %>% rename(from = node, x = x, y = y), by = "from") %>%
  inner_join(nodes %>% rename(to = node, xend = x, yend = y), by = "to") %>%
  mutate(
    dx = xend - x,
    # 弧线方向：右偏正弯，左偏负弯，垂直连线给固定弯曲
    curve_type = case_when(
      dx > 0 ~ "right",
      dx < 0 ~ "left",
      TRUE   ~ "vertical"
    ),
    # 给边标签一个轻微偏移，避免压在线上
    xm = (x + xend) / 2 + ifelse(dx == 0, 0.18, 0),
    ym = (y + yend) / 2 + ifelse(dx == 0, 0, 0.12 * sign(dx))
  )
gene_nodes <- c("NPC2", "TLN1", "TUBA1C", "LRRC32", "PRB2", "SOX9", "RNASE4", "SERPINA3")

# 为每个节点生成可解析标签：基因用 italic()，非基因保持普通
nodes <- nodes %>%
  mutate(
    label_expr = ifelse(
      node %in% gene_nodes,
      paste0("italic('", node, "')"),
      paste0("'", node, "'")
    )
  )


# 4) 绘图
p <- ggplot() +
  # 右弯弧线
  geom_curve(
    data = edges %>% filter(curve_type == "right"),
    aes(x = x, y = y, xend = xend, yend = yend, color = est, linewidth = abs(est)),
    curvature = 0.25,
    alpha = 0.9,
    lineend = "round",
    arrow = arrow(length = unit(0.18, "cm"), type = "closed")
  ) +
  # 左弯弧线
  geom_curve(
    data = edges %>% filter(curve_type == "left"),
    aes(x = x, y = y, xend = xend, yend = yend, color = est, linewidth = abs(est)),
    curvature = -0.25,
    alpha = 0.9,
    lineend = "round",
    arrow = arrow(length = unit(0.18, "cm"), type = "closed")
  ) +
  # 垂直连线也用弧线表示
  geom_curve(
    data = edges %>% filter(curve_type == "vertical"),
    aes(x = x, y = y, xend = xend, yend = yend, color = est, linewidth = abs(est)),
    curvature = 0.20,
    alpha = 0.9,
    lineend = "round",
    arrow = arrow(length = unit(0.18, "cm"), type = "closed")
  ) +
  # 节点：半透明圆，保证文字清晰
  geom_point(
    data = nodes,
    aes(x = x, y = y),
    shape = 21, size = 14, stroke = 0.9,
    fill = alpha("#F2F2F2", 0.45), color = "#4D4D4D"
  )+ geom_label_repel(
    data = edges,
    aes(x = xm, y = ym, label = sprintf("%.2f", est)),
    size = 3.0,
    color = "black",
    fill = "#F7F7F7",      # 统一浅色背景
    alpha = 0.95,          # 轻微透明，避免太“硬”
    label.size = 0.15,
    label.r = unit(0.08, "lines"),
    box.padding = 0.12,
    point.padding = 0.05,
    min.segment.length = 0,
    segment.color = NA,
    segment.size = 0.25,
    force = 1.2,
    max.overlaps = Inf,
    show.legend = FALSE
  )+
  geom_text(
    data = nodes,
    aes(x = x, y = y, label = label_expr),
    parse = TRUE,
    size = 4.3,
    fontface = "bold",
    color = "#1F1F1F"
  ) +
  scale_color_gradient2(
    low = "#2C7BB6", mid = "#BDBDBD", high = "#D7191C", midpoint = 0,
    name = "Standardized coefficient"
  ) +
  scale_linewidth(
    range = c(0.7, 2.8),
    name = "|Standardized coefficient|"
  )+
  coord_cartesian(xlim = c(-5, 5), ylim = c(0.6, 4.4), clip = "off") +
  theme_void(base_size = 13) +
  theme(
    legend.position = "right",
    plot.margin = margin(10, 20, 10, 20)
  )
p
ggsave(
  filename = file.path(output, paste0(labels, "_semplot_ggplot2.png")),
  plot = p, width = 7.5, height = 6, dpi = 300
)
ggsave(
  filename = file.path(output, paste0(labels, "_semplot_ggplot2.pdf")),
  plot = p, width = 7.5, height = 6, dpi = 300
)

all_feature
param <- parameterEstimates(fit)
indirect_effects <- param[grepl("indirect_", param$label), ]
significant_indirect <- indirect_effects[indirect_effects$pvalue < 0.05, ]
write.csv(indirect_effects,file.path(output,paste0(labels,"_indirect_effects.csv")))
write.csv(significant_indirect,file.path(output,paste0(labels,"_significant_indirect.csv")))
