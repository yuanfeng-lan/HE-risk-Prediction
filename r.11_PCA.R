rm(list = ls()); gc()
ORIGINAL_DIR <- "~/data_hdd/R/Yclub/2026.3.4_HE_Cirrhosis"
output <- file.path(ORIGINAL_DIR, "11_PCA")
labels <- "11_PCA_"

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(ORIGINAL_DIR)
library(scatterplot3d)
library(ggplot2)
library(dendextend)
library(tidyr)
library(dplyr)
library(ggdendro)
####GSE139602####
GSE <- "GSE139602"
outdir <- file.path(output,GSE)
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

exp <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE139602_exp.csv"),row.names = 1)
group <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE139602_group.csv"),row.names = 1)
group <- group[colnames(exp),,drop=F]
table(group$characteristics_ch1)
group$group <- ifelse(group$characteristics_ch1=="disease state: Healthy",0,
                      ifelse(group$characteristics_ch1=="disease state: eCLD",1,
                             ifelse(group$characteristics_ch1=="disease state: Compensated Cirrhosis",2,
                                    ifelse(group$characteristics_ch1=="disease state: Decompesated Cirrhosis",3,
                                           ifelse(group$characteristics_ch1=="disease state: Acute-on-chronic liver failure",4,NA)))))
group$type <- ifelse(group$characteristics_ch1=="disease state: Healthy","Healthy",
                     ifelse(group$characteristics_ch1=="disease state: eCLD","eCLD",
                            ifelse(group$characteristics_ch1=="disease state: Compensated Cirrhosis","CC",
                                   ifelse(group$characteristics_ch1=="disease state: Decompesated Cirrhosis","DC",
                                          ifelse(group$characteristics_ch1=="disease state: Acute-on-chronic liver failure","ACLF",NA)))))
group$type <- factor(group$type,levels = c("Healthy","eCLD","CC","DC","ACLF"))
exp <- scale(exp)
exp_data <- as.data.frame(t(exp))
group_data <- group[,2,drop=F]
pca <- prcomp(exp_data, center = FALSE, scale. = FALSE)
scores <- as.data.frame(pca$x[, 1:3])
colnames(scores) <- c("PC1","PC2","PC3")
explained <- round(100 * (pca$sdev^2) / sum(pca$sdev^2), 1)
groups <- group$type
stopifnot(length(groups) == nrow(scores))
pal <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F")  
cols <- pal[as.integer(groups)]
cols <- as.character(cols)

png(file.path(outdir,paste0(labels,GSE,"_pca3d_plot.png")), width = 5, height = 5,units = "in",res = 300)
s3d <- scatterplot3d(scores$PC1, scores$PC2, scores$PC3,
                     color = cols, pch = 19, cex.symbols = 1.2,
                     xlab = paste0("PC1 (", explained[1], "%)"),
                     ylab = paste0("PC2 (", explained[2], "%)"),
                     zlab = paste0("PC3 (", explained[3], "%)"),
                     main = paste0(GSE," PCA 3D "))
op <- par(xpd = NA)
legend("topright",
       legend = levels(groups),
       col = c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F"),
       pch = 19, bty = "n",
       inset = c(-0.12, 0))
par(op)

dev.off()
pdf(file.path(outdir,paste0(labels,GSE,"_pca3d_plot.pdf")), width = 5, height = 5)
s3d <- scatterplot3d(scores$PC1, scores$PC2, scores$PC3,
                     color = cols, pch = 19, cex.symbols = 1.2,
                     xlab = paste0("PC1 (", explained[1], "%)"),
                     ylab = paste0("PC2 (", explained[2], "%)"),
                     zlab = paste0("PC3 (", explained[3], "%)"),
                     main = paste0(GSE," PCA 3D "))
op <- par(xpd = NA)
legend("topright",
       legend = levels(groups),
       col = c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F"),
       pch = 19, bty = "n",
       inset = c(-0.12, 0))
par(op)
dev.off()

#

hc <- hclust(dist(exp_data)) 
known_groups <- group$type
names(known_groups) <- rownames(group)
dend <- as.dendrogram(hc)
group_colors <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F")
sample_colors <- group_colors[known_groups[order.dendrogram(dend)]]
dend_colored <- dend %>%
  set("labels_col", sample_colors) %>%
  set("labels_cex", 0.8) %>%
  set("branches_lwd", 2)
png(file.path(outdir,paste0(labels,GSE,"_hclust_plot.png")), width = 6, height = 5,units = "in",res = 300)
plot(dend_colored, main = paste0(GSE," Cluster Dendrogram"))
unique_groups <- unique(known_groups)
legend("topright", legend = unique_groups, fill = group_colors, 
       title = "Groups", cex = 0.8)
dev.off()
pdf(file.path(outdir,paste0(labels,GSE,"_hclust_plot.pdf")), width = 6, height = 5)
plot(dend_colored, main = paste0(GSE," Cluster Dendrogram"))
unique_groups <- unique(known_groups)
legend("topright", legend = unique_groups, fill = group_colors, 
       title = "Groups", cex = 0.8)
dev.off()

####GSE15654####

GSE <- "GSE15654"
outdir <- file.path(output,GSE)
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

exp = read.csv(file = './00_rawdata/00.rawdata_GSE15654_exp.csv', header = TRUE,row.names = 1) 
group = read.csv(file = './00_rawdata/00.rawdata_GSE15654_group.csv',  header = TRUE,row.names = 1) 
exp <- exp[,rownames(group)]
exp <- scale(exp)
group$type <- factor(group$group,levels = c("Good_prognosis","Poor_prognosis"))

exp_data <- as.data.frame(t(exp))
group_data <- group[,"type",drop=F]
pca <- prcomp(exp_data, center = FALSE, scale. = FALSE)
scores <- as.data.frame(pca$x[, 1:3])
colnames(scores) <- c("PC1","PC2","PC3")
explained <- round(100 * (pca$sdev^2) / sum(pca$sdev^2), 1)
groups <- group$type
stopifnot(length(groups) == nrow(scores))
pal <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F")  
cols <- pal[as.integer(groups)]
cols <- as.character(cols)

png(file.path(outdir,paste0(labels,GSE,"_pca3d_plot.png")), width = 5, height = 5,units = "in",res = 300)
s3d <- scatterplot3d(scores$PC1, scores$PC2, scores$PC3,
                     color = cols, pch = 19, cex.symbols = 1.2,
                     xlab = paste0("PC1 (", explained[1], "%)"),
                     ylab = paste0("PC2 (", explained[2], "%)"),
                     zlab = paste0("PC3 (", explained[3], "%)"),
                     main = paste0(GSE," PCA 3D "))
op <- par(xpd = NA)
legend("topright",
       legend = levels(groups),
       col = c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F"),
       pch = 19, bty = "n",
       inset = c(0, -0.12))
par(op)

dev.off()
pdf(file.path(outdir,paste0(labels,GSE,"_pca3d_plot.pdf")), width = 5, height = 5)
s3d <- scatterplot3d(scores$PC1, scores$PC2, scores$PC3,
                     color = cols, pch = 19, cex.symbols = 1.2,
                     xlab = paste0("PC1 (", explained[1], "%)"),
                     ylab = paste0("PC2 (", explained[2], "%)"),
                     zlab = paste0("PC3 (", explained[3], "%)"),
                     main = paste0(GSE," PCA 3D "))
op <- par(xpd = NA)
legend("topright",
       legend = levels(groups),
       col = c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F"),
       pch = 19, bty = "n",
       inset = c(0, -0.12))
par(op)
dev.off()

#

hc <- hclust(dist(exp_data)) 
known_groups <- group$type
names(known_groups) <- rownames(group)
dend <- as.dendrogram(hc)
group_colors <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F")
sample_colors <- group_colors[known_groups[order.dendrogram(dend)]]
dend_colored <- dend %>%
  set("labels_col", sample_colors) %>%
  set("labels_cex", 0.2) %>%
  set("branches_lwd", 1)
png(file.path(outdir,paste0(labels,GSE,"_hclust_plot.png")), width = 6, height = 5,units = "in",res = 300)
plot(dend_colored, main = paste0(GSE," Cluster Dendrogram"))
unique_groups <- unique(known_groups)
legend("topright", legend = unique_groups, fill = group_colors, 
       title = "Groups", cex = 0.8)
dev.off()
pdf(file.path(outdir,paste0(labels,GSE,"_hclust_plot.pdf")), width = 6, height = 5)
plot(dend_colored, main = paste0(GSE," Cluster Dendrogram"))
unique_groups <- unique(known_groups)
legend("topright", legend = unique_groups, fill = group_colors, 
       title = "Groups", cex = 0.8)
dev.off()

####GSE41919####
GSE <- "GSE41919"
outdir <- file.path(output,GSE)
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

exp = read.csv(file = './00_rawdata/00.rawdata_GSE41919_exp.csv', header = TRUE,row.names = 1) 
group = read.csv(file = './00_rawdata/00.rawdata_GSE41919_group.csv',  header = TRUE,row.names = 1) 
group <- subset(group,group!="non-cirrhotic control")
exp <- exp[,rownames(group)]
exp <- na.omit(exp)
group$type <- factor(group$group,levels = c("cirrhosis with HE","cirrhosis without HE"))
exp <- scale(exp)
exp_data <- as.data.frame(t(exp))
group_data <- group[,"type",drop=F]
pca <- prcomp(exp_data, center = FALSE, scale. = FALSE)
scores <- as.data.frame(pca$x[, 1:3])
colnames(scores) <- c("PC1","PC2","PC3")
explained <- round(100 * (pca$sdev^2) / sum(pca$sdev^2), 1)
groups <- group$type
stopifnot(length(groups) == nrow(scores))
pal <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F")  
cols <- pal[as.integer(groups)]
cols <- as.character(cols)

png(file.path(outdir,paste0(labels,GSE,"_pca3d_plot.png")), width = 5, height = 5,units = "in",res = 300)
s3d <- scatterplot3d(scores$PC1, scores$PC2, scores$PC3,
                     color = cols, pch = 19, cex.symbols = 1.2,
                     xlab = paste0("PC1 (", explained[1], "%)"),
                     ylab = paste0("PC2 (", explained[2], "%)"),
                     zlab = paste0("PC3 (", explained[3], "%)"),
                     main = paste0(GSE," PCA 3D "))
op <- par(xpd = NA)
legend("topright",
       legend = levels(groups),
       col = c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F"),
       pch = 19, bty = "n",
       inset = c(-0.12, -0.14))
par(op)

dev.off()
pdf(file.path(outdir,paste0(labels,GSE,"_pca3d_plot.pdf")), width = 5, height = 5)
s3d <- scatterplot3d(scores$PC1, scores$PC2, scores$PC3,
                     color = cols, pch = 19, cex.symbols = 1.2,
                     xlab = paste0("PC1 (", explained[1], "%)"),
                     ylab = paste0("PC2 (", explained[2], "%)"),
                     zlab = paste0("PC3 (", explained[3], "%)"),
                     main = paste0(GSE," PCA 3D "))
op <- par(xpd = NA)
legend("topright",
       legend = levels(groups),
       col = c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F"),
       pch = 19, bty = "n",
       inset = c(-0.12, -0.14))
par(op)
dev.off()

#

hc <- hclust(dist(exp_data)) 
known_groups <- group$type
names(known_groups) <- rownames(group)
dend <- as.dendrogram(hc)
group_colors <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F")
sample_colors <- group_colors[known_groups[order.dendrogram(dend)]]
dend_colored <- dend %>%
  set("labels_col", sample_colors) %>%
  set("labels_cex", 0.8) %>%
  set("branches_lwd", 2)
png(file.path(outdir,paste0(labels,GSE,"_hclust_plot.png")), width = 6, height = 5,units = "in",res = 300)
plot(dend_colored, main = paste0(GSE," Cluster Dendrogram"))
unique_groups <- unique(known_groups)
legend("topleft", legend = unique_groups, fill = group_colors, 
       title = "Groups", cex = 0.8)
dev.off()
pdf(file.path(outdir,paste0(labels,GSE,"_hclust_plot.pdf")), width = 6, height = 5)
plot(dend_colored, main = paste0(GSE," Cluster Dendrogram"))
unique_groups <- unique(known_groups)
legend("topleft", legend = unique_groups, fill = group_colors, 
       title = "Groups", cex = 0.8)
dev.off()

####GSE57193####

GSE <- "GSE57193"
outdir <- file.path(output,GSE)
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

exp = read.csv(file = './00_rawdata/00.rawdata_GSE57193_exp.csv', header = TRUE,row.names = 1) 
group = read.csv(file = './00_rawdata/00.rawdata_GSE57193_group.csv',  header = TRUE,row.names = 1) 
group <- subset(group,group!="healthy")
exp <- exp[,rownames(group)]
exp <- na.omit(exp)
group$group <- ifelse(group$group=="cirrhosis","cirrhosis without HE","cirrhosis with HE")
group$type <- factor(group$group,levels = c("cirrhosis with HE","cirrhosis without HE"))
exp <- scale(exp)
exp_data <- as.data.frame(t(exp))
group_data <- group[,"type",drop=F]
pca <- prcomp(exp_data, center = FALSE, scale. = FALSE)
scores <- as.data.frame(pca$x[, 1:3])
colnames(scores) <- c("PC1","PC2","PC3")
explained <- round(100 * (pca$sdev^2) / sum(pca$sdev^2), 1)
groups <- group$type
stopifnot(length(groups) == nrow(scores))
pal <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F")  
cols <- pal[as.integer(groups)]
cols <- as.character(cols)

png(file.path(outdir,paste0(labels,GSE,"_pca3d_plot.png")), width = 5, height = 5,units = "in",res = 300)
s3d <- scatterplot3d(scores$PC1, scores$PC2, scores$PC3,
                     color = cols, pch = 19, cex.symbols = 1.2,
                     xlab = paste0("PC1 (", explained[1], "%)"),
                     ylab = paste0("PC2 (", explained[2], "%)"),
                     zlab = paste0("PC3 (", explained[3], "%)"),
                     main = paste0(GSE," PCA 3D "))
op <- par(xpd = NA)
legend("topright",
       legend = levels(groups),
       col = c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F"),
       pch = 19, bty = "n",
       inset = c(-0.12, -0.14))
par(op)

dev.off()
pdf(file.path(outdir,paste0(labels,GSE,"_pca3d_plot.pdf")), width = 5, height = 5)
s3d <- scatterplot3d(scores$PC1, scores$PC2, scores$PC3,
                     color = cols, pch = 19, cex.symbols = 1.2,
                     xlab = paste0("PC1 (", explained[1], "%)"),
                     ylab = paste0("PC2 (", explained[2], "%)"),
                     zlab = paste0("PC3 (", explained[3], "%)"),
                     main = paste0(GSE," PCA 3D "))
op <- par(xpd = NA)
legend("topright",
       legend = levels(groups),
       col = c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F"),
       pch = 19, bty = "n",
       inset = c(-0.12, -0.14))
par(op)
dev.off()

#

hc <- hclust(dist(exp_data)) 
known_groups <- group$type
names(known_groups) <- rownames(group)
dend <- as.dendrogram(hc)
group_colors <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F")
sample_colors <- group_colors[known_groups[order.dendrogram(dend)]]
dend_colored <- dend %>%
  set("labels_col", sample_colors) %>%
  set("labels_cex", 0.8) %>%
  set("branches_lwd", 2)
png(file.path(outdir,paste0(labels,GSE,"_hclust_plot.png")), width = 6, height = 5,units = "in",res = 300)
plot(dend_colored, main = paste0(GSE," Cluster Dendrogram"))
unique_groups <- unique(known_groups)
legend(x = 5.8,y = 42, legend = unique_groups, fill = group_colors, 
       title = "Groups", cex = 0.8)
dev.off()
pdf(file.path(outdir,paste0(labels,GSE,"_hclust_plot.pdf")), width = 6, height = 5)
plot(dend_colored, main = paste0(GSE," Cluster Dendrogram"))
unique_groups <- unique(known_groups)
legend(x = 5.8,y = 42, legend = unique_groups, fill = group_colors, 
       title = "Groups", cex = 0.8)
dev.off()

