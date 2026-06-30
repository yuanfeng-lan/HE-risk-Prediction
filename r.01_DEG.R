rm(list = ls()); gc()
ORIGINAL_DIR <- "~/data_hdd/R/Yclub/2026.3.4_HE_Cirrhosis"
output <- file.path(ORIGINAL_DIR, "01_DEG")
labels <- "00_rawdata_"
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(ORIGINAL_DIR)

library(limma)
library(ggplot2)
library(ggrepel)
library("ggvenn")
####GSE41919####

exp = read.csv(file = './00_rawdata/00.rawdata_GSE41919_exp.csv', header = TRUE,row.names = 1) 
group = read.csv(file = './00_rawdata/00.rawdata_GSE41919_group.csv',  header = TRUE,row.names = 1) 
group <- subset(group,group!="non-cirrhotic control")
exp <- exp[,rownames(group)]
library(limma)
design <- model.matrix(~0+factor(group$group))
row.names(design) <- rownames(group)
colnames(design)=c('with_HE','without_HE')
contrast.matrix<-makeContrasts("with_HE-without_HE",levels=design)
set.seed(123)
##step1
fit <- lmFit(exp,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2) 
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput)
write.csv(nrDEG, file = "./01_DEG/01.DEG_GSE41919_with_HE-without_HE_all.csv")
nrDEG_1 <- subset(nrDEG,abs(logFC)>0.5&P.Value<0.05)
gene01 <- rownames(nrDEG_1)
#valcano_plot
nrDEG$change = ifelse(nrDEG$P.Value < 0.05 & abs(nrDEG$logFC) >= 0.5, 
                      ifelse(nrDEG$logFC> 0 ,'Up','Down'),
                      'NS')
table(nrDEG$change)
write.csv(nrDEG, file = "./01_DEG/01.DEG_GSE41919_with_HE-without_HE_all_2.csv")

nrDEG$Label = ""   
nrDEG <- nrDEG[order(nrDEG$P.Value), ]   
nrDEG$Gene <- rownames(nrDEG)
#up.genes <- head(nrDEG$Gene[which(nrDEG$change == "Up")], 5)
#down.genes <- head(nrDEG$Gene[which(nrDEG$change == "Down")], 5)
#nrDEG.top5.genes <- c(as.character(up.genes), as.character(down.genes))
#nrDEG$Label[match(nrDEG.top5.genes, nrDEG$Gene)] <- nrDEG.top5.genes
p1 <- ggplot(
  nrDEG, 
  aes(x = logFC, 
      y = -log10(P.Value), 
      colour=change)) +
  geom_point(alpha=0.4, size=2) +
  scale_color_manual(values=c("#4DBBD5", "#d2dae2","#E64B35"))+
  geom_vline(xintercept=0,lty=4,col="black",lwd=0.8) +
  geom_vline(xintercept=0.5,lty=2,col="grey50",lwd=0.4) +
  geom_vline(xintercept=-0.5,lty=2,col="grey50",lwd=0.4) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(title = "GSE41919",
       # subtitle = "Cirrhotic Patients with vs. without HE",
       x = expression(log[2]~FC),  # x轴标签，log2中的2为下标
       y = expression(-log[10]~italic(p))  )+
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  
    plot.subtitle = element_text(hjust = 0.5, size = 14), 
    axis.title.x = element_text(size = 12), 
    axis.title.y = element_text(size = 12),  
    axis.text = element_text(size = 12),  
    legend.position="right", 
    legend.text = element_text(size = 12),
    legend.title = element_blank()
  )+geom_text_repel(data = nrDEG, aes(label = Label),
                    size = 3.5,                           
                    box.padding = unit(0.5, "lines"),  
                    point.padding = unit(0.8, "lines"), 
                    segment.color = "black",            
                    show.legend = FALSE,  
                    max.overlaps = 10000) 
ggsave("./01_DEG/01.DEG_GSE41919_DEG_valcano.png",plot = p1,width = 5, height = 5, units = "in")
ggsave("./01_DEG/01.DEG_GSE41919_DEG_valcano.pdf",plot = p1,width = 5, height = 5, units = "in")

####GSE57193####
exp = read.csv(file = './00_rawdata/00.rawdata_GSE57193_exp.csv', header = TRUE,row.names = 1) 
group = read.csv(file = './00_rawdata/00.rawdata_GSE57193_group.csv',  header = TRUE,row.names = 1) 
group <- subset(group,group!="healthy")
exp <- exp[,rownames(group)]
design <- model.matrix(~0+factor(group$group))
row.names(design) <- rownames(group)
colnames(design)=c('cirrhosis','cirrhosis_with_HE')
contrast.matrix<-makeContrasts("cirrhosis_with_HE-cirrhosis",levels=design)
set.seed(123)
##step1
fit <- lmFit(exp,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2) 
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput)
write.csv(nrDEG, file = "./01_DEG/01.DEG_GSE57193_cirrhosis_with_HE-cirrhosis_all.csv")
nrDEG_1 <- subset(nrDEG,abs(logFC)>0.5&P.Value<0.05)
gene01 <- rownames(nrDEG_1)
nrDEG$change = ifelse(nrDEG$P.Value < 0.05 & abs(nrDEG$logFC) >= 0.5, 
                      ifelse(nrDEG$logFC> 0 ,'Up','Down'),
                      'NS')

table(nrDEG$change)
write.csv(nrDEG, file = "./01_DEG/01.DEG_GSE57193_cirrhosis_with_HE-cirrhosis_all_2.csv")

nrDEG$Label = ""  
nrDEG <- nrDEG[order(nrDEG$P.Value), ]   
nrDEG$Gene <- rownames(nrDEG)

#up.genes <- head(nrDEG$Gene[which(nrDEG$change == "Up")], 5)
#down.genes <- head(nrDEG$Gene[which(nrDEG$change == "Down")], 5)
#nrDEG.top5.genes <- c(as.character(up.genes), as.character(down.genes))
#nrDEG$Label[match(nrDEG.top5.genes, nrDEG$Gene)] <- nrDEG.top5.genes
p1 <- ggplot(
  nrDEG, 
  aes(x = logFC, 
      y = -log10(P.Value), 
      colour=change)) +
  geom_point(alpha=0.4, size=2) +
  scale_color_manual(values=c("#4DBBD5", "#d2dae2","#E64B35"))+
  geom_vline(xintercept=0,lty=4,col="black",lwd=0.8) +
  geom_vline(xintercept=0.5,lty=2,col="grey50",lwd=0.4) +
  geom_vline(xintercept=-0.5,lty=2,col="grey50",lwd=0.4) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(title = "GSE57193",
       # subtitle = "Cirrhotic Patients with vs. without HE",
       x = expression(log[2]~FC),  # x轴标签，log2中的2为下标
       y = expression(-log[10]~italic(p))  )+
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  
    plot.subtitle = element_text(hjust = 0.5, size = 14), 
    axis.title.x = element_text(size = 12), 
    axis.title.y = element_text(size = 12),  
    axis.text = element_text(size = 12),  
    legend.position="right", 
    legend.text = element_text(size = 12),
    legend.title = element_blank()
  )+geom_text_repel(data = nrDEG, aes(label = Label),
                    size = 3.5,                           
                    box.padding = unit(0.5, "lines"),  
                    point.padding = unit(0.8, "lines"), 
                    segment.color = "black",            
                    show.legend = FALSE,  
                    max.overlaps = 10000) 
ggsave("./01_DEG/01.DEG_GSE57193_DEG_valcano.png",plot = p1,width = 5, height = 5, units = "in")
ggsave("./01_DEG/01.DEG_GSE57193_DEG_valcano.pdf",plot = p1,width = 5, height = 5, units = "in")

####GSE15654####
GSE15654_exp = read.csv(file = './00_rawdata/00.rawdata_GSE15654_exp.csv', header = TRUE,row.names = 1) 
GSE15654_group = read.csv(file = './00_rawdata/00.rawdata_GSE15654_group.csv',  header = TRUE,row.names = 1) 
GSE15654_exp <- GSE15654_exp[,rownames(GSE15654_group)]
library(limma)
GSE15654_exp <- scale(GSE15654_exp)
design <- model.matrix(~0+factor(GSE15654_group$group))
row.names(design) <- rownames(GSE15654_group)
colnames(design)=c('Good_prognosis','Poor_prognosis')
contrast.matrix<-makeContrasts("Poor_prognosis-Good_prognosis",levels=design)
set.seed(123)
##step1
fit <- lmFit(GSE15654_exp,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2) 
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput)
write.csv(nrDEG, file = "./01_DEG/01.DEG_GSE15654_Poor_prognosis-Good_prognosis_all.csv")
nrDEG_1 <- subset(nrDEG,abs(logFC)>0.5&P.Value<0.05)
gene01 <- rownames(nrDEG_1)
nrDEG$change = ifelse(nrDEG$P.Value < 0.05 & abs(nrDEG$logFC) >= 0.5, 
                      ifelse(nrDEG$logFC> 0 ,'Up','Down'),
                      'NS')

table(nrDEG$change)
write.csv(nrDEG, file = "./01_DEG/01.DEG_GSE15654_Poor_prognosis-Good_prognosis_all_2.csv")

nrDEG$Label = ""   
nrDEG <- nrDEG[order(nrDEG$P.Value), ]   
nrDEG$Gene <- rownames(nrDEG)
#up.genes <- head(nrDEG$Gene[which(nrDEG$change == "Up")], 5)
#down.genes <- head(nrDEG$Gene[which(nrDEG$change == "Down")], 5)
#nrDEG.top5.genes <- c(as.character(up.genes), as.character(down.genes))
#nrDEG$Label[match(nrDEG.top5.genes, nrDEG$Gene)] <- nrDEG.top5.genes
p1 <- ggplot(
  nrDEG, 
  aes(x = logFC, 
      y = -log10(P.Value), 
      colour=change)) +
  geom_point(alpha=0.4, size=2) +
  scale_color_manual(values=c("#4DBBD5", "#d2dae2","#E64B35"))+
  geom_vline(xintercept=0,lty=4,col="black",lwd=0.8) +
  geom_vline(xintercept=0.5,lty=2,col="grey50",lwd=0.4) +
  geom_vline(xintercept=-0.5,lty=2,col="grey50",lwd=0.4) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(title = "GSE15654",
       # subtitle = "Cirrhotic Patients with vs. without HE",
       x = expression(log[2]~FC),  # x轴标签，log2中的2为下标
       y = expression(-log[10]~italic(p))  )+
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  
    plot.subtitle = element_text(hjust = 0.5, size = 14), 
    axis.title.x = element_text(size = 12), 
    axis.title.y = element_text(size = 12),  
    axis.text = element_text(size = 12),  
    legend.position="right", 
    legend.text = element_text(size = 12),
    legend.title = element_blank()
  )+geom_text_repel(data = nrDEG, aes(label = Label),
                    size = 3.5,                           
                    box.padding = unit(0.5, "lines"),  
                    point.padding = unit(0.8, "lines"), 
                    segment.color = "black",            
                    show.legend = FALSE,  
                    max.overlaps = 10000)
ggsave("./01_DEG/01.DEG_GSE15654_DEG_valcano.png",plot = p1,width = 5, height = 5, units = "in")
ggsave("./01_DEG/01.DEG_GSE15654_DEG_valcano.pdf",plot = p1,width = 5, height = 5, units = "in")

#venn_plot


