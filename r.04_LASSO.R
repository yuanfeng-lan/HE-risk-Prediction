rm(list = ls()); gc()
ORIGINAL_DIR <- "~/data_hdd/R/Yclub/2026.3.4_HE_Cirrhosis"
output <- file.path(ORIGINAL_DIR, "04_LASSO")
labels <- "04_LASSO_"
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(ORIGINAL_DIR)
library(glmnet)
library(ggplot2)
library(reshape2)
library(ggsci)
library(survival); library(survminer)
library(ggtext)
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
x <- as.matrix(t(exp))
y <- group$group
set.seed(123)
cvfit=cv.glmnet(x, y, family="gaussian",nlambda=100,alpha=1)
fit <- glmnet(x, y, family = "gaussian",nlambda=100, alpha=1)
min <- cvfit$lambda.min
min <- round(min, digits = 4)
min
x <- coef(fit) 
tmp <- as.data.frame(as.matrix(x)) 
tmp <- tmp[-1,,drop=F]
tmp$coef <- row.names(tmp) 
tmp <- reshape::melt(tmp, id = "coef") 
tmp$variable <- as.numeric(gsub("s", "", tmp$variable)) 
tmp$coef <- gsub('_','-',tmp$coef) 
tmp$lambda <- fit$lambda[tmp$variable+1] 
tmp$norm <- apply(abs(x[-1,]), 2, sum)[tmp$variable+1]
x <- as.matrix(t(exp))
fit2 <- glmnet(x=x, y=y, alpha = 1,family = "gaussian", lambda=cvfit$lambda.min)
fit2$beta
choose_gene=rownames(fit2$beta)[as.numeric(fit2$beta)!=0]

# 修改图例标签的颜色和字体
p <- ggplot(tmp, aes(log(lambda), value, color = coef)) + 
  geom_vline(xintercept = log(cvfit$lambda.min),
             size = 0.8, color = 'grey60',
             alpha = 0.8, linetype = 2) +
  geom_line(size = 1) + 
  xlab(expression(paste( lambda," (log scale)"))) + 
  ylab('Coefficients') + 
  theme_bw(base_rect_size = 2) + 
  scale_color_manual(
    name = "Coefficient",
    values = c(pal_npg()(9), # 限制调色板为 9 个颜色
               pal_d3()(7)), # 限制调色板为 7 个颜色
    breaks = unique(tmp$coef),
    labels = function(x) {
      # 将 choose_gene 高亮为红色，其他基因名改为斜体
      sapply(x, function(gene) {
        if (gene %in% choose_gene) {
          paste0("<span style='color:red;'><i>", gene, "</i></span>")
        } else {
          paste0("<i>", gene, "</i>")
        }
      })
    }
  ) + 
  scale_x_continuous(expand = c(0.01, 0.01)) + 
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 15, color = 'black'),
    axis.text = element_text(size = 12, color = 'black'),
    legend.title = element_blank(),
    legend.text = element_markdown(size = 14), # 使用 element_markdown 支持 HTML 格式
    legend.position = "right",
    legend.direction = "vertical",
    legend.key.height = grid::unit(6, "mm"),
    plot.margin = ggplot2::margin(t = 2, b = 2, l = 2, r = 2, unit = "pt")
  ) + 
  annotate('text', x = -3.0, y = -0.8,
           label = expression(paste("Optimal ", lambda, " = 0.0386")),
           color = 'black', size = 6)+ 
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 2)))
p
ggsave(file.path(output, paste0(labels, GSE, "result01.pdf")), p, w = 6, h = 6)
ggsave(file.path(output, paste0(labels, GSE, "result01.png")), p, w = 6, h = 6)
xx <- data.frame(lambda=cvfit[["lambda"]],
                 cvm=cvfit[["cvm"]],
                 cvsd=cvfit[["cvsd"]], 
                 cvup=cvfit[["cvup"]],
                 cvlo=cvfit[["cvlo"]],
                 nozezo=cvfit[["nzero"]]) 
xx$ll<- log(xx$lambda) 
xx$NZERO<- paste0(xx$nozezo,' vars')
p2 <- ggplot(xx,aes(ll,cvm,color=NZERO))+ 
  geom_errorbar(aes(x=ll,ymin=cvlo,ymax=cvup),
                width=0.05,size=1)+ 
  geom_vline(xintercept = xx$ll[which.min(xx$cvm)],
             size=0.8,color='grey60',alpha=0.8,
             linetype=2)+ 
  geom_point(size=2)+ 
  xlab(expression(paste( lambda," (log scale)")))+
  ylab('Partial Likelihood Deviance')+ 
  theme_bw(base_rect_size = 1.5)+ 
  scale_color_manual(values= c(pal_npg()(10),
                               pal_d3()(10),
                               pal_lancet()(10),
                               pal_aaas()(10),
                               pal_jco()(10)))+ 
  scale_x_continuous(expand = c(0.02,0.02))+ 
  scale_y_continuous(expand = c(0.02,0.02))+ 
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size=15,
                                  color='black'), 
        axis.text = element_text(size=12,
                                 color='black'), 
        legend.title = element_blank(), 
        legend.text = element_text(size=12,
                                   color='black'), 
        legend.position = 'none')+ 
  annotate('text',x= -4.0,y=1.8,
           label=expression(paste("Optimal ", lambda," = 0.0386")),
           color='black',size = 6)+ 
  guides(col=guide_legend(ncol = 6))
ggsave(file.path(output, paste0(labels,GSE,"result02.pdf")),p2,w=6,h=6)
ggsave(file.path(output, paste0(labels,GSE,"result02.png")),p2,w=6,h=6)

lambda <- cvfit$lambda.min
x <- as.matrix(t(exp))
fit2 <- glmnet(x=x, y=y, alpha = 1,family = "gaussian", lambda=cvfit$lambda.min)
fit2$beta
choose_gene=rownames(fit2$beta)[as.numeric(fit2$beta)!=0]
write.csv(choose_gene, file.path(output,paste0(labels,GSE,"result.csv") ))


#GSE15654
GSE <- "GSE15654_"
topgene <- read.csv(file.path(ORIGINAL_DIR,"02_Mfuzz","HE_GSE15654_ALL_gene_DEG.csv"))
topgene <- topgene$x
exp <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE15654_exp.csv"),row.names = 1)
group <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE15654_group.csv"),row.names = 1)
surv <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE15654_sur.csv"),row.names = 1)
surv <- surv[colnames(exp),,drop=F]
exp <- exp[topgene,]
x <- as.matrix(t(exp))
y <- Surv(surv$days_to_death, surv$death)
set.seed(123)
cvfit=cv.glmnet(x, y, family="cox",nlambda=100,alpha=1)
fit <- glmnet(x, y, family = "cox",nlambda=100, alpha=1)
min <- cvfit$lambda.min
min <- round(min, digits = 4)
min
x <- coef(fit) 
tmp <- as.data.frame(as.matrix(x)) 
tmp$coef <- row.names(tmp) 
tmp <- reshape::melt(tmp, id = "coef") 
tmp$variable <- as.numeric(gsub("s", "", tmp$variable)) 
tmp$coef <- gsub('_','-',tmp$coef) 
tmp$lambda <- fit$lambda[tmp$variable+1] 
tmp$norm <- apply(abs(x[-1,]), 2, sum)[tmp$variable+1]
x <- as.matrix(t(exp))
fit2 <- glmnet(x=x, y=y, alpha = 1,family = "cox", lambda=cvfit$lambda.min)
fit2$beta
choose_gene=rownames(fit2$beta)[as.numeric(fit2$beta)!=0]
p <- ggplot(tmp, aes(log(lambda), value, color = coef)) + 
  geom_vline(xintercept = log(cvfit$lambda.min),
             size = 0.8, color = 'grey60',
             alpha = 0.8, linetype = 2) +
  geom_line(size = 1) + 
  xlab(expression(paste( lambda," (log scale)"))) + 
  ylab('Coefficients') + 
  theme_bw(base_rect_size = 2) + 
  scale_color_manual(
    name = "Coefficient",
    values = c(pal_npg()(9), # 限制调色板为 9 个颜色
               pal_d3()(7)), # 限制调色板为 7 个颜色
    breaks = unique(tmp$coef),
    labels = function(x) {
      # 将 choose_gene 高亮为红色，其他基因名改为斜体
      sapply(x, function(gene) {
        if (gene %in% choose_gene) {
          paste0("<span style='color:red;'><i>", gene, "</i></span>")
        } else {
          paste0("<i>", gene, "</i>")
        }
      })
    }
  ) + 
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 15, color = 'black'),
    axis.text = element_text(size = 12, color = 'black'),
    legend.title = element_blank(),
    legend.text = element_markdown(size = 14), # 使用 element_markdown 支持 HTML 格式
    legend.position = "right",
    legend.direction = "vertical",
    legend.key.height = grid::unit(6, "mm"),
    plot.margin = ggplot2::margin(t = 2, b = 2, l = 2, r = 2, unit = "pt")
  ) + 
  annotate('text',x = -4.0, y = 0.0008,
           label = expression(paste("Optimal ", lambda, " = 0.0159")),
           color = 'black', size = 6)+ 
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 2)))
p
ggsave(file.path(output, paste0(labels,GSE,"result01.pdf")),p,w=6,h=6)
ggsave(file.path(output, paste0(labels,GSE,"result01.png")),p,w=6,h=6)
xx <- data.frame(lambda=cvfit[["lambda"]],
                 cvm=cvfit[["cvm"]],
                 cvsd=cvfit[["cvsd"]], 
                 cvup=cvfit[["cvup"]],
                 cvlo=cvfit[["cvlo"]],
                 nozezo=cvfit[["nzero"]]) 
xx$ll<- log(xx$lambda) 
xx$NZERO<- paste0(xx$nozezo,' vars')
p2 <- ggplot(xx,aes(ll,cvm,color=NZERO))+ 
  geom_errorbar(aes(x=ll,ymin=cvlo,ymax=cvup),
                width=0.05,size=1)+ 
  geom_vline(xintercept = xx$ll[which.min(xx$cvm)],
             size=0.8,color='grey60',alpha=0.8,
             linetype=2)+ 
  geom_point(size=2)+ 
  xlab(expression(paste( lambda," (log scale)")))+
  ylab('Partial Likelihood Deviance')+ 
  theme_bw(base_rect_size = 1.5)+ 
  scale_color_manual(values= c(pal_npg()(10),
                               pal_d3()(10),
                               pal_lancet()(10),
                               pal_aaas()(10),
                               pal_jco()(10)))+ 
  #scale_x_continuous(expand = c(0.02,0.02))+ 
  #scale_y_continuous(expand = c(0.02,0.02))+ 
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size=15,
                                  color='black'), 
        axis.text = element_text(size=12,
                                 color='black'), 
        legend.title = element_blank(), 
        legend.text = element_text(size=12,
                                   color='black'), 
        legend.position = 'none')+ 
  annotate('text',x= -4.0,y=3.8,
           label= expression(paste("Optimal ", lambda, " = 0.0159")),
           color='black',size = 6)+ 
  guides(col=guide_legend(ncol = 6))
ggsave(file.path(output, paste0(labels,GSE,"result02.pdf")),p2,w=6,h=6)
ggsave(file.path(output, paste0(labels,GSE,"result02.png")),p2,w=6,h=6)

lambda <- cvfit$lambda.min
x <- as.matrix(t(exp))
fit2 <- glmnet(x=x, y=y, alpha = 1,family = "cox", lambda=cvfit$lambda.min)
fit2$beta
choose_gene=rownames(fit2$beta)[as.numeric(fit2$beta)!=0]
write.csv(choose_gene, file.path(output,paste0(labels,GSE,"result.csv") ))

