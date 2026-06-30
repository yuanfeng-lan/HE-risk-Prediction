rm(list = ls()); gc()
ORIGINAL_DIR <- "~/data_hdd/R/Yclub/2026.3.4_HE_Cirrhosis"
output <- file.path(ORIGINAL_DIR, "02_Mfuzz")
labels <- "02_Mfuzz_"
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(ORIGINAL_DIR)
library(tidyverse)
library(stringr)
library("Mfuzz")

exp_data <- read.csv("./00_rawdata/00.rawdata_GSE139602_exp.csv",row.names = 1)
group <- read.csv("./00_rawdata/00.rawdata_GSE139602_group.csv",row.names = 1)
group$group <- sub("^[^:]*:\\s?", "", group$characteristics_ch1)
exp_df <- exp_data[,rownames(group)]
exp_df <- exp_df %>% t() %>% as.data.frame()
exp_df$group <- group$group
exp_df$group <- factor(exp_df$group,levels = c("Healthy","eCLD","Compensated Cirrhosis"
                                               ,"Decompesated Cirrhosis","Acute-on-chronic liver failure"))
sample1<-aggregate(exp_df[,1:length(rownames(exp_data))],by=list(exp_df$group),mean,na.rm= TRUE)
row.names(sample1)<-sample1[,1]
sample1<-data.frame(t(sample1[,-1]))
colnames(sample1) <- c("Healthy","eCLD","CC","DC","ACLF")
sample1<-as.matrix(sample1)
sample1<- ExpressionSet(assayData = sample1)
sample1 <- filter.NA(sample1, thres = 0.25)
sample1 <- fill.NA(sample1, mode = 'mean')
sample1 <- filter.std(sample1, min.std = 0)
sample1 <- standardise(sample1)

set.seed(123)
cluster_num <- 12
sample1_cluster <- mfuzz(sample1, c = cluster_num, m = mestimate(sample1))

png(file.path(output,"02_Mfuzz_results_plot.png"),width = 12,height = 8,res = 300,units = "in")
mfuzz.plot2(
  sample1,               
  cl = sample1_cluster,   
  mfrow = c(3, 4),        
  time.labels = colnames(sample1),  
  centre = TRUE,          
  x11 = F,               
  col = c("#E64B35"),
  xlab = "",          
  ylab = "Expression",    
  cex.main = 1.2,         
  cex.lab = 1.2,            
  cex.axis = 1.0          
)
dev.off()

pdf(file.path(output,"02_Mfuzz_results_plot.pdf"),width = 12,height = 8)
mfuzz.plot2(
  sample1,               
  cl = sample1_cluster,   
  mfrow = c(3, 4),        
  time.labels = colnames(sample1),  
  centre = TRUE,          
  x11 = F,               
  col = c("#E64B35"),
  xlab = "",          
  ylab = "Expression",    
  cex.main = 1.2,         
  cex.lab = 1.2,            
  cex.axis = 1.0          
)
dev.off()

mfu_output <- file.path(output,"mfuzz")
if (!dir.exists(mfu_output)) {
  dir.create(mfu_output, recursive = TRUE)
}
for(i in 1:cluster_num){
  potname<-names(sample1_cluster$cluster[unname(sample1_cluster$cluster)==i])
  write.csv(sample1_cluster[[4]][potname,i],file.path(output,paste0("mfuzz/mfuzz_",i,".csv")))
}

table(sample1_cluster$cluster)

write_cluster <- function(i){
  df3 <- read.csv(file.path(mfu_output,paste0("mfuzz_",i,".csv")))
  df3 <- df3[order(df3$x,decreasing = T),]
  colnames(df3) <- c("x","Coefficient")
  write.csv(df3,file.path(output,paste0("mfuzz_",i,"_DEG.csv")))
}
write_cluster("2")
write_cluster("8")


#ggven
GSE139602_cluster2 <- read.csv(file.path("02_Mfuzz","mfuzz_2_DEG.csv"),row.names = 1)
GSE139602_cluster8 <- read.csv(file.path("02_Mfuzz","mfuzz_8_DEG.csv"),row.names = 1)
GSE139602_gene_all_02 <- rbind(GSE139602_cluster2,GSE139602_cluster8)
GSE139602_cluster2_UP <- subset(GSE139602_cluster2,Coefficient>=0.4)
GSE139602_cluster8_DOWN <- subset(GSE139602_cluster8,Coefficient>=0.4)

GSE139602_gene_all <- rbind(GSE139602_cluster2_UP,GSE139602_cluster8_DOWN)
write.csv(GSE139602_cluster2_UP,file.path(output,"GSE139602_cluster2_UP_DEG.csv"))
write.csv(GSE139602_cluster8_DOWN,file.path(output,"GSE139602_cluster8_DOWN_DEG.csv"))
write.csv(GSE139602_gene_all,file.path(output,"GSE139602_gene_all_DEG.csv"))

GSE41919_DEG <- read.csv(file.path("01_DEG","01.DEG_GSE41919_with_HE-without_HE_all.csv"))
GSE57193_DEG <- read.csv(file.path("01_DEG","01.DEG_GSE57193_cirrhosis_with_HE-cirrhosis_all.csv"))
GSE41919_DEG_UP <- subset(GSE41919_DEG,P.Value<=0.05&logFC>=0.5)
GSE41919_DEG_DOWN <- subset(GSE41919_DEG,P.Value<=0.05&logFC<=-0.5)
GSE57193_DEG_UP <- subset(GSE57193_DEG,P.Value<=0.05&logFC>=0.5)
GSE57193_DEG_DOWN <- subset(GSE57193_DEG,P.Value<=0.05&logFC<=-0.5)
GSE41919_DEG_all <- rbind(GSE41919_DEG_UP,GSE41919_DEG_DOWN)
GSE57193_DEG_all <- rbind(GSE57193_DEG_UP,GSE57193_DEG_DOWN)

venn_plot <- function(df01,df02,vec){
  parts_vec <- strsplit(vec, "_")[[1]]
  if ("X" %in% names(df01)) {
    gene01 <- df01$X
  } else if ("x" %in% names(df01)) {
    gene01 <- df01$x
  } else {
    stop("")
  }
  if ("X" %in% names(df02)) {
    gene02 <- df02$X
  } else if ("x" %in% names(df02)) {
    gene02 <- df02$x
  } else {
    stop("")
  }
  venn <- list(list1=gene01,
               list2=gene02)
  names(venn) <- c(paste0(parts_vec[1]," ",parts_vec[3]),paste0(parts_vec[2]," ",parts_vec[3]))
  p <- ggvenn(
    venn,  
    c(paste0(parts_vec[1]," ",parts_vec[3]),paste0(parts_vec[2]," ",parts_vec[3])),  
    text_size = 6.5,  
    fill_color = c("#00A087", "#3C5488"),  
    fill_alpha = 0.7,  
    stroke_color = "black", 
    stroke_size = 1  
  )
  
  ggsave(file.path(output, paste0(vec,"_venn_diagram.pdf")),p,w=6,h=5)
  ggsave(file.path(output, paste0(vec,"_venn_diagram.png")),p,w=6,h=5)
  genes <- intersect(gene01,gene02)
  write.csv(genes,file.path(output,paste0(vec,"_gene.csv")))
  write.csv(genes,file.path(output,paste0(vec,"_gene_DEG.csv")))
}
venn_plot(GSE41919_DEG_UP,GSE57193_DEG_UP,"GSE41919_GSE57193_UP")
venn_plot(GSE41919_DEG_DOWN,GSE57193_DEG_DOWN,"GSE41919_GSE57193_DOWN")
venn_plot(GSE41919_DEG_all,GSE57193_DEG_all,"GSE41919_GSE57193_ALL")

GSE41919_GSE57193_UP_gene <- read.csv(file.path(output,"GSE41919_GSE57193_UP_gene.csv"),row.names = 1)
GSE41919_GSE57193_DOWN_gene <- read.csv(file.path(output,"GSE41919_GSE57193_DOWN_gene.csv"),row.names = 1)
GSE41919_GSE57193_all_gene <- read.csv(file.path(output,"GSE41919_GSE57193_ALL_gene.csv"),row.names = 1)

venn_plot(GSE41919_GSE57193_UP_gene,GSE139602_cluster2_UP,"HE_GSE139602_UP")
venn_plot(GSE41919_GSE57193_DOWN_gene,GSE139602_cluster8_DOWN,"HE_GSE139602_DOWN")
venn_plot(GSE41919_GSE57193_all_gene,GSE139602_gene_all,"HE_GSE139602_ALL")

GSE41919_DEG <- read.csv(file.path("01_DEG","01.DEG_GSE41919_with_HE-without_HE_all.csv"))
GSE57193_DEG <- read.csv(file.path("01_DEG","01.DEG_GSE57193_cirrhosis_with_HE-cirrhosis_all.csv"))
GSE41919_DEG_UP <- subset(GSE41919_DEG,P.Value<=0.05&logFC>=0.5)
GSE41919_DEG_DOWN <- subset(GSE41919_DEG,P.Value<=0.05&logFC<=-0.5)
GSE57193_DEG_UP <- subset(GSE57193_DEG,P.Value<=0.05&logFC>=0.5)
GSE57193_DEG_DOWN <- subset(GSE57193_DEG,P.Value<=0.05&logFC<=-0.5)
GSE41919_DEG_all <- rbind(GSE41919_DEG_UP,GSE41919_DEG_DOWN)
GSE57193_DEG_all <- rbind(GSE57193_DEG_UP,GSE57193_DEG_DOWN)
venn_plot(GSE41919_DEG_UP,GSE57193_DEG_UP,"GSE41919_GSE57193_UP")
venn_plot(GSE41919_DEG_DOWN,GSE57193_DEG_DOWN,"GSE41919_GSE57193_DOWN")
venn_plot(GSE41919_DEG_all,GSE57193_DEG_all,"GSE41919_GSE57193_ALL")

GSE41919_GSE57193_UP_gene <- read.csv(file.path(output,"GSE41919_GSE57193_UP_gene.csv"),row.names = 1)
GSE41919_GSE57193_DOWN_gene <- read.csv(file.path(output,"GSE41919_GSE57193_DOWN_gene.csv"),row.names = 1)
GSE41919_GSE57193_all_gene <- read.csv(file.path(output,"GSE41919_GSE57193_ALL_gene.csv"),row.names = 1)


GSE15654_DEG <- read.csv(file.path("01_DEG","01.DEG_GSE15654_Poor_prognosis-Good_prognosis_all.csv"))
GSE15654_DEG_UP <- subset(GSE15654_DEG,P.Value<=0.05&logFC>=0.5)
GSE15654_DEG_DOWN <- subset(GSE15654_DEG,P.Value<=0.05&logFC<=-0.5)
GSE15654_DEG_all_gene <- rbind(GSE15654_DEG_UP,GSE15654_DEG_DOWN)

venn_plot(GSE41919_GSE57193_UP_gene,GSE15654_DEG_UP,"HE_GSE15654_UP")
venn_plot(GSE41919_GSE57193_DOWN_gene,GSE15654_DEG_DOWN,"HE_GSE15654_DOWN")
venn_plot(GSE41919_GSE57193_all_gene,GSE15654_DEG_all_gene,"HE_GSE15654_ALL")



