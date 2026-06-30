rm(list = ls()); gc()
ORIGINAL_DIR <- "~/data_hdd/R/Yclub/2026.3.4_HE_Cirrhosis"
output <- file.path(ORIGINAL_DIR, "00_rawdata")
labels <- "00_rawdata_"

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(output)

library(GEOquery)
library(tidyverse)
####GSE41919####

GSE41919_gset = getGEO('GSE41919', destdir=".", AnnotGPL = F, getGPL = F)
GSE41919_gset[[1]]
GSE41919_pdata <- pData(GSE41919_gset[[1]])
colnames(GSE41919_pdata)
GSE41919_group <- GSE41919_pdata[,11,drop=F]
table(GSE41919_group$characteristics_ch1.1)
GSE41919_group$group <- ifelse(GSE41919_group$characteristics_ch1.1=="disease state: cirrhosis with HE",
                               "cirrhosis with HE",ifelse(GSE41919_group$characteristics_ch1.1=="disease state: cirrhosis without HE","cirrhosis without HE","non-cirrhotic control"))
table(GSE41919_group$group)
GSE41919_exp <- exprs(GSE41919_gset[[1]])
index01 = GSE41919_gset[[1]]@annotation
gpl14550 <- getGEO(index01, destdir = ".")
annotation_table <- Table(gpl14550)
ann <- annotation_table[,c(1,7)]
ann_output <- file.path(output,"annotation_data")
if (!dir.exists(ann_output)) {
  dir.create(ann_output, recursive = TRUE)
}
write.csv(ann,file.path(ann_output,"gpl_14550_ann.csv"))
GSE41919_exp <- as.data.frame(GSE41919_exp)
GSE41919_exp <- GSE41919_exp %>% mutate(ID=rownames(GSE41919_exp))
GSE41919_exp <- GSE41919_exp %>% inner_join(ann,by="ID") 
GSE41919_exp <- GSE41919_exp[!duplicated(GSE41919_exp$GENE_SYMBOL),]
rownames(GSE41919_exp) <- GSE41919_exp$GENE_SYMBOL
GSE41919_exp <- GSE41919_exp[-1,-(20:21)]

write.csv(GSE41919_exp,"00.rawdata_GSE41919_exp.csv")
write.csv(GSE41919_group,"00.rawdata_GSE41919_group.csv")
write.csv(GSE41919_pdata,"00.rawdata_GSE41919_pdata.csv")

####GSE139602####
GSE139602_gset = getGEO('GSE139602', destdir=".", AnnotGPL = F, getGPL = F)
GSE139602_gset[[1]]
GSE139602_pdata <- pData(GSE139602_gset[[1]])
GSE139602_group <- GSE139602_pdata[,10,drop=F]
table(GSE139602_group$characteristics_ch1)
GSE139602_exp <- exprs(GSE139602_gset[[1]])
index01 = GSE139602_gset[[1]]@annotation

GSE139602_exp <- as.data.frame(GSE139602_exp)
gpl13667 <- getGEO(index01, destdir = ".")
annotation_table <- Table(gpl13667)
ann <- annotation_table[,c(1,15)]
GSE139602_exp <- as.data.frame(GSE139602_exp)
GSE139602_exp <- GSE139602_exp %>% mutate(ID=rownames(GSE139602_exp))
GSE139602_exp <- GSE139602_exp %>% inner_join(ann,by="ID") 

GSE139602_exp <- GSE139602_exp %>%
  separate(`Gene Symbol`, into = c("First_Gene"), sep = "///", extra = "drop")
GSE139602_exp <- GSE139602_exp[!duplicated(GSE139602_exp$First_Gene),]

rownames(GSE139602_exp) <- GSE139602_exp$First_Gene
GSE139602_exp <- GSE139602_exp[-1,-(40:41)]

write.csv(GSE139602_exp,file.path(output,"00.rawdata_GSE139602_exp.csv"))
write.csv(GSE139602_group,file.path(output,"00.rawdata_GSE139602_group.csv"))
write.csv(GSE139602_pdata,file.path(output,"00.rawdata_GSE139602_pdata.csv"))

####GSE57193####
GSE57193_gset = getGEO('GSE57193', destdir=".", AnnotGPL = F, getGPL = F)
GSE57193_gset[[1]]
GSE57193_pdata <- pData(GSE57193_gset[[1]])
GSE57193_group <- GSE57193_pdata[,33,drop = F]
GSE57193_group$group <- ifelse(GSE57193_group$`diagnosis:ch1`=="liver cirrhosis","cirrhosis",
                               ifelse(GSE57193_group$`diagnosis:ch1`=="liver cirrhosis and hepatic encephalopathy","cirrhosis with HE","healthy"))
GSE57193_exp <- exprs(GSE57193_gset[[1]])
index01 = GSE57193_gset[[1]]@annotation
gpl14550 <- getGEO(index01, destdir = ".")
annotation_table <- Table(gpl14550)
ann <- annotation_table[,c(1,7)]
GSE57193_exp <- as.data.frame(GSE57193_exp)
GSE57193_exp <- GSE57193_exp %>% mutate(ID=rownames(GSE57193_exp))
GSE57193_exp <- GSE57193_exp %>% inner_join(ann,by="ID") 
GSE57193_exp <- GSE57193_exp[!duplicated(GSE57193_exp$GENE_SYMBOL),]
rownames(GSE57193_exp) <- GSE57193_exp$GENE_SYMBOL
GSE57193_exp <- GSE57193_exp[-1,-(13:14)]
write.csv(GSE57193_exp,"00.rawdata_GSE57193_exp.csv")
write.csv(GSE57193_group,"00.rawdata_GSE57193_group.csv")
write.csv(GSE57193_pdata,"00.rawdata_GSE57193_pdata.csv")

####GSE15654####
GSE15654_gset = getGEO('GSE15654', destdir=".", AnnotGPL = F, getGPL = F)
GSE15654_gset[[1]]
GSE15654_pdata <- pData(GSE15654_gset[[1]])
GSE15654_sur <- GSE15654_pdata[,c(17:24)]
colnames(GSE15654_sur) <- c("days_to_death","death","days_to_decomp","decomp",
                            "days_to_child","child","days_to_hcc","hcc")
GSE15654_sur$days_to_death <- sub(".*: \\s*", "", GSE15654_sur$days_to_death)

GSE15654_sur <- lapply(GSE15654_sur, function(x) {
  sub(".*: \\s*", "", x)
})

GSE15654_sur <- as.data.frame(GSE15654_sur)
rownames(GSE15654_sur) <- rownames(GSE15654_pdata)
GSE15654_exp <- exprs(GSE15654_gset[[1]])
index01 = GSE15654_gset[[1]]@annotation
gpl8432 <- getGEO(index01, destdir = ".")
annotation_table <- Table(gpl8432)
ann <- annotation_table[,c(1,12)]
GSE15654_exp <- as.data.frame(GSE15654_exp)
GSE15654_exp <- GSE15654_exp %>% mutate(ID=rownames(GSE15654_exp))
GSE15654_exp <- GSE15654_exp %>% inner_join(ann,by="ID") 
GSE15654_exp <- GSE15654_exp[!duplicated(GSE15654_exp$Symbol),]
rownames(GSE15654_exp) <- GSE15654_exp$Symbol
GSE15654_exp <- GSE15654_exp[-1,-(217:218)]
GSE15654_group <- GSE15654_pdata[,11,drop = F]
GSE15654_group$group <- ifelse(GSE15654_group$characteristics_ch1.1=="prediction: Poor prognosis",
                               "Poor_prognosis","Good_prognosis")


write.csv(GSE15654_exp,"00.rawdata_GSE15654_exp.csv")
write.csv(GSE15654_group,"00.rawdata_GSE15654_group.csv")
write.csv(GSE15654_pdata,"00.rawdata_GSE15654_pdata.csv")
write.csv(GSE15654_sur,"00.rawdata_GSE15654_sur.csv")

