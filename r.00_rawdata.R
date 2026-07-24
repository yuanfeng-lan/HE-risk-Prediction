rm(list = ls()); gc()
ORIGINAL_DIR <- ""
output <- file.path(ORIGINAL_DIR, "00_rawdata")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(output)

library(GEOquery)
library(tidyverse)
library(dplyr)
library(tidyr)

# Helper function to process expression data
process_expression <- function(gset, ann, id_col, symbol_col, keep_cols = NULL) {
  exp <- as.data.frame(exprs(gset))
  exp <- exp %>% mutate(ID = rownames(exp))
  exp <- exp %>% inner_join(ann, by = "ID")
  exp <- exp[!duplicated(exp[[symbol_col]]), ]
  rownames(exp) <- exp[[symbol_col]]
  if (!is.null(keep_cols)) {
    exp <- exp[, keep_cols, drop = FALSE]
  }
  return(exp)
}

# Helper function to get annotation
get_annotation <- function(gpl_id, cols) {
  gpl <- getGEO(gpl_id, destdir = ".")
  ann <- Table(gpl)[, cols]
  return(ann)
}

#### GSE41919 ####
GSE41919_gset <- getGEO('GSE41919', destdir=".", AnnotGPL = F, getGPL = F)
GSE41919_pdata <- pData(GSE41919_gset[[1]])
GSE41919_group <- GSE41919_pdata[, 11, drop = F]
GSE41919_group$group <- ifelse(GSE41919_group$characteristics_ch1.1 == "disease state: cirrhosis with HE",
                               "cirrhosis with HE", 
                               ifelse(GSE41919_group$characteristics_ch1.1 == "disease state: cirrhosis without HE",
                                      "cirrhosis without HE", "non-cirrhotic control"))

ann <- get_annotation(GSE41919_gset[[1]]@annotation, c(1, 7))
GSE41919_exp <- process_expression(GSE41919_gset[[1]], ann, "ID", "GENE_SYMBOL", 
                                   keep_cols = -c(ncol(ann)+1, ncol(ann)+2))

write.csv(GSE41919_exp, "00.rawdata_GSE41919_exp.csv")
write.csv(GSE41919_group, "00.rawdata_GSE41919_group.csv")
write.csv(GSE41919_pdata, "00.rawdata_GSE41919_pdata.csv")

#### GSE139602 ####
GSE139602_gset <- getGEO('GSE139602', destdir=".", AnnotGPL = F, getGPL = F)
GSE139602_pdata <- pData(GSE139602_gset[[1]])
GSE139602_group <- GSE139602_pdata[, 10, drop = F]

ann <- get_annotation(GSE139602_gset[[1]]@annotation, c(1, 15))
GSE139602_exp <- as.data.frame(exprs(GSE139602_gset[[1]]))
GSE139602_exp <- GSE139602_exp %>% mutate(ID = rownames(GSE139602_exp))
GSE139602_exp <- GSE139602_exp %>% inner_join(ann, by = "ID")
GSE139602_exp <- GSE139602_exp %>%
  separate(`Gene Symbol`, into = "First_Gene", sep = "///", extra = "drop")
GSE139602_exp <- GSE139602_exp[!duplicated(GSE139602_exp$First_Gene), ]
rownames(GSE139602_exp) <- GSE139602_exp$First_Gene
GSE139602_exp <- GSE139602_exp[, -c(1, ncol(ann)+1, ncol(ann)+2)]

write.csv(GSE139602_exp, "00.rawdata_GSE139602_exp.csv")
write.csv(GSE139602_group, "00.rawdata_GSE139602_group.csv")
write.csv(GSE139602_pdata, "00.rawdata_GSE139602_pdata.csv")

#### GSE57193 ####
GSE57193_gset <- getGEO('GSE57193', destdir=".", AnnotGPL = F, getGPL = F)
GSE57193_pdata <- pData(GSE57193_gset[[1]])
GSE57193_group <- GSE57193_pdata[, 33, drop = F]
GSE57193_group$group <- ifelse(GSE57193_group$`diagnosis:ch1` == "liver cirrhosis", "cirrhosis",
                               ifelse(GSE57193_group$`diagnosis:ch1` == "liver cirrhosis and hepatic encephalopathy", 
                                      "cirrhosis with HE", "healthy"))

ann <- get_annotation(GSE57193_gset[[1]]@annotation, c(1, 7))
GSE57193_exp <- process_expression(GSE57193_gset[[1]], ann, "ID", "GENE_SYMBOL",
                                   keep_cols = -c(ncol(ann)+1, ncol(ann)+2))

write.csv(GSE57193_exp, "00.rawdata_GSE57193_exp.csv")
write.csv(GSE57193_group, "00.rawdata_GSE57193_group.csv")
write.csv(GSE57193_pdata, "00.rawdata_GSE57193_pdata.csv")

#### GSE15654 ####
GSE15654_gset <- getGEO('GSE15654', destdir=".", AnnotGPL = F, getGPL = F)
GSE15654_pdata <- pData(GSE15654_gset[[1]])

GSE15654_sur <- GSE15654_pdata[, 17:24]
colnames(GSE15654_sur) <- c("days_to_death", "death", "days_to_decomp", "decomp",
                            "days_to_child", "child", "days_to_hcc", "hcc")
GSE15654_sur <- as.data.frame(lapply(GSE15654_sur, function(x) sub(".*: \\s*", "", x)))
rownames(GSE15654_sur) <- rownames(GSE15654_pdata)

ann <- get_annotation(GSE15654_gset[[1]]@annotation, c(1, 12))
GSE15654_exp <- process_expression(GSE15654_gset[[1]], ann, "ID", "Symbol",
                                   keep_cols = -c(ncol(ann)+1, ncol(ann)+2))

GSE15654_group <- GSE15654_pdata[, 11, drop = F]
GSE15654_group$group <- ifelse(GSE15654_group$characteristics_ch1.1 == "prediction: Poor prognosis",
                               "Poor_prognosis", "Good_prognosis")

write.csv(GSE15654_exp, "00.rawdata_GSE15654_exp.csv")
write.csv(GSE15654_group, "00.rawdata_GSE15654_group.csv")
write.csv(GSE15654_pdata, "00.rawdata_GSE15654_pdata.csv")
write.csv(GSE15654_sur, "00.rawdata_GSE15654_sur.csv")
