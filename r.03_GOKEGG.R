rm(list = ls()); gc()
ORIGINAL_DIR <- "~/data_hdd/R/Yclub/2026.3.4_HE_Cirrhosis"
output <- file.path(ORIGINAL_DIR, "03_GOKEGG")
labels <- "03_GOKEGG_"
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(ORIGINAL_DIR)
library(data.table); library(ggplot2); library(stringr); library(scales); library(cowplot)
library(clusterProfiler); library(enrichplot); library(AnnotationDbi); library(org.Hs.eg.db); 
library(tidyverse)
script.dir <- file.path(output)
files <- list.files(script.dir, pattern = "DEG\\.csv$", full.names = TRUE, ignore.case = TRUE)
output_dir <- file.path(output,"enrich_results")
if(!dir.exists(output_dir)) dir.create(output_dir)
calculate_gene_ratio <- function(ratio_str) {
  if (is.character(ratio_str) && grepl("/", ratio_str)) {
    parts <- strsplit(ratio_str, "/")[[1]]
    num <- as.numeric(parts[1])
    den <- as.numeric(parts[2])
    return(ifelse(den != 0, num / den, 0))  
  } else {
    return(as.numeric(ratio_str))  
  }
}
for (fl in files) {
  sample <- tools::file_path_sans_ext(basename(fl))
  message("=：", sample)
  info <- read.csv(fl,sep = "\t")
  if(is.null(info$x)){
    info <- read.csv(fl)
  }
  
  OUT_DIR <- file.path(output_dir,sample)
  dir.create(OUT_DIR, showWarnings = FALSE)
  genes <- info$x
  entrez <- bitr(genes,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db) 
  res.list <- list()
  go <- enrichGO(gene = entrez$ENTREZID, 
                 OrgDb = org.Hs.eg.db,
                 keyType = "ENTREZID", 
                 ont = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 1,
                 qvalueCutoff = 1,
                 readable = TRUE)
  kegg <- enrichKEGG(gene = entrez$ENTREZID,
                     keyType = "kegg",
                     organism = "hsa",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 1)
  go_results <- as.data.frame(go)
  kegg_results <- as.data.frame(kegg)
  kegg_results$ONTOLOGY <- "KEGG"
  col <- c("ONTOLOGY","ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID",
           "Count")
  go_results <- go_results[,col]
  kegg_results <- kegg_results[,col]
  go_results$ONTOLOGY <- paste0("GO",go_results$ONTOLOGY)
  all_results <- rbind(go_results,kegg_results)
  write.csv(all_results,file.path(OUT_DIR,paste0(sample,"_all_results.csv")))
  all_results <- subset(all_results,pvalue<=0.05)
  all_results <- all_results %>%
    mutate(GeneRatio_numeric = sapply(GeneRatio, calculate_gene_ratio))
  filtered_df <- all_results %>%
    group_by(ONTOLOGY) %>%
    arrange(desc(GeneRatio_numeric)) %>%
    slice_head(n = 5) %>%
    ungroup()
  filtered_df <- filtered_df %>%
    mutate(log10_p_value = -log10(pvalue))
  filtered_df <- filtered_df %>%
    mutate(Description_wrapped = str_wrap(Description, width = 40))
  filtered_df <- filtered_df[order(match(filtered_df$ONTOLOGY, unique(filtered_df$ONTOLOGY)), filtered_df$GeneRatio_numeric), ]
  filtered_df$Description_wrapped <- factor(filtered_df$Description_wrapped,levels = filtered_df$Description_wrapped)
  filtered_df$ONTOLOGY <- gsub("GO","GO ",filtered_df$ONTOLOGY)
  p <- ggplot(filtered_df, aes(x = GeneRatio_numeric, y = Description_wrapped, size = Count, color = log10_p_value)) +
    geom_point(alpha = 0.7) +
    scale_color_gradientn(
      name = expression(-log[10](italic(p))), 
      colors = c("#3C5488", "#F1C40F", "#DC0000"), 
      limits = c(min(filtered_df$log10_p_value), max(filtered_df$log10_p_value)), 
      na.value = "grey50" 
    ) +  
    labs(
      x = "GeneRatio",
      y = "Terms",
      title = ""
    ) +
    facet_wrap(~ ONTOLOGY, ncol = 1, scales = "free_y") +  
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
      axis.text.y = element_text(size = 16),  
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      plot.title = element_text(hjust = 0.5),
      strip.text = element_text(size = 15, face = "bold")  
    )
  
  ggsave(file.path(OUT_DIR,paste0(sample,"_enrichment_bubble_plot.png")), p, width = 8, height = 11)
  ggsave(file.path(OUT_DIR,paste0(sample,"_enrichment_bubble_plot.pdf")), p, width = 8, height = 11)
  
}



#GSE15654
GSE15654_DEG <- read.csv(file.path("01_DEG","01.DEG_GSE15654_Poor_prognosis-Good_prognosis_all.csv"))
GSE15654_DEG_UP <- subset(GSE15654_DEG,P.Value<=0.05&logFC>=0.5)
GSE15654_DEG_DOWN <- subset(GSE15654_DEG,P.Value<=0.05&logFC<=-0.5)
GSE15654_DEG_all_gene <- rbind(GSE15654_DEG_UP,GSE15654_DEG_DOWN)
#GSE15654_DEG_UP
sample <- "GSE15654_DEG_UP"
OUT_DIR <- file.path(output_dir,sample)
dir.create(OUT_DIR, showWarnings = FALSE)
genes <- GSE15654_DEG_UP$X
entrez <- bitr(genes,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db) 
res.list <- list()
go <- enrichGO(gene = entrez$ENTREZID, 
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID", 
               ont = "ALL",
               pAdjustMethod = "BH",
               pvalueCutoff = 1,
               qvalueCutoff = 1,
               readable = TRUE)
kegg <- enrichKEGG(gene = entrez$ENTREZID,
                   keyType = "kegg",
                   organism = "hsa",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 1)
go_results <- as.data.frame(go)
kegg_results <- as.data.frame(kegg)
kegg_results$ONTOLOGY <- "KEGG"
col <- c("ONTOLOGY","ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID",
         "Count")
go_results <- go_results[,col]
kegg_results <- kegg_results[,col]
go_results$ONTOLOGY <- paste0("GO",go_results$ONTOLOGY)
all_results <- rbind(go_results,kegg_results)
write.csv(all_results,file.path(OUT_DIR,paste0(sample,"_all_results.csv")))
all_results <- subset(all_results,pvalue<=0.05)
all_results <- all_results %>%
  mutate(GeneRatio_numeric = sapply(GeneRatio, calculate_gene_ratio))
filtered_df <- all_results %>%
  group_by(ONTOLOGY) %>%
  arrange(desc(GeneRatio_numeric)) %>%
  slice_head(n = 5) %>%
  ungroup()
filtered_df <- filtered_df %>%
  mutate(log10_p_value = -log10(pvalue))
filtered_df <- filtered_df %>%
  mutate(Description_wrapped = str_wrap(Description, width = 50))
filtered_df <- filtered_df[order(match(filtered_df$ONTOLOGY, unique(filtered_df$ONTOLOGY)), filtered_df$GeneRatio_numeric), ]
filtered_df$Description_wrapped <- factor(filtered_df$Description_wrapped,levels = filtered_df$Description_wrapped)
filtered_df$ONTOLOGY <- gsub("GO","GO ",filtered_df$ONTOLOGY)
p <- ggplot(filtered_df, aes(x = GeneRatio_numeric, y = Description_wrapped, size = Count, color = log10_p_value)) +
  geom_point(alpha = 0.7) +
  scale_color_gradientn(
    name = expression(-log[10](italic(p))), 
    colors = c("#3C5488", "#F1C40F", "#DC0000"), 
    limits = c(min(filtered_df$log10_p_value), max(filtered_df$log10_p_value)), 
    na.value = "grey50" 
  ) +  
  labs(
    x = "GeneRatio",
    y = "Terms",
    title = ""
  ) +
  facet_wrap(~ ONTOLOGY, ncol = 1, scales = "free_y") +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
    axis.text.y = element_text(size = 16),  
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 15, face = "bold")  
  )

ggsave(file.path(OUT_DIR,paste0(sample,"_enrichment_bubble_plot.png")), p, width = 9, height = 11)
ggsave(file.path(OUT_DIR,paste0(sample,"_enrichment_bubble_plot.pdf")), p, width = 9, height = 11)

#GSE15654_DEG_DOWN
sample <- "GSE15654_DEG_DOWN"
OUT_DIR <- file.path(output_dir,sample)
dir.create(OUT_DIR, showWarnings = FALSE)
genes <- GSE15654_DEG_DOWN$X
entrez <- bitr(genes,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db) 
res.list <- list()
go <- enrichGO(gene = entrez$ENTREZID, 
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID", 
               ont = "ALL",
               pAdjustMethod = "BH",
               pvalueCutoff = 1,
               qvalueCutoff = 1,
               readable = TRUE)
kegg <- enrichKEGG(gene = entrez$ENTREZID,
                   keyType = "kegg",
                   organism = "hsa",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 1)
go_results <- as.data.frame(go)
kegg_results <- as.data.frame(kegg)
kegg_results$ONTOLOGY <- "KEGG"
col <- c("ONTOLOGY","ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID",
         "Count")
go_results <- go_results[,col]
kegg_results <- kegg_results[,col]
go_results$ONTOLOGY <- paste0("GO",go_results$ONTOLOGY)
all_results <- rbind(go_results,kegg_results)
write.csv(all_results,file.path(OUT_DIR,paste0(sample,"_all_results.csv")))
all_results <- subset(all_results,pvalue<=0.05)
all_results <- all_results %>%
  mutate(GeneRatio_numeric = sapply(GeneRatio, calculate_gene_ratio))
filtered_df <- all_results %>%
  group_by(ONTOLOGY) %>%
  arrange(desc(GeneRatio_numeric)) %>%
  slice_head(n = 5) %>%
  ungroup()
filtered_df <- filtered_df %>%
  mutate(log10_p_value = -log10(pvalue))
filtered_df <- filtered_df %>%
  mutate(Description_wrapped = str_wrap(Description, width = 50))
filtered_df <- filtered_df[order(match(filtered_df$ONTOLOGY, unique(filtered_df$ONTOLOGY)), filtered_df$GeneRatio_numeric), ]
filtered_df$Description_wrapped <- factor(filtered_df$Description_wrapped,levels = filtered_df$Description_wrapped)
filtered_df$ONTOLOGY <- gsub("GO","GO ",filtered_df$ONTOLOGY)
p <- ggplot(filtered_df, aes(x = GeneRatio_numeric, y = Description_wrapped, size = Count, color = log10_p_value)) +
  geom_point(alpha = 0.7) +
  scale_color_gradientn(
    name = expression(-log[10](italic(p))), 
    colors = c("#3C5488", "#F1C40F", "#DC0000"), 
    limits = c(min(filtered_df$log10_p_value), max(filtered_df$log10_p_value)), 
    na.value = "grey50" 
  ) +  
  labs(
    x = "GeneRatio",
    y = "Terms",
    title = ""
  ) +
  facet_wrap(~ ONTOLOGY, ncol = 1, scales = "free_y") +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
    axis.text.y = element_text(size = 16),  
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 15, face = "bold")  
  )

ggsave(file.path(OUT_DIR,paste0(sample,"_enrichment_bubble_plot.png")), p, width = 9, height = 11)
ggsave(file.path(OUT_DIR,paste0(sample,"_enrichment_bubble_plot.pdf")), p, width = 9, height = 11)





