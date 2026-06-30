rm(list = ls()); gc()
ORIGINAL_DIR <- "~/data_hdd/R/Yclub/2026.3.4_HE_Cirrhosis"
output <- file.path(ORIGINAL_DIR, "10_GSEA")
labels <- "10_GSEA_"

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(ORIGINAL_DIR)
library(tidyverse)
library(data.table)
library(dplyr)
library(msigdbr)
library(fgsea)
library(ggplot2)
library(ggridges)
library(dplyr)

exp_df <- read.csv(file.path(ORIGINAL_DIR,"09_SEM","09_SEM__exp_corrected.csv"),row.names = 1)
group_merged <- read.csv(file.path(ORIGINAL_DIR,"09_SEM","09_SEM__group_merged.csv"),row.names = 1)
GSE41919_exp <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE41919_exp.csv"),row.names = 1)
GSE41919_group <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE41919_group.csv"),row.names = 1)
GSE41919_sample <- GSE41919_group %>% filter(GSE41919_group$group!="non-cirrhotic control")
GSE41919_sample <- rownames(GSE41919_sample)
GSE41919_exp <- GSE41919_exp[GSE41919_sample,]
GSE41919_group <- GSE41919_group[GSE41919_sample,]

GSE57193_exp <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE57193_exp.csv"),row.names = 1)
GSE57193_group <- read.csv(file.path(ORIGINAL_DIR,"00_rawdata","00.rawdata_GSE57193_group.csv"),row.names = 1)
GSE57193_sample <- GSE57193_group %>% filter(GSE57193_group$group!="healthy")
GSE57193_sample <- rownames(GSE57193_sample)
GSE57193_exp <- GSE57193_exp[GSE57193_sample,]
GSE57193_group <- GSE57193_group[GSE57193_sample,]
c5bp <- tryCatch(
  msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP"),
  error = function(e) NULL
)
save(c5bp,file = file.path(output,"c5bp.Rdata"))

####PRB2####
target_gene <- "PRB2"

expr <- as.matrix(exp_df)
mode(expr) <- "numeric"
target_vec <- expr["PRB2", , drop = TRUE]
cors <- apply(expr, 1, function(x) {
  suppressWarnings(cor(x, target_vec, method = "spearman", use = "pairwise.complete.obs"))
})

ranked <- sort(cors, decreasing = TRUE)

pathways_list <- list()
if(!is.null(c5bp)) {
  c5_sets <- split(c5bp$gene_symbol, c5bp$gs_name)
  names(c5_sets) <- paste0("C5_BP|", names(c5_sets))
  pathways_list <- c(pathways_list, c5_sets)
}
set.seed(123)
fgsea_res <- fgsea(pathways = pathways_list,
                   stats = ranked,
                   minSize = 15,
                   maxSize = 500,
                   nperm = 10000)

fgsea_res <- fgsea_res %>%
  as.data.frame() %>%
  arrange(padj, pval, desc(NES))
safe_fgsea <- fgsea_res
list_cols <- names(safe_fgsea)[sapply(safe_fgsea, is.list)]
for (col in list_cols) {
  safe_fgsea[[col]] <- sapply(safe_fgsea[[col]], function(x) {
    if (is.null(x)) return(NA_character_)
    if (length(x) == 0) return(NA_character_)
    paste(as.character(x), collapse = ";")
  }, USE.NAMES = FALSE)
}
safe_fgsea <- as.data.frame(safe_fgsea, stringsAsFactors = FALSE)
write.csv(safe_fgsea, file = file.path(output, paste0(labels, "PRB2_GSEA_results01.csv")), row.names = FALSE)


sig <- fgsea_res %>% filter(pval < 0.01)
if(nrow(sig) == 0) stop("No pathways with pval < 0.01")

top_pos <- sig %>% arrange(desc(NES)) %>% head(5)
top_neg <- sig %>% arrange(NES) %>% head(5)
selected <- bind_rows(top_pos, top_neg) %>% distinct(pathway, .keep_all = TRUE)
selected_paths <- selected$pathway
display_names <- sub("^C5_BP\\|GOBP_", "", selected_paths)
display_names <- sub("^C5_BP\\|", "", display_names)
display_names <- str_replace_all(str_to_lower(display_names), "_", " ")
names(display_names) <- selected_paths
eps <- 1e-300
neglog10pval <- -log10(selected$pval + eps)
names(neglog10pval) <- selected$pathway
calc_running_es <- function(stats, geneset) {
  N <- length(stats)
  selected_idx <- which(names(stats) %in% geneset)
  Nh <- length(selected_idx)
  if(Nh == 0) return(rep(0, N))
  weights <- abs(stats)^1
  hit <- integer(N); hit[selected_idx] <- 1
  Phit <- cumsum(hit * weights) / sum(weights[hit == 1])
  Pmiss <- cumsum((1 - hit) / (N - Nh))
  Phit - Pmiss
}
df_list <- lapply(selected_paths, function(pw) {
  es <- calc_running_es(ranked, pathways_list[[pw]])
  data.frame(rank = seq_along(ranked),
             ES = es,
             ESpos = es - min(es),
             pathway = pw,
             pathway_display = display_names[pw],
             NES = selected$NES[match(pw, selected$pathway)],
             neglog10pval = neglog10pval[pw],
             stringsAsFactors = FALSE)
})
df <- bind_rows(df_list)
pos_order <- selected %>% arrange(desc(NES)) %>% pull(pathway)
df$pathway_display <- factor(df$pathway_display, levels = display_names[pos_order])
max_chars <- 40
df$pathway_display_wrapped <- str_wrap(df$pathway_display, width = max_chars)
p <- ggplot(df, aes(x = rank, y = pathway_display_wrapped, height = ESpos, fill = neglog10pval)) +
  geom_ridgeline(stat = "identity", scale = 1, colour = "black", size = 0.2) +
  scale_fill_gradient(low = "#4DBBD5", high = "#E64B35", name = expression(-log[10](italic(p) - value))) +
  labs(x = "Gene rank (by correlation)", y = "", title = bquote(italic(.(target_gene))~"GSEA")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 15, face = "bold")
  )

ggsave(file.path(output,paste0(labels,"PRB2_GSEA_results01.png")),p,width = 6.5,height = 5)
ggsave(file.path(output,paste0(labels,"PRB2_GSEA_results01.pdf")),p,width = 6.5,height = 5)


####SERPINA3####
target_gene <- "SERPINA3"

target_vec <- expr[target_gene, , drop = TRUE]
cors <- apply(expr, 1, function(x) {
  suppressWarnings(cor(x, target_vec, method = "spearman", use = "pairwise.complete.obs"))
})

ranked <- sort(cors, decreasing = TRUE)

c5bp <- tryCatch(
  msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP"),
  error = function(e) NULL
)
pathways_list <- list()
if(!is.null(c5bp)) {
  c5_sets <- split(c5bp$gene_symbol, c5bp$gs_name)
  names(c5_sets) <- paste0("C5_BP|", names(c5_sets))
  pathways_list <- c(pathways_list, c5_sets)
}
set.seed(123)
fgsea_res <- fgsea(pathways = pathways_list,
                   stats = ranked,
                   minSize = 15,
                   maxSize = 500,
                   nperm = 10000)
fgsea_res <- fgsea_res %>%
  as.data.frame() %>%
  arrange(padj, pval, desc(NES))
safe_fgsea <- fgsea_res
list_cols <- names(safe_fgsea)[sapply(safe_fgsea, is.list)]
for (col in list_cols) {
  safe_fgsea[[col]] <- sapply(safe_fgsea[[col]], function(x) {
    if (is.null(x)) return(NA_character_)
    if (length(x) == 0) return(NA_character_)
    paste(as.character(x), collapse = ";")
  }, USE.NAMES = FALSE)
}
safe_fgsea <- as.data.frame(safe_fgsea, stringsAsFactors = FALSE)
write.csv(safe_fgsea, file = file.path(output, paste0(labels, target_gene,"_GSEA_results01.csv")), row.names = FALSE)
library(ggplot2)
library(ggridges)
library(dplyr)

sig <- fgsea_res %>% filter(pval < 0.01)
if(nrow(sig) == 0) stop("No pathways with pval < 0.01")

#top_pos <- sig %>% arrange(desc(NES)) %>% head(5)
#top_neg <- sig %>% arrange(NES) %>% head(5)
#selected <- bind_rows(top_pos, top_neg) %>% distinct(pathway, .keep_all = TRUE)
selected <- sig[c(80,75,51,78,100,371,454,178,180,176),]
selected_paths <- selected$pathway
display_names <- sub("^C5_BP\\|GOBP_", "", selected_paths)
display_names <- sub("^C5_BP\\|", "", display_names)
display_names <- str_replace_all(str_to_lower(display_names), "_", " ")
names(display_names) <- selected_paths
eps <- 1e-300
neglog10pval <- -log10(selected$pval + eps)
names(neglog10pval) <- selected$pathway
calc_running_es <- function(stats, geneset) {
  N <- length(stats)
  selected_idx <- which(names(stats) %in% geneset)
  Nh <- length(selected_idx)
  if(Nh == 0) return(rep(0, N))
  weights <- abs(stats)^1
  hit <- integer(N); hit[selected_idx] <- 1
  Phit <- cumsum(hit * weights) / sum(weights[hit == 1])
  Pmiss <- cumsum((1 - hit) / (N - Nh))
  Phit - Pmiss
}
df_list <- lapply(selected_paths, function(pw) {
  es <- calc_running_es(ranked, pathways_list[[pw]])
  data.frame(rank = seq_along(ranked),
             ES = es,
             ESpos = es - min(es),
             pathway = pw,
             pathway_display = display_names[pw],
             NES = selected$NES[match(pw, selected$pathway)],
             neglog10pval = neglog10pval[pw],
             stringsAsFactors = FALSE)
})
df <- bind_rows(df_list)
pos_order <- selected %>% arrange(desc(NES)) %>% pull(pathway)
df$pathway_display <- factor(df$pathway_display, levels = display_names[pos_order])

max_chars <- 40

df$pathway_display_wrapped <- str_wrap(df$pathway_display, width = max_chars)
p <- ggplot(df, aes(x = rank, y = pathway_display_wrapped, height = ESpos, fill = neglog10pval)) +
  geom_ridgeline(stat = "identity", scale = 1, colour = "black", size = 0.2) +
  scale_fill_gradient(low = "#4DBBD5", high = "#E64B35", name = expression(-log[10](italic(p) - value))) +
  labs(x = "Gene rank (by correlation)", y = "", title = bquote(italic(.(target_gene))~"GSEA")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 15, face = "bold")
  )
ggsave(file.path(output,paste0(labels,target_gene,"_GSEA_results01.png")),p,width = 6.5,height = 5)
ggsave(file.path(output,paste0(labels,target_gene,"_GSEA_results01.pdf")),p,width = 6.5,height = 5)
####RNASE4####
target_gene <- "RNASE4"

target_vec <- expr[target_gene, , drop = TRUE]
cors <- apply(expr, 1, function(x) {
  suppressWarnings(cor(x, target_vec, method = "spearman", use = "pairwise.complete.obs"))
})

ranked <- sort(cors, decreasing = TRUE)

c5bp <- tryCatch(
  msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP"),
  error = function(e) NULL
)
pathways_list <- list()
if(!is.null(c5bp)) {
  c5_sets <- split(c5bp$gene_symbol, c5bp$gs_name)
  names(c5_sets) <- paste0("C5_BP|", names(c5_sets))
  pathways_list <- c(pathways_list, c5_sets)
}
set.seed(123)
fgsea_res <- fgsea(pathways = pathways_list,
                   stats = ranked,
                   minSize = 15,
                   maxSize = 500,
                   nperm = 10000)
fgsea_res <- fgsea_res %>%
  as.data.frame() %>%
  arrange(padj, pval, desc(NES))
safe_fgsea <- fgsea_res
list_cols <- names(safe_fgsea)[sapply(safe_fgsea, is.list)]
for (col in list_cols) {
  safe_fgsea[[col]] <- sapply(safe_fgsea[[col]], function(x) {
    if (is.null(x)) return(NA_character_)
    if (length(x) == 0) return(NA_character_)
    paste(as.character(x), collapse = ";")
  }, USE.NAMES = FALSE)
}
safe_fgsea <- as.data.frame(safe_fgsea, stringsAsFactors = FALSE)
write.csv(safe_fgsea, file = file.path(output, paste0(labels, target_gene,"_GSEA_results01.csv")), row.names = FALSE)
library(ggplot2)
library(ggridges)
library(dplyr)

sig <- fgsea_res %>% filter(pval < 0.01)
if(nrow(sig) == 0) stop("No pathways with pval < 0.01")

#top_pos <- sig %>% arrange(desc(NES)) %>% head(5)
#top_neg <- sig %>% arrange(NES) %>% head(5)
#selected <- bind_rows(top_pos, top_neg) %>% distinct(pathway, .keep_all = TRUE)
selected <- sig[c(49,82,73,55,87,230,255,247,155,194),]
selected_paths <- selected$pathway
display_names <- sub("^C5_BP\\|GOBP_", "", selected_paths)
display_names <- sub("^C5_BP\\|", "", display_names)
display_names <- str_replace_all(str_to_lower(display_names), "_", " ")
names(display_names) <- selected_paths
eps <- 1e-300
neglog10pval <- -log10(selected$pval + eps)
names(neglog10pval) <- selected$pathway
calc_running_es <- function(stats, geneset) {
  N <- length(stats)
  selected_idx <- which(names(stats) %in% geneset)
  Nh <- length(selected_idx)
  if(Nh == 0) return(rep(0, N))
  weights <- abs(stats)^1
  hit <- integer(N); hit[selected_idx] <- 1
  Phit <- cumsum(hit * weights) / sum(weights[hit == 1])
  Pmiss <- cumsum((1 - hit) / (N - Nh))
  Phit - Pmiss
}
df_list <- lapply(selected_paths, function(pw) {
  es <- calc_running_es(ranked, pathways_list[[pw]])
  data.frame(rank = seq_along(ranked),
             ES = es,
             ESpos = es - min(es),
             pathway = pw,
             pathway_display = display_names[pw],
             NES = selected$NES[match(pw, selected$pathway)],
             neglog10pval = neglog10pval[pw],
             stringsAsFactors = FALSE)
})
df <- bind_rows(df_list)
pos_order <- selected %>% arrange(desc(NES)) %>% pull(pathway)
df$pathway_display <- factor(df$pathway_display, levels = display_names[pos_order])

max_chars <- 40 

df$pathway_display_wrapped <- str_wrap(df$pathway_display, width = max_chars)
p <- ggplot(df, aes(x = rank, y = pathway_display_wrapped, height = ESpos, fill = neglog10pval)) +
  geom_ridgeline(stat = "identity", scale = 1, colour = "black", size = 0.2) +
  scale_fill_gradient(low = "#4DBBD5", high = "#E64B35", name = expression(-log[10](italic(p) - value))) +
  labs(x = "Gene rank (by correlation)", y = "", title = bquote(italic(.(target_gene))~"GSEA")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 15, face = "bold")
  )
ggsave(file.path(output,paste0(labels,target_gene,"_GSEA_results01.png")),p,width = 6.5,height = 5)
ggsave(file.path(output,paste0(labels,target_gene,"_GSEA_results01.pdf")),p,width = 6.5,height = 5)

####SOX9####
target_gene <- "SOX9"

target_vec <- expr[target_gene, , drop = TRUE]
cors <- apply(expr, 1, function(x) {
  suppressWarnings(cor(x, target_vec, method = "spearman", use = "pairwise.complete.obs"))
})

ranked <- sort(cors, decreasing = TRUE)

c5bp <- tryCatch(
  msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP"),
  error = function(e) NULL
)
pathways_list <- list()
if(!is.null(c5bp)) {
  c5_sets <- split(c5bp$gene_symbol, c5bp$gs_name)
  names(c5_sets) <- paste0("C5_BP|", names(c5_sets))
  pathways_list <- c(pathways_list, c5_sets)
}
set.seed(123)
fgsea_res <- fgsea(pathways = pathways_list,
                   stats = ranked,
                   minSize = 15,
                   maxSize = 500,
                   nperm = 10000)
fgsea_res <- fgsea_res %>%
  as.data.frame() %>%
  arrange(padj, pval, desc(NES))
safe_fgsea <- fgsea_res
list_cols <- names(safe_fgsea)[sapply(safe_fgsea, is.list)]
for (col in list_cols) {
  safe_fgsea[[col]] <- sapply(safe_fgsea[[col]], function(x) {
    if (is.null(x)) return(NA_character_)
    if (length(x) == 0) return(NA_character_)
    paste(as.character(x), collapse = ";")
  }, USE.NAMES = FALSE)
}
safe_fgsea <- as.data.frame(safe_fgsea, stringsAsFactors = FALSE)
write.csv(safe_fgsea, file = file.path(output, paste0(labels, target_gene,"_GSEA_results01.csv")), row.names = FALSE)
library(ggplot2)
library(ggridges)
library(dplyr)

sig <- fgsea_res %>% filter(pval < 0.01)
if(nrow(sig) == 0) stop("No pathways with pval < 0.01")

#top_pos <- sig %>% arrange(desc(NES)) %>% head(5)
#top_neg <- sig %>% arrange(NES) %>% head(5)
#selected <- bind_rows(top_pos, top_neg) %>% distinct(pathway, .keep_all = TRUE)
selected <- sig[c(57,156,68,234,75,31,63,98,15,53),]
selected_paths <- selected$pathway
display_names <- sub("^C5_BP\\|GOBP_", "", selected_paths)
display_names <- sub("^C5_BP\\|", "", display_names)
display_names <- str_replace_all(str_to_lower(display_names), "_", " ")
names(display_names) <- selected_paths
eps <- 1e-300
neglog10pval <- -log10(selected$pval + eps)
names(neglog10pval) <- selected$pathway
calc_running_es <- function(stats, geneset) {
  N <- length(stats)
  selected_idx <- which(names(stats) %in% geneset)
  Nh <- length(selected_idx)
  if(Nh == 0) return(rep(0, N))
  weights <- abs(stats)^1
  hit <- integer(N); hit[selected_idx] <- 1
  Phit <- cumsum(hit * weights) / sum(weights[hit == 1])
  Pmiss <- cumsum((1 - hit) / (N - Nh))
  Phit - Pmiss
}
df_list <- lapply(selected_paths, function(pw) {
  es <- calc_running_es(ranked, pathways_list[[pw]])
  data.frame(rank = seq_along(ranked),
             ES = es,
             ESpos = es - min(es),
             pathway = pw,
             pathway_display = display_names[pw],
             NES = selected$NES[match(pw, selected$pathway)],
             neglog10pval = neglog10pval[pw],
             stringsAsFactors = FALSE)
})
df <- bind_rows(df_list)
pos_order <- selected %>% arrange(desc(NES)) %>% pull(pathway)
df$pathway_display <- factor(df$pathway_display, levels = display_names[pos_order])
max_chars <- 40 
df$pathway_display_wrapped <- str_wrap(df$pathway_display, width = max_chars)
p <- ggplot(df, aes(x = rank, y = pathway_display_wrapped, height = ESpos, fill = neglog10pval)) +
  geom_ridgeline(stat = "identity", scale = 1, colour = "black", size = 0.2) +
  scale_fill_gradient(low = "#4DBBD5", high = "#E64B35", name = expression(-log[10](italic(p) - value))) +
  labs(x = "Gene rank (by correlation)", y = "", title = bquote(italic(.(target_gene))~"GSEA")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 15, face = "bold")
  )
ggsave(file.path(output,paste0(labels,target_gene,"_GSEA_results01.png")),p,width = 6.5,height = 5)
ggsave(file.path(output,paste0(labels,target_gene,"_GSEA_results01.pdf")),p,width = 6.5,height = 5)

####TUBA1C####

target_gene <- "TUBA1C"

target_vec <- expr[target_gene, , drop = TRUE]
cors <- apply(expr, 1, function(x) {
  suppressWarnings(cor(x, target_vec, method = "spearman", use = "pairwise.complete.obs"))
})

ranked <- sort(cors, decreasing = TRUE)

c5bp <- tryCatch(
  msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP"),
  error = function(e) NULL
)
pathways_list <- list()
if(!is.null(c5bp)) {
  c5_sets <- split(c5bp$gene_symbol, c5bp$gs_name)
  names(c5_sets) <- paste0("C5_BP|", names(c5_sets))
  pathways_list <- c(pathways_list, c5_sets)
}
set.seed(123)
fgsea_res <- fgsea(pathways = pathways_list,
                   stats = ranked,
                   minSize = 15,
                   maxSize = 500,
                   nperm = 10000)
fgsea_res <- fgsea_res %>%
  as.data.frame() %>%
  arrange(padj, pval, desc(NES))
safe_fgsea <- fgsea_res
list_cols <- names(safe_fgsea)[sapply(safe_fgsea, is.list)]
for (col in list_cols) {
  safe_fgsea[[col]] <- sapply(safe_fgsea[[col]], function(x) {
    if (is.null(x)) return(NA_character_)
    if (length(x) == 0) return(NA_character_)
    paste(as.character(x), collapse = ";")
  }, USE.NAMES = FALSE)
}
safe_fgsea <- as.data.frame(safe_fgsea, stringsAsFactors = FALSE)
write.csv(safe_fgsea, file = file.path(output, paste0(labels, target_gene,"_GSEA_results01.csv")), row.names = FALSE)
library(ggplot2)
library(ggridges)
library(dplyr)

sig <- fgsea_res %>% filter(pval < 0.01)
if(nrow(sig) == 0) stop("No pathways with pval < 0.01")

top_pos <- sig %>% arrange(desc(NES)) %>% head(5)
top_neg <- sig %>% arrange(NES) %>% head(5)
selected <- bind_rows(top_pos, top_neg) %>% distinct(pathway, .keep_all = TRUE)
selected_paths <- selected$pathway
display_names <- sub("^C5_BP\\|GOBP_", "", selected_paths)
display_names <- sub("^C5_BP\\|", "", display_names)
display_names <- str_replace_all(str_to_lower(display_names), "_", " ")
names(display_names) <- selected_paths
eps <- 1e-300
neglog10pval <- -log10(selected$pval + eps)
names(neglog10pval) <- selected$pathway
calc_running_es <- function(stats, geneset) {
  N <- length(stats)
  selected_idx <- which(names(stats) %in% geneset)
  Nh <- length(selected_idx)
  if(Nh == 0) return(rep(0, N))
  weights <- abs(stats)^1
  hit <- integer(N); hit[selected_idx] <- 1
  Phit <- cumsum(hit * weights) / sum(weights[hit == 1])
  Pmiss <- cumsum((1 - hit) / (N - Nh))
  Phit - Pmiss
}
df_list <- lapply(selected_paths, function(pw) {
  es <- calc_running_es(ranked, pathways_list[[pw]])
  data.frame(rank = seq_along(ranked),
             ES = es,
             ESpos = es - min(es),
             pathway = pw,
             pathway_display = display_names[pw],
             NES = selected$NES[match(pw, selected$pathway)],
             neglog10pval = neglog10pval[pw],
             stringsAsFactors = FALSE)
})
df <- bind_rows(df_list)
pos_order <- selected %>% arrange(desc(NES)) %>% pull(pathway)
df$pathway_display <- factor(df$pathway_display, levels = display_names[pos_order])
max_chars <- 40 
df$pathway_display_wrapped <- str_wrap(df$pathway_display, width = max_chars)
p <- ggplot(df, aes(x = rank, y = pathway_display_wrapped, height = ESpos, fill = neglog10pval)) +
  geom_ridgeline(stat = "identity", scale = 1, colour = "black", size = 0.2) +
  scale_fill_gradient(low = "#4DBBD5", high = "#E64B35", name = expression(-log[10](italic(p) - value))) +
  labs(x = "Gene rank (by correlation)", y = "", title = bquote(italic(.(target_gene))~"GSEA")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 15, face = "bold")
  )

ggsave(file.path(output,paste0(labels,target_gene,"_GSEA_results01.png")),p,width = 6.5,height = 5)
ggsave(file.path(output,paste0(labels,target_gene,"_GSEA_results01.pdf")),p,width = 6.5,height = 5)

####TLN1####

target_gene <- "TLN1"

target_vec <- expr[target_gene, , drop = TRUE]
cors <- apply(expr, 1, function(x) {
  suppressWarnings(cor(x, target_vec, method = "spearman", use = "pairwise.complete.obs"))
})

ranked <- sort(cors, decreasing = TRUE)

c5bp <- tryCatch(
  msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP"),
  error = function(e) NULL
)
pathways_list <- list()
if(!is.null(c5bp)) {
  c5_sets <- split(c5bp$gene_symbol, c5bp$gs_name)
  names(c5_sets) <- paste0("C5_BP|", names(c5_sets))
  pathways_list <- c(pathways_list, c5_sets)
}
set.seed(123)
fgsea_res <- fgsea(pathways = pathways_list,
                   stats = ranked,
                   minSize = 15,
                   maxSize = 500,
                   nperm = 10000)
fgsea_res <- fgsea_res %>%
  as.data.frame() %>%
  arrange(padj, pval, desc(NES))
safe_fgsea <- fgsea_res
list_cols <- names(safe_fgsea)[sapply(safe_fgsea, is.list)]
for (col in list_cols) {
  safe_fgsea[[col]] <- sapply(safe_fgsea[[col]], function(x) {
    if (is.null(x)) return(NA_character_)
    if (length(x) == 0) return(NA_character_)
    paste(as.character(x), collapse = ";")
  }, USE.NAMES = FALSE)
}
safe_fgsea <- as.data.frame(safe_fgsea, stringsAsFactors = FALSE)
write.csv(safe_fgsea, file = file.path(output, paste0(labels, target_gene,"_GSEA_results01.csv")), row.names = FALSE)
library(ggplot2)
library(ggridges)
library(dplyr)

sig <- fgsea_res %>% filter(pval < 0.01)
if(nrow(sig) == 0) stop("No pathways with pval < 0.01")

top_pos <- sig %>% arrange(desc(NES)) %>% head(5)
top_neg <- sig %>% arrange(NES) %>% head(5)
selected <- bind_rows(top_pos, top_neg) %>% distinct(pathway, .keep_all = TRUE)
selected_paths <- selected$pathway
display_names <- sub("^C5_BP\\|GOBP_", "", selected_paths)
display_names <- sub("^C5_BP\\|", "", display_names)
display_names <- str_replace_all(str_to_lower(display_names), "_", " ")
names(display_names) <- selected_paths
eps <- 1e-300
neglog10pval <- -log10(selected$pval + eps)
names(neglog10pval) <- selected$pathway
calc_running_es <- function(stats, geneset) {
  N <- length(stats)
  selected_idx <- which(names(stats) %in% geneset)
  Nh <- length(selected_idx)
  if(Nh == 0) return(rep(0, N))
  weights <- abs(stats)^1
  hit <- integer(N); hit[selected_idx] <- 1
  Phit <- cumsum(hit * weights) / sum(weights[hit == 1])
  Pmiss <- cumsum((1 - hit) / (N - Nh))
  Phit - Pmiss
}
df_list <- lapply(selected_paths, function(pw) {
  es <- calc_running_es(ranked, pathways_list[[pw]])
  data.frame(rank = seq_along(ranked),
             ES = es,
             ESpos = es - min(es),
             pathway = pw,
             pathway_display = display_names[pw],
             NES = selected$NES[match(pw, selected$pathway)],
             neglog10pval = neglog10pval[pw],
             stringsAsFactors = FALSE)
})
df <- bind_rows(df_list)
pos_order <- selected %>% arrange(desc(NES)) %>% pull(pathway)
df$pathway_display <- factor(df$pathway_display, levels = display_names[pos_order])
max_chars <- 40 
df$pathway_display_wrapped <- str_wrap(df$pathway_display, width = max_chars)
p <- ggplot(df, aes(x = rank, y = pathway_display_wrapped, height = ESpos, fill = neglog10pval)) +
  geom_ridgeline(stat = "identity", scale = 1, colour = "black", size = 0.2) +
  scale_fill_gradient(low = "#4DBBD5", high = "#E64B35", name = expression(-log[10](italic(p) - value))) +
  labs(x = "Gene rank (by correlation)", y = "", title = bquote(italic(.(target_gene))~"GSEA")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 15, face = "bold")
  )


ggsave(file.path(output,paste0(labels,target_gene,"_GSEA_results01.png")),p,width = 6.5,height = 5)
ggsave(file.path(output,paste0(labels,target_gene,"_GSEA_results01.pdf")),p,width = 6.5,height = 5)


####LRRC32####
target_gene <- "LRRC32"

target_vec <- expr[target_gene, , drop = TRUE]
cors <- apply(expr, 1, function(x) {
  suppressWarnings(cor(x, target_vec, method = "spearman", use = "pairwise.complete.obs"))
})

ranked <- sort(cors, decreasing = TRUE)

c5bp <- tryCatch(
  msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP"),
  error = function(e) NULL
)
pathways_list <- list()
if(!is.null(c5bp)) {
  c5_sets <- split(c5bp$gene_symbol, c5bp$gs_name)
  names(c5_sets) <- paste0("C5_BP|", names(c5_sets))
  pathways_list <- c(pathways_list, c5_sets)
}
set.seed(123)
fgsea_res <- fgsea(pathways = pathways_list,
                   stats = ranked,
                   minSize = 15,
                   maxSize = 500,
                   nperm = 10000)
fgsea_res <- fgsea_res %>%
  as.data.frame() %>%
  arrange(padj, pval, desc(NES))
safe_fgsea <- fgsea_res
list_cols <- names(safe_fgsea)[sapply(safe_fgsea, is.list)]
for (col in list_cols) {
  safe_fgsea[[col]] <- sapply(safe_fgsea[[col]], function(x) {
    if (is.null(x)) return(NA_character_)
    if (length(x) == 0) return(NA_character_)
    paste(as.character(x), collapse = ";")
  }, USE.NAMES = FALSE)
}
safe_fgsea <- as.data.frame(safe_fgsea, stringsAsFactors = FALSE)
write.csv(safe_fgsea, file = file.path(output, paste0(labels, target_gene,"_GSEA_results01.csv")), row.names = FALSE)
library(ggplot2)
library(ggridges)
library(dplyr)

sig <- fgsea_res %>% filter(pval < 0.01)
if(nrow(sig) == 0) stop("No pathways with pval < 0.01")

top_pos <- sig %>% arrange(desc(NES)) %>% head(5)
top_neg <- sig %>% arrange(NES) %>% head(5)
selected <- bind_rows(top_pos, top_neg) %>% distinct(pathway, .keep_all = TRUE)
selected_paths <- selected$pathway
display_names <- sub("^C5_BP\\|GOBP_", "", selected_paths)
display_names <- sub("^C5_BP\\|", "", display_names)
display_names <- str_replace_all(str_to_lower(display_names), "_", " ")
names(display_names) <- selected_paths
eps <- 1e-300
neglog10pval <- -log10(selected$pval + eps)
names(neglog10pval) <- selected$pathway
calc_running_es <- function(stats, geneset) {
  N <- length(stats)
  selected_idx <- which(names(stats) %in% geneset)
  Nh <- length(selected_idx)
  if(Nh == 0) return(rep(0, N))
  weights <- abs(stats)^1
  hit <- integer(N); hit[selected_idx] <- 1
  Phit <- cumsum(hit * weights) / sum(weights[hit == 1])
  Pmiss <- cumsum((1 - hit) / (N - Nh))
  Phit - Pmiss
}
df_list <- lapply(selected_paths, function(pw) {
  es <- calc_running_es(ranked, pathways_list[[pw]])
  data.frame(rank = seq_along(ranked),
             ES = es,
             ESpos = es - min(es),
             pathway = pw,
             pathway_display = display_names[pw],
             NES = selected$NES[match(pw, selected$pathway)],
             neglog10pval = neglog10pval[pw],
             stringsAsFactors = FALSE)
})
df <- bind_rows(df_list)
pos_order <- selected %>% arrange(desc(NES)) %>% pull(pathway)
df$pathway_display <- factor(df$pathway_display, levels = display_names[pos_order])
max_chars <- 40 
df$pathway_display_wrapped <- str_wrap(df$pathway_display, width = max_chars)
p <- ggplot(df, aes(x = rank, y = pathway_display_wrapped, height = ESpos, fill = neglog10pval)) +
  geom_ridgeline(stat = "identity", scale = 1, colour = "black", size = 0.2) +
  scale_fill_gradient(low = "#4DBBD5", high = "#E64B35", name = expression(-log[10](italic(p) - value))) +
  labs(x = "Gene rank (by correlation)", y = "", title = bquote(italic(.(target_gene))~"GSEA")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 15, face = "bold")
  )


ggsave(file.path(output,paste0(labels,target_gene,"_GSEA_results01.png")),p,width = 6.5,height = 5)
ggsave(file.path(output,paste0(labels,target_gene,"_GSEA_results01.pdf")),p,width = 6.5,height = 5)

####NPC2####
target_gene <- "NPC2"

target_vec <- expr[target_gene, , drop = TRUE]
cors <- apply(expr, 1, function(x) {
  suppressWarnings(cor(x, target_vec, method = "spearman", use = "pairwise.complete.obs"))
})

ranked <- sort(cors, decreasing = TRUE)

c5bp <- tryCatch(
  msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP"),
  error = function(e) NULL
)
pathways_list <- list()
if(!is.null(c5bp)) {
  c5_sets <- split(c5bp$gene_symbol, c5bp$gs_name)
  names(c5_sets) <- paste0("C5_BP|", names(c5_sets))
  pathways_list <- c(pathways_list, c5_sets)
}
set.seed(123)
fgsea_res <- fgsea(pathways = pathways_list,
                   stats = ranked,
                   minSize = 15,
                   maxSize = 500,
                   nperm = 10000)
fgsea_res <- fgsea_res %>%
  as.data.frame() %>%
  arrange(padj, pval, desc(NES))
safe_fgsea <- fgsea_res
list_cols <- names(safe_fgsea)[sapply(safe_fgsea, is.list)]
for (col in list_cols) {
  safe_fgsea[[col]] <- sapply(safe_fgsea[[col]], function(x) {
    if (is.null(x)) return(NA_character_)
    if (length(x) == 0) return(NA_character_)
    paste(as.character(x), collapse = ";")
  }, USE.NAMES = FALSE)
}
safe_fgsea <- as.data.frame(safe_fgsea, stringsAsFactors = FALSE)
write.csv(safe_fgsea, file = file.path(output, paste0(labels, target_gene,"_GSEA_results01.csv")), row.names = FALSE)
library(ggplot2)
library(ggridges)
library(dplyr)

sig <- fgsea_res %>% filter(pval < 0.01)
if(nrow(sig) == 0) stop("No pathways with pval < 0.01")

top_pos <- sig %>% arrange(desc(NES)) %>% head(5)
top_neg <- sig %>% arrange(NES) %>% head(5)
selected <- bind_rows(top_pos, top_neg) %>% distinct(pathway, .keep_all = TRUE)
selected_paths <- selected$pathway
display_names <- sub("^C5_BP\\|GOBP_", "", selected_paths)
display_names <- sub("^C5_BP\\|", "", display_names)
display_names <- str_replace_all(str_to_lower(display_names), "_", " ")
names(display_names) <- selected_paths
eps <- 1e-300
neglog10pval <- -log10(selected$pval + eps)
names(neglog10pval) <- selected$pathway
calc_running_es <- function(stats, geneset) {
  N <- length(stats)
  selected_idx <- which(names(stats) %in% geneset)
  Nh <- length(selected_idx)
  if(Nh == 0) return(rep(0, N))
  weights <- abs(stats)^1
  hit <- integer(N); hit[selected_idx] <- 1
  Phit <- cumsum(hit * weights) / sum(weights[hit == 1])
  Pmiss <- cumsum((1 - hit) / (N - Nh))
  Phit - Pmiss
}
df_list <- lapply(selected_paths, function(pw) {
  es <- calc_running_es(ranked, pathways_list[[pw]])
  data.frame(rank = seq_along(ranked),
             ES = es,
             ESpos = es - min(es),
             pathway = pw,
             pathway_display = display_names[pw],
             NES = selected$NES[match(pw, selected$pathway)],
             neglog10pval = neglog10pval[pw],
             stringsAsFactors = FALSE)
})
df <- bind_rows(df_list)
pos_order <- selected %>% arrange(desc(NES)) %>% pull(pathway)
df$pathway_display <- factor(df$pathway_display, levels = display_names[pos_order])
max_chars <- 40 
df$pathway_display_wrapped <- str_wrap(df$pathway_display, width = max_chars)
p <- ggplot(df, aes(x = rank, y = pathway_display_wrapped, height = ESpos, fill = neglog10pval)) +
  geom_ridgeline(stat = "identity", scale = 1, colour = "black", size = 0.2) +
  scale_fill_gradient(low = "#4DBBD5", high = "#E64B35", name = expression(-log[10](italic(p) - value))) +
  labs(x = "Gene rank (by correlation)", y = "", title = bquote(italic(.(target_gene))~"GSEA")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 15, face = "bold")
  )


ggsave(file.path(output,paste0(labels,target_gene,"_GSEA_results01.png")),p,width = 6.5,height = 5)
ggsave(file.path(output,paste0(labels,target_gene,"_GSEA_results01.pdf")),p,width = 6.5,height = 5)
