rm(list=ls())
setwd("/home/data/dengmin/02.Xu_lab/tingting")
library(Seurat)
library(dplyr)
library(patchwork)
library(devtools)
library(ggplot2)
library(tidydr)

## how to read multiple rds data in Folder Public_metasasis
## 1. read all rds files in the folder
rds_files <- list.files(path = "/home/data/dengmin/02.Xu_lab/tingting/Public_metastasis", 
                         full.names = TRUE)

## 2. read each rds file and store it in a list
rds_list <- lapply(rds_files, readRDS)

## 3. assign names to the list elements based on the file names
rds_names <- gsub("\\.rds$", "", basename(rds_files))
names(rds_list) <- rds_names

## 4. combine the list of Seurat objects into a single Seurat object
total_seurat <- Reduce(function(x, y) merge(x, y), rds_list)
dim(total_seurat)

#####################################
## Convert ensemble id into symbol
library(AnnotationDbi)
library(org.Hs.eg.db)  

# 转换 Ensembl ID 到 Symbol
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(total_seurat),
  keytype = "ENSEMBL",
  column = "SYMBOL"
)
# 转换为数据框
gene_info <- data.frame(
  ensembl_gene_id = names(gene_symbols),
  external_gene_name = unname(gene_symbols)
)

# 创建映射字典（Ensembl ID → Symbol）
id_to_symbol <- gene_info$external_gene_name
names(id_to_symbol) <- gene_info$ensembl_gene_id
# 替换 Seurat 对象的基因名
new_symbols <- id_to_symbol[rownames(total_seurat)]
rownames(total_seurat@assays$RNA@counts) <- new_symbols
rownames(total_seurat@assays$RNA@data) <- new_symbols

# 删除无 Symbol 的基因
total_seurat <- total_seurat[!is.na(new_symbols), ]
dim(total_seurat)

# 检查转换后的基因名
head(rownames(total_seurat), 10)
# 示例输出：SAMD11, NOC2L, HES4, ...
# 检查是否有重复 Symbol
sum(duplicated(rownames(total_seurat)))  # 应为0

## 5.1 import Normal sample
other.seurat <- readRDS("V2.merge.data.rds")

table(other.seurat$subtype)
Idents(other.seurat) <- "subtype"
normal.seurat <- subset(other.seurat, idents = c("Normal"))

combined_seurat <- merge(total_seurat, y = normal.seurat, add.cell.ids = c("Metastasis", "Normal"), project = "Combined")
combined_seurat

## 5.2 save the combined Seurat object
#saveRDS(combined_seurat, file = "/home/data/dengmin/02.Xu_lab/tingting/Public_metastasis/combined_seurat.rds")

## 6. load the combined Seurat object
#combined_seurat <- readRDS("/home/data/dengmin/02.Xu_lab/tingting/Public_metastasis/combined_seurat.rds")

combined_seurat$tissue <- ifelse(combined_seurat$subtype %in% NA, combined_seurat$tissue,
                          combined_seurat$subtype)

combined_seurat$cell.id <- rownames(combined_seurat@meta.data)

combined_seurat@meta.data <- combined_seurat@meta.data[, c("tissue","cell.id")]
dim(combined_seurat@meta.data)

combined_seurat[["percent.mt"]] <- PercentageFeatureSet(combined_seurat, pattern = "^MT-",
                                                        assay = "RNA")

range(combined_seurat$percent.mt)

combined_seurat$nFeature_RNA <- Matrix::colSums(combined_seurat@assays$RNA@counts > 0)
combined_seurat$nCount_RNA <- Matrix::colSums(combined_seurat@assays$RNA@counts)

# Add these metrics to the metadata
combined_seurat <- AddMetaData(combined_seurat, 
                               metadata = combined_seurat$nFeature_RNA, col.name = "nFeature_RNA")
combined_seurat <- AddMetaData(combined_seurat, metadata = combined_seurat$nCount_RNA, 
                               col.name = "nCount_RNA")
combined_seurat <- subset(combined_seurat, nFeature_RNA >= 500)
range(combined_seurat$nFeature_RNA)

combined_seurat$batch <- ifelse(combined_seurat$tissue %in% "Normal","batch2","batch1")
combined_seurat <- AddMetaData(combined_seurat, 
                               metadata = combined_seurat$batch, col.name = "batch")

saveRDS(combined_seurat, file = "combined_seurat.rds")




