rm(list=ls())
setwd("/home/data/dengmin/02.Xu_lab/tingting")
library(Seurat)
library(dplyr)
library(patchwork)
library(devtools)
library(ggplot2)
library(igraph)
library(harmony)
#remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 64 * 1024^3) 

data2 <- readRDS("combined_seurat.rds")

## 2.remove or not remove batch effect
data2[["RNA"]] <- split(data2[["RNA"]], f = data2$batch)
data2
data2 <- NormalizeData(object = data2)
data2 <- FindVariableFeatures(object = data2)
data2 <- ScaleData(object = data2)
data2 <- RunPCA(object = data2, do.print = FALSE)

data2 <- FindNeighbors(data2, dims = 1:30, reduction = "pca")
data2 <- FindClusters(data2, resolution = 0.8, cluster.name = "unintegrated_clusters")
data2 <- RunUMAP(data2, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

DimPlot(data2, reduction = "umap.unintegrated")

## 3.Remove batch effect using harmony and UMAP generation
data2 <- IntegrateLayers(
  object = data2, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

data2 <- FindNeighbors(data2, reduction = "harmony", dims = 1:30)
data2 <- FindClusters(data2, resolution = 0.5, cluster.name = "harmony_clusters")
data2 <- RunUMAP(data2, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

p1 <- DimPlot(data2,reduction = "umap.harmony", label = T)
p2 <- DimPlot(data2, reduction = "umap.unintegrated", label = T)

ggsave("1.Metastasis_Normal.UMAP.unanno-1.pdf", p1, width = 8, height = 7)
ggsave("1.Metastasis_Normal.UMAP.unanno-2.pdf", p2, width = 8, height = 7)


## 4.DE analysis and annotated_clusters
rb.genes <- rownames(data2)[grep("^RP[SL]",rownames(data2))]
data3 = data2[(!(rownames(data2)  %in% rb.genes))]  

library(presto)
DefaultAssay(data2) <- "RNA"
data3 <- JoinLayers(data2)
Idents(data3) <- "harmony_clusters"
de.markers <- FindAllMarkers(data3, only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0.25)
de.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

write.csv(de.markers, "1.unannotaed_metasta_normal_DEGs.csv",quote = F)


################################################################################
################################################################################
## 5.save the annotated data and then do the inferCNV
Idents(data3) <- "harmony_clusters"
data3 <- RenameIdents(object = data3,
                           '0'="SPP1- macrophage",'1'="SPP1+ macrophage",
                           '2'="SPP1- macrophage",'3'="SPP1- macrophage",
                           '4'="SPP1- macrophage")

data3@meta.data$macro.clustering <- Idents(data3)

saveRDS(data3,"public_metastasis_normal_20250403.rds")



