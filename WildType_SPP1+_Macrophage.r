rm(list=ls())
setwd("/home/data/dengmin/02.Xu_lab/Tingting")
library(Seurat)
library(dplyr)
library(patchwork)
library(devtools)
library(ggplot2)
library(tidydr)

macro.obj <- readRDS("mye.wt.tnbc.rds")

View(macro.obj@meta.data)
table(macro.obj$condition)

###################################################################
### Subclustering macrophage
###################################################################
Idents(macro.obj) <- "condition"
macro.obj <- subset(macro.obj, idents = c("Primary"))
dim(macro.obj)

DefaultAssay(macro.obj) <- "RNA"
macro.obj[['RNA']] <- JoinLayers(macro.obj[['RNA']] )
macro.obj[["RNA"]] <- split(macro.obj[["RNA"]], f = macro.obj$orig.ident)
macro.obj

macro.obj <- NormalizeData(macro.obj)
macro.obj <- FindVariableFeatures(macro.obj)
macro.obj <- ScaleData(macro.obj)
macro.obj <- RunPCA(macro.obj)

ElbowPlot(macro.obj,ndims = 50)

macro.obj = macro.obj %>% RunHarmony("orig.ident", plot_convergence = TRUE)
macro.obj = RunUMAP(macro.obj, reduction = "harmony",dims = 1:20)

macro.obj = FindNeighbors(macro.obj, reduction = "harmony", dims = 1:20)             
macro.obj = FindClusters(macro.obj, resolution = 0.2, algorithm=1)

p1 <- DimPlot(macro.obj)+
      labs(x = "UMAP 1", y = "UMAP 2")+
      tidydr::theme_dr(xlength = 0.2, ylength = 0.2, 
                   arrow=arrow(length = unit(0.2, "inches"), type = "closed"))+
      theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2, hjust = 0.03))
ggsave("1.BRCA1-_unannotated_macrophage.pdf", p1)

p2 <- FeaturePlot(macro.obj, features = c("SPP1","FN1","S100A9"))
p <- p1|p2
p

ggsave("2.BACA1-Macrophage_Feature_plots.pdf", p, width = 14, height = 6)


rb.genes <- rownames(macro.obj)[grep("^RP[SL]",rownames(macro.obj))]
macro.obj2 = macro.obj[(!(rownames(macro.obj)  %in% rb.genes))]  

library(presto)
macro.obj2 <- JoinLayers(macro.obj2)
Idents(macro.obj2) <- "RNA_snn_res.0.2"
de.markers.3 <- FindAllMarkers(macro.obj2, only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0.25)
de.markers.3 %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

de.markers.3 <- subset(de.markers.3)
write.csv(de.markers.3,"2.unannotaed_macrophageS_DEGs.csv",quote = F)


Idents(macro.obj2) <- "RNA_snn_res.0.2"
macro.obj2 <- RenameIdents(object = macro.obj2,
                    '0'="SPP1- macrophage",'1'="SPP1+ macrophage",
                    '2'="SPP1- macrophage",'3'="SPP1- macrophage",
                    '4'="SPP1- macrophage")

macro.obj2@meta.data$macro.clustering <- Idents(macro.obj2)

saveRDS(macro.obj2,"BRCA1-_SPP1+_Macrophage_20250328.rds")


library(presto)
macro.obj2 <- JoinLayers(macro.obj2)
Idents(macro.obj2) <- "macro.clustering"
de.markers.4 <- FindAllMarkers(macro.obj2, only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0.25)
de.markers.4 %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(de.markers.4,"BRCA1-.SPP1+_macrophage.DEGs.csv",quote = F)

macro.colors = c('SPP1+ macrophage' ='#dccce2',
                 'SPP1- macrophage' = '#9ec4be')



p12 <- DimPlot(macro.obj2, cols = macro.colors, pt.size = 3)+
    labs(x = "UMAP 1", y = "UMAP 2")+
    tidydr::theme_dr(xlength = 0.2, ylength = 0.2, 
    arrow=arrow(length = unit(0.2, "inches"), type = "closed"))+
    theme(panel.grid = element_blank(),
    axis.title = element_text(face = 2, hjust = 0.03))
p12
ggsave("3.BRCA1-UMAP_SPP1+_macrophage.pdf", p12, width = 6, height = 5)
