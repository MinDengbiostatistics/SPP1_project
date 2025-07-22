## KEGG enrichment for Epithelial cell population 
rm(list=ls())
setwd("F:/9.Others/4.Tingting/GSE161529_rawdata")

## For epithelial cell DE analysis between 'tumor' = 'normal'
library(Seurat)
library(clusterProfiler)

epi.data <- readRDS("epi.wt.er.rds")
Idents(epi.data) <- "subtype"
table(epi.data$subtype)

## Loading inferCNV grouping information
group.df <- read.table("./Epithelial.inferCNV/infercnv.observation_groupings_Epithelial_normal_er.txt")
group.df.2 <- group.df[order(nrow(group.df):1),]
group.df.2 <- group.df.2
View(group.df.2)

group.df.3 <- group.df.2[1:19461,]

group.df.normal <- subset(group.df.3,group.df.3$Annotation.Group=="2")
group.df.jump <- subset(group.df.3,group.df.3$Annotation.Group=="1")

##Extract cell id
epi.data[["CellName"]] <- colnames(epi.data)
cell.ids <- c(rownames(group.df.normal),rownames(group.df.jump))

normal.tumor.data <- subset(epi.data, subset = CellName %in% cell.ids)
dim(normal.tumor.data)

library(dplyr)
pbmc.markers <- FindAllMarkers(normal.tumor.data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

pbmc.markers.2 <- subset(pbmc.markers,pbmc.markers$p_val_adj<0.05)

## KEGG pathway analysis
library(clusterProfiler)

gene <- bitr(pbmc.markers.2$gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = "org.Hs.eg.db")

set.seed(1)
KEGG <- enrichKEGG(gene$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
#View(KEGG@result)
#write.csv(KEGG@result,'KEGG_Enrichment_Epi.csv',quote = F)

library(ggplot2)
library(stringr)

y <- as.data.frame(KEGG@result)
#View(y)
ref.ID <- c("PI3K-Akt signaling pathway","MAPK signaling pathway","Transcriptional misregulation in cancer",
            "Regulation of actin cytoskeleton","IL-17 signaling pathway","Osteoclast differentiation",
            "Osteoclast differentiation","HIF-1 signaling pathway","Prostate cancer","Th1 and Th2 cell differentiation",
            "PD-L1 and PD-1 checkpoint in cancer","Antigen processing and presentation","Toll-like receptor signaling pathway",
            "EGFR tyrosine kinase inhibitor resistance")
kegg.df <- subset(y,y$Description%in%ref.ID)
kegg.df <- subset(kegg.df,kegg.df$pvalue<0.05)
kegg.df$number <- factor(kegg.df$Description,levels = rev(kegg.df$Description))
p.kegg <- ggplot(data = kegg.df, # 绘图使用的数据
                        aes(x = Count,y = reorder(number,Count)))+ #横纵坐标及排序
         geom_point(aes(size = Count,color = pvalue),size=8)+ # 气泡大小及颜色设置
         theme_bw()+ # 去除背景色++
         scale_colour_gradient(low = "red",high = "purple")+ # 设置气泡渐变色
         labs(x = "Gene Number", y = "",title = "Enriched and Overlapped KEGG Pathways in Epithelial-ER+ cells", # 设置坐标轴标题及图标题
              color = expression(p.adjust),size = "Count")+ # 设置图例颜色及大小
         scale_y_discrete(labels = function(x) str_wrap(x,width = 45))+
         theme(axis.title = element_text(size = 15),
               axis.text = element_text(size = 15,face="bold"),
               plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
               legend.title = element_text(size = 13),legend.text = element_text(size = 11),
               plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))
ggsave(p.kegg,filename = "KEGG_dotplot_Epithelial-ER+.pdf",height = 8,width = 10)

