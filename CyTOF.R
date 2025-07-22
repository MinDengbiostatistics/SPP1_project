rm(list=ls())
library(flowCore)
library(readxl)
library(CATALYST)
library(Rcpp)

setwd("/home/data/dengmin/02.Xu_lab/tingting/CyTOF/20250108_CyTOF")

x.1 = list.files('./data/')
x.1 = sapply(x.1, function(x)(x=paste("./data/",x,sep = '')))

md = data.frame(file_name = x.1,
                  sample_id = c('Control-1','Control-2','Control-3',
                               'MG-1','MG-2','MG-3',
                               'PD1-1','PD1-2','PD1-3',
                               'PD1+nhwd-1','PD1+nhwd-2','PD1+nhwd-3',
                               'PD1+nhwd+tas-1','PD1+nhwd+tas-2','PD1+nhwd+tas-3',
                               'PD1+nhwd+tas-4'),
                 condition = c('Control','Control','Control',
                               'MG','MG','MG',
                               'PD1','PD1','PD1',
                               'PD1+nhwd','PD1+nhwd','PD1+nhwd',
                               'PD1+nhwd+tas','PD1+nhwd+tas','PD1+nhwd+tas',
                               'PD1+nhwd+tas'),
                 patient_id = c('p1','p2','p3','p4','p5','p6','p7','p8','p9','p10',
                                'p11','p12','p13','p14','p15','p16'))
md$sample_id <- factor(md$sample_id, levels = md$sample_id)
head(data.frame(md))


fcs_raw <- read.flowSet(x.1, transformation = FALSE, truncate_max_range = FALSE)
fcs_raw

panel = read.csv('Panel.marker.csv')
head(data.frame(panel))

#######################################################################################################################
lineage_markers = data.frame(panel)[,1]
panel_fcs <- pData(parameters(fcs_raw[[1]]))

library(stringr)
panel$fcs_colname = paste0(str_sub(panel$fcs_colname,4,5),str_sub(panel$fcs_colname,1,3),'Di')
panel = panel[panel$fcs_colname %in% panel_fcs$name,]

# 数据转换
## arcsinh transformation and column subsetting
#这里应该是创建一个用于分析的矩阵
fcs <- fsApply(fcs_raw, function(x, cofactor=5){
  colnames(x) <- panel_fcs$name
  expr <- exprs(x)
  expr <- asinh(expr[, panel$fcs_colname] / cofactor)
  exprs(x) <- expr
  x
})
fcs

## 最终得到的是一个纵坐标为抗体（全部的抗体，包括了谱系以及功能标记物）
expr <- fsApply(fcs, exprs)
dim(expr)
#这里进行数据标准化，将全部的数据转变为0到1范围之内
library(matrixStats)
rng <- colQuantiles(expr, probs = c(0.01, 0.99))
expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
expr01[expr01 < 0] <- 0
expr01[expr01 > 1] <- 1

## Generate sample IDs corresponding to each cell in the 'expr' matrix
sample_ids <- rep(md$sample_id, fsApply(fcs_raw, nrow))

## 最终将整个数据框变为3列，第一列为样本名称，第二列为抗体名称，第三列为具体的值
#######################################################################################################################
expr03 = expr
sample_ids <- rep(md$sample_id, fsApply(fcs_raw, nrow))
table(sample_ids)

## Find and skip duplicates
dups <- which(!duplicated(expr03[, panel$fcs_colname]))

## Data subsampling: create indices by sample
inds <- split(1:length(sample_ids), sample_ids)

###############################################
## How many cells to downsample per-sample
tsne_ncells <- pmin(table(sample_ids), 8000)
#tsne_ncells <- pmin(table(sample_ids), 5000)


## Get subsampled indices
set.seed(1234)
tsne_inds <- lapply(names(inds), function(i){
  s <- sample(inds[[i]], tsne_ncells[i], replace = FALSE)
  intersect(s, dups)
})

tsne_inds <- unlist(tsne_inds)
tsne_expr <- expr03[tsne_inds, panel$fcs_colname]

library(stringr)

## Run t-SNE
library(Rtsne)
library(RColorBrewer)

set.seed(1234)
tsne_out <- Rtsne(tsne_expr, check_duplicates = FALSE, pca = T,pca_scale = T, perplexity = 60)
dr <- data.frame(tSNE1 = tsne_out$Y[, 1], tSNE2 = tsne_out$Y[, 2],
                 expr03[tsne_inds, panel$fcs_colname])

dr$sample_id <- sample_ids[tsne_inds]
mm <- match(dr$sample_id, md$sample_id)
dr$condition <- md$condition[mm]

## Clustering
library(FlowSOM)
fsom <- ReadInput(fcs, transform = FALSE, scale = FALSE)
set.seed(1234)
som <- BuildSOM(fsom, colsToUse = lineage_markers)

###################################################################
## Metal-clustering into 15 clusters with ConsensusClusterPlus
library(ConsensusClusterPlus)
codes <- som$map$codes
plot_outdir <- "consensus_plots"

nmc <- 20
mc <- ConsensusClusterPlus(t(codes), maxK = nmc, reps = 1000,
                           pItem = 0.8, pFeature = 1, title = plot_outdir, plot = "png",
                           clusterAlg = "hc", 
                           innerLinkage = "average",
                           finalLinkage = "average",
                           #distance = "minkowski", 
                           distance = "euclidean",
                           seed = 1234)
## Get cluster ids for each cell
code_clustering1 <- mc[[nmc]]$consensusClass
cell_clustering1 <- code_clustering1[som$map$mapping[,1]]

dr$sample_id <- sample_ids[tsne_inds]
mm <- match(dr$sample_id, md$sample_id)
dr$condition <- md$condition[mm]
dr$cell_clustering1 <- factor(cell_clustering1[tsne_inds], levels = 1:nmc)

color_clusters <- c("#DC050C", "#FB8072",
                    "#1965B0", "#addd8e",
                    "#882E72", "#7BAFDE",
                    "#FF7F00", "#ABA300",
                    "#E7298A", "#8494FF",
                    "#00BFC4", "#C77CFF","#FF61CC","green",
                    "#2ca25f","#e34a33","#e6ab02","#ffff33",
                    "#0000FF","#7E6148CC","#91D1C2CC","grey")

## Plot t-SNE colored by clusters
library(ggplot2)
ggp.total <- ggplot(dr, aes(x=tSNE1, y=tSNE2, color=cell_clustering1))+
  geom_point(size=0.8)+theme_bw()+scale_color_manual(values = color_clusters)+
  guides(color=guide_legend(override.aes = list(size = 4), ncol = 2))
ggp.total

ggsave("Cell_8000_TSNE_panel_unannotation.pdf",ggp.total,width = 12, height = 12)


###################################################################################
ggp = ggplot(dr, aes(x = tSNE1, y = tSNE2, color = cell_clustering1))+
  geom_point(size = 0.2)+theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  scale_color_manual(values = color_clusters)+
  guides(color = guide_legend(override.aes = list(size = 5), ncol = 3))
ggp1 <- ggp + facet_wrap(~sample_id)
#ggp1
ggsave("TSNE_unannotation_panel_Bysample.pdf",ggp1,width = 16, height = 12)


##########################################################################################
library(ComplexHeatmap)
library(dplyr)
plot_clustering_heatmap_wrapper <- function(expr,expr01, 
                                            cell_clustering, color_clusters, 
                                            cluster_merging = NULL){
  
  
  # Calculate the median expression
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  expr01_median <- data.frame(expr01, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  
  # Calculate cluster frequencies
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  # Sort the cell clusters with hierarchical clustering
  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_median[, colnames(expr01)])
  rownames(expr_heat) <- expr01_median$cell_clustering
  expr_heat <- expr_heat[,panel$fcs_colname]
  colnames(expr_heat) <- panel$antigen
  

  # Colors for the heatmap
  color_heat <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
  legend_breaks = seq(from = 0, to = 1, by = 0.2)
  labels_row <- paste0(expr01_median$cell_clustering, " (", clustering_prop, "%)")
  
  # Annotation for the original clusters
  annotation_row <- data.frame(Cluster = factor(expr01_median$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  # Annotation for the merged clusters
  if(!is.null(cluster_merging)){
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$Cluster_merging <- cluster_merging$new_cluster 
    color_clusters2 <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters2) <- levels(cluster_merging$new_cluster)
    annotation_colors$Cluster_merging <- color_clusters2
  }
  
  pheatmap(expr_heat, color = color_heat, cluster_cols = FALSE, 
           cluster_rows = cluster_rows, labels_row = labels_row, 
           display_numbers = TRUE, number_color = "black", 
           fontsize = 8, fontsize_number = 6, legend_breaks = legend_breaks,
           annotation_row = annotation_row, annotation_colors = annotation_colors)
}

pdf("Heatmap_panel_unannotation.pdf")
plot_clustering_heatmap_wrapper(expr = expr[, colnames(expr)], 
                                expr01 = expr01[, colnames(expr)], 
                                cell_clustering = cell_clustering1, 
                                color_clusters = color_clusters)
dev.off()


##############################################################################
# annotated cell clusters
# For panel 1
cluster_merging1 <- data.frame("original_cluster"=c(1,2,3,4,5,6,7,8,9,10,
                                                    11,12,13,14,15,16,17,18,19,20),
                               'new_cluster'=c("CD8+PD1+ Tex","CD8T","remove",
                                               "CD8T","CD8T","CD8T",
                                               "CD4T",'remove',"remove",
                                               "remove","CD4T","MDSC",
                                               "remove",'CD4T','MDSC',
                                               "CD4T","CD8T","CD4+PD1+ Tex",
                                               "CD4T","MDSC"))

levels_clusters_merged <- c("CD8+PD1+ Tex","CD8T",
                            "remove","CD4T","MDSC",
                            "CD4+PD1+ Tex")

###########################################################################
cluster_merging1$new_cluster <- factor(cluster_merging1$new_cluster, 
                                       levels = levels_clusters_merged)

## New clustering1m
mm <- match(cell_clustering1, cluster_merging1$original_cluster)
cell_clustering1m <- cluster_merging1$new_cluster[mm]

mm <- match(code_clustering1, cluster_merging1$original_cluster)
code_clustering1m <- cluster_merging1$new_cluster[mm]
dr$cell_clustering1m <- cell_clustering1m[tsne_inds]


#######################################################
dr2 <- subset(dr, !(dr$cell_clustering1m%in% c("remove")))

###############################################################
## For panel
cols = c("MDSC"="#FB8072","CD8+PD1+ Tex"="#1965B0", 
         "CD4+PD1+ Tex"="#e6ab02", "CD8T"="#addd8e","CD4T"="#00bfc4")

library(dplyr)
dr4 <- subset(dr2, select=c("sample_id","condition","cell_clustering1m","tSNE1","tSNE2"))
sample.id.len <- unique(dr4$sample_id)
sample.id.len

dr4.1 <- subset(dr4, dr4$sample_id %in% sample.id.len[1])
dr4.2 <- subset(dr4, dr4$sample_id %in% sample.id.len[2])
dr4.3 <- subset(dr4, dr4$sample_id %in% sample.id.len[3])

dr4.4 <- subset(dr4, dr4$sample_id %in% sample.id.len[4])
dr4.5 <- subset(dr4, dr4$sample_id %in% sample.id.len[5])
dr4.6 <- subset(dr4, dr4$sample_id %in% sample.id.len[6])

dr4.7 <- subset(dr4, dr4$sample_id %in% sample.id.len[7])
dr4.8 <- subset(dr4, dr4$sample_id %in% sample.id.len[8])
dr4.9 <- subset(dr4, dr4$sample_id %in% sample.id.len[9])

dr4.10 <- subset(dr4, dr4$sample_id %in% sample.id.len[10])
dr4.11 <- subset(dr4, dr4$sample_id %in% sample.id.len[11])
dr4.12 <- subset(dr4, dr4$sample_id %in% sample.id.len[12])

dr4.13 <- subset(dr4, dr4$sample_id %in% sample.id.len[13])
dr4.14 <- subset(dr4, dr4$sample_id %in% sample.id.len[14])
dr4.15 <- subset(dr4, dr4$sample_id %in% sample.id.len[15])
dr4.16 <- subset(dr4, dr4$sample_id %in% sample.id.len[16])

## using TSNE function for each sample
func_tsne_plot <- function(data, sample.id, tSNE1, tSNE2){ggplot(data,  aes(x=tSNE1, y=tSNE2, color=cell_clustering1m)) +
                  geom_point(size = 0.6)+theme_bw()+ 
                  theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank())+
                  scale_color_manual(values = cols)+ 
                  ggtitle(sample.id)+
                  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
}

p1 = func_tsne_plot(dr4.1, sample.id.len[1], tSNE1 = dr4.1$tSNE1, tSNE2 = dr4.1$tSNE2)
p2 = func_tsne_plot(dr4.2, sample.id.len[2], tSNE1 = dr4.2$tSNE1, tSNE2 = dr4.2$tSNE2)
p3 = func_tsne_plot(dr4.3, sample.id.len[3], tSNE1 = dr4.3$tSNE1, tSNE2 = dr4.3$tSNE2)

p4 = func_tsne_plot(dr4.4, sample.id.len[4], tSNE1 = dr4.4$tSNE1, tSNE2 = dr4.4$tSNE2)
p5 = func_tsne_plot(dr4.5, sample.id.len[5], tSNE1 = dr4.5$tSNE1, tSNE2 = dr4.5$tSNE2)
p6 = func_tsne_plot(dr4.6, sample.id.len[6], tSNE1 = dr4.6$tSNE1, tSNE2 = dr4.6$tSNE2)

p7 = func_tsne_plot(dr4.7, sample.id.len[7], tSNE1 = dr4.7$tSNE1, tSNE2 = dr4.7$tSNE2)
p8 = func_tsne_plot(dr4.8, sample.id.len[8], tSNE1 = dr4.8$tSNE1, tSNE2 = dr4.8$tSNE2)
p9 = func_tsne_plot(dr4.9, sample.id.len[9], tSNE1 = dr4.9$tSNE1, tSNE2 = dr4.9$tSNE2)

p10 = func_tsne_plot(dr4.10, sample.id.len[10], tSNE1 = dr4.10$tSNE1, tSNE2 = dr4.10$tSNE2)
p11 = func_tsne_plot(dr4.11, sample.id.len[11], tSNE1 = dr4.11$tSNE1, tSNE2 = dr4.11$tSNE2)
p12 = func_tsne_plot(dr4.12, sample.id.len[12], tSNE1 = dr4.12$tSNE1, tSNE2 = dr4.12$tSNE2)

p13 = func_tsne_plot(dr4.13, sample.id.len[13], tSNE1 = dr4.13$tSNE1, tSNE2 = dr4.13$tSNE2)
p14 = func_tsne_plot(dr4.14, sample.id.len[14], tSNE1 = dr4.14$tSNE1, tSNE2 = dr4.14$tSNE2)
p15 = func_tsne_plot(dr4.15, sample.id.len[15], tSNE1 = dr4.15$tSNE1, tSNE2 = dr4.15$tSNE2)
p16 = func_tsne_plot(dr4.16, sample.id.len[16], tSNE1 = dr4.16$tSNE1, tSNE2 = dr4.16$tSNE2)

library(patchwork)
p.total.1 <- (p1|p2|p3|p4)/(p5|p6|p7|p8)/(p9|p10|p11|p12)/(p13|p14|p15|p16)

ggsave("Panel_tSNE_annotated.pdf", p.total.1, width = 24, height = 12)


pdf("Panel_total_tSNE_Total.pdf", width = 8, height = 6)
ggplot(dr2,aes(x = tSNE1, y = tSNE2, color = cell_clustering1m)) +
  geom_point(size = 0.4) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  scale_color_manual(values = cols) +
  guides(color = guide_legend(override.aes = list(size = 4)))
dev.off()

#################################################################
## calculate cell numbers

df.1 = as.data.frame(table(dr4.1$cell_clustering1m))
colnames(df.1)[2] <- as.character(sample.id.len[1])
df.2 = as.data.frame(table(dr4.2$cell_clustering1m))
colnames(df.2)[2] <- as.character(sample.id.len[2])
df.3 = as.data.frame(table(dr4.3$cell_clustering1m))
colnames(df.3)[2] <- as.character(sample.id.len[3])

df.4 = as.data.frame(table(dr4.4$cell_clustering1m))
colnames(df.4)[2] <- as.character(sample.id.len[4])
df.5 = as.data.frame(table(dr4.5$cell_clustering1m))
colnames(df.5)[2] <- as.character(sample.id.len[5])
df.6 = as.data.frame(table(dr4.6$cell_clustering1m))
colnames(df.6)[2] <- as.character(sample.id.len[6])

df.7 = as.data.frame(table(dr4.7$cell_clustering1m))
colnames(df.7)[2] <- as.character(sample.id.len[7])
df.8 = as.data.frame(table(dr4.8$cell_clustering1m))
colnames(df.8)[2] <- as.character(sample.id.len[8])
df.9 = as.data.frame(table(dr4.9$cell_clustering1m))
colnames(df.9)[2] <- as.character(sample.id.len[9])

df.10 = as.data.frame(table(dr4.10$cell_clustering1m))
colnames(df.10)[2] <- as.character(sample.id.len[10])
df.11 = as.data.frame(table(dr4.11$cell_clustering1m))
colnames(df.11)[2] <- as.character(sample.id.len[11])
df.12 = as.data.frame(table(dr4.12$cell_clustering1m))
colnames(df.12)[2] <- as.character(sample.id.len[12])

df.13 = as.data.frame(table(dr4.13$cell_clustering1m))
colnames(df.13)[2] <- as.character(sample.id.len[13])
df.14 = as.data.frame(table(dr4.14$cell_clustering1m))
colnames(df.14)[2] <- as.character(sample.id.len[14])
df.15 = as.data.frame(table(dr4.15$cell_clustering1m))
colnames(df.15)[2] <- as.character(sample.id.len[15])
df.16 = as.data.frame(table(dr4.16$cell_clustering1m))
colnames(df.16)[2] <- as.character(sample.id.len[16])

df.num <- cbind(df.1, df.2, df.3, df.4, df.5, df.6, df.7, df.8, df.9, df.10, df.11, df.12,
                df.13, df.14, df.15, df.16)

write.csv(df.num, "Panel_annotated_cell_num.csv")


