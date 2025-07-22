## InferCNV analysis for Epithelial cells
rm(list=ls())
setwd("/home/data/dengmin/02.Xu_lab/tingting")

library(Seurat)
library(dplyr)
library(patchwork)
library(devtools)
library(ggplot2)
library(infercnv)

big.data <- readRDS("WT_annotate_scRNA.data_20250323.rds")
unique(big.data$subtype)


Idents(big.data) <- "clustering"
epi.data.1 <- subset(big.data, idents = c("Epithelial cell"))
dim(epi.data.1)

big.data.2 <- readRDS("Epithelial.data.rds")
table(big.data.2$subtype)
Idents(big.data.2) <- "subtype"
epi.data.2 <- subset(big.data.2, idents = c("Normal"))
table(epi.data.2$subtype)

merge.data <- merge(epi.data.1, epi.data.2)
table(merge.data$subtype)
dim(merge.data)

## Add meta datacondition to merge .data with triple negative, HER2+, ER+ to normal, and normal to nromal
merge.data$condition <- "Normal"
merge.data$condition[merge.data$subtype == "Triple negative"] <- "tumor"
merge.data$condition[merge.data$subtype == "Triple negative "] <- "tumor"
merge.data$condition[merge.data$subtype == "HER2+"] <- "tumor"
merge.data$condition[merge.data$subtype == "ER+"] <- "tumor"
table(merge.data$condition)

Idents(merge.data) <- "condition"
#merge.data <- subset(merge.data, downsample = 1000)
table(merge.data$condition)
merge.data <- JoinLayers(merge.data)
tumor_vs_normal <- FindMarkers(merge.data, ident.1 = "tumor", ident.2 = "Normal", min.pct = 0.0001, logfc.threshold = 0.05)
tumor_vs_normal <- tumor_vs_normal[order(tumor_vs_normal$p_val),]
head(tumor_vs_normal)
## transitional_tumor vs normal


## Output DE genes list from above using write.csv
write.csv(tumor_vs_normal, "WT_Epi_tumor_vs_normal_DEGs.csv", quote = FALSE, row.names = TRUE)





## Detect DE genes between tumor and normal
merge.data <- JoinLayers(merge.data)
matrix <- as.matrix(GetAssayData(merge.data[["RNA"]], slot = "counts"))
head(matrix)
dim(matrix)


## 2.Cell annotation file
cellAnnotation = subset(merge.data@meta.data, select = condition)
cellAnnotation$cell.id = rownames(cellAnnotation)

cellAnnotation = cellAnnotation[,match(colnames(cellAnnotation),c("cell.id","condition"))]
head(cellAnnotation)
unique(cellAnnotation$condition)
colnames(cellAnnotation) <- NULL
rownames(cellAnnotation) <- NULL
write.table(cellAnnotation,"v7_WT_cellAnnotation.Epithelial.txt",sep = "\t",quote = F,row.names = F)


## construct inferCNV object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=matrix,
                                    annotations_file="v7_WT_cellAnnotation.Epithelial.txt",
                                    delim="\t",
                                    gene_order_file="gene_pos.txt",
                                    ref_group_names=NULL)

## Run inferCNV
options(mc.cores = parallel::detectCores())
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff = 0.1,  
                             out_dir = "v7_WT_out_inferCNV.Epithelial",  # 如果输出文件夹不存在，会自己创建
                             cluster_by_groups = F,   # cluster
                             denoise = T,
                             #cluster_method = "ward.D2",
                             #k_obs_groups = 5,
                             HMM = F, plot_steps = F, num_threads = 90)


############################################
## infer CNV
infercnv_obj = readRDS("./v7_WT_out_inferCNV.Epithelial/run.final.infercnv_obj")

expr <- infercnv_obj@expr.data
dim(expr)
expr[1:4,1:4]

gn <- rownames(expr)
geneFile <- read.table('gene_pos.txt')
length(geneFile$V1)
length(gn) 

sub_geneFile <-  geneFile[match(gn,geneFile$V1),]

identical(rownames(expr),sub_geneFile$V1)
tumor_expr <- as.data.frame(expr)

tumor_expr$chr <- as.factor(sub_geneFile$V2)
table(tumor_expr$chr)
length(tumor_expr$chr)

chr1.len <- length(which(tumor_expr$chr=="chr1"))
chr2.len <- length(which(tumor_expr$chr=="chr2"))
chr3.len <- length(which(tumor_expr$chr=="chr3"))
chr4.len <- length(which(tumor_expr$chr=="chr4"))
chr5.len <- length(which(tumor_expr$chr=="chr5"))
chr6.len <- length(which(tumor_expr$chr=="chr6"))
chr7.len <- length(which(tumor_expr$chr=="chr7"))
chr8.len <- length(which(tumor_expr$chr=="chr8"))
chr9.len <- length(which(tumor_expr$chr=="chr9"))
chr10.len <- length(which(tumor_expr$chr=="chr10"))
chr11.len <- length(which(tumor_expr$chr=="chr11"))
chr12.len <- length(which(tumor_expr$chr=="chr12"))
chr13.len <- length(which(tumor_expr$chr=="chr13"))
chr14.len <- length(which(tumor_expr$chr=="chr14"))
chr15.len <- length(which(tumor_expr$chr=="chr15"))
chr16.len <- length(which(tumor_expr$chr=="chr16"))
chr17.len <- length(which(tumor_expr$chr=="chr17"))
chr18.len <- length(which(tumor_expr$chr=="chr18"))
chr19.len <- length(which(tumor_expr$chr=="chr19"))
chr20.len <- length(which(tumor_expr$chr=="chr20"))
chr21.len <- length(which(tumor_expr$chr=="chr21"))
chr22.len <- length(which(tumor_expr$chr=="chr22"))

chr.order <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
               "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
               "chr21","chr22")

n <- t(tumor_expr[,-ncol(tumor_expr)])

## For column annotation
library(RColorBrewer)

annotation_col = data.frame(
  chr = factor(rep(chr.order,c(chr1.len,chr2.len,chr3.len,chr4.len,chr5.len,
                               chr6.len,chr7.len,chr8.len,chr9.len,chr10.len,
                               chr11.len,chr12.len,chr13.len,chr14.len,chr15.len,
                               chr16.len,chr17.len,chr18.len,chr19.len,chr20.len,
                               chr21.len,chr22.len)))
)
rownames(annotation_col) = colnames(n)
head(annotation_col)

get_group_color_palette <- function () {
  return(colorRampPalette(RColorBrewer::brewer.pal(12, "Set3")))
}
col.color <- get_group_color_palette()(22)

ann_colors = list(
  chr = c(chr1=col.color[1],chr2 =col.color[2],chr3=col.color[3],chr4=col.color[4],
          chr5=col.color[5],chr6=col.color[6],chr7=col.color[7],chr8=col.color[8],
          chr9=col.color[9],chr10=col.color[10],chr11=col.color[11],chr12=col.color[12],
          chr13=col.color[13],chr14=col.color[14],chr15=col.color[15],chr16=col.color[16],
          chr17=col.color[17],chr18=col.color[18],chr19=col.color[19],chr20=col.color[20],
          chr21=col.color[21],chr22=col.color[22])
)
head(ann_colors)

## For row annotation
## load grouping data
group.df <- read.table("./WT_out_inferCNV.Epithelial/infercnv.observation_groupings.txt")
dim(group.df)
colnames(group.df)
head(group.df)

table(group.df$Dendrogram.Color)

#group.df.v2 <- read.table("./v2_out_inferCNV.Epithelial/infercnv_subclusters.observation_groupings.txt")
transitional_tumor.df <- subset(group.df, 
group.df$Dendrogram.Color %in% c("#FFED6F","#EDBC63","#F0EC87","#F0F9B5",
"#CDD796","#F4F4B9","#E1EBA0"))

transitional_tumor.id <- rownames(transitional_tumor.df)
length(transitional_tumor.id)

tumor.df <- group.df[!rownames(group.df) %in% transitional_tumor.id, ]
tumor.id <- rownames(tumor.df)

write.table(transitional_tumor.id,"id.transitional_tumor.Epithelial.txt",sep = "\t",quote = F,row.names = F)
write.table(tumor.id,"id.tumor.Epithelial.txt",sep = "\t",quote = F,row.names = F)

## how to add tumor.id and transitional_tumor.id to merge.data.2 obj
transitional_tumor.id <- read.table("id.transitional_tumor.Epithelial.txt",sep = "\t",header = T)

head(transitional_tumor.id)
tt.id <- transitional_tumor.id$x

tumor.id <- read.table("id.tumor.Epithelial.txt",sep = "\t",header = T)
t.id <- tumor.id$x


tt.merge.data <- subset(merge.data.2, cells = tt.id)
tt.merge.data$condition <- "transitional_tumor"
dim(tt.merge.data)

t.merge.data <- subset(merge.data.2, cells = t.id)
t.merge.data$condition <- "tumor"
dim(t.merge.data)


merge.data.3 <- merge(tt.merge.data,t.merge.data)
table(merge.data.3$condition)


## Select genes expressed in more than 60% of cells in both the `transitional_tumor` and `tumor` conditions.
## Output the gene lists along with their expression percentages and mean values.
## extract the subtype == "Normal" in merge.data.2
normal.merge.data <- subset(merge.data.2, idents = "Normal")
dim(normal.merge.data)
table(normal.merge.data$subtype)

## use the normal.merge.data do the same thing as above
normal_expr <- normal.merge.data@assays$RNA@counts
dim(normal_expr)
normal_expr[1:4,1:4]  

# Calculate the percentage of cells expressing each gene in each condition
transitional_tumor_expr <- rowSums(tt.merge.data@assays$RNA@counts > 0) / ncol(tt.merge.data@assays$RNA@counts) * 100
tumor_expr <- rowSums(t.merge.data@assays$RNA@counts > 0) / ncol(t.merge.data@assays$RNA@counts) * 100
normal_expr <- rowSums(normal.merge.data@assays$RNA@counts > 0) / ncol(normal.merge.data@assays$RNA@counts) * 100

# Filter genes expressed in more than 60% of cells in both conditions
transitional_tumor_genes <- names(transitional_tumor_expr[transitional_tumor_expr > 10])
tumor_genes <- names(tumor_expr[tumor_expr > 10])
normal_genes <- names(normal_expr[normal_expr > 10])

# Calculate mean expression values for the common genes
transitional_tumor_means <- rowMeans(tt.merge.data@assays$RNA@counts[transitional_tumor_genes, ])
tumor_means <- rowMeans(t.merge.data@assays$RNA@counts[tumor_genes, ])
normal_means <- rowMeans(normal.merge.data@assays$RNA@counts[normal_genes, ])

# Create data frames for output
transitional_tumor_output <- data.frame(
  Gene = transitional_tumor_genes,
  Expression_Percentage = transitional_tumor_expr[transitional_tumor_genes],
  Mean_Expression = transitional_tumor_means
)

tumor_output <- data.frame(
  Gene = tumor_genes,
  Expression_Percentage = tumor_expr[tumor_genes],
  Mean_Expression = tumor_means
)

norma_output <- data.frame(
  Gene = normal_genes,
  Expression_Percentage = normal_expr[normal_genes],
  Mean_Expression = normal_means
)
# Write the outputs to files
write.csv(transitional_tumor_output, "transitional_tumor_genes.csv",  quote = FALSE, row.names = FALSE)
write.csv(tumor_output, "tumor_genes.csv", quote = FALSE, row.names = FALSE)
write.csv(norma_output, "normal_genes.csv", quote = FALSE, row.names = FALSE)

## overlap tumor_genes and transitional_tumor_genes
common_genes_1 <- intersect(tumor_genes, transitional_tumor_genes)
length(common_genes_1)

common_genes_2 <- intersect(normal_genes, transitional_tumor_genes)
length(common_genes_2)

## detecte DEGs from tumor vs normal, and transitional_tumor vs normal
## tumor vs normal
t.normal <- merge(t.merge.data,normal.merge.data)
table(t.normal$condition)
Idents(t.normal) <- "condition"

tumor_vs_normal <- FindMarkers(t.normal, ident.1 = "tumor", ident.2 = "Normal", min.pct = 0.0001, logfc.threshold = 0.05)
tumor_vs_normal <- tumor_vs_normal[order(tumor_vs_normal$p_val),]
head(tumor_vs_normal)
## transitional_tumor vs normal
tt.normal <- merge(tt.merge.data,normal.merge.data)
table(tt.normal$condition)
Idents(tt.normal) <- "condition"

transitional_tumor_vs_normal <- FindMarkers(tt.normal, ident.1 = "transitional_tumor", ident.2 = "Normal", min.pct = 0.0001, 
                                            logfc.threshold = 0.05)
transitional_tumor_vs_normal <- transitional_tumor_vs_normal[order(transitional_tumor_vs_normal$p_val),]

## Output DE genes list from above using write.csv
write.csv(tumor_vs_normal, "tumor_vs_normal_DEGs.csv", quote = FALSE, row.names = TRUE)
write.csv(transitional_tumor_vs_normal, "transitional_tumor_vs_normal_DEGs.csv", quote = FALSE, row.names = TRUE)

