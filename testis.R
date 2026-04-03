library(dplyr)
library(ggplot2)
library(Seurat)
library(gridExtra)
library(reshape2)
library(scales)
library(RColorBrewer)
library(stringr)

# https://www.nature.com/articles/s41420-025-02802-4
counts2 <- read.csv('GSE112013_Combined_UMI_table.txt', sep="\t", row.names = 1)

seurat_object2 <- CreateSeuratObject(counts = counts2, 
                             project = "testis2", 
                             min.cells = 0, 
                             min.features = 0)

seurat_object2 = NormalizeData(object = seurat_object2, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object2 <- FindVariableFeatures(seurat_object2, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_object2)
seurat_object2 <- ScaleData(seurat_object2, features = all.genes)
seurat_object2 <- RunPCA(seurat_object2, features = VariableFeatures(object = seurat_object2))

# UMAP
seurat_object2 <- RunUMAP(seurat_object2, dims = 1:12, reduction.key = "UMAP_")
seurat_object2$UMAP_1 <- as.data.frame(Embeddings(object = seurat_object2[["umap"]]))$UMAP_1
seurat_object2$UMAP_2 <- as.data.frame(Embeddings(object = seurat_object2[["umap"]]))$UMAP_2

#Tsne
seurat_object2 <- RunTSNE(seurat_object2, dims.use = 1:20, reduction.use = "pca")

# Clusterization
seurat_object2 <- FindNeighbors(seurat_object2, dims = 1:12)
seurat_object2 <- FindClusters(seurat_object2, resolution = 1)

Idents(object = seurat_object2) <- seurat_object2@meta.data$seurat_clusters

p <- DotPlot(seurat_object2, 
 features =   c('VIM', 'CD14', 'VWF', 'ACTA2', 'DLK1', 'DAZL', 'MAGEA4', 'UTF1', 'ID4', 'FGFR3', 'KIT', 'DMRT1', 'DMRTB1', 'STRA8', 'SYCP3', 'SPO11', 'MLH3', 'ZPBP', 'TNP1', 'PRM2'))

ggsave("markers_dot_plot.png", plot = p, width = 30, height = 30)

seurat_object2@meta.data$celltype <- NA
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 0] = "Sperm3"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 1] = "Myoid cells" #
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 2] = "Sperm2"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 3] = "Early primary S'cytes"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 4] = "Endothelial cells"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 5] = "Round S'tids"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 6] = "Sperm2"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 7] = "Differentiating S'gonia"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 8] = "Myoid cells"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 9] = "Elongated S'tids"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 10] = "Sperm2"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 11] = "SSCs"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 12] = "Late primary S'cytes"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 13] = "Sperm1"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 14] = "Early primary S'cytes"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 15] = "Macrophage"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 16] = "Late primary S'cytes"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 17] = "Elongated S'tids"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 18] = "Macrophage"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 19] = "Late primary S'cytes"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 20] = "Sertoli cells"

saveRDS(seurat_object2, 'testis.rds')
