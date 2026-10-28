library(dplyr)
library(ggplot2)
library(Seurat)
library(gridExtra)
library(reshape2)
library(scales)
library(RColorBrewer)
library(stringr)

# https://www.nature.com/articles/s41420-025-02802-4
counts <- read.csv('GSE112013_Combined_UMI_table.txt', sep="\t", row.names = 1)

seurat_object <- CreateSeuratObject(counts = counts, 
                             project = "testis", 
                             min.cells = 0, 
                             min.features = 0)

seurat_object = NormalizeData(object = seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all.genes)
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))

# UMAP
seurat_object <- RunUMAP(seurat_object, dims = 1:12, reduction.key = "UMAP_")
seurat_object$UMAP_1 <- as.data.frame(Embeddings(object = seurat_object[["umap"]]))$UMAP_1
seurat_object$UMAP_2 <- as.data.frame(Embeddings(object = seurat_object[["umap"]]))$UMAP_2

# Clusterization
seurat_object <- FindNeighbors(seurat_object, dims = 1:12)
seurat_object <- FindClusters(seurat_object, resolution = 1)

Idents(object = seurat_object) <- seurat_object@meta.data$seurat_clusters

seurat_object@meta.data$celltype <- NA
seurat_object@meta.data$celltype[seurat_object@meta.data$seurat_clusters == 0] = "Sperm3"
seurat_object@meta.data$celltype[seurat_object@meta.data$seurat_clusters == 1] = "Myoid cells" 
seurat_object@meta.data$celltype[seurat_object@meta.data$seurat_clusters == 2] = "Sperm2"
seurat_object@meta.data$celltype[seurat_object@meta.data$seurat_clusters == 3] = "Early primary S'cytes"
seurat_object@meta.data$celltype[seurat_object@meta.data$seurat_clusters == 4] = "Endothelial cells"
seurat_object@meta.data$celltype[seurat_object@meta.data$seurat_clusters == 5] = "Round S'tids"
seurat_object@meta.data$celltype[seurat_object@meta.data$seurat_clusters == 6] = "Sperm2"
seurat_object@meta.data$celltype[seurat_object@meta.data$seurat_clusters == 7] = "Differentiating S'gonia"
seurat_object@meta.data$celltype[seurat_object@meta.data$seurat_clusters == 8] = "Myoid cells"
seurat_object@meta.data$celltype[seurat_object@meta.data$seurat_clusters == 9] = "Elongated S'tids"
seurat_object@meta.data$celltype[seurat_object@meta.data$seurat_clusters == 10] = "Sperm2"
seurat_object@meta.data$celltype[seurat_object@meta.data$seurat_clusters == 11] = "SSCs"
seurat_object@meta.data$celltype[seurat_object@meta.data$seurat_clusters == 12] = "Late primary S'cytes"
seurat_object@meta.data$celltype[seurat_object@meta.data$seurat_clusters == 13] = "Sperm1"
seurat_object@meta.data$celltype[seurat_object@meta.data$seurat_clusters == 14] = "Early primary S'cytes"
seurat_object@meta.data$celltype[seurat_object@meta.data$seurat_clusters == 15] = "Macrophage"
seurat_object@meta.data$celltype[seurat_object@meta.data$seurat_clusters == 16] = "Late primary S'cytes"
seurat_object@meta.data$celltype[seurat_object@meta.data$seurat_clusters == 17] = "Elongated S'tids"
seurat_object@meta.data$celltype[seurat_object@meta.data$seurat_clusters == 18] = "Macrophage"
seurat_object@meta.data$celltype[seurat_object@meta.data$seurat_clusters == 19] = "Late primary S'cytes"
seurat_object@meta.data$celltype[seurat_object@meta.data$seurat_clusters == 20] = "Sertoli cells"

saveRDS(seurat_object, 'testis.rds')
