library(dplyr)
library(tidyr)

library(Seurat)

library(ChIPseeker)
library(ChIPpeakAnno)
library(GenomeInfoDb)
library(AnnotationHub)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)

library(GenomicRanges)
library(GenomicFeatures)
library(stringr)
library(reshape2)

library(ggplot2)
library(ggrepel)
library(ggpubr)
library(grid)
library(gridExtra)
library(png)
library(RColorBrewer)

library(biomaRt)
library(scales)


panel_theme <- theme(
    panel.border       = element_blank(),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    legend.title       = element_text(size = 6),
    legend.key.size    = unit(0.2, "cm"),
    legend.text        = element_text(size = 5),
    axis.title.x       = element_text(size = 6),
    axis.text.x        = element_text(size = 6),
    axis.title.y       = element_text(size = 6),
    axis.text.y        = element_text(size = 6),
    plot.title         = element_text(size = 6, face="plain"),
    plot.tag           = element_text(size = 6, face = "plain"),
    legend.spacing.y   = unit(2, "cm")
)


# Single-cell GSE112013 dataset
seurat_object <- readRDS('testis.rds')
Idents(object = seurat_object) <- seurat_object@meta.data$seurat_clusters

# UMAPs with annotation
panel_a <- DimPlot(seurat_object, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, raster=FALSE) + 
panel_theme +
  theme(plot.title = element_blank(),
        legend.position = "none") + 
  guides(fill = guide_legend(ncol = 3)) +
  labs(title = "", tag = "A")


p1 <- DotPlot(seurat_object, 
 features =   c('VIM', 'CD14', 'VWF', 'ACTA2', 'DLK1', 'DAZL', 'MAGEA4', 'UTF1', 'ID4', 'FGFR3', 'KIT', 'DMRT1', 
 	'DMRTB1', 'STRA8', 'SYCP3', 'SPO11', 'MLH3', 'ZPBP', 'TNP1', 'PRM2')) + panel_theme +
  labs(title = "", tag = "B")

seurat_object_promoter <- readRDS('testis_promoter.rds')
seurat_object_genebody <- readRDS('testis_genebody.rds')

# H3K36me3 in promoter, HiMoRNA rho > 0.5
p2 <- FeaturePlot(seurat_object_promoter, paste0("add_module_score", 1), pt.size = 0.5, reduction = "umap") + 
  scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral"))) +
  panel_theme +
  theme(text = element_text(size = 6)) +
  labs(title = "", tag = "C")

# H3K36me3 in promoter, HiMoRNA rho < -0.5
p3 <- FeaturePlot(seurat_object_promoter, paste0("add_module_score", 2), pt.size = 0.5, reduction = "umap") + 
  scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral"))) +
  panel_theme +
  theme(text = element_text(size = 6)) +
  labs(title = "", tag = "")

# H3K36me3 in genebody, HiMoRNA rho > 0.5
p4 <- FeaturePlot(seurat_object_genebody, paste0("add_module_score", 1), pt.size = 0.5, reduction = "umap") + 
  scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral"))) +
  panel_theme +
  theme(text = element_text(size = 6)) +
  labs(title = "", tag = "")

# H3K36me3 in genebody, HiMoRNA rho < -0.5
p5 <- FeaturePlot(seurat_object_genebody, paste0("add_module_score", 2), pt.size = 0.5, reduction = "umap") + 
  scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral"))) +
  panel_theme +
  theme(text = element_text(size = 6)) +
  labs(title = "", tag = "")

#################################

p <- ggarrange(p2, p3, p4, p5, ncol = 4, nrow = 1, legend="none")

figure3_grob <- arrangeGrob(
  panel_a, p1, p,
  ncol = 2, nrow = 5,
  layout_matrix = rbind(c(1, 1), c(1, 1), c(2, 2), c(2, 2), c(3, 3))
)

ggsave("Supplementary_Figure4.png", plot = figure3_grob, width = 16, height = 20, units = "cm")



#p6, p7, p8,