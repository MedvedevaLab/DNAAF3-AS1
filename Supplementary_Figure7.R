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


panel_a <- DotPlot(seurat_object, 
 features =   c('VIM', 'CD14', 'VWF', 'ACTA2', 'DLK1', 'DAZL', 'MAGEA4', 'UTF1', 'ID4', 'FGFR3', 'KIT', 'DMRT1', 
 	'DMRTB1', 'STRA8', 'SYCP3', 'SPO11', 'MLH3', 'ZPBP', 'TNP1', 'PRM2')) + panel_theme +
  labs(title = "", tag = "A")



# UMAPs with annotation
panel_b <- DimPlot(seurat_object, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, raster=FALSE) + 
panel_theme +
  theme(plot.title = element_blank(),
        legend.position = "none") + 
  guides(fill = guide_legend(ncol = 3)) +
  labs(title = "", tag = "B")


Idents(object = seurat_object) <- seurat_object@meta.data$celltype


celltype_levels <- c("SSCs", "Differentiating S'gonia", "Early primary S'cytes",
                     "Late primary S'cytes", "Round S'tids", "Elongated S'tids",
                     "Sperm1", "Sperm2", "Sperm3", "Macrophage",
                     "Endothelial cells", "Myoid cells", "Sertoli cells")


seurat_object$celltype <- factor(seurat_object$celltype, 
                        levels = celltype_levels)

cols <- rev(c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#D37295", "#499894", "#A0CBE8", "#BAB0AC"))

# UMAPs with annotation
# 
panel_c <- DimPlot(seurat_object, reduction = "umap", group.by = "celltype", label = FALSE, repel = TRUE, raster=FALSE) + 
scale_color_manual(values = rev(cols), guide = guide_legend(reverse = TRUE), drop = FALSE)

# 2. compute label positions from the OTHER identity
emb <- as.data.frame(Embeddings(seurat_object, reduction = 'umap'))
colnames(emb)[1:2] <- c("d1", "d2")
emb$grp <- seurat_object@meta.data[rownames(emb), "seurat_clusters"]

centroids <- emb %>%
  group_by(grp) %>%
  summarise(d1 = median(d1), d2 = median(d2), .groups = "drop")

# 3. drop the label text on top
panel_c <- panel_c + geom_text_repel(
      data = centroids,
      aes(d1, d2, label = grp),
      inherit.aes = FALSE, size = 4, fontface = "bold"
    ) + panel_theme +
  theme(plot.title = element_blank(), legend.position = "right") +
  guides(fill = guide_legend(ncol = 3)) +
  scale_color_manual(values = rev(cols), guide = guide_legend(reverse = TRUE), drop = FALSE) +
  labs(title = "", tag = "C")


extract_module_scores <- function(rds_file, source_label) {
  obj <- readRDS(rds_file)
  md  <- obj@meta.data
  data.frame(
    celltype              = as.character(md$celltype),
    `H3K36me3 rho > 0.5`  = md[["add_module_score1"]],
    `H3K36me3 rho < -0.5` = md[["add_module_score2"]],
    source                = source_label,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

promoter_scores <- extract_module_scores("testis_promoter.rds", "Promoter")
genebody_scores <- extract_module_scores("testis_genebody.rds", "Gene body")

module_scores_long <- rbind(promoter_scores, genebody_scores) %>%
  tidyr::pivot_longer(cols = c("H3K36me3 rho > 0.5", "H3K36me3 rho < -0.5"),
                      names_to = "rho_set", values_to = "module_score") %>%
  dplyr::mutate(
    set      = paste(source, rho_set, sep = " | "),
    celltype = factor(celltype, levels = celltype_levels)
  )

# Four sets shown as four rows: testis_promoter score1, testis_promoter score2,
# testis_genebody score1, testis_genebody score2
set_levels <- c("Promoter | H3K36me3 rho > 0.5",
                "Promoter | H3K36me3 rho < -0.5",
                "Gene body | H3K36me3 rho > 0.5",
                "Gene body | H3K36me3 rho < -0.5")
module_scores_long$set <- factor(module_scores_long$set, levels = set_levels)

panel_d <- ggplot(module_scores_long,
    aes(x = celltype, y = module_score, fill = celltype)) +
  geom_violin(scale = "width", linewidth = 0.2, alpha = 0.9) +
  facet_wrap(set ~ ., scales = "free_y", ncol = 2) +
  scale_fill_manual(values = rev(cols), drop = FALSE, guide = "none") +
  theme_bw() +
  panel_theme +
  theme(axis.text.x      = element_text(angle = 45, hjust = 1, size = 5),
        strip.background = element_blank(),
        strip.text.y     = element_text(size = 6, angle = 0)) +
  labs(title = "H3K36me3 HiMoRNA module scores across testis cell types",
       x = NULL, y = "Module score", tag = "D")



figure7_grob <- arrangeGrob(
  panel_a, panel_b, panel_c, panel_d,
  ncol = 2, nrow = 5,
  layout_matrix = rbind(c(1, 1), c(1, 1),  c(2, 3), c(4, 4), c(4, 4))
)

ggsave("Supplementary_Figure7.png", plot = figure7_grob, width = 16, height = 22, units = "cm")






# seurat_object_promoter <- readRDS('testis_promoter.rds')
# seurat_object_genebody <- readRDS('testis_genebody.rds')

# # H3K36me3 in promoter, HiMoRNA rho > 0.5
# p2 <- FeaturePlot(seurat_object_promoter, paste0("add_module_score", 1), pt.size = 0.5, reduction = "umap") + 
#   scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral"))) +
#   panel_theme +
#   theme(text = element_text(size = 6)) +
#   labs(title = "", tag = "C")

# # H3K36me3 in promoter, HiMoRNA rho < -0.5
# p3 <- FeaturePlot(seurat_object_promoter, paste0("add_module_score", 2), pt.size = 0.5, reduction = "umap") + 
#   scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral"))) +
#   panel_theme +
#   theme(text = element_text(size = 6)) +
#   labs(title = "", tag = "")

# # H3K36me3 in genebody, HiMoRNA rho > 0.5
# p4 <- FeaturePlot(seurat_object_genebody, paste0("add_module_score", 1), pt.size = 0.5, reduction = "umap") + 
#   scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral"))) +
#   panel_theme +
#   theme(text = element_text(size = 6)) +
#   labs(title = "", tag = "")

# # H3K36me3 in genebody, HiMoRNA rho < -0.5
# p5 <- FeaturePlot(seurat_object_genebody, paste0("add_module_score", 2), pt.size = 0.5, reduction = "umap") + 
#   scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral"))) +
#   panel_theme +
#   theme(text = element_text(size = 6)) +
#   labs(title = "", tag = "")

# #################################

# p <- ggarrange(p2, p3, p4, p5, ncol = 4, nrow = 1, legend="none")
