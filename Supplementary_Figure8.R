library(dplyr)
library(tidyr)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(gridExtra)


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

# A. Marker gene expression per cluster
panel_a <- DotPlot(seurat_object, 
 features =   c('VIM', 'CD14', 'VWF', 'ACTA2', 'DLK1', 'DAZL', 'MAGEA4', 'UTF1', 'ID4', 'FGFR3', 'KIT', 'DMRT1', 
 	'DMRTB1', 'STRA8', 'SYCP3', 'SPO11', 'MLH3', 'ZPBP', 'TNP1', 'PRM2')) + panel_theme +
  labs(title = "", tag = "A")

# B. UMAPs with cluster annotation
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

celltype_scale <- function(aesthetic = "colour")
  scale_color_viridis_d(option = "turbo", aesthetics = aesthetic,
                        guide = guide_legend(reverse = TRUE), drop = FALSE)

# C. UMAPs with cell type annotation
panel_c <- DimPlot(seurat_object, reduction = "umap", group.by = "celltype", label = FALSE, repel = TRUE, raster=FALSE) +
celltype_scale("colour")

emb <- as.data.frame(Embeddings(seurat_object, reduction = 'umap'))
colnames(emb)[1:2] <- c("d1", "d2")
emb$grp <- seurat_object@meta.data[rownames(emb), "seurat_clusters"]

centroids <- emb %>%
  group_by(grp) %>%
  summarise(d1 = median(d1), d2 = median(d2), .groups = "drop")

panel_c <- panel_c + geom_text_repel(
      data = centroids,
      aes(d1, d2, label = grp),
      inherit.aes = FALSE, size = 4, fontface = "bold"
    ) + panel_theme +
  theme(plot.title = element_blank(), legend.position = "right") +
  guides(fill = guide_legend(ncol = 3)) +
  celltype_scale("colour") +
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

# C. Himorna target genes modules expression
promoter_scores <- extract_module_scores("testis_promoter.rds", "Promoter")
genebody_scores <- extract_module_scores("testis_genebody.rds", "Gene body")

module_scores_long <- rbind(promoter_scores, genebody_scores) %>%
  tidyr::pivot_longer(cols = c("H3K36me3 rho > 0.5", "H3K36me3 rho < -0.5"),
                      names_to = "rho_set", values_to = "module_score") %>%
  dplyr::mutate(
    set      = paste(source, rho_set, sep = " | "),
    celltype = factor(celltype, levels = celltype_levels)
  )

set_levels <- c("Promoter | H3K36me3 rho > 0.5",
                "Promoter | H3K36me3 rho < -0.5",
                "Gene body | H3K36me3 rho > 0.5",
                "Gene body | H3K36me3 rho < -0.5")
module_scores_long$set <- factor(module_scores_long$set, levels = set_levels)

panel_d <- ggplot(module_scores_long,
    aes(x = celltype, y = module_score, fill = celltype)) +
  geom_violin(scale = "width", linewidth = 0.2, alpha = 0.9) +
  facet_wrap(set ~ ., scales = "free_y", ncol = 2) +
  scale_fill_viridis_d(option = "turbo", drop = FALSE, guide = "none") +
  theme_bw() +
  panel_theme +
  theme(axis.text.x      = element_text(angle = 45, hjust = 1, size = 5),
        strip.background = element_blank(),
        strip.text.y     = element_text(size = 6, angle = 0)) +
  labs(title = "H3K36me3 HiMoRNA module scores across testis cell types",
       x = NULL, y = "Module score", tag = "D")

g <- arrangeGrob(
  panel_a, panel_b, panel_c, panel_d,
  ncol = 2, nrow = 6,
  layout_matrix = rbind(c(1, 1), c(1, 1),  c(2, 3), c(4, 4), c(4, 4))
)

ggsave("Supplementary_Figure8.png", plot = g, width = 16, height = 24, units = "cm")