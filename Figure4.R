library(ggplot2)
library(ggpubr)
library(ChIPseeker)
library(png)
library(gridExtra)
library(dplyr)
library(ChIPpeakAnno)
library(grid)
library("GenometriCorr")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(stringr)
library(ggrepel)
library(Seurat)

panel_theme <- theme(
    panel.border       = element_blank(),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    axis.title.x       = element_text(size = 6),
    axis.text.x        = element_text(size = 6, hjust = 0.5),
    legend.title       = element_text(size = 6),
    legend.key.size    = unit(0.2, "cm"),
    legend.text        = element_text(size = 5),
    axis.ticks.x       = element_blank(),
    axis.text.y        = element_text(size = 6),
    axis.title.y       = element_blank(),
    plot.title         = element_text(size = 6),
    plot.tag           = element_text(size = 6, face = "plain"),
    legend.spacing.y   = unit(2, "cm")
)



# Annotate genomic peaks from a data frame of coordinates.
# `cols` selects the chrom/start/end columns (by name or position).
get_annotation <- function(peaks, txdb, cols = c(1, 2, 3)) {
  coords <- peaks[, cols]
  coords <- cbind(coords, "*")
  colnames(coords) <- c("seqnames", "start", "end", "strand")
  coords$start <- as.numeric(coords$start)
  coords$end   <- as.numeric(coords$end)
  ranges <- toGRanges(coords, format = "BED")
  peakAnno <- annotatePeak(ranges, tssRegion = c(-1000, 100), TxDb = txdb, overlap = "all")
  return(peakAnno)
}

gtf_name <- "gencode.v49.annotation.gtf"
txdb <- makeTxDbFromGFF(gtf_name)

# Panel A
img <- readPNG('metascape.png')
g11 <- rasterGrob(img,
  interpolate=TRUE, width = unit(6,"in"), height=unit(1.2,"in"))

panel_a<- ggplot(geom="blank") +
  annotation_custom(g11) +
  theme_bw() +
  panel_theme +
  theme(
      axis.title.x = element_blank(),
      legend.key.size = unit(0.2, "cm"),
      legend.text = element_text(size = 5)) + labs(tag = "A")



diffbind_results <- read.csv('/data/budkina/himorna_chipseq_data/test/CTD-2587H24.5_H3K36me3_differential_peaks.csv', sep = '\t')
diffbind_sig <- diffbind_results[diffbind_results$FDR < 0.05,] # & abs(diffbind_results$Fold) > 1

diffbind_sig$significance <- "Unchanged"
diffbind_sig$significance[diffbind_sig$Fold < 0] <- "Lost"
diffbind_sig$significance[diffbind_sig$Fold > 0] <- "Gained"
diffbind_sig$significance <- factor(diffbind_sig$significance, levels = c( "Gained", "Lost", "Unchanged"))


diffbind_sig_gained    <- diffbind_sig[diffbind_sig$significance == "Gained", ]

diffbind_sig_gained_bed <- diffbind_sig_gained[, c("seqnames", "start", "end", "strand")]

write.table(
    file      = "diffbind_sig_gained.bed",
    diffbind_sig_gained_bed,
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep       = "\t"
)

diffbind_sig_lost    <- diffbind_sig[diffbind_sig$significance == "Lost", ]

diffbind_sig_lost_bed <- diffbind_sig_lost[, c("seqnames", "start", "end", "strand")]

write.table(
    file      = "diffbind_sig_lost.bed",
    diffbind_sig_lost_bed,
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep       = "\t"
)
             

gained_peak_anno       <- get_annotation(diffbind_sig_gained, txdb)
gained_peak_anno_df    <- as.data.frame(gained_peak_anno@anno)
gained_peak_anno_df$geneId <- str_split_fixed(gained_peak_anno_df$geneId, "\\.", 2)[, 1]
gained_peak_anno_df <- gained_peak_anno_df[gained_peak_anno_df$annotation != 'Promoter',]
  

all_deseq2_results <- read.csv("DESeq2_genes_ASO_all.tsv", sep = "\t")
deseq2_aso1     <- dplyr::filter(all_deseq2_results, perturb_id == "ASO_G0267577_01")
deseq2_aso1_sig <- deseq2_aso1[deseq2_aso1$fdr < 0.05, ]
deseq2_aso3     <- dplyr::filter(all_deseq2_results, perturb_id == "ASO_G0267577_03")
deseq2_aso3_sig <- deseq2_aso3[deseq2_aso3$fdr < 0.05, ]




de_to_plot <-  dplyr::filter(deseq2_aso1, geneID %in% gained_peak_anno_df$geneId)


de_to_plot$log_padj <- -log10(de_to_plot$fdr)
de_to_plot$sig <- de_to_plot$fdr < 0.05

reproduction_genes <- c('BCL2L2', 'PSMA1', 'RFX2', 'RPL30', 'VDAC3', 'SELENOF')
de_to_plot$reproduction <- de_to_plot$geneSymbol %in% reproduction_genes

panel_b <- ggplot(de_to_plot, aes(x = log2FC, y = log_padj, label = ifelse(sig, geneSymbol, NA))) +
    geom_point(shape = 1, alpha = 0.4) +
    geom_text_repel(color = "black", size = 2, max.overlaps = 12) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60") +
    theme_bw() +
    panel_theme +
    labs(title = "ASO_01 differential expression gained peaks",
         x = "log2FC", y = "-log10(FDR)", tag = "B")



# G. Dotplot of reproduction genes across ASO_01 and ASO_03
dotplot_df <- dplyr::filter(all_deseq2_results,
    perturb_id %in% c("ASO_G0267577_01", "ASO_G0267577_03"),
    geneSymbol %in% reproduction_genes)
dotplot_df$ASO <- ifelse(dotplot_df$perturb_id == "ASO_G0267577_01", "ASO_01", "ASO_03")
dotplot_df <- dotplot_df %>%
    dplyr::group_by(ASO, geneSymbol) %>%
    dplyr::slice_min(fdr, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
dotplot_df$log_padj <- -log10(dotplot_df$fdr)
dotplot_df$geneSymbol <- factor(dotplot_df$geneSymbol, levels = rev(reproduction_genes))
dotplot_df$sig <- dotplot_df$fdr < 0.05

panel_c <- ggplot(dotplot_df, aes(x = ASO, y = geneSymbol, fill = log2FC, size = log_padj, color = sig)) +
    geom_point(shape = 21, stroke = 0.6) +
    scale_fill_gradient2(low = "#1F78B4", mid = "white", high = "#E31A1C", midpoint = 0) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = NA), guide = "none") +
    scale_size_continuous(range = c(1, 5)) +
    theme_bw() +
    panel_theme +
    labs(title = "Reproduction genes", x = NULL, size = "-log10(FDR)",
         fill = "log2FC", tag = "C")




# Single-cell GSE112013 dataset
seurat_object <- readRDS('testis.rds')
Idents(object = seurat_object) <- seurat_object@meta.data$celltype
seurat_object$celltype <- factor(seurat_object$celltype, 
                        levels = c("SSCs", "Differentiating S'gonia", "Early primary S'cytes",
                          "Late primary S'cytes", "Round S'tids", "Elongated S'tids", "Sperm1",
                          "Sperm2", "Sperm3", "Macrophage", "Endothelial cells", "Myoid cells",
                          "Sertoli cells"))

cols <- rev(c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#D37295", "#499894", "#A0CBE8", "#BAB0AC"))

# UMAPs with annotation
# 
panel_d <- DimPlot(seurat_object, reduction = "umap", group.by = "celltype", label = FALSE, repel = TRUE, raster=FALSE) + 
panel_theme +
  theme(plot.title = element_blank()) + 
  guides(fill = guide_legend(ncol = 3)) +
  scale_color_manual(values = rev(cols), guide = guide_legend(reverse = TRUE), drop = FALSE) +
  labs(title = "", tag = "D")

# ── panel_d: DNAAF3-AS1 and DNAAF3 expression across cell types ───────────────
# Combine the two per-gene violins into one boxplot panel (both genes side by
# side per cell type) and overlay a smoothed trend across cell types (ordered
# along spermatogenesis) for each gene.
d_genes <- c('CTD-2587H24.5', 'DNAAF3', 'RFX2')
d_expr  <- FetchData(seurat_object, vars = c(d_genes, "celltype"))
d_long  <- tidyr::pivot_longer(d_expr, cols = tidyselect::all_of(d_genes),
                               names_to = "gene", values_to = "expression")
d_long$gene         <- factor(d_long$gene, levels = d_genes,
                              labels = c("DNAAF3-AS1", "DNAAF3", 'RFX2'))
# Drop somatic / non-germline cell types from this panel
d_exclude_celltypes <- c("Macrophage", "Endothelial cells", "Myoid cells", "Sertoli cells")
d_long              <- d_long[!(d_long$celltype %in% d_exclude_celltypes), ]
d_long$celltype     <- factor(d_long$celltype,
                              levels = setdiff(levels(seurat_object$celltype),
                                               d_exclude_celltypes))
d_long$celltype_num <- as.numeric(d_long$celltype)

gene_cols <- c("DNAAF3-AS1" = "#4E79A7",
  "DNAAF3" = "#F28E2B",
  'RFX2' = "#E15759")

panel_e <- ggplot(d_long, aes(x = celltype, y = expression, fill = gene)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.2, alpha = 0.7,
               position = position_dodge(width = 0.8)) +
  geom_smooth(aes(x = celltype_num, color = gene, group = gene),
              method = "loess", se = FALSE, linewidth = 0.6) +
  scale_fill_manual(values = gene_cols, name = "") +
  scale_color_manual(values = gene_cols, name = "") +
  theme_bw() +
  panel_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) +
  labs(title = "", x = NULL,
       y = "Expression (log-normalized)", tag = "E")


g <- arrangeGrob(panel_a, panel_b, panel_c, panel_d, panel_e, ncol = 3, nrow = 3, layout_matrix= rbind(c(1,1,1), c(2,2,3), c(4,5,5)))
ggsave(paste0("Figure4.png"), plot = g, width = 18, height = 24, units = "cm")
