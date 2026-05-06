library(ggplot2)
library(ggbreak)
library(ggpubr)
library(ChIPseeker)
library(png)
library(gridExtra)
library(reshape2)
library(AnnotationHub)
library(biomaRt)
library(scales)
library(RColorBrewer)
library(stringr)
library(GenomeInfoDb)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPpeakAnno)
library(png)
library(grid)
library(ggbreak) 

library("org.Hs.eg.db")
library(clusterProfiler)
require("ensembldb")

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

get_annotation <- function(bed, txdb) 
{
  coords <- bed[,c(1,2,3)]
  coords <- cbind(coords, "*")
  colnames(coords) <- c('seqnames', 'start', 'end', 'strand')
  coords$start <- as.numeric(coords$start)
  coords$end <- as.numeric(coords$end)
  ranges <- toGRanges(coords, format='BED')
  peakAnno <- annotatePeak(ranges, tssRegion=c(-1000, 100), TxDb=txdb, overlap="all")
  return(peakAnno)
}

# A. Metagene 


img <- readPNG('Figure2A/Figure2A1.png')
g11 <- rasterGrob(img,
  interpolate=TRUE, width = unit(3.5,"in"), height=unit(2.2,"in"))

p11 <- ggplot(geom="blank") +
  annotation_custom(g11) + 
  theme_bw() + 
  panel_theme + 
  theme(
      axis.title.x = element_blank(),
      legend.key.size = unit(0.2, "cm"),
      legend.text = element_text(size = 5)) + labs(tag = "A")

img <- readPNG('Figure2A/Figure2A2.png')
g12 <- rasterGrob(img,
  interpolate=TRUE, width = unit(3.5,"in"), height=unit(2.2,"in"))

p12 <- ggplot(geom="blank") +
  annotation_custom(g12) + 
  theme_bw() +
  panel_theme +
  theme(
      axis.title.x = element_blank(),
      legend.key.size = unit(0.2, "cm"),
      legend.text = element_text(size = 5)) + labs(tag = "")

# B. Read ChIP-seq peaks
aso_01 <- read.csv('/data/budkina/himorna_chipseq_data/old_processing/all_peaks/macs3_peaks/C-9_/C-9__peaks.broadPeak', sep = '\t', header = F)
aso_03 <- read.csv('/data/budkina/himorna_chipseq_data/old_processing/all_peaks/macs3_peaks/C-10/C-10_peaks.broadPeak', sep = '\t', header = F)

nca_01 <- read.csv('/data/budkina/himorna_chipseq_data/old_processing/all_peaks/macs3_peaks/C-19/C-19_peaks.broadPeak', sep = '\t', header = F)
nca_02 <- read.csv('/data/budkina/himorna_chipseq_data/old_processing/all_peaks/macs3_peaks/C-20/C-20_peaks.broadPeak', sep = '\t', header = F)

peak_num <- data.frame(c('ASO_01', 'ASO_03', 'Control1', 'Control2'), 
  c(dim(aso_01)[1], dim(aso_03)[1], dim(nca_01)[1], dim(nca_02)[1]),
  c('KD','KD','Control','Control'))

colnames(peak_num) <- c('experiment', 'n_peaks', 'type')

p2 <- ggplot(peak_num, aes(x = experiment, y = n_peaks)) +
    theme_bw() +
    geom_bar(position = "dodge", stat = "identity") +
    scale_fill_brewer(palette = "Set1", guide = guide_legend(reverse = TRUE), drop = FALSE) +
    panel_theme +
    guides(fill = guide_legend(ncol = 3)) +
    labs(title = "", x = "Experiment", y = "Number of peaks", tag = "B")

diffbind_results <- read.csv('/data/budkina/himorna_chipseq_data/old_processing/DiffBind_results/macs_broad/CTD-2587H24.5_H3K36me3_differential_peaks.csv', sep = '\t')
diffbind_sig <- diffbind_results[diffbind_results$FDR < 0.05,] # & abs(diffbind_results$Fold) > 1

diffbind_sig$significance <- "Unchanged"
diffbind_sig$significance[diffbind_sig$Fold < 0] <- "Lost"
diffbind_sig$significance[diffbind_sig$Fold > 0] <- "Gained"
diffbind_sig$significance <- factor(diffbind_sig$significance, levels = c( "Gained", "Lost", "Unchanged"))

colors <- c( "Gained" = "#E31A1C", "Lost" =  "#1F78B4", "Unchanged" = "#CCCCCC")

p3 <- ggplot(diffbind_sig, aes(x = Conc, y = Fold, color = significance)) +
  theme_bw() +
  geom_point(alpha = 0.6, size = 0.8) +
  scale_color_manual(values = colors, name = "Peak Status") +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", alpha = 0.5) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +  
  theme(panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(size = 6),
      axis.text.x = element_text(size = 6, vjust = 0.5, hjust = 1),
      legend.title = element_text(size = 6),
      legend="bottom",
      axis.title.y = element_text(size = 6),
      axis.text.y = element_text(size = 6),
      legend.text = element_text(size = 6),
      plot.title = element_text(size = 6),
      legend.key.size = unit(0.2, "cm"),
      plot.tag = element_text(size = 6, face="plain")) + 
  labs(title = "",
        x = "Log concentration",
        y = "Log2FC (KD/Control)",
        tag = "C")


gtf_name <- "gencode.v49.annotation.gtf"
txdb <- makeTxDbFromGFF(gtf_name)

anno_aso_01 <- get_annotation(aso_01, txdb)
anno_aso_03 <- get_annotation(aso_03, txdb)
anno_nca_01 <- get_annotation(nca_01, txdb)
anno_nca_02 <- get_annotation(nca_02, txdb)

annotation_aso_01 <- anno_aso_01@annoStat
annotation_aso_01$Type <- 'ASO_01'
annotation_aso_03 <- anno_aso_03@annoStat
annotation_aso_03$Type <- 'ASO_03'
annotation_nca_01 <- anno_nca_01@annoStat
annotation_nca_01$Type <- 'Control1'
annotation_nca_02 <- anno_nca_02@annoStat
annotation_nca_02$Type <- 'Control2'
annotation <- rbind(annotation_aso_01, annotation_aso_03, annotation_nca_01, annotation_nca_02)
annotation$Feature <- factor(annotation$Feature, levels = c("Promoter",
"5' UTR",
"1st Exon",
"1st Intron",
"Other Exon",
"Other Intron",
"3' UTR",
"Downstream (<=300)",
"Distal Intergenic"
))

annotation$Type <- factor(annotation$Type, levels = c('ASO_01', 'ASO_03', 'Control1', 'Control2'))

cols <- rev(c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F"))
p4 <- ggplot(annotation, aes(x = Type, y = Frequency, fill =  Feature)) + 
  theme_bw() + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=rev(cols), guide = guide_legend(reverse = TRUE), drop = FALSE) +
  theme(panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 6, vjust = 0.5, hjust = 1),
      legend.title = element_text(size = 6),
      legend="bottom",
      axis.text.y = element_text(size = 6),
      legend.text = element_text(size = 6),
      plot.title = element_text(size = 6),
      legend.key.size = unit(0.2, "cm"),
      plot.tag = element_text(size = 6, face="plain")) + 
      labs(title = "H3K36me3 peaks annotation", tag = "D")

# B. peak annotation
# B. MA plot for DiffBind

p1 <- ggarrange(p11, p12, ncol = 2, nrow = 1, common.legend = TRUE, legend="bottom")

##################################################################################
# C. ChIP-seq peak annotation

#########################################################
# diff_peak_anno_df_promoter <- dplyr::filter(diff_peak_anno_df_up, annotation == 'Promoter')
# diff_peak_anno_df_promoter$geneId <-  str_split_fixed(diff_peak_anno_df_promoter$geneId, "\\.", 2)[,1]
# diff_peak_anno_df_promoter_selected <- diff_peak_anno_df_promoter[diff_peak_anno_df_promoter$distanceToTSS !=0,]

# diffbind_sig_up <- diffbind_sig[diffbind_sig$significance == "Gained",]
# sig_for_annotation <- diffbind_sig_up[,c('seqnames', 'start', 'end')]
# diff_peak_anno_up <- get_annotation(sig_for_annotation, txdb)
# diff_peak_anno_df_up <- as.data.frame(diff_peak_anno_up@anno)

# diff_peak_anno_df_not_promoter <- dplyr::filter(diff_peak_anno_df_up, annotation != 'Promoter')
# diff_peak_anno_df_not_promoter$geneId <-  str_split_fixed(diff_peak_anno_df_not_promoter$geneId, "\\.", 2)[,1]

# ego_BP_up <- enrichGO(gene      = diff_peak_anno_df_not_promoter$geneId,
#                     OrgDb         = org.Hs.eg.db,                    # universe      = universe,
#                     keyType       = 'ENSEMBL',
#                     ont           = "BP",
#                     pAdjustMethod = "BH",
#                     qvalueCutoff  = 0.05) 

# p4 <-  dotplot(ego_BP_up, showCategory = 10 ) +
#   theme(panel.border = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.title.x = element_blank(), 
#         text = element_text(size = 6), 
#         legend.text = element_text(size = 6),
#         legend.key.size = unit(0.2, "cm"),
#         legend.title = element_text(size = 6),
#         axis.text.x = element_text(size = 8),
#         axis.text.y = element_text(size = 6),
#         plot.title = element_text(size = 6, face="plain"),
#         axis.title = element_text(size = 6),
#         plot.tag = element_text(size = 8, face="plain")) + 
#   scale_fill_gradientn(colours=c("#E41A1C", '#FFFFFF', '#377EB8'), limits=c(0, 0.05)) +
#   labs(title = "GSEA GO:BP", tag = "D")

g <- arrangeGrob(p1, p2, p3, p4, ncol = 2, nrow = 3, layout_matrix= rbind(c(1,1), c(2,4), c(3,4)))
ggsave(paste0("Figure2.png"), plot = g, width = 18, height = 20, units = "cm")

