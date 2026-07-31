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
himorna_peaks <- read.csv("ENSG00000267577.tsv", sep = "\t")
himorna_peaks <- dplyr::filter(himorna_peaks, hm == "H3K36me3")

himorna_peaks_pos_cor <- himorna_peaks[himorna_peaks$corr > 0, ]
himorna_peaks_neg_cor <- himorna_peaks[himorna_peaks$corr < 0, ]

himorna_peaks_pos_cor_bed <- himorna_peaks_pos_cor[, c("chrom", "start", "end")]
write.table(
    file      = "himorna_peaks_pos_cor.bed",
    himorna_peaks_pos_cor_bed,
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep       = "\t"
)


himorna_peaks_neg_cor_bed <- himorna_peaks_neg_cor[, c("chrom", "start", "end")]

write.table(
    file      = "himorna_peaks_neg_cor.bed",
    himorna_peaks_neg_cor_bed,
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep       = "\t"
)

rownames(himorna_peaks_pos_cor) <- himorna_peaks_pos_cor$X
rownames(himorna_peaks_neg_cor) <- himorna_peaks_neg_cor$X

peak_anno_pos_cor <- get_annotation(himorna_peaks_pos_cor, txdb, c("chrom", "start", "end"))
peak_anno_neg_cor <- get_annotation(himorna_peaks_neg_cor, txdb, c("chrom", "start", "end"))

anno_stat_pos       <- peak_anno_pos_cor@annoStat
anno_stat_pos$Type  <- "rho > 0"

anno_stat_neg       <- peak_anno_neg_cor@annoStat
anno_stat_neg$Type  <- "rho < 0"

anno_stat_combined           <- rbind(anno_stat_pos, anno_stat_neg)
anno_stat_combined$Type      <- factor(anno_stat_combined$Type, levels = c("rho > 0", "rho < 0"))
anno_stat_combined$Frequency <- as.numeric(anno_stat_combined$Frequency)
anno_stat_combined$Feature   <- factor(anno_stat_combined$Feature,
    levels = rev(c("Promoter", "5' UTR", "1st Exon", "1st Intron",
                   "Other Exon", "Other Intron", "3' UTR", "Distal Intergenic")))

anno_stat_combined <- as.data.frame(tidyr::complete(anno_stat_combined, Type, Feature))

write.table(
    file      = "Supplementary_Table2.csv",
    rbind(as.data.frame(peak_anno_pos_cor@anno),as.data.frame(peak_anno_neg_cor@anno)),
    quote     = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    sep       = "\t"
)

panel_a <- ggplot(anno_stat_combined, aes(x = Feature, y = Frequency, fill = Type)) +
    theme_bw() +
    geom_bar(position = "dodge", stat = "identity") +
    scale_fill_brewer(palette = "Set1", guide = guide_legend(reverse = TRUE), drop = FALSE) +
    panel_theme +
    guides(fill = guide_legend(ncol = 3)) +
    labs(title = "HiMoRNA peaks annotation", y = "Frequency (%)",
         fill = "Spearman correlation", tag = "B") +
    coord_flip()


# Panel B

diffbind_results <- read.csv('/data/budkina/himorna_chipseq_data/test/CTD-2587H24.5_H3K36me3_differential_peaks.csv', sep = '\t')
diffbind_sig <- diffbind_results[diffbind_results$FDR < 0.05,] # & abs(diffbind_results$Fold) > 1

diffbind_sig$significance <- "Unchanged"
diffbind_sig$significance[diffbind_sig$Fold < 0] <- "Lost"
diffbind_sig$significance[diffbind_sig$Fold > 0] <- "Gained"
diffbind_sig$significance <- factor(diffbind_sig$significance, levels = c( "Gained", "Lost", "Unchanged"))

colors <- c( "Gained" = "#E31A1C", "Lost" =  "#1F78B4", "Unchanged" = "#CCCCCC")

panel_b <- ggplot(diffbind_sig, aes(x = Conc, y = Fold, color = significance)) +
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
        tag = "B")

# D. Metagene

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
      legend.text = element_text(size = 5)) + labs(tag = "C")

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

panel_c <- ggarrange(p11, p12, ncol = 2, nrow = 1, common.legend = TRUE, legend="bottom")


gained <- diffbind_sig[diffbind_sig$Fold > 0,]
lost <- diffbind_sig[diffbind_sig$Fold < 0,] 
anno_gained <- get_annotation(gained, txdb)
anno_lost <- get_annotation(lost, txdb)

annotation_gained <- anno_gained@annoStat
annotation_gained$Type <- 'Gained'
annotation_lost <- anno_lost@annoStat
annotation_lost$Type <- 'Lost'
annotation <- rbind(annotation_gained, annotation_lost)
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

annotation$Type <- factor(annotation$Type, levels = c('Gained', 'Lost'))

cols <- rev(c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F"))
panel_d <- ggplot(annotation, aes(x = Type, y = Frequency, fill =  Feature)) + 
  theme_bw() + 
  geom_bar(position="stack", stat="identity") + coord_flip()
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



lost_counts <- diffbind_sig %>%#  dplyr::filter(Fold < 0) %>%
  dplyr::count(seqnames, name = "n_lost")
lost_counts$seqnames <- as.character(lost_counts$seqnames)


# Natural chromosome ordering (chr1..chr22, chrX, chrY, ...)
chrom_levels <- unique(lost_counts$seqnames)
chrom_order  <- chrom_levels[order(
    suppressWarnings(as.numeric(gsub("chr", "", chrom_levels))),
    gsub("chr", "", chrom_levels))]

ref_chr <- if ("chr19" %in% chrom_order) "chr19" else "19"

chrom_lengths <- seqlengths(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Number of genes per chromosome (hg38 knownGene annotation)
genes_gr <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene,
                  single.strand.genes.only = TRUE)
genes_per_chrom <- table(as.character(seqnames(genes_gr)))

# Match chromosome names to the hg38 lengths (add "chr" prefix if missing)
lookup <- ifelse(grepl("^chr", lost_counts$seqnames),
                 lost_counts$seqnames, paste0("chr", lost_counts$seqnames))
lost_counts$chrom_length <- chrom_lengths[lookup]
lost_counts$n_genes      <- as.integer(genes_per_chrom[lookup])
lost_counts <- lost_counts[!is.na(lost_counts$chrom_length) &
                           !is.na(lost_counts$n_genes) &
                           lost_counts$n_genes > 0, ]
# Lost peaks per Mb per gene: normalize by both chromosome length and gene count
lost_counts$lost_per_mb  <- lost_counts$n_lost / lost_counts$chrom_length * 1e6 /
                            lost_counts$n_genes


lost_counts$seqnames  <- factor(lost_counts$seqnames, levels = chrom_order)
lost_counts$highlight <- ifelse(as.character(lost_counts$seqnames) == ref_chr,
                                "chr19", "other")

panel_e <- ggplot(lost_counts, aes(x = seqnames, y = lost_per_mb, fill = highlight)) +
  theme_bw() +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("chr19" = "#1F78B4", "other" = "grey70"), guide = "none") +
  panel_theme +
  theme(axis.title.x = element_text(size = 6),
        axis.title.y = element_text(size = 6, angle = 90),
        axis.text.x  = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Lost peaks (Fold<0) per Mb per gene by chromosome",
       x = "Chromosome", y = "Lost peaks per Mb per gene", tag = "E")



genometri_res <- readRDS('genometri_res_lost.rds')

rel_dist_data <- genometri_res$awhole$relative.distances.data

rel_dist_df <- data.frame(x = rel_dist_data)

proj_pval <- genometri_res$awhole$relative.distances.ecdf.deviation.area.p.value

panel_f <- ggplot(rel_dist_df, aes(x = x)) +
  stat_ecdf(aes(color = "Observed"), geom = "step") +
  geom_abline(aes(color = "Expected (y = 2x)"), slope = 2, intercept = 0) +
  scale_color_manual(values = c("Observed" = "#E31A1C",
                                "Expected (y = 2x)" = "black"),
                     breaks = c("Observed", "Expected (y = 2x)")) +
  guides(color = guide_legend(ncol = 2)) +
  annotate("text", x = 0.02, y = 0.95, hjust = 0, size = 2,
           label = paste0("ecdf deviation area p-value = ", proj_pval)) +
  theme_bw() +
  panel_theme +
  labs(title = "Relative distance GenometriCorr",
       x = "Fractional distance", y = "Cumulative fraction",
       color = "", tag = "F") +
  theme(
    legend.key.size    = unit(0.5, "cm"), 
    legend.spacing.y   = unit(1, "cm")
)






g <- arrangeGrob(panel_a, panel_b, panel_c, panel_d, panel_e, panel_f, ncol = 2, nrow = 5, layout_matrix= rbind(c(1,2), c(3,3), c(4,4), c(5,6), c(5,6)))
ggsave(paste0("Figure3.png"), plot = g, width = 18, height = 26, units = "cm")

# pyGenomeTracks --tracks plot_rfx2.ini --region chr19:5979119-6220000 -o Figure3A.png --plotWidth 90 --dpi 130
