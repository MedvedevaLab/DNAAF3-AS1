library(ggplot2)
library(ggpubr)
library(ChIPseeker)
library(png)
library(gridExtra)
library(dplyr)
library(ChIPpeakAnno)
library(grid)
library(GenometriCorr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

cb_blue   <- "#0072B2"
cb_orange <- "#D55E00"
cb_grey   <- "#999999"

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

# Panel A: Himorna peaks genome annotation
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

write.csv(anno_stat_combined, "Figure3A_data.csv", row.names = FALSE)

# Proportion test for the 1st intron
anno_pos_df <- as.data.frame(peak_anno_pos_cor)
anno_neg_df <- as.data.frame(peak_anno_neg_cor)

is_first_intron <- function(x) grepl("intron 1 of", x$annotation, ignore.case = TRUE)

n_first_intron_pos   <- sum(is_first_intron(anno_pos_df))
n_first_intron_neg   <- sum(is_first_intron(anno_neg_df))
n_first_intron_total <- n_first_intron_pos + n_first_intron_neg

prop_test_first_intron <- prop.test(
    x           = n_first_intron_neg,
    n           = n_first_intron_total,
    p           = 0.5,
    alternative = "less"
)

cat(sprintf(
    "1st intron himorna peaks: %d negatively correlated / %d total (%.1f%%)\n",
    n_first_intron_neg, n_first_intron_total,
    100 * n_first_intron_neg / n_first_intron_total))
print(prop_test_first_intron)

# proportion test for other intron 
is_other_intron <- function(x) {
    grepl("intron", x$annotation, ignore.case = TRUE) &
        !grepl("intron 1 of", x$annotation, ignore.case = TRUE)
}

n_other_intron_pos   <- sum(is_other_intron(anno_pos_df))
n_other_intron_neg   <- sum(is_other_intron(anno_neg_df))
n_other_intron_total <- n_other_intron_pos + n_other_intron_neg

prop_test_other_intron <- prop.test(
    x           = n_other_intron_neg,
    n           = n_other_intron_total,
    p           = 0.5,
    alternative = "less"
)

cat(sprintf(
    "Other intron himorna peaks: %d negatively correlated / %d total (%.1f%%)\n",
    n_other_intron_neg, n_other_intron_total,
    100 * n_other_intron_neg / n_other_intron_total))
print(prop_test_other_intron)

# Panel B: DiffBind results
diffbind_results <- read.csv('CTD-2587H24.5_H3K36me3_differential_peaks.csv', sep = '\t')
diffbind_sig <- diffbind_results[diffbind_results$FDR < 0.05,] # & abs(diffbind_results$Fold) > 1

diffbind_sig$significance <- "Unchanged"
diffbind_sig$significance[diffbind_sig$Fold < 0] <- "Lost"
diffbind_sig$significance[diffbind_sig$Fold > 0] <- "Gained"
diffbind_sig$significance <- factor(diffbind_sig$significance, levels = c( "Gained", "Lost", "Unchanged"))

write.csv(diffbind_sig, "Figure3B_data.csv", row.names = FALSE)

# Panel C: read from png

# Panel D: Gained/Lost peaks genome annotation
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

write.csv(annotation, "Figure3D_data.csv", row.names = FALSE)

# Panel E: Lost peaks (Fold<0) per Mb per gene by chromosome
lost_counts <- diffbind_sig %>% dplyr::filter(Fold < 0) %>%
  dplyr::count(seqnames, name = "n_lost")
lost_counts$seqnames <- as.character(lost_counts$seqnames)

chrom_levels <- unique(lost_counts$seqnames)
chrom_order  <- chrom_levels[order(
    suppressWarnings(as.numeric(gsub("chr", "", chrom_levels))),
    gsub("chr", "", chrom_levels))]

ref_chr <- if ("chr19" %in% chrom_order) "chr19" else "19"

chrom_lengths <- seqlengths(TxDb.Hsapiens.UCSC.hg38.knownGene)

genes_gr <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene,
                  single.strand.genes.only = TRUE)
genes_per_chrom <- table(as.character(seqnames(genes_gr)))

lookup <- ifelse(grepl("^chr", lost_counts$seqnames),
                 lost_counts$seqnames, paste0("chr", lost_counts$seqnames))
lost_counts$chrom_length <- chrom_lengths[lookup]
lost_counts$n_genes      <- as.integer(genes_per_chrom[lookup])
lost_counts <- lost_counts[!is.na(lost_counts$chrom_length) &
                           !is.na(lost_counts$n_genes) &
                           lost_counts$n_genes > 0, ]

lost_counts$lost_per_mb  <- lost_counts$n_lost / lost_counts$chrom_length * 1e6 /
                            lost_counts$n_genes


lost_counts$seqnames  <- factor(lost_counts$seqnames, levels = chrom_order)
lost_counts$highlight <- ifelse(as.character(lost_counts$seqnames) == ref_chr,
                                "chr19", "other")

write.csv(lost_counts, "Figure3E_data.csv", row.names = FALSE)


# Plot all data
anno_stat_combined <- read.csv("Figure3A_data.csv")
anno_stat_combined$Type    <- factor(anno_stat_combined$Type, levels = c("rho > 0", "rho < 0"))
anno_stat_combined$Feature <- factor(anno_stat_combined$Feature,
    levels = rev(c("Promoter", "5' UTR", "1st Exon", "1st Intron",
                   "Other Exon", "Other Intron", "3' UTR", "Distal Intergenic")))

panel_a <- ggplot(anno_stat_combined, aes(x = Feature, y = Frequency, fill = Type)) +
    theme_bw() +
    geom_bar(position = "dodge", stat = "identity") +
    scale_fill_manual(values = c("rho > 0" = cb_orange , "rho < 0" = cb_blue ),
                      guide = guide_legend(reverse = TRUE), drop = FALSE) +
    panel_theme +
    guides(fill = guide_legend(ncol = 3)) +
    labs(title = "", y = "Frequency (%)",
         fill = "Spearman", tag = "A") +
    coord_flip()

saveRDS(panel_a, "Figure3A.rds")

diffbind_sig <- read.csv("Figure3B_data.csv")
diffbind_sig$significance <- factor(diffbind_sig$significance, levels = c("Gained", "Lost", "Unchanged"))
colors <- c( "Gained" =  cb_orange , "Lost" = cb_blue, "Unchanged" = cb_grey)

panel_b <- ggplot(diffbind_sig, aes(x = Conc, y = Fold, color = significance)) +
  theme_bw() +
  geom_point(alpha = 0.9, size = 0.8) +
  scale_color_manual(values = colors, name = "Peak Status") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", alpha = 0.5) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +  
  panel_theme + 
  theme(axis.title.y = element_text(size = 6, angle = 90)) +
  labs(title = "",
        x = "Log concentration",
        y = "Log2FC (KD/Control)",
        tag = "B")

saveRDS(panel_b, "Figure3B.rds")

# Panel C: Metagene
img <- readPNG('Figure2A1.png')
g11 <- rasterGrob(img,
  interpolate=TRUE, width = unit(3,"in"), height=unit(1.8,"in"))

p11 <- ggplot(geom="blank") +
  annotation_custom(g11) + 
  theme_bw() + 
  panel_theme + 
  theme(
      axis.title.x = element_blank(),
      legend.key.size = unit(0.2, "cm"),
      legend.text = element_text(size = 5)) + labs(tag = "C")

img <- readPNG('Figure2A2.png')
g12 <- rasterGrob(img,
  interpolate=TRUE, width = unit(3,"in"), height=unit(1.8,"in"))

p12 <- ggplot(geom="blank") +
  annotation_custom(g12) + 
  theme_bw() +
  panel_theme +
  theme(
      axis.title.x = element_blank(),
      legend.key.size = unit(0.2, "cm"),
      legend.text = element_text(size = 5)) + labs(tag = "")

panel_c <- ggarrange(p11, p12, ncol = 2, nrow = 1, common.legend = TRUE, legend="bottom")
saveRDS(panel_c, "Figure3C.rds")

# Panel D
annotation <- read.csv("Figure3D_data.csv")
annotation$Feature <- factor(annotation$Feature, levels = rev(c("Promoter",
"5' UTR",
"1st Exon",
"1st Intron",
"Other Exon",
"Other Intron",
"3' UTR",
"Downstream (<=300)",
"Distal Intergenic"
)))
annotation$Type <- factor(annotation$Type, levels = c('Gained', 'Lost'))

panel_d <- ggplot(annotation, aes(x = Type, y = Frequency, fill =  Feature)) +
  theme_bw() + 
  geom_bar(position="stack", stat="identity") + coord_flip() +
  scale_fill_viridis_d(option = "viridis", guide = guide_legend(reverse = FALSE), drop = FALSE) +
  panel_theme +
  labs(title = "", tag = "D")

saveRDS(panel_d, "Figure3D.rds")

# Panel E
lost_counts <- read.csv("Figure3E_data.csv")
lost_counts$seqnames <- factor(lost_counts$seqnames, levels = chrom_order)

panel_e <- ggplot(lost_counts, aes(x = seqnames, y = lost_per_mb, fill = highlight)) +
  theme_bw() +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("chr19" = cb_orange, "other" = cb_grey), guide = "none") +
  panel_theme +
  theme(axis.title.x = element_text(size = 6),
        axis.title.y = element_text(size = 6, angle = 90),
        axis.text.x  = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "",
       x = "Chromosome", y = "Lost peaks per Mb per gene", tag = "E")

saveRDS(panel_e, "Figure3E.rds")


# Panel F: Genometricorr distribution
genometri_res <- readRDS("genometri_res_trans_lost.rds")

rel_dist_data <- genometri_res$awhole$relative.distances.data
rel_dist_df <- data.frame(x = rel_dist_data)
proj_pval <- genometri_res$awhole$relative.distances.ecdf.deviation.area.p.value

write.csv(rel_dist_df, "Figure3F_data.csv", row.names = FALSE)


rel_dist_df <- read.csv("Figure3F_data.csv")

expected_df <- data.frame(x = c(0, 0.5), y = c(0, 1))

panel_f <- ggplot(rel_dist_df, aes(x = x)) +
  stat_ecdf(aes(color = "Observed"), geom = "step") +
  geom_line(data = expected_df, aes(x = x, y = y, color = "Expected (y = 2x)")) +
  scale_color_manual(values = c("Observed" = cb_blue,
                                "Expected (y = 2x)" = "black"),
                     breaks = c("Observed", "Expected (y = 2x)")) +
  annotate("text", x = 0.02, y = 0.95, hjust = 0, size = 2,
           label = paste0("ecdf deviation area p-value = ", proj_pval)) +
  theme_bw() +
  panel_theme +
  labs(title = "",
       x = "Fractional distance", y = "Cumulative fraction",
       color = "", tag = "F") +
  theme(
   legend.key.size    = unit(0.5, "cm")) +
   guides(color = guide_legend(ncol = 2))

saveRDS(panel_f, "Figure3F.rds")



panel_a <- readRDS("Figure3A.rds")
panel_b <- readRDS("Figure3B.rds")
panel_c <- readRDS("Figure3C.rds")
panel_d <- readRDS("Figure3D.rds")
panel_e <- readRDS("Figure3E.rds")
panel_f <- readRDS("Figure3F.rds")


g <- arrangeGrob(panel_a, panel_b, panel_c, panel_d, panel_e, panel_f, ncol = 2, nrow = 5, layout_matrix= rbind(c(1,2), c(3,3), c(4,4), c(5,6)))
ggsave(paste0("Figure3.png"), plot = g, width = 17, height = 20, units = "cm")
