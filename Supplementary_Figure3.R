library(ggplot2)
library(ggrepel)
library(ggpubr)
library(grid)
library(gridExtra)
library(ChIPseeker)
library(GenomeInfoDb)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPpeakAnno)

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

# Generate GRanges for peaks
get_tss_dist <- function(ctd_himorna, corr_type, label, txdb) 
{
  coords <- ctd_himorna[,c('chrom', 'start', 'end')]
  coords <- cbind(coords, "*")
  colnames(coords) <- c('seqnames', 'start', 'end', 'strand')
  coords$start <- as.numeric(coords$start)
  coords$end <- as.numeric(coords$end)
  ranges <- toGRanges(coords, format='BED')

  peakAnno <- annotatePeak(ranges, tssRegion=c(-1000, 100), TxDb=txdb, overlap="all")

  p <- plotDistToTSS(peakAnno, distanceColumn = "distanceToTSS", xlab = "", ylab = corr_type) +
  theme_bw() +
  panel_theme +
  # scale_fill_brewer(palette = "Set2", guide = guide_legend(reverse = TRUE), drop = FALSE) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6, vjust = 0.5, hjust = 1),
        legend.title = element_text(size = 6),
        legend="bottom",
        axis.text.y = element_text(size = 6),
        plot.title = element_text(size = 6),
        legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 5),
        plot.tag = element_text(size = 6, face="plain")) +
  labs(tag = label, title = "H3K36me3 peaks (%) (5'->3')")

  return(p)
}


gtf_file <- "gencode.v49.annotation.gtf"
txdb     <- makeTxDbFromGFF(gtf_file)

# B. Himorna peaks counts
# Read himorna peaks
ctd_himorna <- read.csv('ENSG00000267577.tsv', sep = '\t')
ctd_himorna$peak_type <- NA
ctd_himorna[ctd_himorna$corr > 0, 'peak_type'] <- 'rho > 0.5'
ctd_himorna[ctd_himorna$corr < 0, 'peak_type'] <- 'rho < -0.5'

# Count peaks h3k36me3
hm_counts <- as.data.frame(table(ctd_himorna$hm, ctd_himorna$peak_type), dnn = c("Histone mark", "Sperman rho"))

hm_counts$Var1   <- factor(hm_counts$Var1,
    levels = rev(c('H3K36me3', 'H3K27me3', 'H3K4me3', 'H3K27ac', 'H3K4me1', 'H3K4me2', 'H3K9ac', 'H3K9me3', 'H3K79me2', 'H4K20me1')))

hm_counts <- as.data.frame(tidyr::complete(hm_counts, Var2, Var1))

panel_a <- ggplot(hm_counts, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw() +
  geom_bar(position="dodge", stat="identity") +
  geom_text(aes(label = ifelse(is.na(Freq) | Freq == 0, "", Freq)),
            position = position_dodge(width = 0.9),
            hjust = -0.2, size = 1.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_brewer(palette = "Set1", guide = guide_legend(reverse = TRUE), drop = FALSE) +#  scale_y_break(c(10, 500)) +  # break the axis between 10 and 1190
  theme(panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position = 'bottom',
  text = element_text(size = 7),
  axis.text.x = element_text(size = 6),
  axis.ticks.x = element_blank(),
  axis.title.y = element_text(size = 6),
  axis.text.y = element_text(size = 6),
  legend.key.size = unit(0.2, "cm"),
  legend.text = element_text(size = 5),
  plot.title = element_text(size = 6, hjust = 0.5),
  plot.tag = element_text(size = 8, face="plain")) +
  guides(fill=guide_legend(ncol = 2)) +
  labs(title = "Peaks correlated with DNAAF3-AS1", x = "", y = "Number of peaks", fill = "", tag = "A") +
  coord_flip()


# C. Himorna peaks distance from TSS
himorna_peaks <- read.csv("ENSG00000267577.tsv", sep = "\t")
himorna_peaks <- dplyr::filter(himorna_peaks, hm == "H3K36me3")

himorna_peaks_pos_cor <- himorna_peaks[himorna_peaks$corr > 0, ]
himorna_peaks_neg_cor <- himorna_peaks[himorna_peaks$corr < 0, ]
rownames(himorna_peaks_pos_cor) <- himorna_peaks_pos_cor$X
rownames(himorna_peaks_neg_cor) <- himorna_peaks_neg_cor$X

p31 <-  get_tss_dist(himorna_peaks_pos_cor, 'rho > 0', "B", txdb) 
p32 <-  get_tss_dist(himorna_peaks_neg_cor, 'rho < 0', "", txdb) 
panel_b <- ggarrange(p31, p32, ncol = 1, nrow = 2, common.legend = TRUE, legend="bottom")


############

fantom <- read.csv('DESeq2_genes_ASO_all.tsv', sep = '\t')
aso01 <- dplyr::filter(fantom, perturb_id == 'ASO_G0267577_01')
aso01 <- aso01[aso01$fdr < 0.05, ]

aso03 <- dplyr::filter(fantom, perturb_id == 'ASO_G0267577_03')
aso03 <- aso03[aso03$fdr < 0.05, ]

fantom_subset <- dplyr::filter(fantom, geneSymbol %in% c('SETD2',
  'NSD1',
  'WHSC1',
  'WHSC1L1',
  'ASH1L',
  'SETMAR',
  'SMYD2',
  'SMYD5',
  'SETD3',
  'SETDB1',
  'KDM2A',
  'KDM2B',
  'KDM4A',
  'KDM4B',
  'KDM4C',
  'KDM4D',
  'DNAH3','KIF27','RAB3IP','AXDND1','CFAP54','DNAH2','CCDC13','CC2D2B'))

fantom_subset <- dplyr::filter(fantom_subset, perturb_id %in% c('ASO_G0267577_01', 'ASO_G0267577_03'))
fantom_subset$sig <- ""
fantom_subset[fantom_subset$fdr < 0.05, 'sig'] <- '*'
fantom_subset[fantom_subset$fdr < 0.001, 'sig'] <- '**'
fantom_subset[fantom_subset$fdr < 0.0001, 'sig'] <- '***'


panel_c <- ggplot(fantom_subset, aes(x = geneSymbol, y = perturb_id, fill = log2FC))+
theme_bw() +
geom_tile() +
scale_fill_distiller(palette = "RdBu", limits = c(-1.5, 1.5)) +
geom_text(data=fantom_subset, aes(x = geneSymbol, y = perturb_id, label = sig, angle = 0), size=2.5) +
scale_x_discrete(drop = FALSE) +
scale_y_discrete(drop = FALSE) +
theme(panel.border = element_blank(), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      text = element_text(size = 6),
      axis.text.x = element_text(angle = 0, size = 6),
      axis.text.y = element_text(angle = 0, size = 6),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 6, face="plain"),
      plot.tag = element_text(size = 8, face="plain"),
      legend.position="right",
      legend.key.size = unit(0.2, "cm"),
      legend.text = element_text(size = 5),
      legend.title=element_text(size = 6)) +
labs(title = "", tag = "C")




g <- arrangeGrob(panel_a, panel_b, panel_c, ncol = 2, nrow =2, layout_matrix= rbind(c(1,2), c(3,3)))
ggsave(paste0("Supplementary_Figure3.png"), plot = g, width = 18, height = 15, units = "cm")
