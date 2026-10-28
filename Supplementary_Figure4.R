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

cb_blue   <- "#0072B2"
cb_orange <- "#D55E00"
cb_green  <- "#009E73"
cb_grey   <- "#999999"

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

# A. Himorna peaks counts
# Read himorna peaks
ctd_himorna <- read.csv('ENSG00000267577.tsv', sep = '\t')
ctd_himorna$peak_type <- NA
ctd_himorna[ctd_himorna$corr > 0, 'peak_type'] <- 'rho > 0.5'
ctd_himorna[ctd_himorna$corr < 0, 'peak_type'] <- 'rho < -0.5'

# Count peaks H3K36me3
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
  scale_fill_manual(values = c("rho > 0.5" = cb_orange, "rho < -0.5" = cb_blue),
                    guide = guide_legend(reverse = TRUE), drop = FALSE) + panel_theme +
  guides(fill=guide_legend(ncol = 2)) +
  labs(title = "Peaks correlated with DNAAF3-AS1", x = "", y = "Number of peaks", fill = "", tag = "A") +
  coord_flip()

# B. Himorna peaks distance from TSS
himorna_peaks <- read.csv("ENSG00000267577.tsv", sep = "\t")
himorna_peaks <- dplyr::filter(himorna_peaks, hm == "H3K36me3")
himorna_peaks_pos_cor <- himorna_peaks[himorna_peaks$corr > 0, ]
himorna_peaks_neg_cor <- himorna_peaks[himorna_peaks$corr < 0, ]
rownames(himorna_peaks_pos_cor) <- himorna_peaks_pos_cor$X
rownames(himorna_peaks_neg_cor) <- himorna_peaks_neg_cor$X

p31 <-  get_tss_dist(himorna_peaks_pos_cor, 'rho > 0', "B", txdb) 
p32 <-  get_tss_dist(himorna_peaks_neg_cor, 'rho < 0', "", txdb) 
panel_b <- ggarrange(p31, p32, ncol = 1, nrow = 2, common.legend = TRUE, legend="bottom")

# C. Differential expression of H3K36me3 modifier
fantom <- read.csv('DESeq2_genes_ASO_all.tsv', sep = '\t')
aso01 <- dplyr::filter(fantom, perturb_id == 'ASO_G0267577_01')
aso01 <- aso01[aso01$fdr < 0.05, ]
aso03 <- dplyr::filter(fantom, perturb_id == 'ASO_G0267577_03')
aso03 <- aso03[aso03$fdr < 0.05, ]

fantom_subset <- dplyr::filter(fantom, geneSymbol %in% c('SETD2',
  'NSD1',
  'WHSC1',    # current symbol: NSD2
  'WHSC1L1',  # current symbol: NSD3
  'NSD2',
  'NSD3',
  'ASH1L',
  'SETMAR',
  'SMYD2',
  'SMYD5',
  'SETD3',
  'KDM2A',
  'KDM2B',
  'KDM4A',
  'KDM4B',
  'KDM4C',
  'KDM4D'))

fantom_subset <- dplyr::filter(fantom_subset, perturb_id %in% c('ASO_G0267577_01', 'ASO_G0267577_03'))
fantom_subset$geneSymbol <- dplyr::recode(fantom_subset$geneSymbol,
                                          WHSC1 = "NSD2", WHSC1L1 = "NSD3")
fantom_subset$sig <- ""
fantom_subset[fantom_subset$fdr < 0.05,  'sig'] <- '*'
fantom_subset[fantom_subset$fdr < 0.01,  'sig'] <- '**'
fantom_subset[fantom_subset$fdr < 0.001, 'sig'] <- '***'

panel_c <- ggplot(fantom_subset, aes(x = geneSymbol, y = perturb_id, fill = log2FC))+
theme_bw() +
geom_tile() +
scale_fill_gradient2(low = cb_blue, mid = "white", high = cb_orange, midpoint = 0, limits = c(-1.5, 1.5)) +
geom_text(data=fantom_subset, aes(x = geneSymbol, y = perturb_id, label = sig, angle = 0), size=2.5) +
scale_x_discrete(drop = FALSE) +
scale_y_discrete(drop = FALSE) +
panel_theme + theme(
    legend.key.size    = unit(0.6, "cm"),
    axis.title.y       = element_blank()
) + labs(title = "", tag = "C")

g <- arrangeGrob(panel_a, panel_b, panel_c, ncol = 2, nrow =2, layout_matrix= rbind(c(1,2), c(3,3)))
ggsave(paste0("Supplementary_Figure4.png"), plot = g, width = 18, height = 15, units = "cm")
