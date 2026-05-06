library(ggplot2)
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
library(grid)
library(ggsankey)
library(tidyr)
library(clusterProfiler)

# ── Shared theme elements ──────────────────────────────────────────────────────
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

################################ Functions ################################

# Annotate genomic peaks from a data frame of coordinates
get_annotation <- function(peak_data, correlation_label, txdb) {
    coords <- peak_data[, c("chrom", "start", "end")]
    coords <- cbind(coords, "*")
    colnames(coords) <- c("seqnames", "start", "end", "strand")
    coords$start <- as.numeric(coords$start)
    coords$end   <- as.numeric(coords$end)
    peak_granges <- toGRanges(coords, format = "BED")
    peak_annotation <- annotatePeak(peak_granges, tssRegion = c(-1000, 100), TxDb = txdb, overlap = "all")
    return(peak_annotation)
}

################################ Main ################################

gtf_file <- "gencode.v49.annotation.gtf"
txdb     <- makeTxDbFromGFF(gtf_file)

# ── Panel A: Gene tracks ───────────────────────────────────────────────────────
# pyGenomeTracks --tracks plot.ini --region chr19:55158653-55178653 -o Figure1A.png --plotWidth 120 --dpi 130
gene_track_img  <- readPNG("Figure1A.png")
gene_track_grob <- rasterGrob(gene_track_img,
    interpolate = TRUE, width = unit(4, "in"), height = unit(1.5, "in"))

panel_a <- ggplot(geom = "blank") +
    annotation_custom(gene_track_grob) +
    theme_bw() +
    panel_theme +
    labs(tag = "A")

# ── Panel B: HiMoRNA peaks annotation ─────────────────────────────────────────
himorna_peaks <- read.csv("ENSG00000267577.tsv", sep = "\t")
himorna_peaks <- dplyr::filter(himorna_peaks, hm == "H3K36me3")

himorna_peaks_pos_cor <- himorna_peaks[himorna_peaks$corr > 0, ]
himorna_peaks_neg_cor <- himorna_peaks[himorna_peaks$corr < 0, ]
rownames(himorna_peaks_pos_cor) <- himorna_peaks_pos_cor$X
rownames(himorna_peaks_neg_cor) <- himorna_peaks_neg_cor$X

peak_anno_pos_cor <- get_annotation(himorna_peaks_pos_cor, "rho > 0", txdb)
peak_anno_neg_cor <- get_annotation(himorna_peaks_neg_cor, "rho < 0", txdb)

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

panel_b <- ggplot(anno_stat_combined, aes(x = Feature, y = Frequency, fill = Type)) +
    theme_bw() +
    geom_bar(position = "dodge", stat = "identity") +
    scale_fill_brewer(palette = "Set1", guide = guide_legend(reverse = TRUE), drop = FALSE) +
    panel_theme +
    guides(fill = guide_legend(ncol = 3)) +
    labs(title = "HiMoRNA peaks annotation", y = "Frequency (%)",
         fill = "Spearman correlation", tag = "B") +
    coord_flip()

# ── Panel C: Differential promoters sankey (ASO_01) ───────────────────────────
dex_results_aso1      <- read.csv("aso1_dex_results_not_sig.csv", sep = "\t")
dex_results_aso1 <- dex_results_aso1[dex_results_aso1$padj < 0.05,]
dex_results_aso1      <- dex_results_aso1[!is.na(dex_results_aso1$log2fold_aso_control), ]
dex_results_aso1_up   <- dplyr::filter(dex_results_aso1, log2fold_aso_control > 0)
dex_results_aso1_down <- dplyr::filter(dex_results_aso1, log2fold_aso_control < 0)

dex_results_aso3      <- read.csv("aso3_dex_results_not_sig.csv", sep = "\t")
dex_results_aso3 <- dex_results_aso3[dex_results_aso3$padj < 0.05,]
dex_results_aso3      <- dex_results_aso3[!is.na(dex_results_aso3$log2fold_aso_control), ]
dex_results_aso3_up   <- dplyr::filter(dex_results_aso3, log2fold_aso_control > 0)
dex_results_aso3_down <- dplyr::filter(dex_results_aso3, log2fold_aso_control < 0)

coords_up <- dex_results_aso1_up[, c("genomicData.seqnames", "genomicData.start", "genomicData.end", "genomicData.strand")]
colnames(coords_up) <- c("seqnames", "start", "end", "strand")
granges_up <- toGRanges(coords_up, format = "BED")

coords_down <- dex_results_aso1_down[, c("genomicData.seqnames", "genomicData.start", "genomicData.end", "genomicData.strand")]
colnames(coords_down) <- c("seqnames", "start", "end", "strand")
granges_down <- toGRanges(coords_down, format = "BED")

peak_anno_up_aso1   <- annotatePeak(granges_up,   tssRegion = c(-1000, 100), TxDb = txdb, overlap = "all")
peak_anno_down_aso1 <- annotatePeak(granges_down, tssRegion = c(-1000, 100), TxDb = txdb, overlap = "all")

anno_df_up        <- peak_anno_up_aso1@anno
anno_df_up$Type   <- "up"
anno_df_down      <- peak_anno_down_aso1@anno
anno_df_down$Type <- "down"

anno_df_combined             <- rbind(as.data.frame(anno_df_up), as.data.frame(anno_df_down))
anno_df_combined$annotation  <- str_split_fixed(anno_df_combined$annotation, " ", 2)[, 1]
anno_df_combined             <- dplyr::filter(anno_df_combined, annotation %in% c("Promoter", "Exon", "Intron"))

sankey_df <- anno_df_combined %>% make_long(annotation, Type)

panel_c <- ggplot(sankey_df, aes(x = x, next_x = next_x, node = node,
                                  next_node = next_node, fill = factor(node), label = node)) +
    geom_sankey() +
    geom_sankey_label(size = 3, color = "white", fill = "gray40") +
    theme_bw() +
    panel_theme +
    theme(axis.title.x = element_blank()) +
    scale_fill_manual(values = c("#4DAF4A", "#984EA3", "#999999", "#E41A1C", "#377EB8")) +
    labs(title = "Promoter annotation ASO_01", fill = "Feature", tag = "C")

# ── Assemble and save figure ───────────────────────────────────────────────────
figure1_grob <- arrangeGrob(
    panel_a, panel_b, panel_c,
    ncol = 2, nrow = 2,
    layout_matrix = rbind(c(1, 3), c(2, 3))
)
ggsave("Figure1.png", plot = figure1_grob, width = 18, height = 15, units = "cm")



# Promoter changes direction

dex_results_aso1 <- read.csv("aso1_dex_results_not_sig.csv", sep = "\t")
dex_results_aso1 <- dex_results_aso1[dex_results_aso1$padj < 0.05,]

dex_results_aso1_up <- dplyr::filter(dex_results_aso1, log2fold_aso_control > 0)
dex_results_aso1_down <- dplyr::filter(dex_results_aso1, log2fold_aso_control < 0)

coords <- dex_results_aso1_up[, c("genomicData.seqnames", "genomicData.start", "genomicData.end", "genomicData.strand")]
colnames(coords) <- c("seqnames", "start", "end", "strand")
granges <- toGRanges(coords, format = "BED")
granges_anno_up   <- annotatePeak(granges,   tssRegion = c(-1000, 100), TxDb = txdb, overlap = "all")
granges_anno_up_df <- as.data.frame(granges_anno_up@anno)
granges_anno_up_df$effect <- 'up'

coords <- dex_results_aso1_down[, c("genomicData.seqnames", "genomicData.start", "genomicData.end", "genomicData.strand")]
colnames(coords) <- c("seqnames", "start", "end", "strand")
granges <- toGRanges(coords, format = "BED")
granges_anno_down   <- annotatePeak(granges,   tssRegion = c(-1000, 100), TxDb = txdb, overlap = "all")
granges_anno_down_df <- as.data.frame(granges_anno_down@anno)
granges_anno_down_df$effect <- 'down'

genes_with_switch <- c()
for(gene in unique(c(granges_anno_up_df$geneId, granges_anno_down_df$geneId)))
{
    granges_anno_up_df_gene <- granges_anno_up_df[granges_anno_up_df$geneId == gene,]
    granges_anno_down_df_gene <- granges_anno_down_df[granges_anno_down_df$geneId == gene,]
    granges_anno_up_df_gene_gene_promoter <- granges_anno_up_df_gene[granges_anno_up_df_gene$annotation == 'Promoter',]
    granges_anno_down_df_gene_gene_promoter <- granges_anno_down_df_gene[granges_anno_down_df_gene$annotation == 'Promoter',]

    if((('+' %in% granges_anno_up_df_gene_gene_promoter$strand) && # both + and - present
       ('-' %in% granges_anno_down_df_gene_gene_promoter$strand)) ||
      (('-' %in% granges_anno_up_df_gene_gene_promoter$strand) &&
       ('+' %in% granges_anno_down_df_gene_gene_promoter$strand)))
    {
        genes_with_switch <- c(genes_with_switch, gene)
    }
}


n_genes_aso01 <- length(unique(dex_results_aso1$groupID))
n_genes_switch_aso01 <- length(unique(genes_with_switch))

dex_results_aso3      <- read.csv("aso3_dex_results_not_sig.csv", sep = "\t")
dex_results_aso3 <- dex_results_aso3[dex_results_aso3$padj < 0.05,]

dex_results_aso3_up <- dplyr::filter(dex_results_aso3, log2fold_aso_control > 0)
dex_results_aso3_down <- dplyr::filter(dex_results_aso3, log2fold_aso_control < 0)

coords <- dex_results_aso3_up[, c("genomicData.seqnames", "genomicData.start", "genomicData.end", "genomicData.strand")]
colnames(coords) <- c("seqnames", "start", "end", "strand")
granges <- toGRanges(coords, format = "BED")
granges_anno_up   <- annotatePeak(granges,   tssRegion = c(-1000, 100), TxDb = txdb, overlap = "all")
granges_anno_up_df <- as.data.frame(granges_anno_up@anno)
granges_anno_up_df$effect <- 'up'

coords <- dex_results_aso3_down[, c("genomicData.seqnames", "genomicData.start", "genomicData.end", "genomicData.strand")]
colnames(coords) <- c("seqnames", "start", "end", "strand")
granges <- toGRanges(coords, format = "BED")
granges_anno_down   <- annotatePeak(granges,   tssRegion = c(-1000, 100), TxDb = txdb, overlap = "all")
granges_anno_down_df <- as.data.frame(granges_anno_down@anno)
granges_anno_down_df$effect <- 'down'

genes_with_switch <- c()
for(gene in unique(c(granges_anno_up_df$geneId, granges_anno_down_df$geneId)))
{
    granges_anno_up_df_gene <- granges_anno_up_df[granges_anno_up_df$geneId == gene,]
    granges_anno_down_df_gene <- granges_anno_down_df[granges_anno_down_df$geneId == gene,]
    granges_anno_up_df_gene_gene_promoter <- granges_anno_up_df_gene[granges_anno_up_df_gene$annotation == 'Promoter',]
    granges_anno_down_df_gene_gene_promoter <- granges_anno_down_df_gene[granges_anno_down_df_gene$annotation == 'Promoter',]

    if((('+' %in% granges_anno_up_df_gene_gene_promoter$strand) && # both + and - present
       ('-' %in% granges_anno_down_df_gene_gene_promoter$strand)) ||
      (('-' %in% granges_anno_up_df_gene_gene_promoter$strand) &&
       ('+' %in% granges_anno_down_df_gene_gene_promoter$strand)))
    {
        genes_with_switch <- c(genes_with_switch, gene)
    }
}

n_genes_aso03 <- length(unique(dex_results_aso3$groupID))
n_genes_switch_aso03 <- length(unique(genes_with_switch))

df <- data.frame(c(n_genes_switch_aso01, n_genes_aso01 - n_genes_switch_aso01, n_genes_switch_aso03, n_genes_aso03 - n_genes_switch_aso03),
    c('ASO_01','ASO_01','ASO_03','ASO_03'),
    c(TRUE, FALSE, TRUE, FALSE))
colnames(df) <- c('n_genes', 'aso', 'bidirectional')

panel_c <- ggplot(df, aes(x = aso, y = n_genes, fill = bidirectional)) +
    theme_bw() +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_brewer(palette = "Set1", guide = guide_legend(reverse = TRUE), drop = FALSE) +
    panel_theme +
    guides(fill = guide_legend(ncol = 3)) +
    labs(title = "", y = "Number of genes",
         fill = "Bidirectional switch", tag = "C")