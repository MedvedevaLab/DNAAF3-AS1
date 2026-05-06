library(ggplot2)
library(magick)
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
# library(ChIPpeakAnno)
library(stringr)
library(GenomeInfoDb)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPpeakAnno)
library(png)
library(grid)
library(ggsankey)
library(ggbreak) 

library(UpSetR)

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

################################ Main ################################

gtf_file <- "gencode.v49.annotation.gtf"
txdb     <- makeTxDbFromGFF(gtf_file)

# ── Panel A: Gene tracks ───────────────────────────────────────────────────────

# pyGenomeTracks --tracks plot_supplementary.ini --region chr7:4524000-4536000 -o SupplementaryFigure1A.png --plotWidth 35 --dpi 300
gene_track_img  <- readPNG("SupplementaryFigure1A.png") # 
gene_track_grob <- rasterGrob(gene_track_img,
    interpolate = TRUE, width = unit(4, "in"), height = unit(1.5, "in"))

panel_a <- ggplot(geom = "blank") +
    annotation_custom(gene_track_grob) +
    theme_bw() +
    panel_theme +
    labs(tag = "A")

# gene_track_img  <- readPNG("alignment.png") #
# gene_track_grob <- rasterGrob(gene_track_img,
#     interpolate = TRUE, width = unit(4, "in"), height = unit(1.5, "in"))
# panel_b <- ggplot(geom = "blank") +
#     annotation_custom(gene_track_grob) +
#     theme_bw() +
#     panel_theme +
#     labs(tag = "B")

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

panel_b <- ggplot(hm_counts, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw() +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_brewer(palette = "Set1", guide = guide_legend(reverse = TRUE), drop = FALSE) +
  scale_y_break(c(10, 500)) +  # break the axis between 10 and 1190
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
  labs(title = "Peaks correlated with DNAAF3-AS1", x = "", y = "Number of peaks", fill = "", tag = "B") +
  coord_flip()


# C. Himorna peaks distance from TSS
himorna_peaks <- read.csv("ENSG00000267577.tsv", sep = "\t")
himorna_peaks <- dplyr::filter(himorna_peaks, hm == "H3K36me3")

himorna_peaks_pos_cor <- himorna_peaks[himorna_peaks$corr > 0, ]
himorna_peaks_neg_cor <- himorna_peaks[himorna_peaks$corr < 0, ]
rownames(himorna_peaks_pos_cor) <- himorna_peaks_pos_cor$X
rownames(himorna_peaks_neg_cor) <- himorna_peaks_neg_cor$X

p31 <-  get_tss_dist(himorna_peaks_pos_cor, 'rho > 0', "C", txdb) 
p32 <-  get_tss_dist(himorna_peaks_neg_cor, 'rho < 0', "", txdb) 
panel_c <- ggarrange(p31, p32, ncol = 1, nrow = 2, common.legend = TRUE, legend="bottom")

dex_results_aso1      <- read.csv("aso1_dex_results.csv", sep = "\t")
dex_results_aso1$log2fold_aso_control <- - dex_results_aso1$log2fold_control_aso
dex_results_aso1      <- dex_results_aso1[!is.na(dex_results_aso1$log2fold_aso_control), ]
dex_results_aso1_up   <- dplyr::filter(dex_results_aso1, log2fold_aso_control > 0)
dex_results_aso1_down <- dplyr::filter(dex_results_aso1, log2fold_aso_control < 0)

dex_results_aso3      <- read.csv("aso3_dex_results.csv", sep = "\t")
dex_results_aso3$log2fold_aso_control <- - dex_results_aso3$log2fold_control_aso
dex_results_aso3      <- dex_results_aso3[!is.na(dex_results_aso3$log2fold_aso_control), ]
dex_results_aso3_up   <- dplyr::filter(dex_results_aso3, log2fold_aso_control > 0)
dex_results_aso3_down <- dplyr::filter(dex_results_aso3, log2fold_aso_control < 0)

df <- data.frame(c('ASO_01', 'ASO_01', 'ASO_03', 'ASO_03'),
c('up', 'down', 'up', 'down'),
c(dim(dex_results_aso1_up)[1], dim(dex_results_aso1_down)[1], dim(dex_results_aso3_up)[1], dim(dex_results_aso3_down)[1]))

colnames(df) <- c('aso', 'promoter', 'N')
df$promoter <- factor(df$promoter, levels = c('up', 'down'))
panel_d <- ggplot(df, aes(x = aso, y = N, fill = promoter)) +
  theme_bw() +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_brewer(palette = "Set1", guide = guide_legend(reverse = TRUE), drop = FALSE) +
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
  labs(x = "", y = "Number of promoters", fill = "", tag = "D")

coords_up <- dex_results_aso1_up[, c("genomicData.seqnames", "genomicData.start", "genomicData.end", "genomicData.strand")]
colnames(coords_up) <- c("seqnames", "start", "end", "strand")
granges_up <- toGRanges(coords_up, format = "BED")

coords_down <- dex_results_aso1_down[, c("genomicData.seqnames", "genomicData.start", "genomicData.end", "genomicData.strand")]
colnames(coords_down) <- c("seqnames", "start", "end", "strand")
granges_down <- toGRanges(coords_down, format = "BED")

peak_anno_up_aso1   <- annotatePeak(granges_up,   tssRegion = c(-1000, 100), TxDb = txdb, overlap = "all")
peak_anno_down_aso1 <- annotatePeak(granges_down, tssRegion = c(-1000, 100), TxDb = txdb, overlap = "all")

anno_df_01_up_stat        <- peak_anno_up_aso1@annoStat
anno_df_01_up_stat$Type   <- "up"
anno_df_01_up_stat$ASO   <- "ASO_01"

anno_df_01_down_stat      <- peak_anno_down_aso1@annoStat
anno_df_01_down_stat$Type <- "down"
anno_df_01_down_stat$ASO   <- "ASO_01"


coords_up <- dex_results_aso3_up[, c("genomicData.seqnames", "genomicData.start", "genomicData.end", "genomicData.strand")]
colnames(coords_up) <- c("seqnames", "start", "end", "strand")
granges_up <- toGRanges(coords_up, format = "BED")

coords_down <- dex_results_aso3_down[, c("genomicData.seqnames", "genomicData.start", "genomicData.end", "genomicData.strand")]
colnames(coords_down) <- c("seqnames", "start", "end", "strand")
granges_down <- toGRanges(coords_down, format = "BED")

peak_anno_up_aso3   <- annotatePeak(granges_up,   tssRegion = c(-1000, 100), TxDb = txdb, overlap = "all")
peak_anno_down_aso3 <- annotatePeak(granges_down, tssRegion = c(-1000, 100), TxDb = txdb, overlap = "all")

anno_df_03_up_stat        <- peak_anno_up_aso3@annoStat
anno_df_03_up_stat$Type   <- "up"
anno_df_03_up_stat$ASO   <- "ASO_03"

anno_df_03_down_stat      <- peak_anno_down_aso3@annoStat
anno_df_03_down_stat$Type <- "down"
anno_df_03_down_stat$ASO   <- "ASO_03"

df_list <- list(anno_df_01_up_stat, anno_df_01_down_stat, anno_df_03_up_stat, anno_df_03_down_stat)
annotation <- do.call(rbind, df_list)
annotation$Type <- factor(annotation$Type, levels = c('up', 'down'))

annotation$ASO <- factor(annotation$ASO, levels = c('ASO_01', 'ASO_03'))
annotation$Feature <- factor(annotation$Feature, levels = c("Promoter", "5' UTR", "1st Exon" , "1st Intron", "Other Exon" , "Other Intron" , "Downstream (<=300)"
))

cols <- rev(c("#999999", "#FF7F00", "#4DAF4A", "#984EA3", "#c4e8c3", "#e4cae8", "#FFFF33" ))
panel_e <- ggplot(annotation, aes(x = Type, y = Frequency, fill =  Feature)) + 
  facet_wrap(~ASO) +
  theme_bw() + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=rev(cols), guide = guide_legend(reverse = TRUE), drop = FALSE) +
  theme(panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 6, vjust = 0.5, hjust = 1),
      legend.title = element_text(size = 6),
      legend = "bottom",
      axis.text.y = element_text(size = 6),
      legend.text = element_text(size = 6),
      plot.title = element_text(size = 6),
      legend.key.size = unit(0.2, "cm"),
      plot.tag = element_text(size = 6, face="plain")) + 
      labs(title = "Differential usage promoters annotation", tag = "E")

anno_df_03_up        <- peak_anno_up_aso3@anno
anno_df_03_up$Type   <- "up"

anno_df_03_down      <- peak_anno_down_aso3@anno
anno_df_03_down$Type <- "down"

anno_df_combined             <- rbind(as.data.frame(anno_df_03_up), as.data.frame(anno_df_03_down))
anno_df_combined$annotation  <- str_split_fixed(anno_df_combined$annotation, " ", 2)[, 1]
anno_df_combined             <- dplyr::filter(anno_df_combined, annotation %in% c("Promoter", "Exon", "Intron"))

sankey_df <- anno_df_combined %>% make_long(annotation, Type)

panel_f <- ggplot(sankey_df, aes(x = x, next_x = next_x, node = node,
                                  next_node = next_node, fill = factor(node), label = node)) +
    geom_sankey() +
    geom_sankey_label(size = 3, color = "white", fill = "gray40") +
    theme_bw() +
    panel_theme +
    theme(axis.title.x       = element_blank(),
          axis.ticks.x       = element_blank(),
          axis.text.x        = element_blank()) +
    guides(fill=guide_legend(ncol = 3)) +
    scale_fill_manual(values = c("#4DAF4A", "#984EA3", "#999999", "#E41A1C", "#377EB8")) +
    labs(title = "Promoter annotation ASO_03", fill = "Feature", tag = "F")

figure1_grob <- arrangeGrob(
    panel_a, panel_b, panel_c, panel_d, panel_e, panel_f,
    ncol = 2, nrow = 4,
    layout_matrix = rbind(c(1, 1, 2), c(3, 3, 4), c(5, 5, 6), c(5, 5, 6))
)
ggsave("Supplementary_Figure1.png", plot = figure1_grob, width = 18, height = 24, units = "cm")

# stat_pos <- peakAnno_up@annoStat
# stat_pos$Type <- 'Up'

# stat_neg <- peakAnno_down@annoStat
# stat_neg$Type <- 'Down'

# annotation <- rbind(stat_pos, stat_neg)

# # annotation <- dplyr::filter(annotation, Feature != 'Promoter')
# annotation$Type <- factor(annotation$Type, levels = c('Down', 'Up'))

# p42 <- ggplot(annotation, aes(x = Type, node = Feature, fill = Feature, value = Frequency)) +
#    geom_sankey_bump(space = 0, 
#                    type = "alluvial", 
#                    color = "transparent") +
#    theme_bw() + 
#    theme(panel.border = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.position = "bottom",
#         axis.title.x = element_blank(),#axis.text.x=element_blank(),
#         axis.text.x = element_text(size = 6, hjust = 0.5),
#         legend.title = element_text(size = 6),
#         legend.key.size = unit(0.2, "cm"),
#         legend.text = element_text(size = 5),
#         axis.ticks.x = element_blank(),
#         axis.text.y = element_text(size = 6),
#         axis.title.y = element_blank(),
#         plot.title = element_text(size = 6),
#         plot.tag = element_text(size = 6, face="plain"),
#         legend.spacing.y = unit(2, 'cm')) +
#   # scale_y_break(c(20, 90)) +  # break the axis between 10 and 1190
#   scale_fill_manual(values = c("#999999", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#F781BF")) +
#   labs(title = "Promoter annotation ASO_03", tag = "")









#intersect himorna and promoter annotation
# himorna

install.packages("ggupset")

library(ggupset)
himorna_pos_anno <- as.data.frame(anno_pos@anno)
himorna_pos_anno$annotation <- str_split_fixed(himorna_pos_anno$annotation, " ", 2)[,1]
himorna_neg_anno <- as.data.frame(anno_neg@anno)
himorna_neg_anno$annotation <- str_split_fixed(himorna_neg_anno$annotation, " ", 2)[,1]

himorna_genebody_poscor <- dplyr::filter(himorna_pos_anno, annotation %in% c('Exon', 'Intron'))
himorna_genebody_poscor$Type <- "Himorna rho > 0"
himorna_genebody_negcor <- dplyr::filter(himorna_neg_anno, annotation %in% c('Exon', 'Intron'))
himorna_genebody_negcor$Type <- "Himorna rho < 0"

# promoters

promoter_up_01_anno <- as.data.frame(peakAnno_up_01@anno)
promoter_up_01_anno$annotation <- str_split_fixed(promoter_up_01_anno$annotation, " ", 2)[,1]
peakAnno_up_01_gene_body <- dplyr::filter(promoter_up_01_anno, annotation %in% c('Exon', 'Intron'))
peakAnno_up_01_gene_body$Type <- "promoters ASO01 up"

promoter_down_01_anno <- as.data.frame(peakAnno_down_01@anno)
promoter_down_01_anno$annotation <- str_split_fixed(promoter_down_01_anno$annotation, " ", 2)[,1]
peakAnno_down_01_gene_body <- dplyr::filter(promoter_down_01_anno, annotation %in% c('Exon', 'Intron'))
peakAnno_down_01_gene_body$Type <- "promoters ASO01 down"

promoter_up_03_anno <- as.data.frame(peakAnno_up_03@anno)
promoter_up_03_anno$annotation <- str_split_fixed(promoter_up_03_anno$annotation, " ", 2)[,1]
peakAnno_up_03_gene_body <- dplyr::filter(promoter_up_03_anno, annotation %in% c('Exon', 'Intron'))
peakAnno_up_03_gene_body$Type <- "promoters ASO03 up"

promoter_down_03_anno <- as.data.frame(peakAnno_down_03@anno)
promoter_down_03_anno$annotation <- str_split_fixed(promoter_down_03_anno$annotation, " ", 2)[,1]
peakAnno_down_03_gene_body <- dplyr::filter(promoter_down_03_anno, annotation %in% c('Exon', 'Intron'))
peakAnno_down_03_gene_body$Type <- "promoters ASO03 down"

sum_df <- do.call(rbind, list(himorna_genebody_poscor, himorna_genebody_negcor, peakAnno_up_01_gene_body, peakAnno_down_01_gene_body))


sum_df <- sum_df[,c('geneId', 'Type')]


sum_df$Type <- factor(sum_df$Type, levels = c("Himorna rho > 0", "Himorna rho < 0", "promoters ASO01 up", "promoters ASO01 down"))
sum_df %>%
  group_by(geneId) %>%
  summarize(Type = list(Type)) %>%
  ggplot(aes(x = Type)) +
    geom_bar() +
    scale_x_upset()

    , peakAnno_up_03_gene_body, peakAnno_down_03_gene_body

####################3
library(ComplexUpset)
library(tidyr)
library(dplyr)

upset_data <- sum_df %>%
  dplyr::distinct(geneId, Type) %>%
  dplyr::mutate(value = TRUE) %>%
  tidyr::pivot_wider(
    names_from = Type,
    values_from = value,
    values_fill = FALSE,
    values_fn = any
  ) %>%
  dplyr::select(-geneId) %>%
  as.data.frame()

gene_types <- colnames(upset_data)

# --- Create UpSet plot as ggplot object ---
p <- upset(
  upset_data,
  intersect = gene_types,
  base_annotations = list(
    "Intersection Size" = intersection_size(
      counts = TRUE,
      text = list(size = 4),
      fill = "#457B9D"
    )
  # ),
  # set_sizes = (
  #   upset_set_size(
  #     fill = "#E63946",
  #     color = "white"
  #   )
  ),
  themes = upset_modify_themes(list(
    "intersections_matrix" = theme(text = element_text(size = 12)),
    "overall_sizes"        = theme(text = element_text(size = 12))
  )),
  name = "Gene Type Intersections",
  width_ratio = 0.3,
  sort_intersections_by = "cardinality"
)

# --- p is a ggplot object --- 
print(p)

# Modify like any ggplot
p + ggtitle("Gene Type UpSet Plot") +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))



length(intersect(himorna_genebody_poscor$geneId, peakAnno_up_01_gene_body$geneId))
length(intersect(himorna_genebody_poscor$geneId, peakAnno_down_01_gene_body$geneId))
length(intersect(himorna_genebody_negcor$geneId, peakAnno_up_01_gene_body$geneId))
length(intersect(himorna_genebody_negcor$geneId, peakAnno_down_01_gene_body$geneId))

length(intersect(himorna_genebody_poscor$geneId, peakAnno_up_03_gene_body$geneId))
length(intersect(himorna_genebody_poscor$geneId, peakAnno_down_03_gene_body$geneId))
length(intersect(himorna_genebody_negcor$geneId, peakAnno_up_03_gene_body$geneId))
length(intersect(himorna_genebody_negcor$geneId, peakAnno_down_03_gene_body$geneId))



length(intersect(intersect(himorna_genebody_poscor$geneId, peakAnno_up_01_gene_body$geneId), peakAnno_up_03_gene_body$geneId ))

library(ggupset)
library(dplyr)

# Load example data (from ggupset package)
data(tidy_movies, package = "ggupset")

# Prepare data with a list-column of Genres
# The tidy_movies dataset is already in a suitable format where 'Genres' is a list
# of character vectors (e.g., c("Action", "Thriller"))

# Create the UpSet plot
ggplot(tidy_movies, aes(x = Genres)) +
  geom_bar(fill = "darkgray") +
  scale_x_upset(order_by = "freq", n_intersections = 12) + # Transform x-axis to UpSet matrix
  theme_bw() +
  labs(title = "UpSet Plot of Movie Genres", y = "Count")
#   coords <- ctd_himorna_negative_cor[,c('chrom', 'start', 'end')]
#   coords <- cbind(coords, "*")
#   colnames(coords) <- c('seqnames', 'start', 'end', 'strand')
#   coords$start <- as.numeric(coords$start)
#   coords$end <- as.numeric(coords$end)
#   ranges <- toGRanges(coords, format='BED')
#   peak_ids <- rownames(coords) 


#   peakAnno <- annotatePeak(ranges, tssRegion=c(-2000, 2000),
#     TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb="org.Hs.eg.db", overlap="all")

#   peakAnno_df <- as.data.frame(peakAnno@anno)

#   peakAnno_df_distal <- dplyr::filter(peakAnno_df, annotation == 'Distal Intergenic' )

# p5 <- ggplot(peakAnno_df_distal, aes(x = distanceToTSS)) + 
#   theme_bw() + 
#   geom_histogram(bins = 50, alpha = 0.5, position = "identity") +
#   theme(panel.border = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.position="bottom",
#         axis.title.x=element_text(size = 12),#axis.text.x=element_blank(),
#         axis.text.x = element_text(size = 12, vjust = 0.5, hjust=1),
#         legend.title = element_text(size = 12),
#         legend.text = element_text(margin = margin(r = 20, unit = "pt"), size = 12),
#         axis.text.y = element_text(size = 12),
#         plot.title = element_text(size = 12),
#         plot.tag = element_text(size = 12, face="plain"),
#         legend.spacing.y = unit(2, 'cm'))  + 
#   scale_x_continuous(limits = c(-100000, 100000), labels = label_comma())+
#   guides(fill=guide_legend(ncol = 3)) +
#   geom_vline(xintercept = 0, color = "red", linetype=3) +
#   labs(title = "Distal Intergenic peaks distanceToTSS", x = "distanceToTSS", tag = "D")
