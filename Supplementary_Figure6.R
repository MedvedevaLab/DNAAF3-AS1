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
library(scales)

library(biomaRt)

mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensembl <- mart

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

# A. Heatmap with contacts density

df <- read.csv('rnachrom_heatmap.csv', header = F, sep =',')
colnames(df) <- c('experiment', 'chromosome', 'density')

df$chromosome <- factor(df$chromosome, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
"chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
"chr19", "chr20", "chr21", "chr22", "chrX", "chrY"))

df$experiment <- factor(df$experiment, levels = c("Exp.ID: 1; iMARGI (AC010327.4)",
"Exp.ID: 2; iMARGI (AC010327.4)",
"Exp.ID: 3; iMARGI (AC010327.4)",
"Exp.ID: 4; iMARGI (AC010327.4)",
"Exp.ID: 5; iMARGI (AC010327.4)",
"Exp.ID: 7; GRID (AC010327.4)",
"Exp.ID: 8; GRID (AC010327.4)",
"Exp.ID: 9; Red-C (AC010327.4)",
"Exp.ID: 10; Red-C (AC010327.4)"))

# df[is.na(df$density), 'density'] = 0
df$density <- log(df$density, 10)
panel_a <- ggplot(df, aes(x = chromosome, y = experiment, fill = density)) +
  theme_bw() +
  geom_tile() +
  scale_fill_distiller(palette = "Reds", direction = 1, na.value = 'white') +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  panel_theme +
  theme(strip.text.x = element_text(size = 6),
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = 6),
        axis.title.y = element_blank()) +
  labs(title = "", x = 'Chromosome', y = 'Experiment', fill = 'log10(Contacts Density)', tag = "A")



# B. Histogram with cluster size
contact_clusters <- readRDS("contact_clusters_hg.rds")

cluster_borders <- contact_clusters %>%
  group_by(seqnames, cluster) %>%
  summarize(cluster_start = min(start), cluster_end = max(end), .groups = "drop")

# Length of each cluster range, one bar per range, ordered by chromosome
chrom_levels <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                  "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                  "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
                  "chr22", "chrX", "chrY")

cluster_borders <- cluster_borders %>%
  mutate(
    length     = cluster_end - cluster_start,
    seqnames   = factor(as.character(seqnames), levels = chrom_levels)
  ) %>%
  filter(!is.na(seqnames)) %>%
  arrange(seqnames, cluster_start)

# fixed x order: ranges sorted by chromosome (then by start within a chromosome)
cluster_borders$range_id <- factor(seq_len(nrow(cluster_borders)),
                                   levels = seq_len(nrow(cluster_borders)))

panel_b <- ggplot(cluster_borders,
                  aes(x = range_id, y = length, fill = seqnames)) +
  theme_bw() +
  geom_col() +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(
      length(levels(cluster_borders$seqnames))), drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  panel_theme +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(ncol = 8)) +
  labs(title = "", x = "Cluster range (sorted by chromosome)",
       y = "Range length (bp)", fill = "Chromosome", tag = "B")


# C, D. GenometriCorr graphical reports, loaded from PNG as raster panels
# (same readPNG/rasterGrob/annotation_custom approach as panel_a in Figure1.R)
read_png_panel <- function(png_file, tag) {
    img  <- readPNG(png_file)
    grob <- rasterGrob(img, interpolate = TRUE)
    ggplot(geom = "blank") +
        annotation_custom(grob) +
        theme_bw() +
        panel_theme +
        labs(tag = tag)
}

panel_c <- read_png_panel("genometri_res_lost.png",       "C")
panel_d <- read_png_panel("genometri_res_trans_lost.png", "D")

g <- arrangeGrob(panel_a, panel_b, panel_c, panel_d, ncol = 2, nrow = 3,
                 layout_matrix = rbind(c(1, 2), c(3, 3), c(4, 4)))
ggsave(paste0("Supplementary_Figure6.png"), plot = g, width = 18, height = 20, units = "cm")

