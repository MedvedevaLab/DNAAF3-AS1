library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(png)
library(RColorBrewer)

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

df$density <- log(df$density, 10)
panel_a <- ggplot(df, aes(x = chromosome, y = experiment, fill = density)) +
  theme_bw() +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", direction = 1, na.value = 'white') +
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

cluster_borders$range_id <- factor(seq_len(nrow(cluster_borders)),
                                   levels = seq_len(nrow(cluster_borders)))

qual_palette <- unique(c(brewer.pal(8, "Dark2"),
                         brewer.pal(9, "Set1"),
                         brewer.pal(8, "Set2"),
                         brewer.pal(12, "Paired")))
chrom_cols <- setNames(qual_palette[seq_along(chrom_levels)], chrom_levels)
chrom_centre <- tapply(as.integer(cluster_borders$range_id),
                       cluster_borders$seqnames,
                       function(v) round(median(v)))
chrom_centre <- chrom_centre[!is.na(chrom_centre)]

panel_b <- ggplot(cluster_borders,
                  aes(x = range_id, y = length, fill = seqnames)) +
  theme_bw() +
  geom_col() +
  scale_fill_manual(values = chrom_cols, drop = FALSE, guide = "none") +
  scale_x_discrete(breaks = as.character(chrom_centre),
                   labels = names(chrom_centre)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  panel_theme +
  theme(axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "", x = "Chromosome",
       y = "Range length (bp)", fill = "Chromosome", tag = "B")


# C, D. GenometriCorr graphical reports
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
ggsave(paste0("Supplementary_Figure7.png"), plot = g, width = 18, height = 20, units = "cm")

