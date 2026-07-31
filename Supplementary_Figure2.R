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
library("ggVennDiagram")

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


df <- read.csv('cell_profiler.csv', sep = '\t', dec = ',')
df <- df[df$KD.pubID == 'CTD-2587H24.5', ]
df <- df[df$perturb_id %in% c("ASO_G0267577_01", "ASO_G0267577_03"),]

params <- c(
  "GrowthRate",
  "Area",
  "Compactness",
  "Eccentricity",
  "Extent",
  "FormFactor",
  "MajorAxisLength",
  "MaxFeretDiameter",
  "MaximumRadius",
  "MedianRadius",
  "MinFeretDiameter",
  "MinorAxisLength",
  "Perimeter",
  "Solidity"
)

mean_cols <- paste0(params, "_mean")
fdr_cols  <- paste0(params, "_fdr")

# means -> long
mean_long <- melt(df[, c("perturb_id", mean_cols)],
                  id.vars = "perturb_id",
                  variable.name = "parameter",
                  value.name = "mean")
mean_long$parameter <- sub("_mean$", "", mean_long$parameter)

# fdr -> long
fdr_long <- melt(df[, c("perturb_id", fdr_cols)],
                 id.vars = "perturb_id",
                 variable.name = "parameter",
                 value.name = "fdr")
fdr_long$parameter <- sub("_fdr$", "", fdr_long$parameter)

# combine
plot_df <- merge(mean_long, fdr_long, by = c("perturb_id", "parameter"))

# keep column order as requested
plot_df$parameter <- factor(plot_df$parameter, levels = params)

# significance stars
plot_df$star <- ifelse(!is.na(plot_df$fdr) & plot_df$fdr < 0.05, "*", "")

panel_a <- ggplot(plot_df, aes(x = parameter, y = perturb_id, fill = mean)) +
  geom_tile(color = "white") +
  geom_text(aes(label = star), size = 6, vjust = 0.75) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1) +
  labs(x = NULL, y = NULL, fill = "mean") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Panel A
img <- readPNG('MDM2.png')
g1 <- rasterGrob(img,
  interpolate=TRUE, width = unit(6,"in"), height=unit(1.2,"in"))

panel_b1<- ggplot(geom="blank") +
  annotation_custom(g1) +
  theme_bw() +
  panel_theme +
  theme(
      axis.title.x = element_blank(),
      legend.key.size = unit(0.2, "cm"),
      legend.text = element_text(size = 5)) + labs(tag = "B")


img <- readPNG('NR2F2.png')
g2 <- rasterGrob(img,
  interpolate=TRUE, width = unit(6,"in"), height=unit(1.2,"in"))

panel_b2<- ggplot(geom="blank") +
  annotation_custom(g2) +
  theme_bw() +
  panel_theme +
  theme(
      axis.title.x = element_blank(),
      legend.key.size = unit(0.2, "cm"),
      legend.text = element_text(size = 5)) + labs(tag = "B")


img <- readPNG('FASTK.png')
g3 <- rasterGrob(img,
  interpolate=TRUE, width = unit(6,"in"), height=unit(1.2,"in"))

panel_b3<- ggplot(geom="blank") +
  annotation_custom(g3) +
  theme_bw() +
  panel_theme +
  theme(
      axis.title.x = element_blank(),
      legend.key.size = unit(0.2, "cm"),
      legend.text = element_text(size = 5)) + labs(tag = "B")


img <- readPNG('CFDP1.png')
g4 <- rasterGrob(img,
  interpolate=TRUE, width = unit(6,"in"), height=unit(1.2,"in"))

panel_b4<- ggplot(geom="blank") +
  annotation_custom(g4) +
  theme_bw() +
  panel_theme +
  theme(
      axis.title.x = element_blank(),
      legend.key.size = unit(0.2, "cm"),
      legend.text = element_text(size = 5)) + labs(tag = "B")

panel_b <- ggarrange(panel_b1, panel_b2, panel_b3, panel_b4, ncol = 1, common.legend = TRUE, legend = "bottom")


figure2_grob <- arrangeGrob(
    panel_a, panel_b1, panel_b2, panel_b3, panel_b4,
    ncol = 1, nrow = 5,
    layout_matrix = rbind(c(1), c(2), c(3), c(4), c(5))
)
ggsave("Supplementary_Figure2.png", plot = figure2_grob, width = 18, height = 24, units = "cm")


