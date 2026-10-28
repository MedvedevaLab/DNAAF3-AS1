library(ggplot2)
library(png)
library(gridExtra)
library(reshape2)
library(RColorBrewer)
library(grid)

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
mean_long <- melt(df[, c("perturb_id", mean_cols)],
                  id.vars = "perturb_id",
                  variable.name = "parameter",
                  value.name = "mean")
mean_long$parameter <- sub("_mean$", "", mean_long$parameter)

fdr_long <- melt(df[, c("perturb_id", fdr_cols)],
                 id.vars = "perturb_id",
                 variable.name = "parameter",
                 value.name = "fdr")
fdr_long$parameter <- sub("_fdr$", "", fdr_long$parameter)
plot_df <- merge(mean_long, fdr_long, by = c("perturb_id", "parameter"))
plot_df$parameter <- factor(plot_df$parameter, levels = params)
plot_df$star <- ifelse(!is.na(plot_df$fdr) & plot_df$fdr < 0.05, "*", "")

panel_a <- ggplot(plot_df, aes(x = parameter, y = perturb_id, fill = mean)) +
  geom_tile(color = "white") +
  geom_text(aes(label = star), size = 6, vjust = 0.75) +
  scale_fill_gradient2(low = cb_blue, mid = "white", high = cb_orange, midpoint = 1) +
  labs(x = NULL, y = NULL, fill = "ASO/Control", tag = 'A') +
  panel_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Panels B-E: gene tracks with differential promoter usage
image_panel <- function(png_file, tag) {
  img  <- readPNG(png_file)
  grob <- rasterGrob(img, interpolate = TRUE,
                     width = unit(1, "npc"), height = unit(1, "npc"))
  ggplot(geom = "blank") +
    annotation_custom(grob) +
    theme_void() +
    theme(plot.tag = element_text(size = 6, face = "plain"),
          plot.margin = margin(0, 0, 0, 0)) +
    labs(tag = tag)
}

g1 <- readPNG('MDM2.png')
panel_b <- image_panel('MDM2.png',  "B")
g2 <- readPNG('NR2F2.png')
panel_c <- image_panel('NR2F2.png', "C")
g3 <- readPNG('FASTK.png')
panel_d <- image_panel('FASTK.png', "D")
g4 <- readPNG('CFDP1.png')
panel_e <- image_panel('CFDP1.png', "E")


rdbu_ramp <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))
rdbu_at <- function(value, lo = -1, hi = 1, n = 256L) {
  pal <- rdbu_ramp(n)
  pal[pmin(n, pmax(1L, round((value - lo) / (hi - lo) * (n - 1L)) + 1L))]
}

promoter_cats <- data.frame(
  label = c("Significant (p < 0.05), coef < 0",
            "Significant (p < 0.05), coef > 0",
            "Not significant, coef < 0",
            "Not significant, coef > 0"),
  value = c(1, -1, 0.3, -0.3),
  stringsAsFactors = FALSE)
promoter_cats$label <- factor(promoter_cats$label, levels = promoter_cats$label)
promoter_cat_colours <- setNames(rdbu_at(promoter_cats$value), promoter_cats$label)

promoter_legend_plot <- ggplot(promoter_cats, aes(x = label, fill = label)) +
  geom_bar() +
  scale_fill_manual(values = promoter_cat_colours, name = "Promoter") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  panel_theme +
  theme(legend.position = "bottom",
        legend.spacing.x = unit(0.4, "cm"),
        legend.text = element_text(size = 5, margin = margin(r = 0.4, unit = "cm")))

get_legend <- function(p) {
  g   <- ggplotGrob(p)
  idx <- which(vapply(g$grobs, function(x) x$name, character(1)) == "guide-box")
  g$grobs[[idx]]
}
promoter_legend <- get_legend(promoter_legend_plot)

figure3_grob <- arrangeGrob(
    panel_a, panel_b, panel_c, panel_d, panel_e, promoter_legend,
    ncol = 2, nrow = 4,
    layout_matrix = rbind(c(1,1), c(2,3), c(4,5), c(6,6)),
    heights = unit.c(unit(1, "null"), unit(1, "null"),
                     unit(1, "null"), unit(1.2, "cm"))
)
ggsave("Supplementary_Figure3.png", plot = figure3_grob, width = 18, height = 24, units = "cm")