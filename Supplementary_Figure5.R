library(ggplot2)
library(dplyr)
library(png)
library(grid)
library(gridExtra)

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
    plot.title         = element_text(size = 6, face = "plain"),
    plot.tag           = element_text(size = 6, face = "plain")
)

# Bigwigs whose metagene clusters we show. Prefixes follow peak_calling.py:
# bw name with spaces -> "_" and no ".bw".
PREFIXES <- c("ASO_G0267577_01_signal", "ASO_G0267577_03_signal")
# Only these clusters are shown as metagene panels.
CLUSTERS <- c("cluster_7", "cluster_8")

# ── One cluster metagene panel: read {prefix}.{cluster}.png as a raster ────────
# (same approach as panel_a in Figure1.R)
read_cluster_panel <- function(prefix, cluster, tag) {
    img  <- readPNG(paste0(prefix, ".", cluster, ".png"))
    grob <- rasterGrob(img, interpolate = TRUE)
    ggplot(geom = "blank") +
        annotation_custom(grob) +
        theme_bw() +
        panel_theme +
        labs(title = paste(prefix, cluster), tag = tag)
}

# ── Cluster-size barplot panel (as in cluster_sizes_barplot.R) ─────────────────
cluster_sizes_panel <- function(prefix, tag) {
    bed <- read.table(paste0(prefix, ".clusters.bed"), sep = "\t",
                      header = FALSE, comment.char = "#",
                      stringsAsFactors = FALSE)
    sizes <- data.frame(cluster = bed[[ncol(bed)]], stringsAsFactors = FALSE) %>%
        count(cluster, name = "size") %>%
        mutate(cluster = factor(cluster,
            levels = unique(cluster[order(as.integer(gsub("\\D", "", cluster)))])))

    ggplot(sizes, aes(x = cluster, y = size, fill = cluster)) +
        geom_col() +
        geom_text(aes(label = size), vjust = -0.3, size = 2) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        theme_bw() +
        panel_theme +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none") +
        labs(x = "Cluster", y = "Number of regions",
             title = paste("Metagene cluster sizes,", prefix), tag = tag)
}

# ── Assemble: for each prefix, cluster_7 / cluster_8 metagenes + size barplot ──
panels <- list()
tags   <- LETTERS
i      <- 1
for (prefix in PREFIXES) {
    panels <- c(panels, list(cluster_sizes_panel(prefix, tags[i])))
    i <- i + 1
    for (cluster in CLUSTERS) {
        panels <- c(panels, list(read_cluster_panel(prefix, cluster, tags[i])))
        i <- i + 1
    }
}

# each prefix contributes one row: size barplot, cluster_7, cluster_8
g <- arrangeGrob(grobs = panels, ncol = 3, nrow = length(PREFIXES))
ggsave("Supplementary_Figure5.png", plot = g,
       width = 18, height = 12, units = "cm")
