library(ggplot2)
library(dplyr)
library(png)
library(grid)
library(gridExtra)

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
    plot.title         = element_text(size = 6, face = "plain"),
    plot.tag           = element_text(size = 6, face = "plain")
)


PREFIXES <- c("ASO_G0267577_01_signal", "ASO_G0267577_03_signal")
CLUSTERS <- c("cluster_7", "cluster_8")

read_cluster_panel <- function(prefix, cluster, tag) {
    img  <- readPNG(paste0(prefix, ".", cluster, ".png"))
    grob <- rasterGrob(img, interpolate = TRUE)
    ggplot(geom = "blank") +
        annotation_custom(grob) +
        theme_bw() +
        panel_theme +
        labs(title = paste(prefix, cluster), tag = tag)
}

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
        scale_fill_viridis_d(option = "viridis") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        theme_bw() +
        panel_theme +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none") +
        labs(x = "Cluster", y = "Number of regions",
             title = paste("Metagene cluster sizes,", prefix), tag = tag)
}

read_png_panel <- function(png_file, tag, title = "") {
    img  <- readPNG(png_file)
    grob <- rasterGrob(img, interpolate = TRUE)
    ggplot(geom = "blank") +
        annotation_custom(grob) +
        theme_bw() +
        panel_theme +
        labs(title = title, tag = tag)
}

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

panel_metascape <- read_png_panel("chipseq_metascape.png", tags[i],
                                  "ChIP-seq peaks Metascape enrichment")
panels <- c(panels, list(panel_metascape))

n_prefix     <- length(PREFIXES)
layout_matrix <- rbind(
    matrix(seq_len(3 * n_prefix), ncol = 3, byrow = TRUE),
    rep(length(panels), 3)
)
g <- arrangeGrob(grobs = panels, layout_matrix = layout_matrix)
ggsave("Supplementary_Figure6.png", plot = g,
       width = 18, height = 18, units = "cm")
