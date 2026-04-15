library(ggplot2)
library(ggrepel)
library(ggpubr)
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
    plot.title         = element_text(size = 6, face="plain"),
    plot.tag           = element_text(size = 6, face = "plain"),
    legend.spacing.y   = unit(2, "cm")
)

# A. Alignment rate
# B. Alignment rate
# C. Fraction of reads in peaks
alignment_stat <- read.csv('alignment_stat.csv', sep = '\t')
alignment_stat$Sample_name <- c('ASO_01_signal', 'ASO_03_signal', 'NCA_01_signal', 'NCA_02_signal', 
	'ASO_01_input', 'ASO_03_input', 'NCA_01_input', 'NCA_02_input')
alignment_stat$Sample_name <- factor(alignment_stat$Sample_name, levels = c('ASO_01_signal', 'ASO_03_signal', 'NCA_01_signal', 'NCA_02_signal', 
	'ASO_01_input', 'ASO_03_input', 'NCA_01_input', 'NCA_02_input'))
p1 <- ggplot(alignment_stat, aes(x = Sample_name, y = Alignment_rate)) +
    theme_bw() +
    geom_bar(position = "dodge", stat = "identity") +
    scale_fill_brewer(palette = "Set1", guide = guide_legend(reverse = TRUE), drop = FALSE) +
    panel_theme +
    theme(axis.text.x = element_text(size = 6, angle = 90)) +
    guides(fill = guide_legend(ncol = 3)) +
    labs(title = "", x = "Experiment", y = "Alignment rate", tag = "A")

p2 <- ggplot(alignment_stat, aes(x = Sample_name, y = Deduplication_rate)) +
    theme_bw() +
    geom_bar(position = "dodge", stat = "identity") +
    scale_fill_brewer(palette = "Set1", guide = guide_legend(reverse = TRUE), drop = FALSE) +
    panel_theme +
    theme(axis.text.x = element_text(size = 6, angle = 90)) +
    guides(fill = guide_legend(ncol = 3)) +
    labs(title = "", x = "Experiment", y = "Deduplication rate", tag = "B")


frip_stat <- read.csv('frip_stat.csv', sep = '\t')

p3 <- ggplot(frip_stat, aes(x = Sample, y = FRIP)) +
    theme_bw() +
    geom_bar(position = "dodge", stat = "identity") +
    scale_fill_brewer(palette = "Set1", guide = guide_legend(reverse = TRUE), drop = FALSE) +
    panel_theme +
    theme(axis.text.x = element_text(size = 6, angle = 90)) +
    guides(fill = guide_legend(ncol = 3)) +
    labs(title = "", x = "Experiment", y = "FRIP", tag = "C")


p <- ggarrange(p1, p2, p3, ncol = 3, nrow = 1, common.legend = F)

load("C-9_.spp.Rdata")

# vertical lines data
vlines_peak    <- data.frame(x = cc.peak$x,                 type = "CC Peak")
vlines_phantom <- data.frame(x = crosscorr$phantom.cc$x,    type = "Phantom Peak")

# build plot
p41 <- ggplot(crosscorr$cross.correlation, aes(x = x, y = y)) +
  geom_line() +
  geom_vline(data = vlines_peak,    aes(xintercept = x), linetype = "dashed", color = "red") +
  geom_vline(data = vlines_phantom, aes(xintercept = x), linetype = "dashed", color = "blue") +
  labs(
    x        = sprintf("strand-shift (%s)", paste(cc.peak$x, collapse = ",")),
    y        = "cross-correlation",
    title    = "ASO_01",
    subtitle = sprintf("NSC=%g, RSC=%g, Qtag=%d",
                       crosscorr$phantom.coeff,
                       crosscorr$rel.phantom.coeff,
                       crosscorr$phantom.quality.tag), tag = "D") +
  theme_classic(base_size = 10) +
  theme(
    plot.title    = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title    = element_text(size  = 9),
    plot.margin   = margin(t = 2, r = 0.5, b = 4, l = 3.5, unit = "mm")
  )


load("C-10.spp.Rdata")

# vertical lines data
vlines_peak    <- data.frame(x = cc.peak$x,                 type = "CC Peak")
vlines_phantom <- data.frame(x = crosscorr$phantom.cc$x,    type = "Phantom Peak")

# build plot
p42 <- ggplot(crosscorr$cross.correlation, aes(x = x, y = y)) +
  geom_line() +
  geom_vline(data = vlines_peak,    aes(xintercept = x), linetype = "dashed", color = "red") +
  geom_vline(data = vlines_phantom, aes(xintercept = x), linetype = "dashed", color = "blue") +
  labs(
    x        = sprintf("strand-shift (%s)", paste(cc.peak$x, collapse = ",")),
    y        = "cross-correlation",
    title    = "ASO_03",
    subtitle = sprintf("NSC=%g, RSC=%g, Qtag=%d",
                       crosscorr$phantom.coeff,
                       crosscorr$rel.phantom.coeff,
                       crosscorr$phantom.quality.tag)
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title    = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title    = element_text(size  = 9),
    plot.margin   = margin(t = 2, r = 0.5, b = 4, l = 3.5, unit = "mm")
  )

load("C-19.spp.Rdata")

# vertical lines data
vlines_peak    <- data.frame(x = cc.peak$x,                 type = "CC Peak")
vlines_phantom <- data.frame(x = crosscorr$phantom.cc$x,    type = "Phantom Peak")

# build plot
p43 <- ggplot(crosscorr$cross.correlation, aes(x = x, y = y)) +
  geom_line() +
  geom_vline(data = vlines_peak,    aes(xintercept = x), linetype = "dashed", color = "red") +
  geom_vline(data = vlines_phantom, aes(xintercept = x), linetype = "dashed", color = "blue") +
  labs(
    x        = sprintf("strand-shift (%s)", paste(cc.peak$x, collapse = ",")),
    y        = "cross-correlation",
    title    = "NCA_01",
    subtitle = sprintf("NSC=%g, RSC=%g, Qtag=%d",
                       crosscorr$phantom.coeff,
                       crosscorr$rel.phantom.coeff,
                       crosscorr$phantom.quality.tag)
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title    = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title    = element_text(size  = 9),
    plot.margin   = margin(t = 2, r = 0.5, b = 4, l = 3.5, unit = "mm")
  )


load("C-20.spp.Rdata")

# vertical lines data
vlines_peak    <- data.frame(x = cc.peak$x,                 type = "CC Peak")
vlines_phantom <- data.frame(x = crosscorr$phantom.cc$x,    type = "Phantom Peak")

# build plot
p44 <- ggplot(crosscorr$cross.correlation, aes(x = x, y = y)) +
  geom_line() +
  geom_vline(data = vlines_peak,    aes(xintercept = x), linetype = "dashed", color = "red") +
  geom_vline(data = vlines_phantom, aes(xintercept = x), linetype = "dashed", color = "blue") +
  labs(
    x        = sprintf("strand-shift (%s)", paste(cc.peak$x, collapse = ",")),
    y        = "cross-correlation",
    title    = "NCA_02",
    subtitle = sprintf("NSC=%g, RSC=%g, Qtag=%d",
                       crosscorr$phantom.coeff,
                       crosscorr$rel.phantom.coeff,
                       crosscorr$phantom.quality.tag)
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title    = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title    = element_text(size  = 9),
    plot.margin   = margin(t = 2, r = 0.5, b = 4, l = 3.5, unit = "mm")
  )



g <- arrangeGrob(p, p41, p42, p43, p44, ncol = 2, nrow = 3, layout_matrix= rbind(c(1,1), c(2,3), c(4,5)))
ggsave(paste0("Supplementary_Figure2.png"), plot = g, width = 18, height = 25, units = "cm")

