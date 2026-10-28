library(ggplot2)
library(gridExtra)
library(ggpubr)

panel_theme <- theme(
    strip.background = element_blank(),
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

# Panel A:
panel_a <- readRDS('decoupler_p53_activity_all_ASOs.rds') + 
	panel_theme + 
	theme(axis.text.y = element_blank(),
	    axis.title.y = element_text(size = 6)) + 
	labs(x = "p53 pathway activity (mlm score)", y = "ASO",
         fill = "",
         title = "PROGENy p53 pathway activity across ASOs",
         tag = "A")

# Panel B: decoupler PROGENy pathway activity ASO_G0267577_03
panel_b <- readRDS('decoupler_ASO_G0267577_03_pathway_activities.rds') + panel_theme + labs(title = "ASO_03 PROGENy pathway activities", tag = "B")

# Panel C: decoupler p53 genes 
panel_c1 <- readRDS('decoupler_ASO_G0267577_01_p53_targets.rds') + panel_theme +  labs(title = "ASO_01 p53 target genes", tag = "C") + theme(legend.position = "bottom")
panel_c2 <- readRDS('decoupler_ASO_G0267577_03_p53_targets.rds') + panel_theme +  labs(title = "ASO_03 p53 target genes", tag = "") + theme(legend.position = "bottom")
panel_c <- ggarrange(panel_c1, panel_c2, ncol = 2, common.legend = TRUE, legend = "bottom")

figure1_grob <- arrangeGrob(
    panel_a, panel_b, panel_c,
    ncol = 2, nrow = 2,
    layout_matrix = rbind(c(1, 2), c(3, 3))
)

ggsave("Supplementary_Figure2.png", plot = figure1_grob, width = 18, height = 18, units = "cm")
