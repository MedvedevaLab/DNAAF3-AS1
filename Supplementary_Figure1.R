library(ggplot2)
library(ggpubr)
library(png)
library(gridExtra)
library(reshape2)
library(dplyr)
library(grid)
library(ggVennDiagram)
library(GenomicFeatures)

cb_blue   <- "#0072B2"
cb_orange <- "#D55E00"
cb_green  <- "#009E73"
cb_grey   <- "#999999"

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

gtf_file <- "gencode.v49.annotation.gtf"
txdb     <- makeTxDbFromGFF(gtf_file)

# Panel A: Gene tracks
# pyGenomeTracks --tracks plot_supplementary.ini --region chr7:4524000-4536000 -o SupplementaryFigure1A.png --plotWidth 50 --dpi 300
gene_track_img  <- readPNG("SupplementaryFigure1A.png") # 
gene_track_grob <- rasterGrob(gene_track_img,
    interpolate = TRUE, width = unit(4, "in"), height = unit(1.5, "in"))

panel_a <- ggplot(geom = "blank") +
    annotation_custom(gene_track_grob) +
    theme_bw() +
    panel_theme +
    labs(tag = "A")

# Panel B: Expression in tissues
# https://www.gtexportal.org/home/downloads/adult-gtex%23qtl
df <- read.table('GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct', sep='\t', header=TRUE, skip=2,
                 check.names=FALSE, stringsAsFactors=FALSE)

result <- melt(as.data.frame(df[grepl("ENSG00000267577|ENSG00000167646", df$Name), ]))

result$value <- log10(result$value+1)
panel_b <- ggplot(result, aes(x = variable , y = value , fill = Description )) +
  theme_bw() + scale_fill_manual(values = c("DNAAF3-AS1" = cb_orange, "DNAAF3" = cb_blue)) +
  geom_bar(position = "dodge", stat = "identity") + panel_theme +
  guides(fill=guide_legend(ncol = 2))+
  theme(axis.title.y = element_blank(),
    axis.text.y = element_text(size = 5)) +
  labs(title = "GTEx11", y = "log10(TPM+1)", fill = "", tag = "B") + coord_flip()

# Panel C: control vs knockdown expression at the antisense locus
ENSEMBLE_GENES_INDEX <- 85706
CONTROLS <- c(188, 204, 216, 256, 264, 270, 280, 276)
ASO01    <- c(240, 226)
ASO03    <- c(259, 195)

matr <- read.csv(
  "expression.tpm.gene.nAnT-iCAGE.ASO.human_dermal_fibroblasts.tsv",
  sep = "\t", row.names = 1
)

matr_control <- matr[ENSEMBLE_GENES_INDEX:nrow(matr), CONTROLS]
matr_01      <- matr[ENSEMBLE_GENES_INDEX:nrow(matr), ASO01]
matr_03      <- matr[ENSEMBLE_GENES_INDEX:nrow(matr), ASO03]

gene_map <- c(
  ENSG00000267577 = "DNAAF3-AS1",
  ENSG00000167646 = "DNAAF3",
  ENSG00000129990 = "SYT5"
)
genes_of_interest <- names(gene_map)

build_long <- function(m, group, aso) {
  sub <- m[genes_of_interest, , drop = FALSE]
  long <- melt(as.matrix(sub), varnames = c("gene_id", "sample"),
               value.name = "count")
  long$group <- group
  long$aso   <- aso
  long
}

antisense_df <- rbind(
  build_long(matr_control, "Control",   "Control"),
  build_long(matr_01,      "Knockdown", "ASO01"),
  build_long(matr_03,      "Knockdown", "ASO03")
)

antisense_df$count <- as.numeric(as.character(antisense_df$count))
antisense_df$gene  <- factor(gene_map[as.character(antisense_df$gene_id)],
                             levels = gene_map)
antisense_df$group <- factor(antisense_df$group,
                             levels = c("Control", "Knockdown"))
antisense_df$aso   <- factor(antisense_df$aso,
                             levels = c("Control", "ASO01", "ASO03"))

panel_c <- ggplot(antisense_df, aes(x = group, y = count)) +
  geom_boxplot(aes(fill = group), outlier.shape = NA, width = 0.6,
               alpha = 0.6) +
  geom_jitter(aes(colour = aso), width = 0.15, size = 1) +
  facet_wrap(~ gene, scales = "free_y") +
  scale_fill_manual(values = c(Control = "grey70", Knockdown = "grey40"),
                    guide = "none") +
  scale_colour_manual(values = c(Control = cb_grey,
                                 ASO01   = cb_blue,
                                 ASO03   = cb_orange),
                      name = "") +
  stat_compare_means(comparisons = list(c("Control", "Knockdown")),
                     method = "wilcox.test", size = 1.5) +
  labs(x = NULL, y = "CAGE expression (TPM)",
       title = " ", tag = "C") +
  theme_bw() +
  panel_theme +
  theme(axis.text.x = element_text(size = 5))

# Mean expression per gene within each ASO group
mean_expr_by_gene_aso <- antisense_df %>%
  dplyr::group_by(gene, aso) %>%
  dplyr::summarise(mean_count = mean(count), .groups = "drop")
print(as.data.frame(mean_expr_by_gene_aso))

fold_01 <- mean(antisense_df[(antisense_df$gene == "DNAAF3") & (antisense_df$aso == "ASO01"), "count"]) / mean(antisense_df[(antisense_df$gene == "DNAAF3") & (antisense_df$aso == "Control"), "count"])
fold_03 <- mean(antisense_df[(antisense_df$gene == "DNAAF3") & (antisense_df$aso == "ASO03"), "count"]) / mean(antisense_df[(antisense_df$gene == "DNAAF3") & (antisense_df$aso == "Control"), "count"])

# Panel D: DE genes venn diagram
fantom <- read.csv('DESeq2_genes_ASO_all.tsv', sep = '\t')
aso01 <- dplyr::filter(fantom, perturb_id == 'ASO_G0267577_01')
aso01 <- aso01[aso01$fdr < 0.05, ]

aso03 <- dplyr::filter(fantom, perturb_id == 'ASO_G0267577_03')
aso03 <- aso03[aso03$fdr < 0.05, ]

x <- list(ASO_G0267577_01=aso01$geneSymbol,
          ASO_G0267577_03=aso03$geneSymbol)
 
panel_d <- ggVennDiagram(x, set_size = 2, label_size = 2, edge_size = 0.01) +
  scale_fill_gradient(low="grey90",high = cb_orange) +
  scale_x_continuous(expand = expansion(mult = 1)) +
  panel_theme +
  theme(legend.position = "none",
        axis.title.x    = element_blank(),
        axis.text.x     = element_blank(),
        axis.title.y    = element_blank(),
        axis.text.y     = element_blank(),
        plot.margin = margin(0, 0, 0, 0)) +
  labs(tag = "D")

figure1_grob <- arrangeGrob(
    panel_a, panel_b, panel_c, panel_d,
    ncol = 2, nrow = 3,
    layout_matrix = rbind(c(1, 2), c(3, 2), c(4, 2))
)

ggsave("Supplementary_Figure1.png", plot = figure1_grob, width = 18, height = 18, units = "cm")
