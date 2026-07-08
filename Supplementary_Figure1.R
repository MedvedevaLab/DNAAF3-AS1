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


################################ Main ################################

gtf_file <- "gencode.v49.annotation.gtf"
txdb     <- makeTxDbFromGFF(gtf_file)

# Panel A: Gene tracks

# pyGenomeTracks --tracks plot_supplementary.ini --region chr7:4524000-4536000 -o SupplementaryFigure1A.png --plotWidth 35 --dpi 300
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
  theme_bw() + scale_fill_manual(values = c("DNAAF3-AS1" = "#E41A1C", "DNAAF3" = "#377EB8")) +
  geom_bar(position = "dodge", stat = "identity") +
  theme(panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position = 'bottom',
  text = element_text(size = 6),
  axis.text.x = element_text(size = 5),
  axis.ticks.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.y = element_text(size = 6),
  legend.key.size = unit(0.2, "cm"),
  legend.text = element_text(size = 5),
  plot.title = element_text(size = 6, hjust = 0.5),
  plot.tag = element_text(size = 8, face="plain")) +
  guides(fill=guide_legend(ncol = 2)) +
  labs(title = "GTEx11", fill = "", tag = "B") + coord_flip()

# Panel C: DE genes venn diagram

fantom <- read.csv('DESeq2_genes_ASO_all.tsv', sep = '\t')
aso01 <- dplyr::filter(fantom, perturb_id == 'ASO_G0267577_01')
aso01 <- aso01[aso01$fdr < 0.05, ]

aso03 <- dplyr::filter(fantom, perturb_id == 'ASO_G0267577_03')
aso03 <- aso03[aso03$fdr < 0.05, ]

x <- list(ASO_G0267577_01=aso01$geneSymbol,
          ASO_G0267577_03=aso03$geneSymbol)
 
panel_c <- ggVennDiagram(x, set_size = 2, label_size = 2, edge_size = 0.01) +
  scale_fill_gradient(low="grey90",high = "#E31A1C") +
  scale_x_continuous(expand = expansion(mult = 1)) +
  panel_theme +
  theme(legend.position = "none",
        axis.title.x    = element_blank(),
        axis.text.x     = element_blank(),
        axis.title.y    = element_blank(),
        axis.text.y     = element_blank(),
        plot.margin = margin(0, 0, 0, 0)) +
  labs(tag = "C")


panel_d <- readRDS('decoupler_ASO_G0267577_03_pathway_activities.rds') + panel_theme + labs(tag = "D")
panel_e1 <- readRDS('decoupler_ASO_G0267577_01_p53_targets.rds') + panel_theme +  labs(title = "ASO_01 p53 target genes", tag = "E") + theme(legend.position = "bottom")
panel_e2 <- readRDS('decoupler_ASO_G0267577_03_p53_targets.rds') + panel_theme +  labs(title = "ASO_03 p53 target genes", tag = "")+ theme(legend.position = "bottom")

panel_e <- ggarrange(panel_e1, panel_e2, ncol = 2, common.legend = TRUE, legend = "bottom")

figure1_grob <- arrangeGrob(
    panel_a, panel_b, panel_c, panel_d, panel_e,
    ncol = 2, nrow = 4,
    layout_matrix = rbind(c(1, 3), c(2, 3), c(4, 3), c(5, 5))
)
ggsave("Supplementary_Figure1.png", plot = figure1_grob, width = 18, height = 24, units = "cm")
