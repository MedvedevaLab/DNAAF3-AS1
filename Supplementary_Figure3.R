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

library(biomaRt)
library(scales)


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


# Single-cell GSE112013 dataset
seurat_object <- readRDS('testis.rds')
Idents(object = seurat_object) <- seurat_object@meta.data$seurat_clusters

p1 <- DotPlot(seurat_object, 
 features =   c('VIM', 'CD14', 'VWF', 'ACTA2', 'DLK1', 'DAZL', 'MAGEA4', 'UTF1', 'ID4', 'FGFR3', 'KIT', 'DMRT1', 
 	'DMRTB1', 'STRA8', 'SYCP3', 'SPO11', 'MLH3', 'ZPBP', 'TNP1', 'PRM2'))


seurat_object_promoter <- readRDS('testis_promoter.rds')
seurat_object_genebody <- readRDS('testis_genebody.rds')

p2 <- FeaturePlot(seurat_object_promoter, paste0("add_module_score", 1), pt.size = 0.5, reduction = "umap") + 
  scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral"))) +
  panel_theme +
  theme(text = element_text(size = 6)) +
  labs(title = "", tag = "B")


p3 <- FeaturePlot(seurat_object_genebody, paste0("add_module_score", 1), pt.size = 0.5, reduction = "umap") + 
  scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral"))) +
  panel_theme +
  theme(text = element_text(size = 6)) +
  labs(title = "", tag = "")


p4 <- FeaturePlot(seurat_object_genebody, paste0("add_module_score", 2), pt.size = 0.5, reduction = "umap") + 
  scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral"))) +
  panel_theme +
  theme(text = element_text(size = 6)) +
  labs(title = "", tag = "")

#################################

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
p5 <- ggplot(df, aes(x = chromosome, y = experiment, fill = density)) +
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
  labs(title = "", x = 'Chromosome', y = 'Experiment', fill = 'log10(Contacts Density)', tag = "C")


#################################

df1 <- read.csv("rnachrom1.csv",sep = ',')
df1$experiment <- 'iMARGI1'

df2 <- read.csv("rnachrom2.csv")
df2$experiment <- 'iMARGI2'

df3 <- read.csv("rnachrom3.csv")
df3$experiment <- 'iMARGI3'

df4 <- read.csv("rnachrom4.csv")
df4$experiment <- 'iMARGI4'

df5 <- read.csv("rnachrom5.csv")
df5$experiment <- 'iMARGI5'

df7 <- read.csv("rnachrom7.csv")
df7$experiment <- 'GRIDseq1'

df8 <- read.csv("rnachrom8.csv")
df8$experiment <- 'GRIDseq2'

df10 <- read.csv("rnachrom10.csv")
df10$experiment <- 'Red-C'

df_all <- do.call(rbind, list(df1, df2, df3, df4, df5, df7, df8, df10))

# Target genes intersection

gtf_name <- "gencode.v49.annotation.gtf"
txdb <- makeTxDbFromGFF(gtf_name)
genes_gr <- genes(txdb)   # returns a GRanges with gene_id in metadata

gr_regions <- GRanges(
  seqnames = df_all$Chromosome,
  ranges   = IRanges(start = df_all$Gene.start,
                     end   = df_all$Gene.end)
)

hits <- findOverlaps(gr_regions, genes_gr, ignore.strand = TRUE)

results_genes <- data.frame(
  region_chrom = as.character(seqnames(gr_regions))[queryHits(hits)],
  region_start = start(gr_regions)[queryHits(hits)],
  region_end   = end(gr_regions)[queryHits(hits)],
  gene_id      = genes_gr$gene_id[subjectHits(hits)],
  gene_chrom   = as.character(seqnames(genes_gr))[subjectHits(hits)],
  gene_start   = start(genes_gr)[subjectHits(hits)],
  gene_end     = end(genes_gr)[subjectHits(hits)],
  gene_strand  = as.character(strand(genes_gr))[subjectHits(hits)],
  stringsAsFactors = FALSE
)

target_genes <- unique(results_genes$gene_id)   
target_genes <- str_split_fixed(target_genes, "\\.", 2)[, 1]


himorna_peaks <- read.csv("ENSG00000267577.tsv", sep = "\t")
himorna_peaks <- dplyr::filter(himorna_peaks, hm == "H3K36me3")
himorna_contacts <- dplyr::filter(himorna_peaks, gene_ens %in% target_genes)

fantom <- read.csv('DESeq2_genes_ASO_all.tsv', sep = '\t')
aso01 <- dplyr::filter(fantom, perturb_id == 'ASO_G0267577_03')
aso01 <- aso01[aso01$fdr < 0.05, ]

aso01_contacts <- dplyr::filter(aso01, geneID %in% target_genes)
aso01_contacts$in_himorna <- aso01_contacts$geneID  %in% himorna_contacts$gene_ens
aso01_contacts$log_padj <-  -log(aso01_contacts$fdr,10)

colors <- c('TRUE' = "#E31A1C", 'FALSE' = "#000000")

p6 <- ggplot(aso01_contacts, aes(x = log2FC,
  y = log_padj,
  label = ifelse(in_himorna, geneSymbol, NA),
  color = in_himorna)) +
  geom_point(shape = 1) + 
  geom_text_repel( color = 'black', alpha = 1,  hjust = 0, vjust = 0, size = 4, max.overlaps = 15) +
  scale_color_manual(values = colors, name = "") +
  theme_bw() +
  panel_theme +
  theme(text = element_text(size = 4)) +
    labs(title = "ASO 03", x = "", y = "", tag = "D")

# pyGenomeTracks --tracks plot_smyd3.ini --region chr1:245739000-246520000 -o Figure3A.png--plotWidth 50 --dpi 130
gene_track_img  <- readPNG("Figure3A.png")
gene_track_grob <- rasterGrob(gene_track_img,
    interpolate = TRUE, width = unit(4, "in"), height = unit(1.5, "in"))

p7 <- ggplot(geom = "blank") +
    annotation_custom(gene_track_grob) +
    theme_bw() +
    panel_theme +
    labs(tag = "")

p <- ggarrange(p2, p3, p4, ncol = 3, nrow = 1, common.legend = FALSE, legend="bottom")

figure3_grob <- arrangeGrob(
  p1, p, p5, p6, p7,
  ncol = 2, nrow = 3,
  layout_matrix = rbind(c(1, 1), c(2, 2), c(3, 4), c(5, 5))
)

ggsave("Supplementary_Figure3.png", plot = figure3_grob, width = 16, height = 26, units = "cm")
