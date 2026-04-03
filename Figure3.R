library(GenomicRanges)
library(GenomicFeatures)
library(stringr)
library(ggplot2)
library(png)
library(ggrepel)
library(dplyr)
library(ggplot2)
library(Seurat)
library(gridExtra)
library(RColorBrewer)

library(ggplot2)
library(ggpubr)
library(ChIPseeker)
library(png)
library(gridExtra)
library(reshape2)
library(AnnotationHub)
library(biomaRt)
library(scales)
library(RColorBrewer)
library(stringr)
library(GenomeInfoDb)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPpeakAnno)
library(grid)
library(ggsankey)
library(tidyr)
library(clusterProfiler)

panel_theme <- theme(
    panel.border       = element_blank(),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    axis.title.x       = element_text(size = 6),
    axis.text.x        = element_text(size = 6, hjust = 0.5),
    legend.title       = element_text(size = 6),
    legend.key.size    = unit(0.2, "cm"),
    legend.text        = element_text(size = 5),
    axis.ticks.x       = element_blank(),
    axis.text.y        = element_text(size = 6),
    axis.title.y       = element_blank(),
    plot.title         = element_text(size = 6),
    plot.tag           = element_text(size = 6, face = "plain"),
    legend.spacing.y   = unit(2, "cm")
)


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
panel_a <- ggplot(df, aes(x = chromosome, y = experiment , fill=density)) +
  theme_bw() +
  geom_tile() +
  scale_fill_distiller(palette = "Reds", direction =1, na.value = 'white') +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 6),
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = 6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 6, face="plain"),
        axis.title = element_text(size = 6),
        plot.tag = element_text(size = 8, face="plain")) +
  labs(title = "", x = 'Chromosome', y = 'Experiment', fill = 'log10(Contacts Density)', tag = "A")


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
# counts <- as.data.frame(table(df_all[,'Target.gene']))
# counts <- as.data.frame(table(df_all[,'Chromosome']))

############### Target genes intersection

gtf_name <- "gencode.v49.annotation.gtf"
txdb <- makeTxDbFromGFF(gtf_name)
genes_gr <- genes(txdb)   # returns a GRanges with gene_id in metadata

# gr_regions <- GRanges(
#   seqnames = result$seqnames,
#   ranges   = IRanges(start = result$start - 10000,
#                      end   = result$end   + 10000)
# )
gr_regions <- GRanges(
  seqnames = df_all$Chromosome,
  ranges   = IRanges(start = df_all$Gene.start,
                     end   = df_all$Gene.end)
)

# ── 4. Find overlaps ─────────────────────────────────────────────────────────
hits <- findOverlaps(gr_regions, genes_gr, ignore.strand = TRUE)

# ── 5. Build results table ───────────────────────────────────────────────────
results_genes <- data.frame(
  region_chrom = as.character(seqnames(gr_regions))[queryHits(hits)],
  region_start = start(gr_regions)[queryHits(hits)],   # extended start
  region_end   = end(gr_regions)[queryHits(hits)],     # extended end
  gene_id      = genes_gr$gene_id[subjectHits(hits)],
  gene_chrom   = as.character(seqnames(genes_gr))[subjectHits(hits)],
  gene_start   = start(genes_gr)[subjectHits(hits)],
  gene_end     = end(genes_gr)[subjectHits(hits)],
  gene_strand  = as.character(strand(genes_gr))[subjectHits(hits)],
  stringsAsFactors = FALSE
)

target_genes <- unique(results_genes$gene_id)   
target_genes <- str_split_fixed(target_genes, "\\.", 2)[, 1]

###############################################################


himorna_peaks <- read.csv("ENSG00000267577.tsv", sep = "\t")
# himorna_peaks <- dplyr::filter(himorna_peaks, hm == "H3K36me3")

himorna_contacts <- dplyr::filter(himorna_peaks, gene_ens %in% target_genes)

fantom <- read.csv('DESeq2_genes_ASO_all.tsv', sep = '\t')
aso01 <- dplyr::filter(fantom, perturb_id == 'ASO_G0267577_01')
aso01 <- aso01[aso01$fdr < 0.05, ]


aso03 <- dplyr::filter(fantom, perturb_id == 'ASO_G0267577_03')
aso03 <- aso03[aso03$fdr < 0.05, ]

intersect(aso01$geneID, himorna_peaks$gene_ens)
intersect(aso03$geneID, himorna_peaks$gene_ens)

intersect(aso01$geneSymbol, df_all[,'Target.gene'])
intersect(aso03$geneSymbol, df_all[,'Target.gene'])


aso01_contacts <- dplyr::filter(aso01, geneID %in% target_genes)
aso03_contacts <- dplyr::filter(aso03, geneID %in% target_genes)

aso01_contacts$in_himorna <- aso01_contacts$geneID  %in% himorna_contacts$gene_ens
aso01_contacts$log_padj <-  -log(aso01_contacts$fdr,10)

colors <- c('TRUE' = "#E31A1C", 'FALSE' = "#000000")

panel_b1 <- ggplot(aso01_contacts, aes(x = log2FC,
	y = log_padj,
	label = ifelse(in_himorna, geneSymbol, NA),
	color = in_himorna)) +
  geom_point(shape = 1) + 
  geom_text_repel( color = 'black', alpha = 1,  hjust = 0, vjust = 0, size = 4, max.overlaps = 15) +
  scale_color_manual(values = colors, name = "") +
  theme_bw() +
  theme(panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    text = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x=element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position="bottom") +
    labs(title = "ASO 01", x = "", y = "", tag = "B")

aso03_contacts$in_himorna <- aso03_contacts$geneID  %in% himorna_contacts$gene_ens
aso03_contacts$log_padj <-  -log(aso03_contacts$fdr,10)

colors <- c('TRUE' = "#E31A1C", 'FALSE' = "#000000")

panel_b2 <- ggplot(aso03_contacts, aes(x = log2FC,
	y = log_padj,
	label = ifelse(in_himorna, geneSymbol, NA),
	color = in_himorna)) +
  geom_point(shape = 1) + 
  geom_text_repel( color = 'black', alpha = 1,  hjust = 0, vjust = 0, size = 4, max.overlaps = 15) +
  scale_color_manual(values = colors, name = "") +
  theme_bw() +
  theme(panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    text = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x=element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position="bottom") +
    labs(title = "ASO 03 ", x = "", y = "", tag = "")




# Annotate genomic peaks from a data frame of coordinates
get_annotation <- function(peak_data, correlation_label, txdb) {
    coords <- peak_data[, c("chrom", "start", "end")]
    coords <- cbind(coords, "*")
    colnames(coords) <- c("seqnames", "start", "end", "strand")
    coords$start <- as.numeric(coords$start)
    coords$end   <- as.numeric(coords$end)
    peak_granges <- toGRanges(coords, format = "BED")
    peak_annotation <- annotatePeak(peak_granges, tssRegion = c(-1000, 100), TxDb = txdb, overlap = "all")
    return(peak_annotation)
}

gtf_file <- "gencode.v49.annotation.gtf"
txdb     <- makeTxDbFromGFF(gtf_file)
library(biomaRt)

mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensembl <- mart


library(dplyr)
library(ggplot2)
library(Seurat)

seurat_object2 <- readRDS('testis.rds')


seurat_object2@meta.data$celltype <- NA
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 15] = "Macrophage"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 18] = "Macrophage"

seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 0] = "Sperm3"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 1] = "Myoid cells" #
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 2] = "Sperm2"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 3] = "Early primary S'cytes"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 4] = "Endothelial cells"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 5] = "Round S'tids"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 6] = "Sperm2"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 7] = "Differentiating S'gonia"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 8] = "Myoid cells"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 9] = "Elongated S'tids"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 10] = "Sperm2"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 11] = "SSCs"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 12] = "Late primary S'cytes"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 13] = "Sperm1"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 14] = "Early primary S'cytes"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 16] = "Late primary S'cytes"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 17] = "Elongated S'tids"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 19] = "Late primary S'cytes"
seurat_object2@meta.data$celltype[seurat_object2@meta.data$seurat_clusters == 20] = "Sertoli cells"



Idents(object = seurat_object2) <- seurat_object2@meta.data$celltype
seurat_object2$celltype <- factor(seurat_object2$celltype, 
                        levels = c("SSCs", "Differentiating S'gonia", "Early primary S'cytes", "Late primary S'cytes", "Round S'tids", "Elongated S'tids", "Sperm1", "Sperm2", "Sperm3",
"Macrophage", "Endothelial cells", "Myoid cells", "Sertoli cells"))

# UMAPs with annotation
panel_c <- DimPlot(seurat_object2, reduction = "umap", group.by = "celltype", label = FALSE, label.size = 6, repel = TRUE, raster=FALSE) + 
theme(legend.position = "bottom", legend.text=element_text(size=8))+
  labs(tag = "C")

panel_d1 <- VlnPlot(
  seurat_object2,
  features = c('CTD-2587H24.5'),
  alpha = 0, log = TRUE,
  sort = 'decreasing') + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 6),
        strip.background = element_blank(),
        legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6),
        plot.title = element_text(size = 6, face="plain"),
        plot.tag = element_text(size = 8, face="plain")) +
  labs(title = "DNAAF3-AS1", tag = "D")

panel_d2 <- VlnPlot(
  seurat_object2,
  features = c('DNAAF3'),
  alpha = 0, log = TRUE,
  sort = 'decreasing') + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 6),
        strip.background = element_blank(),
        legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6),
        plot.title = element_text(size = 6, face="plain"),
        plot.tag = element_text(size = 8, face="plain")) +
  labs(title = "DNAAF3", tag = "")

himorna_peaks <- read.csv("ENSG00000267577.tsv", sep = "\t")
himorna_peaks <- dplyr::filter(himorna_peaks, hm == "H3K36me3")

himorna_peaks_pos_cor <- himorna_peaks[himorna_peaks$corr > 0, ]
himorna_peaks_neg_cor <- himorna_peaks[himorna_peaks$corr < 0, ]
rownames(himorna_peaks_pos_cor) <- himorna_peaks_pos_cor$X
rownames(himorna_peaks_neg_cor) <- himorna_peaks_neg_cor$X


peak_anno_pos_cor <- get_annotation(himorna_peaks_pos_cor, "rho > 0", txdb)
peak_anno_neg_cor <- get_annotation(himorna_peaks_neg_cor, "rho < 0", txdb)

peak_anno_pos_cor_promoter <- dplyr::filter(as.data.frame(peak_anno_pos_cor@anno), annotation == "Promoter")
peak_anno_pos_cor_distal <- dplyr::filter(as.data.frame(peak_anno_pos_cor@anno), annotation == "Distal Intergenic")
peak_anno_pos_cor_other <- dplyr::filter(as.data.frame(peak_anno_pos_cor@anno), !(annotation %in% c("Promoter", "Distal Intergenic")))
peak_anno_neg_cor_promoter <- dplyr::filter(as.data.frame(peak_anno_neg_cor@anno), annotation == "Promoter")
peak_anno_neg_cor_distal <- dplyr::filter(as.data.frame(peak_anno_neg_cor@anno), annotation == "Distal Intergenic")
peak_anno_neg_cor_other <- dplyr::filter(as.data.frame(peak_anno_neg_cor@anno), !(annotation %in% c("Promoter", "Distal Intergenic")))

library(stringr)

peak_anno_pos_cor_promoter$geneId <-  str_split_fixed(peak_anno_pos_cor_promoter$geneId , "\\.", 2)[, 1]
peak_anno_pos_cor_distal$geneId <-  str_split_fixed(peak_anno_pos_cor_distal$geneId , "\\.", 2)[, 1]
peak_anno_pos_cor_other$geneId <-  str_split_fixed(peak_anno_pos_cor_other$geneId , "\\.", 2)[, 1]
peak_anno_neg_cor_promoter$geneId <-  str_split_fixed(peak_anno_neg_cor_promoter$geneId , "\\.", 2)[, 1]
peak_anno_neg_cor_distal$geneId <-  str_split_fixed(peak_anno_neg_cor_distal$geneId , "\\.", 2)[, 1]
peak_anno_neg_cor_other$geneId <-  str_split_fixed(peak_anno_neg_cor_other$geneId , "\\.", 2)[, 1]

promoter_names_up = getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), filters = 'ensembl_gene_id', 
    values = peak_anno_pos_cor_promoter$geneId, mart = ensembl, useCache = TRUE)

distal_names_up = getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), filters = 'ensembl_gene_id', 
    values = peak_anno_pos_cor_distal$geneId, mart = ensembl, useCache = TRUE)

other_names_up = getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), filters = 'ensembl_gene_id', 
    values = peak_anno_pos_cor_other$geneId, mart = ensembl, useCache = TRUE)


promoter_names_down = getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), filters = 'ensembl_gene_id', 
    values = peak_anno_neg_cor_promoter$geneId, mart = ensembl, useCache = TRUE)

distal_names_down = getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), filters = 'ensembl_gene_id', 
    values = peak_anno_neg_cor_distal$geneId, mart = ensembl, useCache = TRUE)

other_names_down = getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), filters = 'ensembl_gene_id', 
    values = peak_anno_neg_cor_other$geneId, mart = ensembl, useCache = TRUE)


markers <- list()
markers$pos <- unique(promoter_names_up$hgnc_symbol)
markers$neg <- unique(promoter_names_down$hgnc_symbol)

seurat_object2 <- AddModuleScore(seurat_object2, features = markers, name = "add_module_score")

panel_e1 <- FeaturePlot(seurat_object2, paste0("add_module_score",1), pt.size = 0.5, reduction = "umap") + scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral")))  +
 theme(axis.title.x=element_blank(),#axis.text.x=element_blank(),
        text = element_text(size = 6), #, face = "bold"
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        plot.title = element_text(size = 6, face="plain"),
        axis.title = element_text(size = 6),
        plot.tag = element_text(size = 8, face="plain")) 

panel_e2 <- FeaturePlot(seurat_object2, paste0("add_module_score",2), pt.size = 0.5, reduction = "umap") + scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral")))  +
 theme(axis.title.x=element_blank(),#axis.text.x=element_blank(),
        text = element_text(size = 6), #, face = "bold"
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        plot.title = element_text(size = 6, face="plain"),
        axis.title = element_text(size = 6),
        plot.tag = element_text(size = 8, face="plain")) 

figure1_grob <- arrangeGrob(
  panel_a, panel_b1,
  panel_c, panel_d1,
  panel_e2, panel_d2,
  ncol = 2, nrow = 3,
  layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 6))
)

ggsave("Figure3.png", plot = figure1_grob, width = 30, height = 35, units = "cm")

