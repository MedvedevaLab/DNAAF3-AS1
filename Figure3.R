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

# pyGenomeTracks --tracks plot_smyd3.ini --region chr1:245749342-246507312 -o Figure3A.png --plotWidth 120 --dpi 130

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
Idents(object = seurat_object) <- seurat_object@meta.data$celltype
seurat_object$celltype <- factor(seurat_object$celltype, 
                        levels = c("SSCs", "Differentiating S'gonia", "Early primary S'cytes",
                          "Late primary S'cytes", "Round S'tids", "Elongated S'tids", "Sperm1",
                          "Sperm2", "Sperm3", "Macrophage", "Endothelial cells", "Myoid cells",
                          "Sertoli cells"))

cols <- rev(c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#D37295", "#499894", "#A0CBE8", "#BAB0AC"))

# UMAPs with annotation
panel_a <- DimPlot(seurat_object, reduction = "umap", group.by = "celltype", label = FALSE, repel = TRUE, raster=FALSE) + 
panel_theme +
  theme(plot.title = element_blank()) + 
  guides(fill = guide_legend(ncol = 3)) +
  scale_color_manual(values = rev(cols), guide = guide_legend(reverse = TRUE), drop = FALSE) +
  labs(title = "", tag = "A")

panel_b1 <- VlnPlot(
  seurat_object, group.by = "celltype",
  features = c('CTD-2587H24.5'),
  alpha = 0, log = TRUE) + 
  panel_theme +
  theme(strip.text.x = element_text(size = 6),
        strip.background = element_blank(),
        legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_fill_manual(values = rev(cols), guide = guide_legend(reverse = TRUE), drop = FALSE) +
  labs(title = "DNAAF3-AS1", tag = "")

panel_b2 <- VlnPlot(
  seurat_object, group.by = "celltype",
  features = c('DNAAF3'),
  alpha = 0, log = TRUE) + 
  panel_theme +
  theme(strip.text.x = element_text(size = 6),
        strip.background = element_blank(),
        legend.position="none",
        axis.title.x = element_blank()) +
  scale_fill_manual(values = rev(cols), guide = guide_legend(reverse = TRUE), drop = FALSE) +
  labs(title = "DNAAF3", tag = "")

p <- ggarrange(panel_a, panel_b1, panel_b2, ncol = 1, nrow = 3, common.legend = TRUE, legend="bottom")
ggsave("test.png", plot = p, width = 10, height = 20, units = "cm")

seurat_object <- readRDS('testis_promoter.rds')

panel_c <- FeaturePlot(seurat_object, paste0("add_module_score", 2), pt.size = 0.5, reduction = "umap") + 
  scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral"))) +
  panel_theme +
  theme(text = element_text(size = 6)) +
  labs(title = "", tag = "C")

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
panel_d <- ggplot(df, aes(x = chromosome, y = experiment, fill = density)) +
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
  labs(title = "", x = 'Chromosome', y = 'Experiment', fill = 'log10(Contacts Density)', tag = "D")

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

###############################################################

himorna_peaks <- read.csv("ENSG00000267577.tsv", sep = "\t")
himorna_peaks <- dplyr::filter(himorna_peaks, hm == "H3K36me3")
himorna_contacts <- dplyr::filter(himorna_peaks, gene_ens %in% target_genes)


fantom <- read.csv('DESeq2_genes_ASO_all.tsv', sep = '\t')
aso01 <- dplyr::filter(fantom, perturb_id == 'ASO_G0267577_01')
aso01 <- aso01[aso01$fdr < 0.05, ]

aso01_contacts <- dplyr::filter(aso01, geneID %in% target_genes)
aso01_contacts$in_himorna <- aso01_contacts$geneID  %in% himorna_contacts$gene_ens
aso01_contacts$log_padj <-  -log(aso01_contacts$fdr,10)

colors <- c('TRUE' = "#E31A1C", 'FALSE' = "#000000")

panel_e <- ggplot(aso01_contacts, aes(x = log2FC,
  y = log_padj,
  label = ifelse(in_himorna, geneSymbol, NA),
  color = in_himorna)) +
  geom_point(shape = 1) + 
  geom_text_repel( color = 'black', alpha = 1,  hjust = 0, vjust = 0, size = 4, max.overlaps = 15) +
  scale_color_manual(values = colors, name = "") +
  theme_bw() +
  panel_theme +
  theme(text = element_text(size = 4)) +
    labs(title = "ASO 01", x = "", y = "", tag = "E")


###############################################################
# Heatmap
dex_results_aso1      <- read.csv("aso1_dex_results_not_sig.csv", sep = "\t")
dex_results_aso1      <- dex_results_aso1[!is.na(dex_results_aso1$log2fold_aso_control), ]
dex_results_aso1$groupID  <- str_split_fixed(dex_results_aso1$groupID , "\\.", 2)[, 1]

dex_results_aso3      <- read.csv("aso3_dex_results_not_sig.csv", sep = "\t")
dex_results_aso3      <- dex_results_aso3[!is.na(dex_results_aso3$log2fold_aso_control), ]
dex_results_aso3$groupID  <- str_split_fixed(dex_results_aso3$groupID , "\\.", 2)[, 1]

# Inner join on groupID + featureID; keep annotation from aso1, fall back to aso3
merged_dex <- inner_join(dex_results_aso1, dex_results_aso3,
                         by = c("groupID", "featureID"),
                         suffix = c("", ".y"))

# Keep features with padj < 0.05 in at least one experiment
merged_dex <- merged_dex %>%
  filter(!is.na(padj) | !is.na(padj.y)) %>%
  filter(replace_na(padj, 1) < 0.05 | replace_na(padj.y, 1) < 0.05)

# Keep only target genes and the three annotation categories
merged_dex <- merged_dex %>%
  filter(groupID %in% target_genes)

# Annotate promoters
coords <- merged_dex[, c("genomicData.seqnames", "genomicData.start", "genomicData.end", "genomicData.strand")]
colnames(coords) <- c("seqnames", "start", "end", "strand")
granges<- toGRanges(coords, format = "BED")
peak_anno_granges_aso1   <- annotatePeak(granges, tssRegion = c(-1000, 100), TxDb = txdb, overlap = "all")
anno_df <- as.data.frame(peak_anno_granges_aso1@anno)
anno_df$annotation <- str_split_fixed(anno_df$annotation, " ", 2)[, 1]
anno_df <- dplyr::filter(anno_df, annotation %in% c("Promoter", "Exon", "Intron"))
rownames(anno_df) <- paste0(anno_df$seqnames, anno_df$start - 1, anno_df$end, anno_df$strand)

merged_dex$annotation <- anno_df[paste0(merged_dex$genomicData.seqnames, 
  merged_dex$genomicData.start,
  merged_dex$genomicData.end,
  merged_dex$genomicData.strand),'annotation']

merged_dex <- merged_dex %>%
  filter(annotation %in% c("Exon", "Intron"))

# Row labels: groupID@featureID

mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensembl <- mart

gene_names = getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), filters = 'ensembl_gene_id', 
    values = merged_dex$groupID, mart = ensembl, useCache = TRUE)
rownames(gene_names) <- gene_names$ensembl_gene_id

merged_dex$gene_name <- gene_names[merged_dex$groupID,'hgnc_symbol']
merged_dex[is.na(merged_dex$gene_name), 'gene_name'] = merged_dex$groupID



merged_dex$row_label <- paste0(merged_dex$gene_name, "@", merged_dex$featureID)
rownames(merged_dex) <- merged_dex$row_label
df1 <- merged_dex[ ,c('row_label', 'padj', 'log2fold_aso_control', 'annotation')]

df1$condition <- 'ASO_01'
colnames(df1) <- c('row_label', 'padj', 'log2fold_aso_control', 'annotation', 'condition')
df2 <- merged_dex[ ,c('row_label', 'padj.y', 'log2fold_aso_control.y', 'annotation')]
df2$condition <- 'ASO_03'
colnames(df2) <- c('row_label', 'padj', 'log2fold_aso_control', 'annotation', 'condition')

heatmap_long <- rbind(df1, df2)
# Factor levels
heatmap_long$annotation <- factor(heatmap_long$annotation, levels = c("Exon", "Intron"))
heatmap_long$condition  <- factor(heatmap_long$condition,  levels = c("ASO_01", "ASO_03"))

# # Cluster rows within each annotation group and concatenate the ordered labels
# clust_order <- lapply(c("Exon", "Intron"), function(ann) {
#   mat <- merged_dex %>%
#     filter(annotation == ann) %>%
#     dplyr::select(row_label, aso, aso.y, control, control.y) %>%
#     tibble::column_to_rownames("row_label") %>%
#     as.matrix()
#   if (nrow(mat) < 2) return(rownames(mat))
#   rownames(mat)[hclust(dist(mat, method = "euclidean"), method = "complete")$order]
# })
# row_order <- rev(unlist(clust_order))
# heatmap_long$row_label <- factor(heatmap_long$row_label, levels = row_order)

# heatmap_long <- as.data.frame(heatmap_long)
# heatmap_long$value <- scale(heatmap_long$value)


panel_f <- ggplot(heatmap_long, aes(x = condition, y = row_label, color = log2fold_aso_control, size = -log10(padj))) +
  geom_point() +
  facet_grid(annotation ~ ., scales = "free_y", space = "free_y") +
  scale_color_gradientn(
    colors = c('darkblue', "#1f77b4", "white", "#d62728", 'darkred'),
    values = rescale(c(-14, -2, 0, 2, 4))) +
  theme_bw() +
  panel_theme +
  theme(
    strip.background = element_blank(),
    strip.text.y     = element_text(size = 6, angle = 0),
    axis.text.y      = element_text(size = 4),
    legend.title     = element_text(size = 6),
    legend.text      = element_text(size = 5))+
    labs(title = "", x = "", y = "", tag = "F")


figure3_grob <- arrangeGrob(
  p, panel_c, panel_d, panel_e, panel_f, #  panel_e1,
  ncol = 2, nrow = 3,
  layout_matrix = rbind(c(1, 3), c(1, 4), c(2, 5))
)

ggsave("Figure3.png", plot = figure3_grob, width = 16, height = 26, units = "cm")
