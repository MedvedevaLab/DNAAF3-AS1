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
library(biomaRt)


library(dplyr)
library(ggplot2)
library(Seurat)
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

mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensembl <- mart

himorna_peaks <- read.csv("ENSG00000267577.tsv", sep = "\t")
himorna_peaks <- dplyr::filter(himorna_peaks, hm == "H3K36me3")

himorna_peaks_pos_cor <- himorna_peaks[himorna_peaks$corr > 0, ]
himorna_peaks_neg_cor <- himorna_peaks[himorna_peaks$corr < 0, ]
rownames(himorna_peaks_pos_cor) <- himorna_peaks_pos_cor$X
rownames(himorna_peaks_neg_cor) <- himorna_peaks_neg_cor$X

peak_anno_pos_cor <- get_annotation(himorna_peaks_pos_cor, "rho > 0", txdb)
peak_anno_pos_cor_anno <- as.data.frame(peak_anno_pos_cor@anno)
peak_anno_neg_cor <- get_annotation(himorna_peaks_neg_cor, "rho < 0", txdb)
peak_anno_neg_cor_anno <- as.data.frame(peak_anno_neg_cor@anno)

substrings <- c("Exon", "Intron")

peak_anno_pos_cor_promoter <- peak_anno_pos_cor_anno[peak_anno_pos_cor_anno$annotation == "Promoter",]
peak_anno_pos_cor_distal <- peak_anno_pos_cor_anno[peak_anno_pos_cor_anno$annotation == "Distal Intergenic",]
peak_anno_pos_cor_genebody <- peak_anno_pos_cor_anno[grepl(paste(substrings, collapse = "|"), peak_anno_pos_cor_anno$annotation),]

peak_anno_neg_cor_promoter <-  peak_anno_neg_cor_anno[peak_anno_neg_cor_anno$annotation == "Promoter",]
peak_anno_neg_cor_distal <- peak_anno_neg_cor_anno[peak_anno_neg_cor_anno$annotation == "Distal Intergenic",]
peak_anno_neg_cor_genebody <- peak_anno_neg_cor_anno[grepl(paste(substrings, collapse = "|"), peak_anno_neg_cor_anno$annotation),]

library(stringr)

peak_anno_pos_cor_promoter$geneId <-  str_split_fixed(peak_anno_pos_cor_promoter$geneId , "\\.", 2)[, 1]
peak_anno_pos_cor_distal$geneId <-  str_split_fixed(peak_anno_pos_cor_distal$geneId , "\\.", 2)[, 1]
peak_anno_pos_cor_genebody$geneId <-  str_split_fixed(peak_anno_pos_cor_genebody$geneId , "\\.", 2)[, 1]
peak_anno_neg_cor_promoter$geneId <-  str_split_fixed(peak_anno_neg_cor_promoter$geneId , "\\.", 2)[, 1]
peak_anno_neg_cor_distal$geneId <-  str_split_fixed(peak_anno_neg_cor_distal$geneId , "\\.", 2)[, 1]
peak_anno_neg_cor_genebody$geneId <-  str_split_fixed(peak_anno_neg_cor_genebody$geneId , "\\.", 2)[, 1]

promoter_names_up = getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), filters = 'ensembl_gene_id', 
    values = peak_anno_pos_cor_promoter$geneId, mart = ensembl, useCache = TRUE)

distal_names_up = getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), filters = 'ensembl_gene_id', 
    values = peak_anno_pos_cor_distal$geneId, mart = ensembl, useCache = TRUE)

genebody_names_up = getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), filters = 'ensembl_gene_id', 
    values = peak_anno_pos_cor_genebody$geneId, mart = ensembl, useCache = TRUE)


promoter_names_down = getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), filters = 'ensembl_gene_id', 
    values = peak_anno_neg_cor_promoter$geneId, mart = ensembl, useCache = TRUE)

distal_names_down = getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), filters = 'ensembl_gene_id', 
    values = peak_anno_neg_cor_distal$geneId, mart = ensembl, useCache = TRUE)

genebody_names_down = getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), filters = 'ensembl_gene_id', 
    values = peak_anno_neg_cor_genebody$geneId, mart = ensembl, useCache = TRUE)


# Single-cell GSE112013 dataset
seurat_object2 <- readRDS('testis.rds')
Idents(object = seurat_object2) <- seurat_object2@meta.data$celltype
seurat_object2$celltype <- factor(seurat_object2$celltype, 
                        levels = c("SSCs", "Differentiating S'gonia", "Early primary S'cytes",
                          "Late primary S'cytes", "Round S'tids", "Elongated S'tids", "Sperm1",
                          "Sperm2", "Sperm3", "Macrophage", "Endothelial cells", "Myoid cells",
                          "Sertoli cells"))


markers <- list()
markers$pos <- unique(promoter_names_up$hgnc_symbol)
markers$neg <- unique(promoter_names_down$hgnc_symbol)

seurat_object2 <- AddModuleScore(seurat_object2, features = markers, name = "add_module_score")
saveRDS(seurat_object2, 'testis_promoter.rds')

markers <- list()
markers$pos <- unique(genebody_names_up$hgnc_symbol)
markers$neg <- unique(genebody_names_down$hgnc_symbol)

seurat_object2 <- AddModuleScore(seurat_object2, features = markers, name = "add_module_score")
saveRDS(seurat_object2, 'testis_genebody.rds')