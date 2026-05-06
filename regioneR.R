library(regioneR)

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


dim(df_all[df_all$Density..up.50.kb..up.5.kb.!=0,])
dim(df_all[df_all$Density..up.5.kb..gene.body.!=0,])
dim(df_all[df_all$Density..gene.body.!=0,])
dim(df_all[df_all$Density..gene.body..down.5.kb.!=0,])
dim(df_all[df_all$Density..down.5.kb..down.50.kb.!=0,])

df_all_gene_body <- df_all[df_all$Density..gene.body.!=0,]

# Target genes intersection

gtf_name <- "gencode.v49.annotation.gtf"
txdb <- makeTxDbFromGFF(gtf_name)
genes_gr <- genes(txdb)   # returns a GRanges with gene_id in metadata

gr_regions <- GRanges(
  seqnames = df_all$Chromosome,
  ranges   = IRanges(start = df_all$Gene.start,
                     end   = df_all$Gene.end)
)


# library(GenomicRanges)

# # Assuming your object is called 'gr'
# start(gr_regions) <- start(gr_regions) - 5000
# end(gr_regions)   <- end(gr_regions)   + 5000


# library(GenomicRanges)
# library(GenomicFeatures)
# gtf_name <- "gencode.v49.annotation.gtf"
# txdb <- makeTxDbFromGFF(gtf_name)
# genes_gr <- genes(txdb)   # returns a GRanges with gene_id in metadata
# hits <- findOverlaps(gr_regions, genes_gr, ignore.strand = TRUE)


himorna_peaks <- read.csv("ENSG00000267577.tsv", sep = "\t")
himorna_peaks <- dplyr::filter(himorna_peaks, hm == "H3K36me3")

himorna_peaks_pos_cor <- himorna_peaks[himorna_peaks$corr > 0, ]
himorna_peaks_neg_cor <- himorna_peaks[himorna_peaks$corr < 0, ]
rownames(himorna_peaks_pos_cor) <- himorna_peaks_pos_cor$X
rownames(himorna_peaks_neg_cor) <- himorna_peaks_neg_cor$X

coords <- himorna_peaks_pos_cor[, c("chrom", "start", "end")]
coords <- cbind(coords, "*")
colnames(coords) <- c("seqnames", "start", "end", "strand")
coords$start <- as.numeric(coords$start)
coords$end   <- as.numeric(coords$end)
peak_granges <- toGRanges(coords, format = "BED")

pt.reg <- permTest(A=gr_regions, B=peak_granges, ntimes=50, randomize.function=resampleRegions, evaluate.function=numOverlaps, alternative = 'auto')
overlapPermTest(A=gr_regions, B=peak_granges, ntimes=500)
hits <- findOverlaps(gr_regions, peak_granges, ignore.strand = TRUE)

gr_regions[queryHits(hits)]

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
