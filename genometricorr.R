library(ggplot2)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(GenometriCorr)
library(ChIPpeakAnno)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

contact_clusters <- readRDS("contact_clusters_hg.rds")

cluster_borders <- contact_clusters %>%
  group_by(seqnames, cluster) %>%
  summarize(cluster_start = min(start), cluster_end = max(end), .groups = "drop")

write.table(
    file      = "SupplementaryTable5.csv",
    cluster_borders,
    quote     = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    sep       = "\t"
)

gr_regions <- GRanges(
  seqnames = cluster_borders$seqnames,
  ranges   = IRanges(start = cluster_borders$cluster_start,
                     end   = cluster_borders$cluster_end)
)

diffbind_results <- read.csv('CTD-2587H24.5_H3K36me3_differential_peaks.csv', sep = '\t')
diffbind_sig <- diffbind_results[diffbind_results$FDR < 0.05,] # & abs(diffbind_results$Fold) > 1

diffbind_sig$significance <- "Unchanged"
diffbind_sig$significance[diffbind_sig$Fold < 0] <- "Lost"
diffbind_sig$significance[diffbind_sig$Fold > 0] <- "Gained"
diffbind_sig$significance <- factor(diffbind_sig$significance, levels = c( "Gained", "Lost", "Unchanged"))

# For all lost peaks
all_peaks <- dplyr::filter(diffbind_sig, Fold < 0)

peak_granges<- GRanges(
  seqnames = all_peaks$seqnames,
  ranges   = IRanges(start = all_peaks$start,
                     end   = all_peaks$end)
)

genometri_res <- GenometriCorrelation(peak_granges, gr_regions, 
  awhole.only = TRUE, 
  permut.number = 500, 
  keep.distributions = TRUE, 
  showProgressBar = TRUE,
  alternative = "attraction")

graphical.report(genometri_res, pdffile = "genometri_res_lost.pdf", show.all = TRUE)
saveRDS(genometri_res, "genometri_res_lost.rds")


# For all gained peaks
all_peaks <- dplyr::filter(diffbind_sig, Fold > 0)

peak_granges<- GRanges(
  seqnames = all_peaks$seqnames,
  ranges   = IRanges(start = all_peaks$start,
                     end   = all_peaks$end)
)

genometri_res <- GenometriCorrelation(peak_granges, gr_regions, 
  awhole.only = TRUE, 
  permut.number = 500, 
  keep.distributions = TRUE, 
  showProgressBar = TRUE,
  alternative = "attraction")

graphical.report(genometri_res, pdffile = "genometri_res_gained.pdf", show.all = TRUE)
saveRDS(genometri_res, "genometri_res_gained.rds")

# # For all lost peaks without chr19
cluster_borders_trans <- cluster_borders[cluster_borders$seqnames != "chr19", ]
gr_regions_trans <- GRanges(
  seqnames = cluster_borders_trans$seqnames,
  ranges   = IRanges(start = cluster_borders_trans$cluster_start,
                     end   = cluster_borders_trans$cluster_end)
)


all_peaks_trans  <- all_peaks[all_peaks$seqnames != "chr19", ]


peak_granges_trans<- GRanges(
  seqnames = all_peaks_trans$seqnames,
  ranges   = IRanges(start = all_peaks_trans$start,
                     end   = all_peaks_trans$end)
)

genometri_res_trans <- GenometriCorrelation(peak_granges_trans, gr_regions_trans, 
  awhole.only = TRUE, 
  permut.number = 500, 
  keep.distributions = TRUE, 
  showProgressBar = TRUE,
  alternative = "attraction")

graphical.report(genometri_res_trans, pdffile = "genometri_res_trans_lost.pdf", show.all = TRUE)
saveRDS(genometri_res_trans, "genometri_res_trans_lost.rds")

