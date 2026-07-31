library(dplyr)
library(tidyr)
library(GenomicRanges)

# FANTOM6 RADICL-seq contacts for AC010327.4

BEDPE_COLS <- c('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2',
                'name', 'score', 'strand1', 'strand2', 'qual1', 'qual2', 'anno1')

radicl_files <- c(
  "hiPSC_1_ACAGTG_DNAAF3_AS1.bedpe",
  "hiPSC_1_ACTTGA_DNAAF3_AS1.bedpe",
  "hiPSC_1_CAGATC_DNAAF3_AS1.bedpe",
  "hiPSC_1_GCCAAT_DNAAF3_AS1.bedpe",
  "hiPSC_2_CTTGTA_DNAAF3_AS1.bedpe",
  "hiPSC_2_GATCAG_DNAAF3_AS1.bedpe",
  "hiPSC_2_GGCTAC_DNAAF3_AS1.bedpe",
  "hiPSC_2_TAGCTT_DNAAF3_AS1.bedpe"
)

radicl_dfs <- lapply(radicl_files, function(f) {
  df <- read.csv(f, sep = '\t', header = FALSE)
  colnames(df) <- BEDPE_COLS
  df
})

radicl_chromatin_coords <- do.call(rbind, lapply(radicl_dfs, `[`, c("chrom2", "start2", "end2")))
radicl_chromatin_coords$chrom2 <- paste0("chr", radicl_chromatin_coords$chrom2)

radicl_gr <- GRanges(
  seqnames = radicl_chromatin_coords$chrom2,
  ranges   = IRanges(start = as.numeric(radicl_chromatin_coords$start2) - 100000,
                     end   = as.numeric(radicl_chromatin_coords$end2)   + 100000)
)

# RNAChrom2 database contacts for AC010327.4

rnachrom_experiments <- list(
  list(file = "rnachrom1.csv",  experiment = "iMARGI1"),
  list(file = "rnachrom2.csv",  experiment = "iMARGI2"),
  list(file = "rnachrom3.csv",  experiment = "iMARGI3"),
  list(file = "rnachrom4.csv",  experiment = "iMARGI4"),
  list(file = "rnachrom5.csv",  experiment = "iMARGI5"),
  list(file = "rnachrom7.csv",  experiment = "GRIDseq1"),
  list(file = "rnachrom8.csv",  experiment = "GRIDseq2"),
  list(file = "rnachrom10.csv", experiment = "Red-C")
)

rnachrom_contacts <- do.call(rbind, lapply(rnachrom_experiments, function(exp) {
  df <- read.csv(exp$file)
  df$experiment <- exp$experiment
  df
}))

rnachrom_gr <- GRanges(
  seqnames = rnachrom_contacts$Chromosome,
  ranges   = IRanges(start = rnachrom_contacts$Gene.start - 10000,
                     end   = rnachrom_contacts$Gene.end   + 10000)
)

combined_contacts <- rbind(#cbind(as.data.frame(radicl_gr),   source = 'radicl',   experiment = 'radicl'),
  cbind(as.data.frame(rnachrom_gr), source = 'rnachrom', experiment = rnachrom_contacts$experiment)
)

# Clusterisation
contact_clusters_list <- list()
for (chr in unique(combined_contacts$seqnames)) {
  chr_contacts <- combined_contacts[combined_contacts$seqnames == chr, ]
  chr_gr <- GRanges(chr, IRanges(start = chr_contacts$start, end = chr_contacts$end))
  chr_gr$source     <- chr_contacts$source
  chr_gr$experiment <- chr_contacts$experiment
  merged_regions    <- reduce(chr_gr, min.gapwidth = 10000)
  overlap_hits      <- findOverlaps(chr_gr, merged_regions)
  chr_gr$cluster    <- subjectHits(overlap_hits)
  cluster_df        <- as.data.frame(chr_gr)
  cluster_df$cluster <- paste0(chr, "_", cluster_df$cluster)
  contact_clusters_list[[chr]] <- cluster_df
}

contact_clusters <- do.call(rbind, contact_clusters_list)

saveRDS(contact_clusters, "contact_clusters_hg.rds")
