library(DiffBind)
library(ggplot2)

BiocParallel::register(BiocParallel::MulticoreParam(workers = 32, progressbar = TRUE))

samples <- c(
'CTD-2587H24.5_H3K36me3',
'FGD5-AS1_H3K9ac',
'EMX2OS_H3K9ac',
'AC005592.2_H3K9me3',
'RP11-398K22.12_H3K27ac',
'JPX_H3K9me3',
'MAPKAPK5-AS1_H3K27me3',
'MAPKAPK5-AS1_H3K4me3')

for(sample in samples)
{
  # Step 1: Load metadata (ensure 'Condition' column exists in sample.csv)
  db <- dba(sampleSheet = paste0(sample, ".csv"))

  # Step 2: Initial correlation heatmap (quality check)
  # pdf(paste0(sample, "_DiffBind_initial_heatmap.pdf"), width=8, height=8)
  # dba.plotHeatmap(db)
  # dev.off()

  # Step 3: Count reads (use full peaks for broad marks)
  db <- dba.count(db, 
                  summits = FALSE,
                  minOverlap = 1)     # Your estimated fragment length , fragmentSize = 200

  # Step 4: Define contrast (Control vs KD)
  db <- dba.contrast(db, design = ~Condition, contrast = c('Condition', 'KD', 'Control'), #                   categories = DBA_CONDITION, 
                     minMembers = 2)

  db <- dba.normalize(db, method = DBA_ALL_METHODS) #DBA_ALL_METHODS  , background=TRUE , normalize=DBA_NORM_NATIVE, background=TRUE

  # Step 5: Perform differential analysis
  db <- dba.analyze(db, method = DBA_ALL_METHODS)

  # Step 6: MA plot with normalized counts
  # pdf(paste0(sample, "_DiffBind_MA_plot.pdf"), width=8, height=6)
  # dba.plotMA(db, bNormalized = TRUE)
  # dev.off()

  # Step 7: Heatmap of differential binding sites
  # pdf(paste0(sample, "_DiffBind_differential_heatmap.pdf"), width=8, height=10)
  # dba.plotHeatmap(db, 
  #                 contrast = 1, 
  #                 correlations = FALSE, 
  #                 colScheme = "RdBu")  # Blue-White-Red
  # dev.off()

  # Step 8: Extract significant peaks (FDR < 0.05, |log2FC| > 1)
  diff_peaks_deseq <- dba.report(db, 
                           th = 0.05, 
                           # fold = 0.5,
                           bUsePval = FALSE,
                           bAll = TRUE,
                           method = DBA_DESEQ2)
  if(is.null(diff_peaks_deseq))
  {
    next
  }
  diff_peaks_edger <- dba.report(db, 
                           th = 0.05, 
                           # fold = 1,
                           bUsePval = FALSE,
                           bAll = TRUE,
                           method = DBA_EDGER)


  write.table(diff_peaks_deseq, file = paste0(sample, "_differential_peaks.csv"), quote=F, sep="\t", row.names=T, col.names=T)

  gr <- diff_peaks_deseq[diff_peaks_deseq$Fold > 0, ]
  df <- data.frame(seqnames=seqnames(gr),
    starts=start(gr)-1,
    ends=end(gr),
    names=c(rep(".", length(gr))),
    scores=c(rep(".", length(gr))),
    strands=strand(gr))

  write.table(df, file = paste0(sample, "_differential_peaks_up.bed"), quote=F, sep="\t", row.names=F, col.names=F)

  gr <- diff_peaks_deseq[diff_peaks_deseq$Fold < 0, ]
  df <- data.frame(seqnames=seqnames(gr),
    starts=start(gr)-1,
    ends=end(gr),
    names=c(rep(".", length(gr))),
    scores=c(rep(".", length(gr))),
    strands=strand(gr))

  write.table(df, file = paste0(sample, "_differential_peaks_down.bed"), quote=F, sep="\t", row.names=F, col.names=F)
}