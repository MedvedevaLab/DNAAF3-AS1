library(DiffBind)
library(ggplot2)

BiocParallel::register(BiocParallel::MulticoreParam(workers = 32, progressbar = TRUE))

db <- dba(sampleSheet = "CTD-2587H24.5_H3K36me3.csv")

# Count reads (use full peaks for broad marks)
db <- dba.count(db, 
                summits = FALSE,
                minOverlap = 1)

# Define contrast (Control vs KD)
db <- dba.contrast(db, design = ~Condition, contrast = c('Condition', 'KD', 'Control'),
                   minMembers = 2)

db <- dba.normalize(db, method = DBA_ALL_METHODS)

# Perform differential analysis
db <- dba.analyze(db, method = DBA_ALL_METHODS)

diff_peaks_deseq <- dba.report(db, 
                         th = 0.05, 
                         bUsePval = FALSE,
                         bAll = TRUE,
                         method = DBA_DESEQ2)
if(is.null(diff_peaks_deseq))
{
  next
}
diff_peaks_edger <- dba.report(db, 
                         th = 0.05, 
                         bUsePval = FALSE,
                         bAll = TRUE,
                         method = DBA_EDGER)

write.table(diff_peaks_deseq, file = "CTD-2587H24.5_H3K36me3_differential_peaks.csv", quote=F, sep="\t", row.names=T, col.names=T)

gr <- diff_peaks_deseq[diff_peaks_deseq$Fold > 0, ]
df <- data.frame(seqnames=seqnames(gr),
  starts=start(gr)-1,
  ends=end(gr),
  names=c(rep(".", length(gr))),
  scores=c(rep(".", length(gr))),
  strands=strand(gr))

write.table(df, file = "CTD-2587H24.5_H3K36me3_differential_peaks_up.bed", quote=F, sep="\t", row.names=F, col.names=F)

gr <- diff_peaks_deseq[diff_peaks_deseq$Fold < 0, ]
df <- data.frame(seqnames=seqnames(gr),
  starts=start(gr)-1,
  ends=end(gr),
  names=c(rep(".", length(gr))),
  scores=c(rep(".", length(gr))),
  strands=strand(gr))

write.table(df, file = "CTD-2587H24.5_H3K36me3_differential_peaks_down.bed", quote=F, sep="\t", row.names=F, col.names=F)