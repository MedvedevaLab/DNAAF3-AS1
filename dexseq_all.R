library(GenomicFeatures)
library(stringr)
library(ggplot2)
library(ChIPseeker)
library(ChIPpeakAnno)
library(Rsamtools)
library(GenomicAlignments)
library(GenomeInfoDb)
library(DEXSeq)
library(dplyr)

get_dex_results <- function(promoter_regions, bam_files) {
    summarized_experiment <- summarizeOverlaps(
        promoter_regions, bam_files,
        singleEnd = TRUE,
        fragments = FALSE,
        ignore.strand = FALSE
    )

    colData(summarized_experiment)$condition <- factor(c(
        "aso", "aso",
        "control", "control", "control", "control",
        "control", "control", "control", "control"
    ))

    colData(summarized_experiment)$libType <- factor(rep("single-end", 10))

    dexseq_dataset <- DEXSeqDataSetFromSE(
        summarized_experiment,
        design = ~ sample + exon + condition:exon
    )

    dexseq_dataset <- estimateSizeFactors(dexseq_dataset)
    dexseq_dataset <- estimateDispersions(dexseq_dataset)
    dexseq_dataset <- testForDEU(dexseq_dataset)
    dexseq_dataset <- estimateExonFoldChanges(dexseq_dataset, fitExpToVar = "condition")

    dexseq_results     <- DEXSeqResults(dexseq_dataset)
    dexseq_results_df  <- as.data.frame(dexseq_results)

    per_gene_qvalues <- perGeneQValue(dexseq_results)
    deu_genes        <- names(per_gene_qvalues[per_gene_qvalues < 1.1])
    print(length(deu_genes))

    dexseq_results_deu_genes <- dplyr::filter(dexseq_results_df, groupID %in% deu_genes)
    sig_exon_results         <- dplyr::filter(dexseq_results_deu_genes, padj < 1.1)

    return(sig_exon_results)
}

get_annotation <- function(fantom_tsv, gtf)
{
    matr <- read.csv(fantom_tsv, sep='\t', row.names = 1)
    coords <- as.data.frame(str_split_fixed(rownames(matr) , "_", 4))
    colnames(coords) <- c('seqnames', 'start', 'end', 'strand')
    ranges <- toGRanges(coords, format='BED')
    txdb <- makeTxDbFromGFF(gtf)
    peakAnno <- annotatePeak(ranges, tssRegion=c(-1000, 100), TxDb=txdb, overlap="all")
    return(peakAnno)
}

# ── Input files ────────────────────────────────────────────────────────────────
fantom_promoter_tsv <- "/data/budkina/CHASERR_pulldown/expression.count.promoter.nAnT-iCAGE.ASO.human_dermal_fibroblasts.fixed.tsv"
gtf_file            <- "gencode.v49.annotation.gtf"

# ── Build promoter GRanges from FANTOM peak annotation ────────────────────────
peak_annotation  <- get_annotation(fantom_promoter_tsv, gtf_file)
promoter_regions <- peak_annotation@anno

excluded_annotations <- c("5' UTR", "3' UTR", "Downstream (<=300)", "Distal Intergenic", "Downstream")
promoter_regions <- promoter_regions[!(promoter_regions$annotation %in% excluded_annotations)]

promoter_regions$gene_id   <- promoter_regions$geneId
promoter_regions$tx_id     <- match(promoter_regions$gene_id, unique(promoter_regions$gene_id))
promoter_regions$tx_name   <- promoter_regions$gene_id
promoter_regions$exon_name <- paste(
    promoter_regions$geneId,
    seqnames(promoter_regions),
    ranges(promoter_regions),
    strand(promoter_regions)
)
promoter_regions$exon_id   <- seq_along(promoter_regions)
names(promoter_regions)    <- promoter_regions$exon_name

# Rank each promoter within its gene to assign exonic part numbers
promoter_regions_ranked <- data.frame(promoter_regions) %>%
    group_by(gene_id) %>%
    mutate(occurrence = row_number()) %>%
    ungroup()

promoter_regions$exon_rank    <- promoter_regions_ranked$occurrence
promoter_regions$exonic_part  <- promoter_regions$exon_rank
names(promoter_regions)       <- sprintf("%s:E%0.3d", promoter_regions$gene_id, promoter_regions$exonic_part)
seqlevelsStyle(promoter_regions) <- "UCSC"

# ── Columns to retain in output tables ────────────────────────────────────────
output_columns <- c(
    "groupID", "featureID", "exonBaseMean", "dispersion",
    "stat", "pvalue", "padj", "aso", "control",
    "log2fold_control_aso",
    "genomicData.seqnames", "genomicData.start", "genomicData.end",
    "genomicData.width", "genomicData.strand"
)

# ── ASO 1 ───────────────────────────────────────────────────────────
bam_files_aso1 <- BamFileList(c(
    "/data/budkina/ctd/fantom/CNhi10491_150515_SN549_0169_AC6YUKACXX_NoIndex_L007_R1_001.F6-001-RNA-003-03G07-216C.ACC.10003.bam",
    "/data/budkina/ctd/fantom/CNhi10492_150515_SN549_0169_AC6YUKACXX_NoIndex_L008_R1_001.F6-001-RNA-003-03H07-75B7.GCT.10003.bam",
    "/data/budkina/ctd/fantom/CNhi10486_150515_SN549_0169_AC6YUKACXX_NoIndex_L002_R1_001.F6-001-RNA-003-04E03-FFF3.ACG.10003.bam",
    "/data/budkina/ctd/fantom/CNhi10488_150515_SN549_0169_AC6YUKACXX_NoIndex_L004_R1_001.F6-001-RNA-003-03E09-8DBE.ACG.10003.bam",
    "/data/budkina/ctd/fantom/CNhi10489_150515_SN549_0169_AC6YUKACXX_NoIndex_L005_R1_001.F6-001-RNA-003-01E03-4FC3.GCT.10003.bam",
    "/data/budkina/ctd/fantom/CNhi10494_150515_SN549_0170_BC6YUVACXX_NoIndex_L003_R1_001.F6-001-RNA-003-02E03-8F86.GCT.10003.bam",
    "/data/budkina/ctd/fantom/CNhi10495_150515_SN549_0170_BC6YUVACXX_NoIndex_L004_R1_001.F6-001-RNA-003-04E09-F9BA.TAC.10003.bam",
    "/data/budkina/ctd/fantom/CNhi10496_150515_SN549_0170_BC6YUVACXX_NoIndex_L005_R1_001.F6-001-RNA-003-02E09-3AEE.GCG.10003.bam",
    "/data/budkina/ctd/fantom/CNhi10497_150515_SN549_0170_BC6YUVACXX_NoIndex_L006_R1_001.F6-001-RNA-003-01E09-21E9.TAC.10003.bam",
    "/data/budkina/ctd/fantom/CNhi10497_150515_SN549_0170_BC6YUVACXX_NoIndex_L006_R1_001.F6-001-RNA-003-03E03-93F8.ATG.10003.bam"
))

sig_deu_aso1 <- get_dex_results(promoter_regions, bam_files_aso1)
results_table_aso1 <- sig_deu_aso1[, output_columns]
results_table_aso1[, "genomicData.start"] <- results_table_aso1[, "genomicData.start"] - 1
results_table_aso1[, "genomicData.end"] <- results_table_aso1[, "genomicData.end"]

results_table_aso1 <- results_table_aso1 %>%
  mutate(featureID = sub("^E", "P", featureID))
 
results_table_aso1$log2fold_aso_control <- -results_table_aso1$log2fold_control_aso

rownames(results_table_aso1) <- paste0(results_table_aso1$groupID, ":", results_table_aso1$featureID)

results_table_aso1$groupID <- str_split_fixed(results_table_aso1$groupID, "\\.", 2)[, 1]

write.table(
    file      = "aso1_dex_results_not_sig.csv",
    results_table_aso1,
    quote     = FALSE,
    row.names = TRUE,
    col.names = TRUE,
    sep       = "\t"
)


results_table_aso1      <- read.csv("aso1_dex_results_not_sig.csv", sep = "\t")

smyd_aso1 <- results_table_aso1[results_table_aso1$groupID == 'ENSG00000185420',
    c("genomicData.seqnames", "genomicData.start", "genomicData.end", "featureID", "padj" , "genomicData.strand", "log2fold_aso_control")]
smyd_aso1$featureID <- ifelse(smyd_aso1$padj < 0.05, paste0(smyd_aso1$featureID, "*"), smyd_aso1$featureID)
smyd_aso1$padj <- as.numeric(smyd_aso1$padj > 0.05)
smyd_aso1_up <- smyd_aso1[smyd_aso1$log2fold_aso_control > 0,]
smyd_aso1_up$log2fold_aso_control <- NULL
smyd_aso1_down <- smyd_aso1[smyd_aso1$log2fold_aso_control < 0,]
smyd_aso1_down$log2fold_aso_control <- NULL

write.table(
    file      = "smyd3_aso_01_up.bed6",
    smyd_aso1_up,
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep       = "\t"
)

write.table(
    file      = "smyd3_aso_01_down.bed6",
    smyd_aso1_down,
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep       = "\t"
)

# ── ASO 3 ───────────────────────────────────────────────────────────
bam_files_aso3 <- BamFileList(c(
    "/data/budkina/ctd/fantom/CNhi10487_150515_SN549_0169_AC6YUKACXX_NoIndex_L003_R1_001.F6-001-RNA-003-03F07-416D.ACC.10003.bam",
    "/data/budkina/ctd/fantom/CNhi10495_150515_SN549_0170_BC6YUVACXX_NoIndex_L004_R1_001.F6-001-RNA-003-03E07-93DC.ACG.10003.bam",
    "/data/budkina/ctd/fantom/CNhi10486_150515_SN549_0169_AC6YUKACXX_NoIndex_L002_R1_001.F6-001-RNA-003-04E03-FFF3.ACG.10003.bam",
    "/data/budkina/ctd/fantom/CNhi10488_150515_SN549_0169_AC6YUKACXX_NoIndex_L004_R1_001.F6-001-RNA-003-03E09-8DBE.ACG.10003.bam",
    "/data/budkina/ctd/fantom/CNhi10489_150515_SN549_0169_AC6YUKACXX_NoIndex_L005_R1_001.F6-001-RNA-003-01E03-4FC3.GCT.10003.bam",
    "/data/budkina/ctd/fantom/CNhi10494_150515_SN549_0170_BC6YUVACXX_NoIndex_L003_R1_001.F6-001-RNA-003-02E03-8F86.GCT.10003.bam",
    "/data/budkina/ctd/fantom/CNhi10495_150515_SN549_0170_BC6YUVACXX_NoIndex_L004_R1_001.F6-001-RNA-003-04E09-F9BA.TAC.10003.bam",
    "/data/budkina/ctd/fantom/CNhi10496_150515_SN549_0170_BC6YUVACXX_NoIndex_L005_R1_001.F6-001-RNA-003-02E09-3AEE.GCG.10003.bam",
    "/data/budkina/ctd/fantom/CNhi10497_150515_SN549_0170_BC6YUVACXX_NoIndex_L006_R1_001.F6-001-RNA-003-01E09-21E9.TAC.10003.bam",
    "/data/budkina/ctd/fantom/CNhi10497_150515_SN549_0170_BC6YUVACXX_NoIndex_L006_R1_001.F6-001-RNA-003-03E03-93F8.ATG.10003.bam"
))

sig_deu_aso3 <- get_dex_results(promoter_regions, bam_files_aso3)
results_table_aso3 <- sig_deu_aso3[, output_columns]
results_table_aso3[, "genomicData.start"] <- results_table_aso3[, "genomicData.start"] - 1
results_table_aso3[, "genomicData.end"] <- results_table_aso3[, "genomicData.end"]

results_table_aso3 <- results_table_aso3 %>%
  mutate(featureID = sub("^E", "P", featureID))

results_table_aso3$log2fold_aso_control <- - results_table_aso3$log2fold_control_aso
results_table_aso3$groupID <- str_split_fixed(results_table_aso3$groupID, "\\.", 2)[, 1]

rownames(results_table_aso3) <- paste0(results_table_aso3$groupID, ":", results_table_aso3$featureID)

write.table(
    file      = "aso3_dex_results_not_sig.csv",
    results_table_aso3,
    quote     = FALSE,
    row.names = TRUE,
    col.names = TRUE,
    sep       = "\t"
)


results_table_aso3      <- read.csv("aso3_dex_results_not_sig.csv", sep = "\t")

smyd_aso3 <- results_table_aso3[results_table_aso3$groupID == 'ENSG00000185420',
    c("genomicData.seqnames", "genomicData.start", "genomicData.end", "featureID", "padj" , "genomicData.strand", "log2fold_aso_control")]
smyd_aso3$featureID <- ifelse(smyd_aso3$padj < 0.05, paste0(smyd_aso3$featureID, "*"), smyd_aso3$featureID)
smyd_aso3$padj <- as.numeric(smyd_aso3$padj > 0.05)
smyd_aso3_up <- smyd_aso3[smyd_aso3$log2fold_aso_control > 0,]
smyd_aso3_up$log2fold_aso_control <- NULL
smyd_aso3_down <- smyd_aso3[smyd_aso3$log2fold_aso_control < 0,]
smyd_aso3_down$log2fold_aso_control <- NULL

write.table(
    file      = "smyd3_aso_03_up.bed6",
    smyd_aso3_up,
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep       = "\t"
)

write.table(
    file      = "smyd3_aso_03_down.bed6",
    smyd_aso3_down,
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep       = "\t"
)
