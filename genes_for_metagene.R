library(GenomicFeatures)
library(ChIPpeakAnno)
library(biomaRt)
library(stringr)
library(AnnotationHub)

# Read reference
txdb <- makeTxDbFromGFF("gencode.v49.annotation.gtf")
gr <- toGRanges(txdb, feature = "gene")
df <- data.frame(
  seqnames = seqnames(gr),
  starts   = start(gr) - 1,
  ends     = end(gr),
  names    = rep(".", length(gr)),
  scores   = rep(".", length(gr)),
  strands  = strand(gr)
)
rownames(df) <- names(gr)

ENSEMBLE_GENES_INDEX <- 85706
CONTROLS <- c(188, 204, 216, 256, 264, 270, 280, 276)

classify_expression <- function(expr_mean) {
  max_expr <- max(expr_mean$expr)
  list(
    high   = rownames(dplyr::filter(expr_mean, expr >= 0.5 * max_expr)),
    medium = rownames(dplyr::filter(expr_mean, expr >= 0.25 * max_expr & expr < 0.5 * max_expr)),
    low    = rownames(dplyr::filter(expr_mean, expr < 0.25 * max_expr)),
    zero   = rownames(dplyr::filter(expr_mean, expr == 0))
  )
}

matr <- read.csv(
  "expression.count.gene.nAnT-iCAGE.ASO.human_dermal_fibroblasts.tsv",
  sep = "\t", row.names = 1
)

matr_control <- matr[ENSEMBLE_GENES_INDEX:nrow(matr), CONTROLS]

ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

gene_lengths <- getBM(
  attributes = c("ensembl_gene_id", "ensembl_transcript_id", "transcript_length", "cds_length"),
  filters    = "ensembl_gene_id",
  values     = rownames(matr_control),
  mart       = ensembl,
  useCache   = FALSE
)
write.table(gene_lengths, "gene_lengths.csv", quote = FALSE, sep = "\t", row.names = TRUE)

gene_canonical_transcript <- getBM(
  attributes = c("ensembl_gene_id", "ensembl_transcript_id", "transcript_is_canonical"),
  filters    = "ensembl_gene_id",
  values     = rownames(matr_control),
  mart       = ensembl,
  useCache   = FALSE
)
gene_canonical_transcript <- gene_canonical_transcript[
  !is.na(gene_canonical_transcript$transcript_is_canonical),
]

gene_lengths <- merge(gene_canonical_transcript, gene_lengths,
                      by = c("ensembl_gene_id", "ensembl_transcript_id"))
gene_lengths <- gene_lengths[!is.na(gene_lengths$transcript_length), ]

# Compute TPM
x <- sapply(colnames(matr_control), function(col) {
  as.numeric(matr_control[gene_lengths$ensembl_gene_id, col]) / gene_lengths$transcript_length
})
matr_control_tpm <- t(t(x) * 1e6 / colSums(x))

# Mean expression (log10-transformed)
expr_mean <- data.frame(
  expr      = rowMeans(matr_control_tpm),
  row.names = gene_lengths$ensembl_gene_id
)
expr_mean <- log10(expr_mean + 1)

# Classify genes by expression tier and write outputs
rownames(df) <- str_split_fixed(rownames(df), "\\.", 2)[, 1]
tiers <- classify_expression(expr_mean)

saveRDS(tiers, 'expr_tiers.rds')
tiers <- readRDS('expr_tiers.rds')
for (tier in names(tiers)) {
  genes <- tiers[[tier]]
  genes <- setdiff(genes, c('ENSG00000108010'))
  write.csv(as.data.frame(genes), paste0("Figure2A/ctd_", tier, ".txt"), quote = FALSE, row.names = FALSE)
  bed <- df[intersect(genes, rownames(df)), ]
  write.table(bed, file = paste0("Figure2A/ctd_", tier, ".bed"), quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
}
