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

# wget https://fantom.gsc.riken.jp/6/datafiles/Core_FANTOM6/RELEASE_latest/expression_tables/Set1_human_dermal_fibroblasts/CAGE/expression.count.gene.nAnT-iCAGE.ASO.human_dermal_fibroblasts.tsv.bz2

matr <- read.csv(
  "expression.tpm.gene.nAnT-iCAGE.ASO.human_dermal_fibroblasts.tsv",
  sep = "\t", row.names = 1
)

matr_control_tpm <- matr[ENSEMBLE_GENES_INDEX:nrow(matr), CONTROLS]

# Mean expression (log10-transformed)
expr_mean <- data.frame(
  expr      = rowMeans(data.matrix(matr_control_tpm)),
  row.names = rownames(matr_control_tpm)
)
expr_mean <- log10(expr_mean + 1)

# Classify genes by expression tier and write outputs
tiers <- classify_expression(expr_mean)

saveRDS(tiers, 'expr_tiers.rds')
tiers <- readRDS('expr_tiers.rds')
for (tier in names(tiers)) {
  genes <- tiers[[tier]]
  genes <- setdiff(genes, c('ENSG00000108010')) # contains unfiltered pcr duplicate
  write.csv(as.data.frame(genes), paste0("ctd_", tier, ".txt"), quote = FALSE, row.names = FALSE)
  bed <- df[intersect(genes, rownames(df)), ]
  write.table(bed, file = paste0("ctd_", tier, ".bed"), quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
}
