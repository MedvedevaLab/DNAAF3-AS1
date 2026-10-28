library(ggplot2)
library(ggpubr)
library(stringr)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(BiocParallel)
library(pheatmap)

BiocParallel::register(BiocParallel::MulticoreParam(workers = 1, progressbar = TRUE))

gsea_sources <- c("GO:BP")
run_gsea <- function(logfc_vector, gene_names, sources = gsea_sources,
                     pval_cutoff = 0.05) {
    ranked <- sort(setNames(logfc_vector, gene_names), decreasing = TRUE)
    results <- list()
    for (src in sources) {
        message("Running GSEA: ", src)
        results[[src]] <- switch(src,
            "GO:BP" = gseGO(
                gene          = ranked,
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENSEMBL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = pval_cutoff
            )
        )
    }
    results
}

write_gsea_table <- function(results, file) {
    write.table(
        file      = file,
        as.data.frame(results[["GO:BP"]]),
        quote     = FALSE,
        row.names = TRUE,
        col.names = TRUE,
        sep       = "\t"
    )
}

all_deseq2_results <- read.csv("DESeq2_genes_ASO_all.tsv", sep = "\t")

deseq2_aso1 <- dplyr::filter(all_deseq2_results, perturb_id == "ASO_G0267577_01")
deseq2_aso3 <- dplyr::filter(all_deseq2_results, perturb_id == "ASO_G0267577_03")

gsea_results_aso1 <- run_gsea(deseq2_aso1$log2FC, deseq2_aso1$geneID)
write_gsea_table(gsea_results_aso1, "SupplementaryTable1E.csv")

gsea_results_aso3 <- run_gsea(deseq2_aso3$log2FC, deseq2_aso3$geneID)
write_gsea_table(gsea_results_aso3, "SupplementaryTable1F.csv")

