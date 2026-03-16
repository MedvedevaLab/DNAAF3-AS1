library(ggplot2)
library(ggpubr)
library(stringr)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(BiocParallel)

BiocParallel::register(BiocParallel::MulticoreParam(workers = 1, progressbar = TRUE))

# Prepare a named, sorted log2FC vector for GSEA
prepare_ranked_list <- function(logfc_vector, gene_names) {
    names(logfc_vector) <- gene_names
    sort(logfc_vector, decreasing = TRUE)
}

# Run GO-based GSEA for a single ontology
run_gsea_go <- function(ranked_gene_list, ontology = "BP", pval_cutoff = 0.05) {
    gsea_result <- gseGO(
        gene          = ranked_gene_list,
        OrgDb         = org.Hs.eg.db,
        keyType       = "ENSEMBL",
        ont           = ontology,
        pAdjustMethod = "BH",
        pvalueCutoff  = pval_cutoff
    )
    gsea_result
}

# Run MSigDB-based GSEA for a given category/subcategory
run_gsea_msigdb <- function(ranked_gene_list, category = "H",
                            subcategory = NULL, pval_cutoff = 0.05) {
    msigdb_gene_sets <- msigdbr(
        species     = "Homo sapiens",
        category    = category,
        subcategory = subcategory
    )
    term_to_gene <- msigdb_gene_sets[, c("gs_name", "ensembl_gene")]

    gsea_result <- GSEA(
        geneList      = ranked_gene_list,
        TERM2GENE     = term_to_gene,
        pAdjustMethod = "BH",
        pvalueCutoff  = pval_cutoff
    )
    gsea_result
}

#' Dispatch GSEA to the appropriate database based on source label
run_gsea_for_source <- function(logfc_vector, gene_names, gsea_source, pval_cutoff = 0.05) {
    ranked_gene_list <- prepare_ranked_list(logfc_vector, gene_names)
    gsea_result <- switch(gsea_source,
        "GO:BP" = run_gsea_go(ranked_gene_list, ontology = "BP", pval_cutoff = pval_cutoff),
        "GO:MF" = run_gsea_go(ranked_gene_list, ontology = "MF", pval_cutoff = pval_cutoff),
        "GO:CC" = run_gsea_go(ranked_gene_list, ontology = "CC", pval_cutoff = pval_cutoff),
        "MSigDB_Hallmark"       = run_gsea_msigdb(ranked_gene_list, category = "H",  subcategory = NULL,           pval_cutoff = pval_cutoff),
        "MSigDB_C2_Reactome"    = run_gsea_msigdb(ranked_gene_list, category = "C2", subcategory = "CP:REACTOME",  pval_cutoff = pval_cutoff),
        "MSigDB_C6_Oncogenic"   = run_gsea_msigdb(ranked_gene_list, category = "C6", subcategory = NULL,           pval_cutoff = pval_cutoff),
        "MSigDB_C7_Immunologic" = run_gsea_msigdb(ranked_gene_list, category = "C7", subcategory = NULL,           pval_cutoff = pval_cutoff),
        stop(paste("Unknown GSEA source:", gsea_source))
    )
    gsea_result
}

# ── Load DESeq2 results ────────────────────────────────────────────────────────
all_deseq2_results <- read.csv("DESeq2_genes_ASO_all.tsv", sep = "\t")

deseq2_aso1 <- dplyr::filter(all_deseq2_results, perturb_id == "ASO_G0267577_01")
deseq2_aso3 <- dplyr::filter(all_deseq2_results, perturb_id == "ASO_G0267577_03")

deseq2_aso1_sig <- deseq2_aso1[deseq2_aso1$fdr < 0.05, ]
deseq2_aso3_sig <- deseq2_aso3[deseq2_aso3$fdr < 0.05, ]

# ── GSEA sources to test ───────────────────────────────────────────────────────
gsea_sources <- c(
    "GO:BP",
    "GO:MF",
    "MSigDB_Hallmark",
    "MSigDB_C2_Reactome",
    "MSigDB_C6_Oncogenic",
    "MSigDB_C7_Immunologic"
)

# ── ASO 1 ───────────────────────────────────────────────────────────
ranked_logfc_aso1 <- deseq2_aso1_sig$log2FC
names(ranked_logfc_aso1) <- deseq2_aso1_sig$geneID

gsea_results_aso1 <- list()
for (gsea_source in gsea_sources) {
    message("Running GSEA: ", gsea_source)
    gsea_results_aso1[[gsea_source]] <- run_gsea_for_source(
        ranked_logfc_aso1, names(ranked_logfc_aso1), gsea_source
    )
    write.table(
        file      = paste0(gsea_source, "_aso01.csv"),
        as.data.frame(gsea_results_aso1[[gsea_source]]),
        quote     = FALSE,
        row.names = TRUE,
        col.names = TRUE,
        sep       = "\t"
    )
}
saveRDS(gsea_results_aso1, "GSEA_aso01.rds")

# ── ASO 3 ───────────────────────────────────────────────────────────
ranked_logfc_aso3 <- deseq2_aso3_sig$log2FC
names(ranked_logfc_aso3) <- deseq2_aso3_sig$geneID

gsea_results_aso3 <- list()
for (gsea_source in gsea_sources) {
    message("Running GSEA: ", gsea_source)
    gsea_results_aso3[[gsea_source]] <- run_gsea_for_source(
        ranked_logfc_aso3, names(ranked_logfc_aso3), gsea_source
    )
    write.table(
        file      = paste0(gsea_source, "_aso03.csv"),
        as.data.frame(gsea_results_aso3[[gsea_source]]),
        quote     = FALSE,
        row.names = TRUE,
        col.names = TRUE,
        sep       = "\t"
    )
}
saveRDS(gsea_results_aso3, "GSEA_aso03.rds")
