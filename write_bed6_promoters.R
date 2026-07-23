# Coloured BED6 promoter tracks for pyGenomeTracks, as in dexseq_all.R, but
# built from the FANTOM differential promoter usage tables instead of a local
# DEXSeq run.
library(stringr)
library(dplyr)

# ── Input files (differential promoter usage results per ASO) ─────────────────
files <- c("01" = "CTD2587H24_ASO_G0267577_01.tsv",
           "03" = "CTD2587H24_ASO_G0267577_03.tsv")

# ── Genes of interest (file prefix -> Ensembl gene id) ────────────────────────
genes_of_interest <- c(
    smyd3 = "ENSG00000185420",
    ext1  = "ENSG00000182197",
    rfx23 = "ENSG00000087903",
    setd2 = "ENSG00000181555",
    kdm4a = "ENSG00000066135",
    kdm4b = "ENSG00000089094",
    kdm4c = "ENSG00000107077",
    ash1l = "ENSG00000116539",
    nsd3  = "ENSG00000147548",
    nsd1  = "ENSG00000165671",
    setd3 = "ENSG00000183576"
)

# ── Read one ASO table and expand prmtrID into BED coordinates ────────────────
# prmtrID is chr_start_end_strand with a BED-style (0-based) start, so the
# coordinates go into the BED6 unchanged.
read_dpu <- function(path) {
    df <- read.csv(path, sep = "\t")

    coords <- as.data.frame(str_split_fixed(df$prmtrID, "_", 4))
    colnames(coords) <- c("seqnames", "start", "end", "strand")

    df$seqnames <- coords$seqnames
    df$start    <- as.numeric(coords$start)
    df$end      <- as.numeric(coords$end)
    df$strand   <- coords$strand

    # unversioned gene id, and a promoter number per gene in coordinate order
    df$groupID <- str_split_fixed(df$geneID, "\\.", 2)[, 1]
    df <- df %>%
        group_by(groupID) %>%
        arrange(start, .by_group = TRUE) %>%
        mutate(featureID = sprintf("P%0.3d", row_number())) %>%
        ungroup() %>%
        as.data.frame()

    df
}

# ── Extract one gene's promoters and write a coloured BED6 track ──────────────
write_gene_bed6 <- function(results_table, gene_id, out_file) {
    gene_df <- results_table[results_table$groupID == gene_id,
        c("seqnames", "start", "end", "featureID",
          "prmtrPv", "strand", "prmtrCoef")]

    gene_df <- gene_df[!is.na(gene_df$prmtrPv) & !is.na(gene_df$prmtrCoef), ]
    gene_df <- gene_df[gene_df$prmtrPv < 0.05, ] # !!!!!
    gene_df$featureID <- ifelse(gene_df$prmtrPv < 0.05,
                                paste0(gene_df$featureID, "*"), gene_df$featureID)

    # score column drives the track colour: saturated for significant
    # promoters, faded otherwise; sign follows the ASO-vs-control coefficient
    gene_df$plot_color <- ifelse(gene_df$prmtrPv <  0.05 & gene_df$prmtrCoef > 0,  1,
                          ifelse(gene_df$prmtrPv <  0.05 & gene_df$prmtrCoef < 0, -1,
                          ifelse(gene_df$prmtrPv >= 0.05 & gene_df$prmtrCoef > 0,  0.3,
                          ifelse(gene_df$prmtrPv >= 0.05 & gene_df$prmtrCoef < 0, -0.3, NA))))

    gene_df <- gene_df[, c("seqnames", "start", "end",
                           "featureID", "plot_color", "strand")]

    write.table(
        file      = out_file,
        gene_df,
        quote     = FALSE,
        row.names = FALSE,
        col.names = FALSE,
        sep       = "\t"
    )
}

# ── Write BED6 tracks for every gene of interest for one ASO ──────────────────
write_all_genes <- function(results_table, aso_suffix) {
    for (gene_prefix in names(genes_of_interest)) {
        write_gene_bed6(
            results_table,
            genes_of_interest[[gene_prefix]],
            sprintf("%s_aso_%s.bed6", gene_prefix, aso_suffix)
        )
    }
}

for (aso_suffix in names(files)) {
    write_all_genes(read_dpu(files[[aso_suffix]]), aso_suffix)
}
