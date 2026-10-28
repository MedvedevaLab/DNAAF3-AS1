library(stringr)
library(dplyr)

files <- c("01" = "CTD2587H24_ASO_G0267577_01.tsv",
           "03" = "CTD2587H24_ASO_G0267577_03.tsv")

genes_of_interest <- c(
    mdm2  = "ENSG00000135679",
    nr2f2 = "ENSG00000185551",
    fastk = "ENSG00000164896",
    cfdp1 = "ENSG00000153774"
)

read_dpu <- function(path) {
    df <- read.csv(path, sep = "\t")

    coords <- as.data.frame(str_split_fixed(df$prmtrID, "_", 4))
    colnames(coords) <- c("seqnames", "start", "end", "strand")

    df$seqnames <- coords$seqnames
    df$start    <- as.numeric(coords$start)
    df$end      <- as.numeric(coords$end)
    df$strand   <- coords$strand
    df$groupID <- str_split_fixed(df$geneID, "\\.", 2)[, 1]
    df <- df %>%
        group_by(groupID) %>%
        arrange(start, .by_group = TRUE) %>%
        mutate(featureID = sprintf("P%0.3d", row_number())) %>%
        ungroup() %>%
        as.data.frame()
    df
}

write_gene_bed6 <- function(results_table, gene_id, out_file) {
    gene_df <- results_table[results_table$groupID == gene_id,
        c("seqnames", "start", "end", "featureID",
          "prmtrPv", "strand", "prmtrCoef")]
    gene_df <- gene_df[!is.na(gene_df$prmtrPv) & !is.na(gene_df$prmtrCoef), ]
    gene_df$featureID <- ifelse(gene_df$prmtrPv < 0.05,
                                paste0(gene_df$featureID, "*"), gene_df$featureID)
    gene_df$plot_color <- ifelse(gene_df$prmtrPv <  0.05 & gene_df$prmtrCoef < 0,  1,
                          ifelse(gene_df$prmtrPv <  0.05 & gene_df$prmtrCoef > 0, -1,
                          ifelse(gene_df$prmtrPv >= 0.05 & gene_df$prmtrCoef < 0,  0.3,
                          ifelse(gene_df$prmtrPv >= 0.05 & gene_df$prmtrCoef > 0, -0.3, NA))))
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

# pyGenomeTracks --tracks plot_mdm2.ini  --region chr12:68805000-68848000  -o MDM2.png  --plotWidth 15 --dpi 300
# pyGenomeTracks --tracks plot_nr2f2.ini --region chr15:96323000-96342500  -o NR2F2.png --plotWidth 15 --dpi 300
# pyGenomeTracks --tracks plot_fastk.ini --region chr7:151074000-151083000 -o FASTK.png --plotWidth 15 --dpi 300
# pyGenomeTracks --tracks plot_cfdp1.ini --region chr16:75291000-75436000  -o CFDP1.png --plotWidth 15 --dpi 300