library(ggplot2)
library(ggpubr)
library(ChIPseeker)
library(png)
library(gridExtra)
library(reshape2)
library(dplyr)
library(ChIPpeakAnno)
library(grid)
library(clusterProfiler)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BiocParallel)

BiocParallel::register(BiocParallel::MulticoreParam(workers = 1, progressbar = TRUE))

cb_blue   <- "#0072B2"  # blue
cb_orange <- "#D55E00"  # vermillion
cb_grey   <- "#999999"  # neutral grey

panel_theme <- theme(
    panel.border       = element_blank(),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    axis.title.x       = element_text(size = 6),
    axis.text.x        = element_text(size = 6, hjust = 0.5),
    legend.title       = element_text(size = 6),
    legend.key.size    = unit(0.2, "cm"),
    legend.text        = element_text(size = 5),
    axis.ticks.x       = element_blank(),
    axis.text.y        = element_text(size = 6),
    axis.title.y       = element_blank(),
    plot.title         = element_text(size = 6),
    plot.tag           = element_text(size = 6, face = "plain"),
    legend.spacing.y   = unit(2, "cm")
)

# Annotate genomic peaks from a data frame of coordinates
get_annotation <- function(peak_data, correlation_label, txdb) {
    coords <- peak_data[, c("chrom", "start", "end")]
    coords <- cbind(coords, "*")
    colnames(coords) <- c("seqnames", "start", "end", "strand")
    coords$start <- as.numeric(coords$start)
    coords$end   <- as.numeric(coords$end)
    peak_granges <- toGRanges(coords, format = "BED")
    peak_annotation <- annotatePeak(peak_granges, tssRegion = c(-1000, 100), TxDb = txdb, overlap = "all")
    return(peak_annotation)
}

make_ranked <- function(de) {
    de <- de[!is.na(de$log2FC) & !is.na(de$geneID), ]
    de$gene <- sub("\\..*$", "", de$geneID)
    de <- de[order(-abs(de$log2FC)), ]
    de <- de[!duplicated(de$gene), ]
    sort(setNames(de$log2FC, de$gene), decreasing = TRUE)
}

gsea_go_bp <- function(ranked) {
    gseGO(
        geneList      = ranked,
        OrgDb         = org.Hs.eg.db,
        keyType       = "ENSEMBL",
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        verbose       = FALSE
    )
}

# dotplot of GSEA terms (x = NES), ordered by NES, coloured by adjusted p-value
gsea_dotplot <- function(d, title, tag) {
    d$Description <- factor(d$Description, levels = d$Description[order(d$NES)])
    ggplot(d, aes(x = NES, y = Description, size = setSize, color = p.adjust)) +
        geom_point() +
        geom_vline(xintercept = 0, linetype = "dashed",
                   linewidth = 0.2, colour = "grey60") +
        theme_bw() + panel_theme +
        labs(title = title, tag = tag, y = NULL) + theme(legend.key.size = unit(0.6, "cm"))
}

gtf_file <- "gencode.v49.annotation.gtf"
txdb     <- makeTxDbFromGFF(gtf_file)

# Panel A: Gene tracks
# pyGenomeTracks --tracks plot.ini --region chr19:55158653-55178653 -o Figure1A.png --plotWidth 60 --dpi 300
gene_track_img  <- readPNG("Figure1A.png")
gene_track_grob <- rasterGrob(gene_track_img,
    interpolate = TRUE, width = unit(4, "in"), height = unit(1.5, "in"))

panel_a <- ggplot(geom = "blank") +
    annotation_custom(gene_track_grob) +
    theme_bw() +
    panel_theme +
    labs(tag = "A")

# Panel B: Expression in tissues
# https://www.gtexportal.org/home/downloads/adult-gtex%23qtl
df <- read.table('GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct', sep='\t', header=TRUE, skip=2,
                 check.names=FALSE, stringsAsFactors=FALSE)

result <- melt(as.data.frame(df[grepl("ENSG00000267577|ENSG00000167646", df$Name), ]))
result <- result[order(result$Description, result$value),]
result_show <- dplyr::filter(result, variable %in% tail(result[order(result$Description, result$value),])[,'variable'])
result_show$value <- log10(result_show$value+1)

panel_b <- ggplot(result_show, aes(x = variable , y = value , fill = Description )) +
  theme_bw() + scale_fill_manual(values = c("DNAAF3-AS1" = cb_orange, "DNAAF3" = cb_blue)) +
  geom_bar(position = "dodge", stat = "identity") + panel_theme +
  labs(title = "", fill = "", tag = "B") + coord_flip()

# Panel C: Number of significant DE genes per ASO in FANTOM6, highlight ASO_G0267577_01 and ASO_G0267577_03
all_deseq2_results <- read.csv("DESeq2_genes_ASO_all.tsv", sep = "\t")

highlight_asos <- c("ASO_G0267577_01", "ASO_G0267577_03")

sig_counts <- all_deseq2_results %>%
    dplyr::filter(fdr < 0.05) %>%
    dplyr::count(perturb_id, name = "n_sig") %>%
    dplyr::arrange(n_sig)
sig_counts$perturb_id <- factor(sig_counts$perturb_id, levels = sig_counts$perturb_id)
sig_counts$highlight  <- ifelse(sig_counts$perturb_id %in% highlight_asos,
                                as.character(sig_counts$perturb_id), "other")

x_labels <- ifelse(levels(sig_counts$perturb_id) %in% highlight_asos,
                   levels(sig_counts$perturb_id), "")

panel_c <- ggplot(sig_counts, aes(x = perturb_id, y = n_sig, fill = highlight)) +
    theme_bw() +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("ASO_G0267577_01" = cb_blue,
                                 "ASO_G0267577_03" = cb_orange,
                                 "other"           = cb_grey),
                      breaks = highlight_asos) +
    scale_x_discrete(labels = x_labels) +
    panel_theme +
    theme(axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5)) +
    labs(title = "",
         y = "Number of DE genes", x = NULL, fill = NULL, tag = "C")

# Panel D: PROGENY pathway activities in DE genes
panel_d = readRDS('decoupler_ASO_G0267577_01_pathway_activities.rds') + panel_theme  + labs(title = "ASO_01 PROGENy pathway activities", tag = "D")

# Panel E,F: GSEA GO:BP
deseq2_aso1     <- dplyr::filter(all_deseq2_results, perturb_id == "ASO_G0267577_01")
deseq2_aso1_sig <- deseq2_aso1[deseq2_aso1$fdr < 0.05, ]

write.table(
    file      = "SupplementaryTable1A.csv",
    deseq2_aso1_sig, 
    quote     = FALSE,
    row.names = TRUE,
    col.names = TRUE,
    sep       = "\t"
)

deseq2_aso3     <- dplyr::filter(all_deseq2_results, perturb_id == "ASO_G0267577_03")
deseq2_aso3_sig <- deseq2_aso3[deseq2_aso3$fdr < 0.05, ]

write.table(
    file      = "SupplementaryTable1B.csv",
    deseq2_aso3_sig, 
    quote     = FALSE,
    row.names = TRUE,
    col.names = TRUE,
    sep       = "\t"
)

# GSEA GO:BP on the ASO_01 log2FC ranking (all tested genes)
ranked_aso1 <- make_ranked(deseq2_aso1)
gsea_aso1   <- gsea_go_bp(ranked_aso1)
gsea_df     <- as.data.frame(gsea_aso1)

write.table(
    file      = "SupplementaryTable1E.csv",
    gsea_df,
    quote     = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    sep       = "\t"
)

# GSEA GO:BP on the ASO_01 log2FC ranking (all tested genes)
ranked_aso3 <- make_ranked(deseq2_aso3)
gsea_aso3   <- gsea_go_bp(ranked_aso3)
gsea_df     <- as.data.frame(gsea_aso3)

n_top   <- 10
pos      <- gsea_df[gsea_df$NES > 0, ]
neg      <- gsea_df[gsea_df$NES < 0, ]
top_pos  <- head(pos[order(pos$qvalue), ], n_top)
top_neg  <- head(neg[order(neg$qvalue), ], n_top)

write.table(
    file      = "SupplementaryTable1F.csv",
    gsea_df,
    quote     = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    sep       = "\t"
)

panel_e <- gsea_dotplot(top_pos, "NES > 0", "E")
panel_f <- gsea_dotplot(top_neg, "NES < 0", "F")

# Shared legend for panels E and F
ef_data     <- rbind(panel_e$data[, c("setSize", "p.adjust")],
                     panel_f$data[, c("setSize", "p.adjust")])
padj_range  <- range(ef_data$p.adjust)
size_range  <- range(ef_data$setSize)


shared_scales <- list(
    scale_color_viridis_c(name = "p.adjust", limits = padj_range,
                          option = "viridis", direction = -1,
                          guide = guide_colorbar(reverse = TRUE)),
    scale_size_continuous(name = "setSize", limits = size_range,
                          range = c(0.5, 3))
)
panel_e <- panel_e + shared_scales
panel_f <- panel_f + shared_scales

ef_legend <- ggpubr::get_legend(panel_e)
panel_ef  <- arrangeGrob(
    arrangeGrob(panel_e + theme(legend.position = "none"),
                panel_f + theme(legend.position = "none"),
                ncol = 2),
    ef_legend,
    ncol    = 1,
    heights = c(1, 0.2)
)

figure1_grob <- arrangeGrob(
    panel_a, panel_b, panel_c, panel_d, panel_ef,
    layout_matrix = rbind(
        c(1, 2),
        c(3, 4),
        c(5, 5)
    ),
    heights = c(1, 1, 1)
)

ggsave("Figure1.png", plot = figure1_grob, width = 18, height = 22, units = "cm")