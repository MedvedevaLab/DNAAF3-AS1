library(ggplot2)
library(ggpubr)
library(png)
library(gridExtra)
library(reshape2)
library(dplyr)
library(grid)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GenomicFeatures)
library(BiocParallel)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

BiocParallel::register(BiocParallel::MulticoreParam(workers = 1, progressbar = TRUE))

cb_blue   <- "#0072B2"
cb_orange <- "#D55E00"
cb_grey   <- "#999999"

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

# strip Ensembl version
strip_version <- function(ids) sub("\\..*$", "", ids)

# GO:BP over-representation analysis on a gene set against a universe
enrich_go_bp <- function(genes, universe) {
    enrichGO(
        gene          = strip_version(genes),
        universe      = strip_version(universe),
        OrgDb         = org.Hs.eg.db,
        keyType       = "ENSEMBL",
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.05
    )
}

# dotplot of GO terms ordered by gene ratio, coloured by adjusted p-value
go_dotplot <- function(d, title, tag) {
    d$GeneRatio <- sapply(strsplit(as.character(d$GeneRatio), "/"),
                          function(x) as.numeric(x[1]) / as.numeric(x[2]))
    wrapped <- vapply(as.character(d$Description),
                      function(s) paste(strwrap(s, width = 40), collapse = "\n"),
                      character(1))
    d$Description <- factor(wrapped, levels = wrapped[order(d$GeneRatio)])
    ggplot(d, aes(x = GeneRatio, y = Description, size = Count, color = p.adjust)) +
        geom_point() +
        theme_bw() + panel_theme +
        labs(title = title, tag = tag, y = NULL) +
        theme(legend.key.size = unit(0.6, "cm"),
              plot.tag.position = c(0, 1),
              plot.tag = element_text(size = 6, face = "plain", hjust = 0, vjust = 1))
}

narrow_panel <- function(g, panel_width = unit(2.2, "cm")) {
    gt <- ggplot2::ggplotGrob(g)
    panel_cols <- unique(gt$layout$l[grepl("^panel", gt$layout$name)])
    gt$widths[panel_cols] <- panel_width
    gt
}

gtf_file <- "gencode.v49.annotation.gtf"
txdb     <- makeTxDbFromGFF(gtf_file)

# Panel A: Gene tracks
# pyGenomeTracks --tracks plot_new.ini --trackLabelHAlign right --region chr19:55158653-55181000 -o Figure1A.png --plotWidth 45 --dpi 300
gene_track_img  <- readPNG("Figure1A.png")
gene_track_grob <- rasterGrob(gene_track_img,
    interpolate = TRUE, width = unit(3.8, "in"), height = unit(1.5, "in"))

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
result_show$Description = factor(result_show$Description, c("DNAAF3-AS1", "DNAAF3"))
panel_b <- ggplot(result_show, aes(x = variable , y = value , fill = Description )) +
  theme_bw() + scale_fill_manual(values = c("DNAAF3-AS1" = cb_orange, "DNAAF3" = cb_blue)) +
  geom_bar(position = position_dodge2(reverse = TRUE), stat = "identity") + panel_theme +
  theme(axis.title.x = element_text(size = 6)) +
  labs(title = "", y = "log10(TPM+1)", fill = "", tag = "B") + coord_flip()

# Panel C: Number of significant DE genes per ASO in FANTOM6, highlight ASO_G0267577_01 and ASO_G0267577_03
# https://fantom.gsc.riken.jp/6/datafiles/Core_FANTOM6/RELEASE_latest/analysis/DEGs/01_combined/DESeq2_genes_ASO_all.tsv.bz2

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
    theme(axis.title.x       = element_text(size = 6),
        axis.title.y       = element_text(size = 5, angle = 90)) +
    theme(axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5)) +
    labs(title = "",
        x = "ASO",
        y = "Number of DE genes", fill = NULL, tag = "C")

# Panel D: PROGENY pathway activities in DE genes
panel_d = readRDS('decoupler_ASO_G0267577_01_pathway_activities.rds') + panel_theme  + labs(title = "ASO_01 PROGENy pathway activities", tag = "D")

# Panel E,F: GO:BP over-representation, up- and down-regulated genes separately
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

# GO:BP over-representation on the ASO_01 significant genes (up/down vs all tested)
aso1_universe <- unique(deseq2_aso1$geneID)
ego_aso1_up   <- as.data.frame(enrich_go_bp(deseq2_aso1_sig$geneID[deseq2_aso1_sig$log2FC > 0], aso1_universe))
ego_aso1_up$de   <- "up"
ego_aso1_down <- as.data.frame(enrich_go_bp(deseq2_aso1_sig$geneID[deseq2_aso1_sig$log2FC < 0], aso1_universe))
ego_aso1_down$de <- "down"


# GO:BP over-representation on the ASO_03 significant genes (up/down vs all tested)
aso3_universe <- unique(deseq2_aso3$geneID)
ego_aso3_up   <- as.data.frame(enrich_go_bp(deseq2_aso3_sig$geneID[deseq2_aso3_sig$log2FC > 0], aso3_universe))
ego_aso3_up$de   <- "up"
ego_aso3_down <- as.data.frame(enrich_go_bp(deseq2_aso3_sig$geneID[deseq2_aso3_sig$log2FC < 0], aso3_universe))
ego_aso3_down$de <- "down"

n_top    <- 10
top_up   <- head(ego_aso3_up[order(ego_aso3_up$qvalue), ], n_top)
top_down <- head(ego_aso3_down[order(ego_aso3_down$qvalue), ], n_top)

panel_e <- go_dotplot(top_up, "", "E")
panel_f <- go_dotplot(top_down, "", "F")

# Shared legend for panels E and F
ef_data     <- rbind(panel_e$data[, c("Count", "p.adjust")],
                     panel_f$data[, c("Count", "p.adjust")])
padj_range  <- range(ef_data$p.adjust)
size_range  <- range(ef_data$Count)


shared_scales <- list(
    scale_color_viridis_c(name = "p.adjust", limits = padj_range,
                          option = "viridis", direction = -1,
                          guide = guide_colorbar(reverse = TRUE)),
    scale_size_continuous(name = "Count", limits = size_range,
                          range = c(0.5, 3))
)
panel_e <- panel_e + shared_scales + theme(axis.text.x = element_text(size = 5))

panel_f <- panel_f + shared_scales + theme(axis.text.x = element_text(size = 5))

ef_legend <- ggpubr::get_legend(panel_e)
panel_ef  <- arrangeGrob(
    arrangeGrob(narrow_panel(panel_e + theme(legend.position = "none")),
                narrow_panel(panel_f + theme(legend.position = "none")),
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
    heights = c(1, 1, 1.2)
)

ggsave("Figure1.png", plot = figure1_grob, width = 18, height = 22, units = "cm")