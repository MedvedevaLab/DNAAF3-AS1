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

enrich_go_bp <- function(genes, universe) {
    enrichGO(
        gene          = genes,
        universe      = universe,
        OrgDb         = org.Hs.eg.db,
        keyType       = "ENSEMBL",
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.05
    )
}

gtf_file <- "gencode.v49.annotation.gtf"
txdb     <- makeTxDbFromGFF(gtf_file)

# Panel A: Gene tracks
# pyGenomeTracks --tracks plot.ini --region chr19:55158653-55178653 -o Figure1A.png --plotWidth 120 --dpi 130
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
  theme_bw() + scale_fill_manual(values = c("DNAAF3-AS1" = "#E41A1C", "DNAAF3" = "#377EB8")) +
  geom_bar(position = "dodge", stat = "identity") +
  theme(panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position = 'bottom',
  text = element_text(size = 6),
  axis.text.x = element_text(size = 5),
  axis.ticks.x = element_blank(),
  axis.title.y = element_text(size = 6),
  axis.text.y = element_text(size = 6),
  legend.key.size = unit(0.2, "cm"),
  legend.text = element_text(size = 5),
  plot.title = element_text(size = 6, hjust = 0.5),
  plot.tag = element_text(size = 8, face="plain")) +
  guides(fill=guide_legend(ncol = 2)) +
  labs(title = "GTEx11", fill = "", tag = "B") + coord_flip()

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
    scale_fill_manual(values = c("ASO_G0267577_01" = "#E41A1C",
                                 "ASO_G0267577_03" = "#377EB8",
                                 "other"           = "grey70"),
                      breaks = highlight_asos) +
    scale_x_discrete(labels = x_labels) +
    panel_theme +
    theme(axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5)) +
    labs(title = "Significant DE genes per ASO (FDR < 0.05)",
         y = "Number of DE genes", x = NULL, fill = NULL, tag = "C")

# Panel D: PROGENY pathway activities in DE genes

panel_d = readRDS('decoupler_ASO_G0267577_01_pathway_activities.rds') + panel_theme + labs(tag = "D")

# Panel E,F: GO:BP enrichment
deseq2_aso1     <- dplyr::filter(all_deseq2_results, perturb_id == "ASO_G0267577_01")
deseq2_aso1_sig <- deseq2_aso1[deseq2_aso1$fdr < 0.05, ]
deseq2_aso3     <- dplyr::filter(all_deseq2_results, perturb_id == "ASO_G0267577_03")
deseq2_aso3_sig <- deseq2_aso3[deseq2_aso3$fdr < 0.05, ]

deseq2_aso1_up   <- deseq2_aso1_sig$geneID[deseq2_aso1_sig$log2FC > 0]
deseq2_aso1_down <- deseq2_aso1_sig$geneID[deseq2_aso1_sig$log2FC < 0]
deseq2_aso3_up   <- deseq2_aso3_sig$geneID[deseq2_aso3_sig$log2FC > 0]
deseq2_aso3_down <- deseq2_aso3_sig$geneID[deseq2_aso3_sig$log2FC < 0]

aso1_universe    <- unique(deseq2_aso1$geneID)
aso3_universe    <- unique(deseq2_aso3$geneID)

ego_aso1_up   <- enrich_go_bp(deseq2_aso1_up,   aso1_universe)
ego_aso1_down <- enrich_go_bp(deseq2_aso1_down, aso1_universe)

panel_e <- dotplot(ego_aso1_up, showCategory = 10) +
    theme_bw() + panel_theme +
    labs(title = "GO:BP – ASO_01 up DE genes", tag = "E")

panel_f <- dotplot(ego_aso1_down, showCategory = 7) +
    theme_bw() + panel_theme +
    labs(title = "GO:BP – ASO_01 down DE genes", tag = "F")

# ── Shared legend for panels E and F ──────────────────────────────────────────
# The two dotplots have different colour (p.adjust) and size (Count) scales, so
# harmonise their limits to a common range, then draw a single shared legend.
ef_data     <- rbind(panel_e$data[, c("Count", "p.adjust")],
                     panel_f$data[, c("Count", "p.adjust")])
padj_range  <- range(ef_data$p.adjust)
count_range <- range(ef_data$Count)

shared_scales <- list(
    scale_color_gradient(name = "p.adjust", limits = padj_range,
                         low = "red", high = "blue",
                         guide = guide_colorbar(reverse = TRUE)),
    scale_size_continuous(name = "Count", limits = count_range,
                          range = c(0.5, 3))
)
panel_e <- panel_e + shared_scales
panel_f <- panel_f + shared_scales

# Extract one common legend and remove the per-panel legends.
ef_legend <- ggpubr::get_legend(panel_e)
panel_ef  <- arrangeGrob(
    arrangeGrob(panel_e + theme(legend.position = "none"),
                panel_f + theme(legend.position = "none"),
                ncol = 2),
    ef_legend,
    ncol    = 1,
    heights = c(1, 0.2)
)

# ── Assemble and save figure ───────────────────────────────────────────────────
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