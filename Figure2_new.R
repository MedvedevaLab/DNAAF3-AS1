library(ggplot2)
library(dplyr)
library(stringr)
library(GenomicRanges)
library(GenomicFeatures)
library(biomaRt)
library(ggrepel)
library(ggpubr)
library(grid)

library(png)
library(gridExtra)
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
    axis.title.y       = element_text(size = 6),
    plot.title         = element_text(size = 6),
    plot.tag           = element_text(size = 6, face = "plain")
)

files <- c(ASO_G0267577_01 = "CTD2587H24_ASO_G0267577_01.tsv",
           ASO_G0267577_03 = "CTD2587H24_ASO_G0267577_03.tsv")

strip_version <- function(x) sub("\\..*$", "", x)


# ── Promoters inside gene bodies (promoters_in_genebody.R) ────────────────────
# Gene bodies (full locus, exons + introns), keyed by unversioned ensembl id
gtf_file <- "gencode.v49.annotation.gtf"
txdb        <- makeTxDbFromGFF(gtf_file)
gene_bodies <- genes(txdb)
names(gene_bodies) <- strip_version(names(gene_bodies))

# Exact CDS coordinates from biomaRt, cached between runs
all_genes <- unique(strip_version(
  unlist(lapply(files, function(f) read.csv(f, sep = "\t")$geneID))))
all_genes <- all_genes[grepl("^ENSG", all_genes)]

cds_cache <- "biomart_cds_coords.rds"
if (file.exists(cds_cache)) {
  cds_raw <- readRDS(cds_cache)
} else {
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

  # genomic_coding_start/end are per-exon, so query in chunks and rbind
  chunks <- split(all_genes, ceiling(seq_along(all_genes) / 500))
  cds_raw <- do.call(rbind, lapply(chunks, function(g) {
    getBM(
      attributes = c("ensembl_gene_id", "ensembl_transcript_id",
                     "chromosome_name", "strand",
                     "genomic_coding_start", "genomic_coding_end"),
      filters    = "ensembl_gene_id",
      values     = g,
      mart       = mart
    )
  }))
  saveRDS(cds_raw, cds_cache)
}

# drop non-coding transcripts / exons with no CDS
cds_raw <- cds_raw[!is.na(cds_raw$genomic_coding_start) &
                   !is.na(cds_raw$genomic_coding_end), ]
cds_raw <- cds_raw[!grepl("_|\\.", cds_raw$chromosome_name), ]

# per gene: full span of the coding region, over all coding transcripts
cds <- data.frame(
  ensembl_gene_id = names(tapply(cds_raw$genomic_coding_start,
                                 cds_raw$ensembl_gene_id, min)),
  cdsMin = as.numeric(tapply(cds_raw$genomic_coding_start,
                             cds_raw$ensembl_gene_id, min)),
  cdsMax = as.numeric(tapply(cds_raw$genomic_coding_end,
                             cds_raw$ensembl_gene_id, max)),
  cdsChr = paste0("chr", tapply(cds_raw$chromosome_name,
                                cds_raw$ensembl_gene_id, function(x) x[1])),
  cdsStrand = as.numeric(tapply(cds_raw$strand,
                                cds_raw$ensembl_gene_id, function(x) x[1])),
  stringsAsFactors = FALSE
)

# rows whose promoter falls inside the gene body of its own gene, annotated
# with the exact CDS span of that gene
promoters_in_genebody <- function(df) {
  coords <- as.data.frame(str_split_fixed(df$prmtrID, "_", 4))
  colnames(coords) <- c("seqnames", "start", "end", "strand")
  prmtr <- GRanges(seqnames = coords$seqnames,
                   ranges   = IRanges(start = as.numeric(coords$start) + 1,
                                      end   = as.numeric(coords$end)),
                   strand   = coords$strand)

  # NA when the gene id is absent from gencode (e.g. FANTOM-CAT CATG* models)
  df$geneID_clean <- strip_version(df$geneID)
  idx <- match(df$geneID_clean, names(gene_bodies))
  ok  <- !is.na(idx)

  bodyChr   <- rep(NA_character_, nrow(df))
  bodyStart <- rep(NA_integer_,   nrow(df))
  bodyEnd   <- rep(NA_integer_,   nrow(df))

  matched <- gene_bodies[idx[ok]]
  bodyChr[ok]   <- as.character(seqnames(matched))
  bodyStart[ok] <- start(matched)
  bodyEnd[ok]   <- end(matched)

  hit <- ok &
    as.character(seqnames(prmtr)) == bodyChr &
    start(prmtr) <= bodyEnd &
    end(prmtr)   >= bodyStart
  hit[is.na(hit)] <- FALSE

  # exact CDS span from biomaRt
  ci <- match(df$geneID_clean, cds$ensembl_gene_id)
  df$cdsChr <- cds$cdsChr[ci]
  df$cdsMin <- cds$cdsMin[ci]
  df$cdsMax <- cds$cdsMax[ci]

  # is the gene coding at all (present in the CDS table)?
  df$isCoding <- !is.na(ci)

  # does the promoter overlap the coding region span (ATG .. stop)?
  df$inCDS <- df$isCoding &
    as.character(seqnames(prmtr)) == df$cdsChr &
    start(prmtr) <= df$cdsMax &
    end(prmtr)   >= df$cdsMin
  df$inCDS[is.na(df$inCDS)] <- FALSE

  df$geneID_clean <- NULL
  df[hit, ]
}


# ── Differentially expressed genes per ASO (DESeq2) ────────────────────────────
fantom <- read.csv("DESeq2_genes_ASO_all.tsv", sep = "\t")

de_genes <- lapply(names(files), function(aso) {
  d <- dplyr::filter(fantom, perturb_id == aso, !is.na(fdr), fdr < 0.05)
  unique(strip_version(d$geneID))
})
names(de_genes) <- names(files)


# ── Counts per ASO ────────────────────────────────────────────────────────────
count_categories <- function(aso) {
  df <- read.csv(files[[aso]], sep = "\t")
  df <- df[!is.na(df$geneFDR) & !is.na(df$prmtrPv) & !is.na(df$prmtrCoef), ]

  # genes with a significant gene-level effect and >= 1 significant promoter
  sig <- df[df$geneFDR < 0.05 & df$prmtrPv < 0.05, ]
  genes_sig <- unique(strip_version(sig$geneID))

  # ... of those, the ones also differentially expressed in the same ASO
  genes_de <- intersect(genes_sig, de_genes[[aso]])

  # ... genes with >= 1 significant up promoter located inside the gene body
  # (promoters_in_genebody.R)
  up <- promoters_in_genebody(sig[sig$prmtrCoef > 0, ])
  genes_body <- unique(strip_version(up$geneID))

  data.frame(
    aso      = aso,
    category = c("Genes with significant promoters",
                 "... also differentially expressed",
                 "... with an intragenic up promoter"),
    n        = c(length(genes_sig), length(genes_de), length(genes_body)),
    stringsAsFactors = FALSE
  )
}

counts <- do.call(rbind, lapply(names(files), count_categories))
counts$category <- factor(counts$category,
                          levels = c("Genes with significant promoters",
                                     "... also differentially expressed",
                                     "... with an intragenic up promoter"))
print(counts)


# ── Panel A ───────────────────────────────────────────────────────────────────
panelA <- ggplot(counts, aes(x = aso, y = n, fill = category)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = n),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 2) +
  scale_fill_manual(values = c("Genes with significant promoters"   = "#4C72B0",
                               "... also differentially expressed"  = "#DD8452",
                               "... with an intragenic up promoter" = "#55A868")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Number of genes", fill = NULL, tag = "A") +
  theme_classic() +
  panel_theme +
  guides(fill = guide_legend(ncol = 1))


# ── Panels B and C: promoter volcanoes ────────────────────────────────────────
# All promoters of the genes with a significant promoter switch (geneFDR < 0.05
# and >= 1 promoter with prmtrPv < 0.05). Promoters that are activated inside
# the coding region — significant, prmtrCoef > 0, inside the gene body and
# overlapping the CDS span (promoters_in_genebody.R) — are highlighted.
volcano_data <- function(aso) {
  df <- read.csv(files[[aso]], sep = "\t")
  df <- df[!is.na(df$geneFDR) & !is.na(df$prmtrPv) & !is.na(df$prmtrCoef), ]

  sig       <- df[df$geneFDR < 0.05 & df$prmtrPv < 0.05, ]
  genes_sig <- unique(sig$geneID)

  # promoters activated in the coding region
  body    <- promoters_in_genebody(sig[sig$prmtrCoef > 0, ])
  in_cds  <- body$prmtrID[body$inCDS]

  df <- df[df$geneID %in% genes_sig, ]
  df$class <- ifelse(df$prmtrID %in% in_cds, "Activated in coding region",
              ifelse(df$prmtrPv < 0.05,      "Significant promoter",
                                             "Not significant"))
  df$class <- factor(df$class, levels = c("Activated in coding region",
                                          "Significant promoter",
                                          "Not significant"))
  df$aso <- aso
  df
}

n_label <- 10   # labels per ranking, per panel

volcano_panel <- function(aso, tag) {
  d <- volcano_data(aso)

  # label the coding-region activations with the largest coefficient and,
  # separately, the ones with the smallest p-value
  hits    <- d[d$class == "Activated in coding region", ]
  hits    <- hits[!duplicated(hits$geneName), ]
  by_coef <- head(hits[order(-hits$prmtrCoef), ], n_label)
  by_pv   <- head(hits[order(hits$prmtrPv), ],    n_label)
  lab     <- unique(rbind(by_coef, by_pv))

  ggplot(d, aes(x = prmtrCoef, y = -log10(prmtrPv))) +
    geom_point(data = subset(d, class == "Not significant"),
               colour = "grey80", size = 0.4) +
    geom_point(data = subset(d, class == "Significant promoter"),
               colour = "#4C72B0", size = 0.5) +
    geom_point(data = subset(d, class == "Activated in coding region"),
               colour = "#C44E52", size = 0.8) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               linewidth = 0.2, colour = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed",
               linewidth = 0.2, colour = "grey40") +
    ggrepel::geom_text_repel(data = lab, aes(label = geneName),
                             size = 1.7, min.segment.length = 0,
                             segment.size = 0.2, max.overlaps = Inf,
                             box.padding = 0.2) +
    labs(x = "Promoter coefficient", y = expression(-log[10] ~ italic(p)),
         title = aso, tag = tag) +
    theme_classic() +
    panel_theme +
    theme(legend.position = "none")
}

panelB <- volcano_panel("ASO_G0267577_01", "B")
panelC <- volcano_panel("ASO_G0267577_03", "C")






# Panel A
img <- readPNG('Figure2B.png')
g11 <- rasterGrob(img,
  interpolate=TRUE, width = unit(3,"in"), height=unit(1.2,"in"))

panelD<- ggplot(geom="blank") +
  annotation_custom(g11) +
  theme_bw() +
  panel_theme +
  theme(
      axis.title.x = element_blank(),
      legend.key.size = unit(0.2, "cm"),
      legend.text = element_text(size = 5)) + labs(tag = "D")




g <- arrangeGrob(panelA, panelB, panelC, panelD, ncol = 2, nrow = 3, layout_matrix= rbind(c(1,4), c(2,3), c(2,3)))
ggsave(paste0("Figure2_new.png"), plot = g, width = 18, height = 20, units = "cm")


library(clusterProfiler)
library(org.Hs.eg.db)



enrich_go_bp <- function(genes, universe) {
    enrichGO(
        gene          = genes,
        universe      = universe,
        OrgDb         = org.Hs.eg.db,
        keyType       = "ENSEMBL",
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.1,
        qvalueCutoff  = 0.1
    )
}



test_01 <- res[(res$inCDS) & (res$sample == 'CTD2587H24_ASO_G0267577_01'),]
test_03 <- res[(res$inCDS) & (res$sample == 'CTD2587H24_ASO_G0267577_03'),]

aso1_universe    <- unique(fantom$geneID)

ego_aso1_down <- enrich_go_bp(intersect(test_01$geneID, test_03$geneID), aso1_universe)

panel_e <- dotplot(ego_aso1_down, showCategory = 10) +
    theme_bw() + panel_theme +
    labs(title = "GO:BP – ASO_01 up DE genes", tag = "E")