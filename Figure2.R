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

category_levels <- c("Genes with promoter switch",
                 "... also differentially expressed",
                 "... with an intragenic promoter")

strip_version <- function(x) sub("\\..*$", "", x)

gtf_file <- "gencode.v49.annotation.gtf"
txdb        <- makeTxDbFromGFF(gtf_file)

# All genes from annotatiom
gene_bodies <- genes(txdb)
names(gene_bodies) <- strip_version(names(gene_bodies))

# Exact CDS coordinates from biomaRt
all_genes_list <- lapply(files, function(f) read.csv(f, sep = "\t")$geneID)
all_genes <- unique(strip_version(unlist(all_genes_list)))
all_genes <- all_genes[grepl("^ENSG", all_genes)]

cds_cache <- "biomart_cds_coords.rds"
if (file.exists(cds_cache)) {
  cds_raw <- readRDS(cds_cache)
} else {
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

  # query per-exon genomic_coding_start/end
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

promoters_in_genebody <- function(df) {
  coords <- as.data.frame(str_split_fixed(df$prmtrID, "_", 4))
  colnames(coords) <- c("seqnames", "start", "end", "strand")
  prmtr <- GRanges(seqnames = coords$seqnames,
                   ranges   = IRanges(start = as.numeric(coords$start) + 1,
                                      end   = as.numeric(coords$end)),
                   strand   = coords$strand)

  # Exclude FANTOM-CAT CATG*
  df$geneID_clean <- strip_version(df$geneID)
  idx <- match(df$geneID_clean, names(gene_bodies))
  ok  <- !is.na(idx)

  # Add CDS span
  ci <- match(df$geneID_clean, cds$ensembl_gene_id)
  df$cdsChr <- cds$cdsChr[ci]
  df$cdsMin <- cds$cdsMin[ci]
  df$cdsMax <- cds$cdsMax[ci]

  # Is the gene present in the CDS table
  df$isCoding <- !is.na(ci)

  # Does the promoter overlap the coding region
  df$inCDS <- df$isCoding &
    as.character(seqnames(prmtr)) == df$cdsChr &
    start(prmtr) <= df$cdsMax &
    end(prmtr)   >= df$cdsMin

  df$inCDS[is.na(df$inCDS)] <- FALSE

  df$geneID_clean <- NULL
  df[ok,]
}

#  Differentially expressed genes per ASO
fantom <- read.csv("DESeq2_genes_ASO_all.tsv", sep = "\t")

de_genes <- lapply(names(files), function(aso) {
  d <- dplyr::filter(fantom, perturb_id == aso, !is.na(fdr), fdr < 0.05)
  unique(strip_version(d$geneID))
})
names(de_genes) <- names(files)


#   aso_data[[aso]]$df   – Promoter table
#   aso_data[[aso]]$sig  – Significant promoter switch (geneFDR < 0.05)
#   aso_data[[aso]]$body – those promoters that fall inside the gene body, annotated with the CDS span (inCDS)
aso_data <- lapply(names(files), function(aso) {
  df <- read.csv(files[[aso]], sep = "\t")
  df <- df[!is.na(df$geneFDR) & !is.na(df$prmtrPv) & !is.na(df$prmtrCoef), ]

  sig  <- df[df$geneFDR < 0.05, ]
  body <- promoters_in_genebody(sig)
  body$aso    <- aso
  body$sample <- paste0("CTD2587H24_", aso) 

  list(df = df, sig = sig, body = body)
})
names(aso_data) <- names(files)


write.table(
    file      = "SupplementaryTable2C.csv",
    aso_data$ASO_G0267577_01$body,
    quote     = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    sep       = "\t"
)

write.table(
    file      = "SupplementaryTable2D.csv",
    aso_data$ASO_G0267577_03$body,
    quote     = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    sep       = "\t"
)

res <- do.call(rbind, lapply(aso_data, function(x) x$body))

df1 <- aso_data$ASO_G0267577_01$body[aso_data$ASO_G0267577_01$body$inCDS,]
df3 <- aso_data$ASO_G0267577_03$body[aso_data$ASO_G0267577_03$body$inCDS,]

intersect(head( df1[order(-abs(df1$prmtrCoef)),'geneName'], 30),head( df3[order(-abs(df3$prmtrCoef)),'geneName'],30))

count_categories <- function(aso) {
  d <- aso_data[[aso]]

  genes_sig <- unique(strip_version(d$sig$geneID))
  genes_de  <- intersect(genes_sig, de_genes[[aso]])
  genes_body <- unique(strip_version(d$body$geneID[d$body$inCDS]))

  data.frame(
    aso      = aso,
    category = category_levels,
    n        = c(length(genes_sig), length(genes_de), length(genes_body)),
    stringsAsFactors = FALSE
  )
}

counts <- do.call(rbind, lapply(names(files), count_categories))
counts$category <- factor(counts$category,
                          levels = category_levels)
print(counts)

# Panel A 
panelA <- ggplot(counts, aes(x = aso, y = n, fill = category)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = n),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 2) +
  scale_fill_manual(values = c("Genes with promoter switch" = "#4C72B0",
    "... also differentially expressed" = "#DD8452",
    "... with an intragenic promoter" = "#55A868")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Number of genes", fill = NULL, tag = "A") +
  theme_classic() +
  panel_theme +
  guides(fill = guide_legend(ncol = 1))


# Panel B: metascape analysis
img <- readPNG('Figure2B.png')
g11 <- rasterGrob(img,
  interpolate=TRUE, width = unit(3,"in"), height=unit(1.2,"in"))

panelB<- ggplot(geom="blank") +
  annotation_custom(g11) +
  theme_bw() +
  panel_theme +
  theme(
      axis.title.x = element_blank(),
      legend.key.size = unit(0.2, "cm"),
      legend.text = element_text(size = 5)) + labs(tag = "B")

# Panels C and D: promoter volcanoes 
volcano_data <- function(aso) {
  d         <- aso_data[[aso]]
  df        <- d$df
  genes_sig <- unique(d$sig$geneID)

  # promoters activated in the coding region
  in_cds  <- d$body$prmtrID[d$body$inCDS]

  df <- df[df$geneID %in% genes_sig, ]
  df$class <- ifelse((df$prmtrID %in% in_cds) & (df$prmtrPv < 0.05), "Activated in coding region",
              ifelse(df$prmtrPv < 0.05,      "Significant promoter",
                                             "Not significant"))
  df$class <- factor(df$class, levels = c("Activated in coding region",
                                          "Significant promoter",
                                          "Not significant"))
  df$aso <- aso
  df
}

n_label <- 20   # labels per ranking, per panel

volcano_panel <- function(aso, tag) {
  d <- volcano_data(aso)

  # label the coding-region activations with the largest coefficient and,
  # separately, the ones with the smallest p-value
  hits    <- d[d$class == "Activated in coding region", ]
  # hits    <- hits[!duplicated(hits$geneName), ]
  by_coef <- head(hits[order(-abs(hits$prmtrCoef)), ], n_label)
  by_pv   <- head(hits[order(hits$prmtrPv), ],    n_label)
  lab     <- unique(rbind(by_coef, by_pv))
  print(head(hits[order(-hits$prmtrCoef), ]))
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

panelC <- volcano_panel("ASO_G0267577_01", "C")
panelD <- volcano_panel("ASO_G0267577_03", "D")


g <- arrangeGrob(panelA, panelB, panelC, panelD, ncol = 2, nrow = 3, layout_matrix= rbind(c(1,2), c(3,4), c(3,4)))
ggsave(paste0("Figure2.png"), plot = g, width = 18, height = 20, units = "cm")
