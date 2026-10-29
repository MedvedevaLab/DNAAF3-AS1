library(ggplot2)
library(dplyr)
library(stringr)
library(GenomicRanges)
library(GenomicFeatures)
library(ggrepel)
library(grid)
library(png)
library(gridExtra)

cb_blue   <- "#0072B2"
cb_orange <- "#D55E00"
cb_green  <- "#009E73"
cb_grey   <- "#999999"
okabe_ito <- c(cb_blue, cb_orange, cb_green)

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

category_levels <- c("Genes with differential promoter usage",
                 "... also differentially expressed",
                 "... with a promoter in a coding region")

strip_version <- function(x) sub("\\..*$", "", x)

txdb <- makeTxDbFromGFF("gencode.v49.annotation.gtf")

gene_bodies <- genes(txdb)
names(gene_bodies) <- strip_version(names(gene_bodies))

# Exact CDS coordinates from the GENCODE annotation
cds_by_gene <- cdsBy(txdb, by = "gene")
cds_span    <- unlist(range(cds_by_gene))

standard_chroms <- paste0("chr", c(1:22, "X", "Y"))

cds <- data.frame(
  ensembl_gene_id = strip_version(names(cds_span)),
  cdsMin    = start(cds_span),
  cdsMax    = end(cds_span),
  cdsChr    = as.character(seqnames(cds_span)),
  cdsStrand = ifelse(as.character(strand(cds_span)) == "-", -1L, 1L),
  stringsAsFactors = FALSE
)
cds <- cds[cds$cdsChr %in% standard_chroms, ]
cds <- cds[!duplicated(cds$ensembl_gene_id), ]

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

fantom <- read.csv("DESeq2_genes_ASO_all.tsv", sep = "\t")

de_genes <- lapply(names(files), function(aso) {
  d <- dplyr::filter(fantom, perturb_id == aso, !is.na(fdr), fdr < 0.05)
  unique(strip_version(d$geneID))
})
names(de_genes) <- names(files)

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


a <- aso_data$ASO_G0267577_01$body
b <- aso_data$ASO_G0267577_03$body

intersect(a[a$inCDS,'geneID'],b[b$inCDS,'geneID'])

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

# Panel A
# Number of genes with differential promoter usage
panel_a <- ggplot(counts, aes(x = aso, y = n, fill = category)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = n),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 2) +
  scale_fill_manual(values = okabe_ito) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Number of genes", fill = NULL, tag = "A") +
  theme_classic() +
  panel_theme +
  guides(fill = guide_legend(ncol = 1))

# Panel B
# Metascape analysis for genes with differential promoter usage (fdr < 0.05) in cds
img <- readPNG('Figure2B.png')
g11 <- rasterGrob(img,
  interpolate=TRUE, width = unit(3.4,"in"), height=unit(1.4,"in"))

panel_b<- ggplot(geom="blank") +
  annotation_custom(g11) +
  theme_bw() +
  panel_theme +
  theme(
      axis.title.x = element_blank(),
      legend.key.size = unit(0.2, "cm"),
      legend.text = element_text(size = 5)) + labs(tag = "B")

# Panels C and D
# Promoter volcanoes
volcano_data <- function(aso) {
  d         <- aso_data[[aso]]
  df        <- d$df
  genes_sig <- unique(d$sig$geneID)
  in_cds  <- d$body$prmtrID[d$body$inCDS]
  df <- df[df$geneID %in% genes_sig, ]
  df$class <- ifelse((df$prmtrID %in% in_cds) & (df$prmtrPv < 0.05), "prmtrPv < 0.05, activated in coding region",
              ifelse(df$prmtrPv < 0.05,      "prmtrPv < 0.05",
                                             "Not significant"))
  df$class <- factor(df$class, levels = c("prmtrPv < 0.05, activated in coding region",
                                          "prmtrPv < 0.05",
                                          "Not significant"))
  df$aso <- aso
  df
}

n_label <- 20

volcano_panel <- function(aso, tag) {
  d <- volcano_data(aso)
  hits    <- d[d$class %in% c("prmtrPv < 0.05","prmtrPv < 0.05, activated in coding region"), ]
  by_coef <- head(hits[order(-abs(hits$prmtrCoef)), ], n_label)
  by_pv   <- head(hits[order(hits$prmtrPv), ],    n_label)
  lab     <- unique(rbind(by_coef, by_pv))

  ggplot(d, aes(x = prmtrCoef, y = -log10(prmtrPv))) +
    geom_point(data = subset(d, class == "Not significant"),
               aes(colour = class), size = 0.4) +
    geom_point(data = subset(d, class == "prmtrPv < 0.05"),
               aes(colour = class), size = 0.5) +
    geom_point(data = subset(d, class == "prmtrPv < 0.05, activated in coding region"),
               aes(colour = class), size = 0.8) +
    scale_colour_manual(values = c("prmtrPv < 0.05, activated in coding region" = cb_orange,
                                   "prmtrPv < 0.05"        = cb_blue,
                                   "Not significant"             = "grey80"),
                        drop = FALSE) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               linewidth = 0.2, colour = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed",
               linewidth = 0.2, colour = "grey40") +
    ggrepel::geom_text_repel(data = lab, aes(label = geneName),
                             size = 1.7, min.segment.length = 0,
                             segment.size = 0.2, max.overlaps = Inf,
                             box.padding = 0.2) +
    labs(x = "promoter-level logFC/gene-level logFC", y = expression(-log[10] ~ italic(p)),
         title = aso, tag = tag, colour = NULL) +
    theme_classic() +
    panel_theme +
    guides(colour = guide_legend(ncol = 1,
                                 override.aes = list(size = 1.5)))
  }

get_legend <- function(p) {
  g   <- ggplotGrob(p)
  idx <- which(vapply(g$grobs, function(x) x$name, character(1)) == "guide-box")
  g$grobs[[idx]]
}

panel_c <- volcano_panel("ASO_G0267577_01", "C")
panel_d <- volcano_panel("ASO_G0267577_03", "D")

legend_cd <- get_legend(panel_c)
panelC_nl <- panel_c + theme(legend.position = "none")
panelD_nl <- panel_d + theme(legend.position = "none")

panel_cd <- arrangeGrob(
  arrangeGrob(panelC_nl, panelD_nl, ncol = 2),
  legend_cd,
  ncol    = 1,
  heights = unit.c(unit(1, "null"), unit(1, "cm"))
)

g <- arrangeGrob(panel_a, panel_b, panel_cd, ncol = 2, nrow = 3, layout_matrix= rbind(c(1,2), c(3,3), c(3,3)))
ggsave("Figure2.png", plot = g, width = 18, height = 20, units = "cm")
