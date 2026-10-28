library(ggplot2)
library(ggpubr)
library(png)
library(gridExtra)
library(dplyr)
library(grid)
library(stringr)
library(ggrepel)
library(Seurat)

cb_blue   <- "#0072B2"
cb_orange <- "#D55E00"

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

# Panel A
# Single-cell GSE112013 dataset
seurat_object <- readRDS('testis.rds')
Idents(object = seurat_object) <- seurat_object@meta.data$celltype
seurat_object$celltype <- factor(seurat_object$celltype, 
                        levels = c("SSCs", "Differentiating S'gonia", "Early primary S'cytes",
                          "Late primary S'cytes", "Round S'tids", "Elongated S'tids", "Sperm1",
                          "Sperm2", "Sperm3", "Macrophage", "Endothelial cells", "Myoid cells",
                          "Sertoli cells"))

panel_a <- DimPlot(seurat_object, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 2, repel = TRUE, raster=FALSE) +
panel_theme +
  guides(fill = guide_legend(ncol = 3)) +
  scale_color_viridis_d(option = "turbo", guide = guide_legend(reverse = TRUE), drop = FALSE) +
  labs(title = "", tag = "A") +
  theme(plot.title = element_blank(),
    axis.title.y = element_text(size = 6, angle = 90))

# Panel B expression across cell types

d_genes <- c('CTD-2587H24.5', 'DNAAF3')
d_expr  <- FetchData(seurat_object, vars = c(d_genes, "celltype"))
d_long  <- tidyr::pivot_longer(d_expr, cols = tidyselect::all_of(d_genes),
                               names_to = "gene", values_to = "expression")
d_long$gene         <- factor(d_long$gene, levels = d_genes,
                              labels = c("DNAAF3-AS1", "DNAAF3"))
d_exclude_celltypes <- c("Macrophage", "Endothelial cells", "Myoid cells", "Sertoli cells")
d_long              <- d_long[!(d_long$celltype %in% d_exclude_celltypes), ]
d_long$celltype     <- factor(d_long$celltype,
                              levels = setdiff(levels(seurat_object$celltype),
                                               d_exclude_celltypes))
d_long$celltype_num <- as.numeric(d_long$celltype)

gene_cols <- c("DNAAF3-AS1" = cb_orange,
  "DNAAF3" = cb_blue)

panel_b <- ggplot(d_long, aes(x = celltype, y = expression, fill = gene)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.2, alpha = 0.7,
               position = position_dodge(width = 0.8)) +
  geom_smooth(aes(x = celltype_num, color = gene, group = gene),
              method = "loess", se = FALSE, linewidth = 0.6) +
  scale_fill_manual(name = "", values = gene_cols) +
  scale_color_manual(name = "", values = gene_cols) +
  theme_bw() +
  panel_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5), 
    axis.title.y = element_text(size = 6, angle = 90)) +
  labs(title = "", x = NULL,
       y = "Expression (log-normalized)", tag = "B")

# Fraction of cells with non-zero CTD-2587H24.5 (DNAAF3-AS1) expression per cell type
frac_expr <- FetchData(seurat_object, vars = c("CTD-2587H24.5", "celltype"))
colnames(frac_expr)[1] <- "expression"
frac_nonzero <- frac_expr %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(
    n_cells          = dplyr::n(),
    n_nonzero        = sum(expression > 0),
    fraction_nonzero = sum(expression > 0) / dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(celltype)

write.csv(frac_nonzero, "Supplementary_TableX.csv", row.names = FALSE)

# Panel C: correlation of each himorna module score with CTD-2587H24.5 (DNAAF3-AS1) expression across all cells
gene_candidates <- c("CTD-2587H24.5", "DNAAF3-AS1", "ENSG00000267577")

target_expression <- function(obj) {
  hit <- gene_candidates[gene_candidates %in% rownames(obj)][1]
  if (is.na(hit)) stop("CTD-2587H24.5 / DNAAF3-AS1 not found in object features")
  as.numeric(FetchData(obj, vars = hit, assay = "RNA", layer = "data")[[1]])
}

correlate_scores <- function(rds_file, source_label) {
  obj  <- readRDS(rds_file)
  md   <- obj@meta.data
  expr <- target_expression(obj)

  keep <- !(md$celltype %in% d_exclude_celltypes)
  md   <- md[keep, ]
  expr <- expr[keep]

  score_cols <- c("H3K36me3 rho > 0.5"  = "add_module_score1",
                  "H3K36me3 rho < -0.5" = "add_module_score2")

  do.call(rbind, lapply(names(score_cols), function(rho_set) {
    score <- md[[score_cols[[rho_set]]]]
    ct    <- cor.test(score, expr, method = "spearman")
    rho   <- unname(ct$estimate)
    n     <- sum(complete.cases(score, expr))
    se_z  <- sqrt((1 + rho^2 / 2) / (n - 3))
    z     <- atanh(rho)
    data.frame(
      set    = paste(source_label, rho_set, sep = " | "),
      rho    = rho,
      ci_low = tanh(z - qnorm(0.975) * se_z),
      ci_high = tanh(z + qnorm(0.975) * se_z),
      pval   = ct$p.value,
      n      = n,
      stringsAsFactors = FALSE
    )
  }))
}

cor_df <- rbind(
  correlate_scores("testis_promoter.rds", "Promoter"),
  correlate_scores("testis_genebody.rds", "Gene body")
)

cor_df$set     <- factor(cor_df$set, levels = c("Promoter | H3K36me3 rho > 0.5",
  "Promoter | H3K36me3 rho < -0.5",
  "Gene body | H3K36me3 rho > 0.5",
  "Gene body | H3K36me3 rho < -0.5"))

cor_df$p_label <- ifelse(cor_df$pval < 2.2e-16, "p < 2.2e-16",
                         paste0("p = ", signif(cor_df$pval, 2)))

panel_c <- ggplot(cor_df, aes(x = set, y = rho, fill = set)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2, linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey40") +
  scale_fill_viridis_d(guide = "none") +
  scale_x_discrete(labels = function(x) str_wrap(x, 16)) +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  theme_bw() +
  panel_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5), 
    axis.title.y = element_text(size = 6, angle = 90)) +
  labs(title = "",
       x = NULL, y = "Spearman rho", tag = "C")


g <- arrangeGrob(panel_a, panel_b, panel_c,  ncol = 2, nrow = 2, layout_matrix= rbind(c(1, 1), c(2,3)))
ggsave(paste0("Figure4.png"), plot = g, width = 18, height = 20, units = "cm")



