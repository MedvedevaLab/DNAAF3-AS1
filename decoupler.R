library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

cb_blue   <- "#0072B2"
cb_orange <- "#D55E00"

# Bar plot of pathway activities
plot_pathways <- function(pathway_df, aso_id) {
    p_bar <- ggplot(pathway_df, aes(x = score, y = source, fill = score)) +
        geom_col() +
        geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
        scale_fill_gradient2(low = cb_blue, mid = "grey95",
                             high = cb_orange, midpoint = 0) +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        labs(x = "Pathway activity (mlm score)", y = "",
             fill = "Activity",
             title = paste0("PROGENy pathway activities, ", aso_id))

    bar_file <- paste0("decoupler_", aso_id, "_pathway_activities.png")
    ggsave(bar_file, p_bar, width = 6, height = 5, dpi = 300)
    saveRDS(p_bar, paste0("decoupler_", aso_id, "_pathway_activities.rds"))
}

# Target-gene scatter for the most responsive pathways
plot_targets <- function(deg, top_pathway, net, aso_id) {
    df <- net %>%
        filter(source == top_pathway) %>%
        inner_join(rownames_to_column(deg, "target"), by = "target") %>%
        mutate(direction = ifelse(weight * stat > 0, "Concordant", "Discordant"))

    p_targets <- ggplot(df, aes(x = weight, y = stat, color = direction)) +
        geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
        geom_vline(xintercept = 0, color = "grey60", linetype = "dashed") +
        geom_point(alpha = 0.7) +
        ggrepel::geom_text_repel(
            data = df %>% arrange(desc(abs(weight * stat))) %>% head(15),
            aes(label = target), size = 2.5, max.overlaps = 20) +
        scale_color_manual(values = c(Concordant = cb_orange,
                                      Discordant = cb_blue)) +
        theme_bw() +
        labs(x = paste0(top_pathway, " weight"),
             y = "Gene statistic (signed -log10 FDR)",
             color = "",
             title = paste0(top_pathway, " target genes — ", aso_id))

    tgt_file <- paste0("decoupler_", aso_id, "_", top_pathway, "_targets.png")
    ggsave(tgt_file, p_targets, width = 6, height = 5, dpi = 300)
    saveRDS(p_targets, paste0("decoupler_", aso_id, "_", top_pathway, "_targets.rds"))
}

run_aso <- function(deseq, net, aso_id, save_plots = FALSE) {
    message("Running decoupleR for ", aso_id)
    deg <- deseq %>%
        filter(perturb_id == aso_id) %>%
        arrange(desc(abs(stat))) %>%
        distinct(geneSymbol, .keep_all = TRUE) %>%
        column_to_rownames("geneSymbol")

    mat <- as.matrix(deg[, "stat", drop = FALSE])
    colnames(mat) <- aso_id
    acts <- run_mlm(mat = mat, net = net,
                    .source = "source", .target = "target", .mor = "weight",
                    minsize = 5)

    pathway_df <- acts %>%
        filter(statistic == "mlm") %>%
        arrange(score) %>%
        mutate(source = factor(source, levels = source))

    saveRDS(pathway_df, paste0("decoupler_", aso_id, "_pathway_df.rds"))

    if (save_plots) {
        plot_pathways(pathway_df, aso_id)

        top_pathway <- pathway_df$source[which.max(pathway_df$score)]
        top_pathway <- as.character(top_pathway)
        plot_targets(deg, top_pathway, net, aso_id)

        top_pathway <- pathway_df$source[which.min(pathway_df$score)]
        top_pathway <- as.character(top_pathway)
        plot_targets(deg, top_pathway, net, aso_id)
    }

    invisible(pathway_df)
}

# Load DE results 
deseq <- read.csv("data//DESeq2_genes_ASO_all.tsv", sep = "\t")

deseq <- deseq %>%
    filter(!is.na(geneSymbol), geneSymbol != "",
           !is.na(log2FC), !is.na(fdr)) %>%
    mutate(stat = sign(log2FC) * -log10(pmax(fdr, .Machine$double.xmin)))

net <- get_progeny(organism = "human", top = 500)
plot_asos <- c("ASO_G0267577_01", "ASO_G0267577_03")
all_asos <- unique(deseq$perturb_id)

for (aso_id in all_asos) {
    run_aso(deseq, net, aso_id, save_plots = aso_id %in% plot_asos)
}
pathway_list <- setNames(
    lapply(all_asos, function(aso_id)
        readRDS(paste0("decoupler_", aso_id, "_pathway_df.rds"))),
    all_asos)

# Barplot of p53 pathway activity across all ASOs
p53_df <- bind_rows(lapply(names(pathway_list), function(aso_id) {
    pathway_list[[aso_id]] %>%
        mutate(source = as.character(source)) %>%
        filter(source == "p53") %>%
        transmute(aso_id = aso_id, score = score)
}))

# Sort ASOs by p53 activity and flag the two highlighted ASOs
highlight_cols <- c("ASO_G0267577_01" = cb_blue,
                    "ASO_G0267577_03" = cb_orange)

p53_df <- p53_df %>%
    arrange(score) %>%
    mutate(aso_id = factor(aso_id, levels = aso_id),
           highlight = ifelse(as.character(aso_id) %in% names(highlight_cols),
                              as.character(aso_id), "Other"))

fill_values <- c(highlight_cols, Other = "grey70")

p_p53 <- ggplot(p53_df, aes(x = score, y = aso_id, fill = highlight)) +
    geom_col() +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = fill_values,
                      breaks = names(highlight_cols)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(x = "p53 pathway activity (mlm score)", y = "",
         fill = "",
         title = "PROGENy p53 pathway activity across ASOs")

ggsave("decoupler_p53_activity_all_ASOs.png", p_p53,
       width = 6, height = max(4, 0.2 * nrow(p53_df)), dpi = 300)
saveRDS(p_p53, "decoupler_p53_activity_all_ASOs.rds")

# Quantile of the highlighted ASOs within the p53 score distribution
p53_ecdf <- ecdf(p53_df$score)
for (aso_id in names(highlight_cols)) {
    aso_score <- p53_df$score[as.character(p53_df$aso_id) == aso_id]
    q <- p53_ecdf(aso_score)
    message(sprintf("%s: p53 score = %.3f, quantile = %.3f (%.1f%%)",
                    aso_id, aso_score, q, 100 * q))
}
