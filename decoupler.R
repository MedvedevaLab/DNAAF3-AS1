library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

# Bar plot of pathway activities
plot_pathways <- function(pathway_df, aso_id) {
    p_bar <- ggplot(pathway_df, aes(x = score, y = source, fill = score)) +
        geom_col() +
        geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
        scale_fill_gradient2(low = "#377EB8", mid = "whitesmoke",
                             high = "#E41A1C", midpoint = 0) +
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
plot_targets <- function(top_pathway, net, aso_id) {
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
        scale_color_manual(values = c(Concordant = "#E41A1C",
                                      Discordant = "#377EB8")) +
        theme_bw() +
        labs(x = paste0(top_pathway, " weight"),
             y = "Gene statistic (signed -log10 FDR)",
             color = "",
             title = paste0(top_pathway, " target genes — ", aso_id))

    tgt_file <- paste0("decoupler_", aso_id, "_", top_pathway, "_targets.png")
    ggsave(tgt_file, p_targets, width = 6, height = 5, dpi = 300)
    saveRDS(p_targets, paste0("decoupler_", aso_id, "_", top_pathway, "_targets.rds"))
}

run_aso <- function(deseq, net, aso_id) {
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

    plot_pathways(pathway_df, aso_id)

    bar_df <- acts %>%
      filter(statistic == "mlm") %>%
      arrange(score) %>%
      mutate(source = factor(source, levels = source))
    
    top_pathway <- bar_df$source[which.max(bar_df$score)]
    top_pathway <- as.character(top_pathway)
    plot_targets(top_pathway, net, aso_id)

    top_pathway <- bar_df$source[which.min(bar_df$score)]
    top_pathway <- as.character(top_pathway)
    plot_targets(top_pathway, net, aso_id)

}

# Load DE results 
deseq <- read.csv("DESeq2_genes_ASO_all.tsv", sep = "\t")

deseq <- deseq %>%
    filter(!is.na(geneSymbol), geneSymbol != "",
           !is.na(log2FC), !is.na(fdr)) %>%
    mutate(stat = sign(log2FC) * -log10(pmax(fdr, .Machine$double.xmin)))

net <- get_progeny(organism = "human", top = 500)

run_aso(deseq, net, "ASO_G0267577_01")
run_aso(deseq, net, "ASO_G0267577_03")
