library(reshape2)
library(ggplot2)

df <- read.table('GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct', sep='\t', header=TRUE, skip=2,
                 check.names=FALSE, stringsAsFactors=FALSE)

result <- melt(as.data.frame(df[grepl("ENSG00000267577|ENSG00000167646", df$Name), ]))

result$value <- log10(result$value+1)
panel_a <- ggplot(result, aes(x = variable , y = value , fill = Description )) +
  theme_bw() +
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
  labs(title = "GTEx11", fill = "", tag = "") + coord_flip()

ggsave("Supplementary_Figure3.png", plot = panel_a, width = 18, height = 18, units = "cm")
