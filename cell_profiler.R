library(ggplot2)
library(reshape2)

df <- read.csv('cell_profiler.csv', sep = '\t', dec = ',')
df <- df[df$KD.pubID == 'CTD-2587H24.5', ]
df <- df[df$perturb_id %in% c("ASO_G0267577_01", "ASO_G0267577_03"),]

params <- c(
  "GrowthRate",
  "Area",
  "Compactness",
  "Eccentricity",
  "Extent",
  "FormFactor",
  "MajorAxisLength",
  "MaxFeretDiameter",
  "MaximumRadius",
  "MedianRadius",
  "MinFeretDiameter",
  "MinorAxisLength",
  "Perimeter",
  "Solidity"
)

mean_cols <- paste0(params, "_mean")
fdr_cols  <- paste0(params, "_fdr")

# means -> long
mean_long <- melt(df[, c("perturb_id", mean_cols)],
                  id.vars = "perturb_id",
                  variable.name = "parameter",
                  value.name = "mean")
mean_long$parameter <- sub("_mean$", "", mean_long$parameter)

# fdr -> long
fdr_long <- melt(df[, c("perturb_id", fdr_cols)],
                 id.vars = "perturb_id",
                 variable.name = "parameter",
                 value.name = "fdr")
fdr_long$parameter <- sub("_fdr$", "", fdr_long$parameter)

# combine
plot_df <- merge(mean_long, fdr_long, by = c("perturb_id", "parameter"))

# keep column order as requested
plot_df$parameter <- factor(plot_df$parameter, levels = params)

# significance stars
plot_df$star <- ifelse(!is.na(plot_df$fdr) & plot_df$fdr < 0.05, "*", "")

p <- ggplot(plot_df, aes(x = parameter, y = perturb_id, fill = mean)) +
  geom_tile(color = "white") +
  geom_text(aes(label = star), size = 6, vjust = 0.75) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1) +
  labs(x = NULL, y = NULL, fill = "mean") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)
