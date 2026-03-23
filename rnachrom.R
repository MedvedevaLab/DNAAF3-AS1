df1 <- read.csv("rnachrom1.csv",sep = ',')
df1$experiment <- 'iMARGI1'

df2 <- read.csv("rnachrom2.csv")
df2$experiment <- 'iMARGI2'

df3 <- read.csv("rnachrom3.csv")
df3$experiment <- 'iMARGI3'

df4 <- read.csv("rnachrom4.csv")
df4$experiment <- 'iMARGI4'

df5 <- read.csv("rnachrom5.csv")
df5$experiment <- 'iMARGI5'

df7 <- read.csv("rnachrom7.csv")
df7$experiment <- 'GRIDseq1'

df8 <- read.csv("rnachrom8.csv")
df8$experiment <- 'GRIDseq2'


df_all <- do.call(rbind, list(df1, df2, df3, df4, df5, df7, df8))
counts <- as.data.frame(table(df_all[,'Target.gene']))
counts <- as.data.frame(table(df_all[,'Chromosome']))
head(counts[order(-counts$Freq),], 70)

library(GenomicRanges)

clusters_list <- list()
for (chr in unique(df_all[,'Chromosome'] ))
{
	chr_df <- df_all[df_all[,'Chromosome'] == chr,]
	gr <- GRanges(chr, IRanges(start = chr_df[,'Gene.start'],
	                                end = chr_df[,'Gene.end']))

	# reduce() merges overlapping/nearby intervals into clusters
	reduced <- reduce(gr, min.gapwidth = 10000)

	# map original intervals back to clusters
	hits <- findOverlaps(gr, reduced)
	gr$cluster <- subjectHits(hits)
	cluster_df <- as.data.frame(gr)

	cluster_df$cluster <- paste0(chr, "_", cluster_df$cluster)

	clusters_list[[chr]] <- cluster_df
}
result <- do.call(rbind,clusters_list )


for(cluster in unique(result$cluster))
{
	cluster_regions <- dplyr::filter(df_all, Gene.start %in% result[result$cluster == cluster,'start'])
	if(length(unique(cluster_regions$experiment))>1)
	{
		print(cluster)
	}
}
table(result$cluster)