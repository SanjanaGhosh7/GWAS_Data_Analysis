library(ggplot2)
library(readr)

#Lollipop Plot for top 10 gene frequency

data <- read_csv("Top_10_genes.csv", show_col_types = FALSE)

ggplot(data, aes(x = reorder(Gene, Count), y = Count)) +
  geom_segment(aes(x = Gene, xend = Gene, y = 0, yend = Count), color = "black") +
  geom_point(color = "darkblue", size = 4) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 10 Genes with highest frequency associated with Multiple Traits",
       x = "Gene",
       y = "Association Count")

ggsave("gene_frequency_plot.pdf", width = 10, height = 8)


#Stacked bar plot for Top 10 genes  

library(dplyr)

df <- read_tsv("ChrID_and_mutations.tsv", show_col_types = FALSE)

df$SNP_Count <- sapply(strsplit(df$`SNP-Risk Allele`, ","), length)

gene_chr_snp <- df %>%
  group_by(Gene, `Chromosome ID`) %>%
  summarise(Total_SNPs = sum(SNP_Count, na.rm = TRUE)) %>%
  ungroup()

chrom_colors <- c(
  "2" = "#ff7f0e", "9" = "#d62728", "11" = "#2ca02c", "15" = "#9467bd",
  "17" = "#f7b6d2", "19" = "#17becf")


ggplot(gene_chr_snp, aes(x = reorder(Gene, -Total_SNPs), y = Total_SNPs, fill = as.factor(`Chromosome ID`))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = chrom_colors) +
  coord_flip() +
  theme_minimal() +
  labs(title = "SNP Count per Gene (Colored by Chromosome)",
       x = "Gene",
       y = "Number of SNPs",
       fill = "Chromosome ID")
ggsave("ChrID_and_mutations_plot.pdf", width = 8, height = 6)




#Log scale bar plot for Odds ratio

df <- read_csv("top_10_OR.csv", show_col_types = FALSE)


ggplot(df, aes(x = reorder(MAPPED_GENE, `OR or BETA`), y = `OR or BETA`)) +
  geom_bar(stat = "identity", fill = "dodgerblue") +
  scale_y_log10() +
  theme_minimal() +
  labs(title = "Top Genes by Odds Ratio (log scale)",
       x = "Gene",
       y = "Odds Ratio (log10 scale)")
ggsave("top_10_OR_plot.pdf", width = 8, height = 6)







