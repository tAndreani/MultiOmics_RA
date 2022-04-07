# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: Heatmap with dots on ROS related genes
# Data: Annotated DESeq2 result

# Load library ------------------------------------------------------------

library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)

# Set plot parameter ------------------------------------------------------

theme_set(theme_bw())
theme_update(
  text = element_text(family = "Helvetica", size = 8),
  axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
  axis.title.x = element_blank(),
  plot.title = element_text(size = 8),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 8),
  legend.key.size = unit(4, "mm"),
  legend.position = "none"
)

# Load data ---------------------------------------------------------------

input <- fread("./results/Annotated_deseq_wald_test_normoxia.csv")
View(input)

# Subset and tidy table ---------------------------------------------------

# Define gene list
geneList <- c(
  "MMP9", "IL6", "IL1B", "MPO", "NCF1", "NCF2", "NCF4",
  "HIF1A", "SLC2A3", "SLC2A6", "HK1", "HK2", "HK3", "ACOD1"
)

# Or: Gene list for selecting Taqman primers
geneList <- c(
  "ACOD1", "NCF1", "SHPK", "NFKB1", "RELA", "HIF1A", "SLC2A6", "HK2", "PFKL",
  "ALDOC", "ENO3", "PGK1", "SDHA", "SDHB", "SDHC", "SDHD"
)

# Tidy table
dataTable <- input %>%
  filter(symbol %in% geneList) %>%
  mutate(
    minusLog10padj = -log10(padj),
    FDR_Cutoff = case_when(padj < 0.05 & abs(log2FoldChange) > log2(1.5) ~ "Yes", TRUE ~ "No"),
    comparison = "LPS"
  )
View(dataTable)

# Define y axis order
sorted_level <- geneList

# Dot plot
basic_plot <- dataTable %>%
  ggplot(aes(x = comparison, y = symbol)) +
  geom_point(aes(fill = log2FoldChange, color = FDR_Cutoff, size = minusLog10padj),
    pch = 21, stroke = 1.5
  ) +
  labs(
    title = "LPS human macrophages",
    x = NULL, y = NULL,
    size = "-log10(FDR)",
    color = "Significant",
    fill = "log2(FoldChange)"
  )

## Note: The code below was changed to adapt to Taqman primer plotting
basic_plot <- dataTable %>%
  ggplot(aes(
    x = comparison, y = symbol,
    label = round(log2FoldChange, digits = 2)
  )) +
  geom_point(aes(fill = log2FoldChange, color = FDR_Cutoff, size = minusLog10padj),
    pch = 21, stroke = 1.5
  ) +
  geom_text(nudge_x = 0.2) +
  labs(
    title = "LPS human macrophages",
    subtitle = "Labeled with log2(FoldChange)",
    x = NULL, y = NULL,
    size = "-log10(FDR)",
    color = "Significant",
    fill = "log2(FoldChange)"
  )

p <- basic_plot +
  scale_y_discrete(limits = rev(sorted_level)) +
  scale_color_manual(values = c("gray", "#A50026")) +
  scale_fill_gradientn(
    colors = rev(brewer.pal(11, "Spectral")[c(1, 2, 3, 4, 6, 10)]),
    breaks = seq(-2, 10, 2),
    limits = c(-2, 11)
  )

ggsave("./plot/dotPlot_gene_macrophage_v3.tif",
  plot = p, device = "tiff",
  units = "mm", width = 45, height = 87.5, dpi = 300
)
