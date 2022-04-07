# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: Heatmap with dots on ROS related genes
# Data: Annotated DESeq2 result

# Load library ------------------------------------------------------------

library(data.table)
library(purrr)
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
  legend.key.size = unit(4, "mm")
)

# Load data ---------------------------------------------------------------

# List annotated DESeq2 result csv file path
csvPath <- list.files("./result/Annotated_DESeq2_Wald_Test_Result", full.names = TRUE)
csvPath

# Import multiple dataframes and bind into one
input <- map_df(csvPath, read.csv)
View(input)

unique(input$Timepoint)

# Subset and tidy table ---------------------------------------------------

# Define gene list
geneList <- c(
  "MMP9", "IL6", "IL1B", "MPO", "NCF1", "NCF2", "NCF4",
  "HIF1A", "SLC2A3", "SLC2A6", "HK1", "HK2", "HK3", "ACOD1"
)

# Tidy table
dataTable <- input %>%
  mutate(
    symbol = str_to_upper(symbol),
    minusLog10padj = -log10(padj),
    Timepoint = str_replace_all(Timepoint, "_", " "),
    FDR_Cutoff = case_when(padj < 0.05 & abs(log2FoldChange) > log2(1.5) ~ "Yes", TRUE ~ "No")
  ) %>%
  filter(symbol %in% geneList)
View(dataTable)

# Order timepoints
dataTable$Timepoint <- factor(dataTable$Timepoint, levels = c(
  "2 days", "2 weeks", "3 weeks", "4 weeks",
  "6 weeks", "7 weeks", "8 weeks", "10 weeks"
))

# Define y axis order
sorted_level <- geneList

# Dot plot
basic_plot <- dataTable %>%
  ggplot(aes(x = Timepoint, y = symbol)) +
  geom_point(aes(fill = log2FoldChange, color = FDR_Cutoff, size = minusLog10padj),
    pch = 21, stroke = 1.5
  ) +
  labs(
    title = "CIA mouse paws",
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
  ) +
  scale_size_continuous(limits = c(0, 200), breaks = seq(5, 40, 10))

ggsave("./plot/GSEA/dotPlot_gene_paw_v3.tif",
  plot = p, device = "tiff",
  units = "mm", width = 90, height = 92, dpi = 300
)
