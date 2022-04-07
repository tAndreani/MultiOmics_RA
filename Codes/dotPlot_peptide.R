# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: Heatmap with dots with peptides
# Data: FC and q-value results

# Load library ------------------------------------------------------------

library(data.table)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(stringr)
library(purrr)

# Set plot parameter ------------------------------------------------------

theme_set(theme_bw())
theme_update(
  text = element_text(family = "Helvetica", size = 7),
  axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
  plot.title = element_text(size = 7)
)

# Load data ---------------------------------------------------------------

# FC and q-value results
input <- fread("./data/IPA_Spotfire_Input/Metabolon_log2FC_qValue_notRounded_Spotfire_20210708.csv")
View(input)

unique(input$Timepoint)

# Tidy table --------------------------------------------------------------

# Extract significant peptides
peptideList <- input %>%
  filter(FDR < 0.05, abs(log2FoldChange) > log2(1.5), `Super Pathway` == "Peptide") %>%
  group_by(`Biochemical Name`) %>%
  summarise(n = n()) %>%
  filter(n > 1) %>%
  pull(`Biochemical Name`) 
peptideList

# Tidy table
dataTable <- input %>%
  filter(`Biochemical Name` %in% peptideList) %>%
  mutate(
    minusLog10padj = -log10(FDR),
    FDR_Cutoff = case_when(FDR < 0.05 & abs(log2FoldChange) > log2(1.5) ~ "Yes", TRUE ~ "No"),
    `Biochemical Name` = str_replace_all(`Biochemical Name`, "\\*", ""))
View(dataTable)

# Order timepoints
dataTable$Timepoint <- factor(dataTable$Timepoint, levels = c(
  "2 days", "2 weeks", "3 weeks", "4 weeks",
  "6 weeks", "7 weeks", "8 weeks", "10 weeks"
))

# Define sorted levels
sorted_level <- rev(unique(dataTable$`Biochemical Name`))

# Dot plot ----------------------------------------------------------------

p <- dataTable %>%
  ggplot(aes(x = Timepoint, y = `Biochemical Name`)) +
  geom_point(aes(fill = log2FoldChange, color = FDR_Cutoff, size = minusLog10padj),
    pch = 21, stroke = 1, show.legend = FALSE
  ) +
  labs(
    title = "Peptides significant at >= 2 time points",
    x = NULL, y = NULL,
    size = "-log10(FDR)",
    color = "Significant",
    fill = "log2(FoldChange)"
  ) +
  scale_y_discrete(limits = sorted_level) +
  scale_color_manual(values = c("gray", "#A50026")) +
  scale_fill_gradientn(
    colors = rev(brewer.pal(10, "Spectral")[c(2, 3, 4, 6, 8, 9, 10)]),
    breaks = c(-4, -2, 0, 2, 4),
    limits = c(-5.5, 5.5)
  ) +
  scale_size_continuous(limits = c(0, 10), breaks = c(2, 4, 6, 8))

ggsave("./plot/dotPlot/dotPlot_peptide_v2.tif",
  plot = p, device = "tiff",
  units = "mm", width = 80, height = 85, dpi = 300
)
