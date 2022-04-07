# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: Heatmap with dots with energy and related cofactors
# Data: FC and q-value results

# Load library ------------------------------------------------------------

library(data.table)
library(dplyr)
library(stringr)
library(purrr)

library(ggplot2)
library(RColorBrewer)

# Set plot parameter ------------------------------------------------------

theme_set(theme_bw())
theme_update(
  text = element_text(family = "Helvetica", size = 7),
  axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
  plot.title = element_text(size = 7),
  legend.title = element_text(size = 7),
  legend.text = element_text(size = 7),
  legend.key.size = unit(2.5, "mm"),
  legend.spacing = unit(0, "mm"),
  legend.position = "right",
  legend.direction = "vertical"
)

# Load data ---------------------------------------------------------------

# FC and q-value results
input <- fread("./data/IPA_Spotfire_Input/Metabolon_log2FC_qValue_notRounded_Spotfire_20210708.csv")
View(input)

unique(input$Timepoint)

# Tidy table --------------------------------------------------------------

# Extract significant energy, cofactors and vitamins
carbList <- input %>%
  filter(
    FDR < 0.05, abs(log2FoldChange) > log2(1.5),
    `Super Pathway` %in% c("Energy", "Cofactors and Vitamins")
  ) %>%
  group_by(`Biochemical Name`) %>%
  summarise(n = n()) %>%
  filter(n > 1) %>%
  pull(`Biochemical Name`)

carbList

# Define metabolites to include
sorted_metabolite <- c(
  "nicotinamide", "1-methylnicotinamide", "N1-Methyl-2-pyridone-5-carboxamide",
  "nicotinamide N-oxide", "nicotinate ribonucleoside", "trigonelline (N'-methylnicotinate)", "quinolinate",
  "itaconate", "isocitrate", "malate", "tricarballylate"
)

# Tidy table
dataTable <- input %>%
  filter(`Biochemical Name` %in% sorted_metabolite) %>%
  mutate(
    minusLog10padj = -log10(FDR),
    FDR_Cutoff = case_when(FDR < 0.05 & abs(log2FoldChange) > log2(1.5) ~ "Yes", TRUE ~ "No"),
    `Biochemical Name` = str_replace_all(`Biochemical Name`, "\\*", "")
  )
View(dataTable)

# Order timepoints
dataTable$Timepoint <- factor(dataTable$Timepoint, levels = c(
  "2 days", "2 weeks", "3 weeks", "4 weeks",
  "6 weeks", "7 weeks", "8 weeks", "10 weeks"
))

# Dot plot ----------------------------------------------------------------

p <- dataTable %>%
  ggplot(aes(x = Timepoint, y = `Biochemical Name`)) +
  geom_point(aes(fill = log2FoldChange, color = FDR_Cutoff, size = minusLog10padj),
    pch = 21, stroke = 1
  ) +
  labs(
    title = "Energy and related cofactors significant at >= 2 time points",
    x = NULL, y = NULL,
    size = "-log10(FDR)",
    color = "Significant",
    fill = "log2(FoldChange)"
  ) +
  scale_y_discrete(limits = rev(sorted_metabolite)) +
  scale_color_manual(values = c("gray", "#A50026")) +
  scale_fill_gradientn(
    colors = rev(brewer.pal(10, "Spectral")[c(2, 3, 4, 6, 8, 9, 10)]),
    breaks = c(-4, -2, 0, 2, 4),
    limits = c(-5.5, 5.5)
  ) +
  scale_size_continuous(limits = c(0, 10), breaks = c(2, 4, 6, 8))

ggsave("./plot/dotPlot/dotPlot_energy_cofactors.tif",
  plot = p, device = "tiff",
  units = "mm", width = 112, height = 72, dpi = 300
)
