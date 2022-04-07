# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: Lineplot and barplot on proportion of significant metabolites
# Data: log2FC and qValue table from Metabolon

# Load library ------------------------------------------------------------

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)

# Set plot parameter ------------------------------------------------------

theme_set(theme_bw())
theme_update(
  text = element_text(family = "Helvetica", size = 7),
  axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
  axis.title.x = element_blank(),
  plot.title = element_text(size = 7),
  legend.title = element_blank(),
  legend.key.size = unit(4, "mm")
  #legend.position = "top"
)

# Load data ---------------------------------------------------------------

fc_qVal <- fread("./data/IPA_Spotfire_Input/Metabolon_log2FC_qValue_notRounded_Spotfire_20210708.csv")
View(fc_qVal)

colnames(fc_qVal)

# Count total metabolite per group ----------------------------------------

abs_total <- fc_qVal %>%
  group_by(`Super Pathway`, Timepoint) %>%
  summarise(totalCount = n())
View(abs_total)

# Count significant meatbolites -------------------------------------------

abs_gp_up <- fc_qVal %>%
  filter(FDR < 0.05, log2FoldChange > log2(1.5)) %>%
  group_by(`Super Pathway`, Timepoint) %>%
  summarise(upCount = n())

abs_gp_down <- fc_qVal %>%
  filter(FDR < 0.05, log2FoldChange < log2(1 / 1.5)) %>%
  group_by(`Super Pathway`, Timepoint) %>%
  summarise(downCount = n())

# Merge and calculate proportion ------------------------------------------

mergeTb <- abs_total %>%
  left_join(abs_gp_up, by = c("Super Pathway", "Timepoint")) %>%
  left_join(abs_gp_down, by = c("Super Pathway", "Timepoint")) %>%
  mutate(upCount = replace_na(upCount, 0), downCount = replace_na(downCount, 0)) %>%
  mutate(Up = upCount/totalCount * 100,
         Down = downCount/totalCount * 100) %>%
  pivot_longer(cols = 6:7, names_to = "Direction", values_to = "Proportion")
View(mergeTb)

# Lineplot ----------------------------------------------------------------

# Define orders
mergeTb$`Super Pathway` <- factor(mergeTb$`Super Pathway`,
  levels = c(
    "Amino Acid", "Peptide",
    "Carbohydrate", "Energy",
    "Lipid", "Nucleotide",
    "Cofactors and Vitamins", "Xenobiotics",
    "Partially Characterized Molecules", "N/A"
  )
)

mergeTb$Timepoint <- factor(mergeTb$Timepoint,
  levels = c(
    "2 days", "2 weeks", "3 weeks", "4 weeks",
    "6 weeks", "7 weeks", "8 weeks", "10 weeks"
  )
)

mergeTb$Direction <- factor(mergeTb$Direction, levels = c("Up", "Down"))

# Line plot
linePlot <- mergeTb %>%
  ggplot(aes(x = Timepoint, y = Proportion, group = Direction)) +
  geom_line(aes(color = Direction), size = 1, show.legend = FALSE) +
  geom_point(aes(color = Direction), size = 1.5) +
  labs(
    y = "% significant metabolites per metabolite group",
    caption = "FDR < 0.05 and |log2FoldChange| > 0.58"
  ) +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(vars(`Super Pathway`), nrow = 5)

ggsave("./plot/linePlot/linePlot_proportion_metabolite.tif",
  plot = linePlot, device = "tiff",
  units = "mm", width = 89, height = 150, dpi = 300
)

# Barplot -----------------------------------------------------------------

# Re-calculate proportion
mergeTb <- abs_total %>%
  left_join(abs_gp_up, by = c("Super Pathway", "Timepoint")) %>%
  left_join(abs_gp_down, by = c("Super Pathway", "Timepoint")) %>%
  mutate(upCount = replace_na(upCount, 0), downCount = replace_na(downCount, 0)) %>%
  mutate(Proportion = (upCount + downCount)/totalCount * 100)

# Define orders
mergeTb$Timepoint <- factor(mergeTb$Timepoint,
                            levels = c(
                              "2 days", "2 weeks", "3 weeks", "4 weeks",
                              "6 weeks", "7 weeks", "8 weeks", "10 weeks"
                            )
)

pathLabel <- unique(fc_qVal$`Super Pathway`)

# Normalized percentage
norPlot <- mergeTb %>%
  ggplot(aes(x = Timepoint, y = Proportion, fill = `Super Pathway`)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Normalized fraction of significant metabolites per metabolite group",
    y = "Normalized fraction",
    caption = "FDR < 0.05 and |log2FoldChange| > 0.58"
  ) +
  scale_fill_manual(values = brewer.pal(10, "Paired"), breaks = pathLabel) 

ggsave("./plot/barPlot/barPlot_proportion_metabolite_normalized.tif",
       plot = norPlot, device = "tiff",
       units = "mm", width = 100, height = 65, dpi = 300
)
