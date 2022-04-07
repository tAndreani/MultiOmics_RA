# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: Heatmap with dots on GSEA adjusted p values
# Data: GSEA result on GO biological process pathways

# Load library ------------------------------------------------------------

library(data.table)
library(purrr)
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)

# Set plot parameter ------------------------------------------------------

theme_set(theme_bw())
theme_update(
  text = element_text(family = "Helvetica", size = 7),
  axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
  plot.title = element_text(size = 7),
  legend.title = element_text(size = 7),
  legend.text = element_text(size = 7),
  legend.key.size = unit(2, "mm"),
  legend.spacing = unit(0, "mm"),
  legend.position = "bottom",
  legend.direction = "vertical"
)

# Load data ---------------------------------------------------------------

# Define file pattern
csvPattern <- "fgsea_go_bp_*"

# List paths
csvPath <- list.files("./result/fgseaResults", full.names = TRUE, pattern = csvPattern)
csvPath

# Import multiple dataframes and bind into one
gseaRes <- map_df(csvPath, read.csv)
View(gseaRes)

length(unique(gseaRes$pathway)) # 4750 unique GO BP pathways
unique(gseaRes$Timepoint)

# Tidy data ---------------------------------------------------------------

# Pathway for manuscript
pathwayList <- c(
  "GO:0006096 glycolytic process",
  "GO:0006099 tricarboxylic acid cycle",
  "GO:0006119 oxidative phosphorylation",
  "GO:0022904 respiratory electron transport chain",
  "GO:0006809 nitric oxide biosynthetic process",
  "GO:0072593 reactive oxygen species metabolic process",
  "GO:0002224 toll-like receptor signaling pathway",
  "GO:0042119 neutrophil activation",
  "GO:0030225 macrophage differentiation",
  "GO:0050900 leukocyte migration",
  "GO:0072676 lymphocyte migration"
)

# Tidy table
dataTable <- gseaRes %>%
  filter(pathway %in% pathwayList) %>%
  mutate(
    pathwayName = str_to_sentence(substring(pathway, 12)),
    minusLog10padj = -log10(padj),
    Timepoint = str_replace_all(Timepoint, "_", " "),
    Direction = case_when(NES < 0 ~ "Depleted", TRUE ~ "Enriched"),
    Significance = case_when(padj < 0.05 ~ "Yes", TRUE ~ "No")
  )

# Order timepoints
dataTable$Timepoint <- factor(dataTable$Timepoint, levels = c(
  "2 days", "2 weeks", "3 weeks", "4 weeks",
  "6 weeks", "7 weeks", "8 weeks", "10 weeks"
))

# Rank pathway names
sorted_level <- str_to_sentence(substring(pathwayList, 12))

# Heatmap on -log10 padj --------------------------------------------------

# Check color palette
show_col(rev(brewer.pal(11, "RdYlBu")))

# Dot plot
p <- dataTable %>%
  ggplot(aes(x = Timepoint, y = pathwayName)) +
  geom_point(aes(fill = Direction, color = Significance, size = minusLog10padj),
    pch = 21, stroke = 1
  ) +
  labs(x = NULL, y = NULL, size = "-log10(FDR)", color = "FDR < 0.05", title = "Gene Ontology biological process") +
  scale_y_discrete(limits = rev(sorted_level)) +
  scale_fill_manual(values = c("#74ADD1", "#F46D43")) +
  scale_color_manual(values = c("gray", "#A50026")) +
  scale_size_continuous(limits = c(0, 10), breaks = c(2, 4, 6, 8))

ggsave("./plot/GSEA/heatmap_go_bp_v2.tif",
  plot = p, device = "tiff",
  units = "mm", width = 90, height = 110, dpi = 300
)
