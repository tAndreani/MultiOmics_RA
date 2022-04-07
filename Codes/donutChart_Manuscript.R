# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Experiment: Untargeted metabolomics on EDTA plasma from CIA mice
# This R script generates donut charts on metabolite groups

# Load library ------------------------------------------------------------

library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Load data ---------------------------------------------------------------

fc_qVal <- fread("./data/IPA_Spotfire_Input/Metabolon_log2FC_qValue_notRounded_Spotfire_20210708.csv")
View(fc_qVal)

colnames(fc_qVal)

# Set theme ---------------------------------------------------------------

theme_set(theme_void())
theme_update(
  text = element_text(family = "Helvetica", size = 7),
  axis.text.x = element_blank(),
  plot.title = element_text(size = 7),
  legend.title = element_blank(),
  legend.text = element_text(size = 7),
  legend.key.size = unit(4, "mm")
)

# Super pathway fractions in whole dataset --------------------------------

# Calculate fractions and define boundaries
df_fraction <- fc_qVal %>%
  group_by(`Super Pathway`) %>%
  summarise(count = n()) %>%
  mutate(fraction = count / sum(count)) %>%
  mutate(ymax = cumsum(fraction)) %>%
  mutate(ymin = c(0, head(ymax, n = -1))) %>%
  mutate(percentage = paste0(round(fraction * 100, digits = 1), "%")) %>%
  mutate(labelPosition = (ymax + ymin)/2)
df_fraction

# Donut chart for super pathway fractions
pathLabel <- unique(fc_qVal$`Super Pathway`)

p <- df_fraction %>%
  ggplot(aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = `Super Pathway`)) +
  geom_rect(color = "#999999") +
  coord_polar(theta = "y") +
  xlim(c(2, 4)) + # Donut thickness
  scale_fill_manual(values = brewer.pal(10, "Paired"), breaks = pathLabel) +
  labs(title = "Total number of metabolites: 929") +
  geom_text(x = 3.5, aes(y = labelPosition, label = percentage), size = 2.2)

ggsave("./plot/donutChart/donutChart_metabolite_class.tif",
       plot = p, device = "tiff", bg = "white",
       units = "mm", width = 100, height = 100, dpi = 300
)
