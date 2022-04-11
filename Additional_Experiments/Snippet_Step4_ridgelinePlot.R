# Introduction ------------------------------------------------------------

# This R snippet visualizes marker intensities
# Data: Expression data output from sce object
# Configuration: >= 8 Cores, >= 64 Gb RAM
# Output: Ridgeline plots and density plots facet by marker, color by condition

# Load packages -----------------------------------------------------------

install.packages("styler")

library(tidyverse)
library(ggridges)
library(data.table)

# Distributions of marker intensities: FlowSOM cell populations --------

# Load expression data of flowSOM clusters
setDTthreads(threads = 8)
df <- fread("./clustering/clusterExpression_wo_6weeks_meta24_20210131.csv", header = T, nThread = 8)
glimpse(df)

# Change to long format
df <- df %>%
  pivot_longer(
    cols = 1:12,
    names_to = "antigen",
    values_to = "expression",
    values_drop_na = TRUE
  )

df$cluster_id <- as.factor(df$cluster_id)

# Divide the df into 3 parts for plotting
v1 <- c("CD45", "CD3", "CD19", "CD11b")
v2 <- c("CD45R", "CD138", "CD4", "CD8")
v3 <- c("CD62L", "CD44", "CXCR5", "PD-1")

dfAntigen <- data.frame(v1, v2, v3)
dfAntigen

# Ridgeline plot
for (i in seq(1:3)) {
  print(paste0("Start ridge plot on markers: ", unlist(dfAntigen[i])))

  ridgePlot <- df %>%
    filter(antigen %in% unlist(dfAntigen[i])) %>%
    ggplot(aes(x = expression, y = cluster_id, fill = condition, color = condition)) +
    geom_density_ridges(scale = 4, alpha = 0.2) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_ridges() +
    facet_grid(. ~ antigen, scales = "free")

  plotName <- paste0("./plot/clusteringPlot_wo_6weeks_typeMarkers/meta12_markerIntensity_v", i, "_20210202.pdf")

  # Save plot as pdf
  ggsave(plotName, plot = ridgePlot, width = 40, height = 24, units = "cm")

  print(paste0(plotName, " is done"))
}

# Distributions of marker intensities: Annotated clusters -----------------

# List cluster expression csv files
exprFile <- list.files(path = "./clustering", pattern = "markerExpression_", full = TRUE)
exprFile

# Loop to create ridge plot on different timepoints
for (i in seq(1, length(exprFile))) {

  # Extract timepoint
  tp <- str_split_fixed(exprFile[i], "_", 5)[, 4]

  # Load data frame
  df <- fread(exprFile[i], header = T, nThread = 8)

  # Change to long format
  df <- df %>%
    select(sample_id, condition, patient_id, cluster_id, CD44, CD62L, CXCR5, `PD-1`) %>%
    pivot_longer(
      cols = 5:8,
      names_to = "antigen",
      values_to = "expression",
      values_drop_na = TRUE
    ) %>%
    filter(cluster_id %in% c("CD4 T cells", "CD8 T cells", "B cells"))
  
  # Rank
  df$condition <- factor(df$condition, levels = c("Ctrl", "CIA"))

  # Ridgeline plot
  ridgePlot <- ggplot(df, aes(x = expression, y = cluster_id, fill = condition, color = condition)) +
    geom_density_ridges(scale = 1, alpha = 0.2) +
    labs(title = tp) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    coord_cartesian(clip = "off") +
    theme_ridges() +
    facet_grid(cols = vars(antigen), scales = "free_x")

  # Save plot
  pName <- str_c("./plot/clusteringPlot_wo_6weeks_typeMarkers/", tp, "_markerIntensity_20210202.png")

  ggsave(pName, plot = ridgePlot, width = 20, height = 20, units = "cm")

  print(paste0(pName, " is done."))

  rm(df)
  rm(ridgePlot)
  rm(pName)
}
