# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: PCA colored with scores and time points
# Data: DESeq2 object dds

# Load library ------------------------------------------------------------

library(DESeq2)
library(limma)
library(dplyr)
library(stringr)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(data.table)
library(gridExtra)

# Set plot parameter ------------------------------------------------------

theme_set(theme_bw())
theme_update(
  text = element_text(family = "Helvetica", size = 7),
  axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
  plot.title = element_text(size = 7),
  legend.title = element_text(size = 7),
  legend.key.size = unit(4, "mm")
)

# PCA: Score --------------------------------------------------------------

# Load vsd object
vsd <- readRDS("./data/paw_RNAseq/DESeq2_preFilter_vst_batchRm_object.rds")

# Select only the 500 genes showing the highest variance
ntop <- 500

Pvars <- rowVars(assay(vsd))
select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, length(Pvars)))]

# Run PCA
PCA <- prcomp(t(assay(vsd)[select, ]), scale = F)
percentVar <- round(100 * PCA$sdev^2 / sum(PCA$sdev^2), 1)

# Prepare data frame containing PCs
dataGG <- data.frame(
  PC1 = PCA$x[, 1], PC2 = PCA$x[, 2],
  PC3 = PCA$x[, 3], PC4 = PCA$x[, 4],
  mouseID = colData(vsd)$mouseID,
  Treatment = colData(vsd)$Treatment,
  Timepoint = colData(vsd)$Timepoint,
  Batch = colData(vsd)$Batch,
  Score = colData(vsd)$Score
)

# Create data frame to extract convex hull points
pca_hull <-
  dataGG %>%
  select(PC1, PC2, Treatment) %>%
  group_by(Treatment) %>%
  slice(chull(PC1, PC2))

# PCA with score
p_score <- dataGG %>% ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Score, shape = Treatment), size = 3) +
  geom_polygon(data = pca_hull, aes(fill = Treatment), alpha = 0.2, show.legend = TRUE) +
  labs(
    x = paste0("PC1: ", round(percentVar[1], 0), "% variance"),
    y = paste0("PC2: ", round(percentVar[2], 0), "% variance")
  ) +
  scale_shape_manual(values = c(1, 16)) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) +
  labs(color = "Paw score", title = "PCA on paw transcriptome (color by score)")

ggsave("./plot/DESeq2/PCA_score_polygon.tif",
  plot = p_score, device = "tiff",
  units = "mm", width = 90, height = 75, dpi = 300
)

# PCA: Timepoint ----------------------------------------------------------

# Load match table for timepoints
d <- fread("./data/metadata/matchTable_absolutDay_Timepoint.csv")
d

d$Timepoint <- str_replace_all(unique(d$Timepoint), " ", "_")

# Create data frame to extract convex hull points
pca_hull <-
  dataGG %>%
  filter(Treatment == "CIA") %>%
  left_join(d[,2:3], by = "Timepoint") %>%
  select(PC1, PC2, Timepoint, stage) %>%
  group_by(stage) %>%
  slice(chull(PC1, PC2))

# PCA with time
p_tp <- dataGG %>% ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Timepoint, shape = Treatment), size = 3) +
  geom_polygon(data = pca_hull, aes(fill = Treatment), alpha = 0.2, show.legend = TRUE) +
  labs(
    x = paste0("PC1: ", round(percentVar[1], 0), "% variance"),
    y = paste0("PC2: ", round(percentVar[2], 0), "% variance"),
    color = "Time", title = "PCA on paw transcriptome (color by time)"
  ) +
  scale_shape_manual(values = c(1, 16)) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) +
  scale_color_brewer(palette = "Paired", 
                     labels = str_replace_all(unique(dataGG$Timepoint), "_", " ")) 

ggsave("./plot/DESeq2/PCA_timepoint_polygon.tif",
       plot = p_tp, device = "tiff",
       units = "mm", width = 90, height = 75, dpi = 300
)

