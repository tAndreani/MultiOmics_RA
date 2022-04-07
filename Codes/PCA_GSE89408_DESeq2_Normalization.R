# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: PCA on vst transformed expression data
# Data: DESeq2 object dds

# Load library ------------------------------------------------------------

library(data.table)
library(DESeq2)
library(limma)
library(dplyr)
library(stringr)

library(ggplot2)
library(gplots)
library(RColorBrewer)

library(pheatmap)

library(gridExtra)
library(vsn)
library(hexbin) # Need to be installed manually

library(BiocParallel)
register(MulticoreParam(4))

# Set plot parameter ------------------------------------------------------

theme_set(theme_bw())
theme_update(
  text = element_text(family = "Helvetica", size = 7),
  axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
  plot.title = element_text(size = 7),
  legend.title = element_text(size = 7),
  legend.key.size = unit(4, "mm")
)

# Load dds object ---------------------------------------------------------

dds <- readRDS(file = "./data/DESeq2_dds_early_disease.rds")

# Pre-filtering dataset ---------------------------------------------------

nrow(dds) # Before filtering: 60699 genes

# Keep only rows that CPM values above 1 in at least 5% of samples
n <- ncol(dds) * 0.05
keep <- rowSums(fpm(dds, robust = FALSE) > 1) >= n
table(keep)

dds <- dds[keep, ]
nrow(dds) # After filtering: 23053 genes

# Transformation ----------------------------------------------------------

# Variance stabilizing transformation (VST)
vsd <- vst(dds, blind = FALSE)

# Save vsd object
saveRDS(vsd, file = "./data/DESeq2_preFilter_vst_object.rds")

# PCA: Condition ----------------------------------------------------------

# Load vsd object
vsd <- readRDS("./data/DESeq2_preFilter_vst_object.rds")

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
  condition = colData(vsd)$condition
)
head(dataGG)

# Create data frame to extract convex hull points
pca_hull <-
  dataGG %>%
  select(PC1, PC2, condition) %>%
  group_by(condition) %>%
  slice(chull(PC1, PC2))

# PCA plot
pcaPlot <- dataGG %>% ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = condition), size = 2, alpha = 0.8) +
  geom_polygon(data = pca_hull, aes(fill = condition), alpha = 0.2, show.legend = TRUE) +
  labs(
    color = NULL, fill = NULL, title = "PCA on human synovium transcriptome",
    x = paste0("PC1: ", round(percentVar[1], 0), "% variance"),
    y = paste0("PC2: ", round(percentVar[2], 0), "% variance")
  ) +
  scale_color_brewer(
    palette = "Set2",
    breaks = c("Healthy", "Arthralgia", "UA", "Early_RA"),
    labels = c("Healthy", "Arthralgia", "UA", "Early RA")
  ) +
  scale_fill_brewer(
    palette = "Set2",
    breaks = c("Healthy", "Arthralgia", "UA", "Early_RA"),
    labels = c("Healthy", "Arthralgia", "UA", "Early RA")
  )

ggsave("./plot/PCA_condition.tif",
  plot = pcaPlot, device = "tiff",
  units = "mm", width = 90, height = 70, dpi = 300
)
