# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: PCA on vst transformed expression data
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

library(vsn)
library(hexbin) # Need to be installed manually

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

dds <- readRDS(file = "./data/DESeq2_dds_macrophage_normoxia.rds")

# Pre-filtering dataset ---------------------------------------------------

nrow(dds) # Before filtering: 58825 genes

# Keep only rows that CPM values above 1 in at least 50% of samples
n <- ncol(dds) * 0.5
keep <- rowSums(fpm(dds, robust = FALSE) > 1) >= n
table(keep)

dds <- dds[keep, ]
nrow(dds) # After filtering: 13311 genes

# Transformation ----------------------------------------------------------

# Variance stabilizing transformation (VST)
vsd <- vst(dds, blind = FALSE)

# Save vsd object
saveRDS(vsd, file = "./data/DESeq2_preFilter_vst_object_normoxia.rds")

# PCA: Condition ----------------------------------------------------------

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
  Treatment = colData(vsd)$condition,
  Donor = colData(vsd)$donor
) %>% mutate(
  Donor = str_replace_all(Donor, "donor_", "Donor "),
  Treatment = str_replace_all(Treatment, "unstim", "Unstimulated")
)
head(dataGG)

# Create data frame to extract convex hull points
pca_hull <-
  dataGG %>%
  select(PC1, PC2, Treatment) %>%
  group_by(Treatment) %>%
  slice(chull(PC1, PC2))

# PCA plot
pcaPlot <- dataGG %>% ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Treatment, shape = Donor), size = 3, alpha = 0.8) +
  geom_polygon(data = pca_hull, aes(fill = Treatment), alpha = 0.2, show.legend = TRUE) +
  labs(
    title = "PCA on human macrophage transcriptome",
    color = NULL, fill = NULL, shape = NULL,
    x = paste0("PC1: ", round(percentVar[1], 0), "% variance"),
    y = paste0("PC2: ", round(percentVar[2], 0), "% variance")
  )

ggsave("./plot/PCA_macrophage_treatment.tif",
  plot = pcaPlot, device = "tiff",
  units = "mm", width = 90, height = 65, dpi = 300
)
