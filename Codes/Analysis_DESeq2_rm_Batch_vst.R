# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: PCA, batch-corrected and vst transformed expression data
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

library(BiocParallel)
register(MulticoreParam(4))

# Load dds object ---------------------------------------------------------

dds <- readRDS(file = "./data/paw_RNAseq/DESeq2_dds.rds")

# Pre-filtering dataset ---------------------------------------------------

nrow(dds) # Before filtering: 53646 genes

# Keep only rows that CPM values above 1 in at least 20% of samples
n <- ncol(dds) * 0.2
keep <- rowSums(fpm(dds, robust = FALSE) > 1) >= n
table(keep)

dds <- dds[keep, ]
nrow(dds) # After filtering: 16803 genes

# Transformation ----------------------------------------------------------

# Raw data
mat_raw <- assay(dds)

# Variance stabilizing transformation (VST)
vsd <- vst(dds, blind = FALSE)

## Plot SD of each gene against the mean
vst_before <- meanSdPlot(mat_raw, ranks = FALSE, bins = 100)$gg +
  ggtitle("Before VST")

vst_after <- meanSdPlot(assay(vsd), ranks = FALSE, bins = 100)$gg +
  ggtitle("After VST")

vst_fig1 <- grid.arrange(vst_before, vst_after, ncol = 2)

ggsave(
  filename = "./plot/DESeq2/meanSdPlot_vst.png", vst_fig1,
  width = 12, height = 6
)

## Plot count distribution
vst_before2 <- mat_raw[, 1:20] %>%
  as_tibble() %>%
  mutate(geneID = rownames(.)) %>%
  pivot_longer(-geneID, names_to = "name", values_to = "value") %>%
  ggplot(aes(x = name, y = value)) +
  geom_violin() +
  geom_point(size = .5, color = "gray") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = NULL, y = "Gene count", title = "Before VST")

vst_after2 <- assay(vsd)[, 1:20] %>%
  as_tibble() %>%
  mutate(geneID = rownames(.)) %>%
  pivot_longer(-geneID, names_to = "name", values_to = "value") %>%
  ggplot(aes(x = name, y = value)) +
  geom_violin() +
  geom_point(size = .5, color = "gray") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = NULL, y = "Gene count", title = "After VST")

vst_fig2 <- grid.arrange(vst_before2, vst_after2, ncol = 2)

ggsave(
  filename = "./plot/DESeq2/countDistribution_vst.png", vst_fig2,
  width = 8, height = 4
)

# Remove batch effect -----------------------------------------------------

# Before removing batch effect: PCA
pca_rm1 <- plotPCA(vsd, "Batch") + 
  ggtitle("Before batch effect removal") +
  labs(color = "Batch") +
  theme_bw()

# Remove batch effect: PCA
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$Batch)

pca_rm2 <- plotPCA(vsd, "Batch") + 
  ggtitle("After batch effect removal") +
  labs(color = "Batch") +
  theme_bw()

rm_fig <- grid.arrange(pca_rm1, pca_rm2, ncol = 2)

# Save grid figure
ggsave(filename = "./plot/DESeq2/PCA_batchRemoval.png", rm_fig, width = 10, height = 6)

# Save a single PCA plot after batch removal
ggsave(filename = "./plot/DESeq2/PCA_batch.png", pca_rm2, width = 6, height = 6)

# Sample distances --------------------------------------------------------

# Calculate Euclidean distance between samples
distsRL <- dist(t(assay(vsd)))
mat <- as.matrix(distsRL)

# Check for batch effect
colnames(mat) <- colData(vsd)$Batch
rownames(mat) <- paste0(colData(vsd)$Treatment, "_", colData(vsd)$Timepoint, "_", colData(vsd)$mouseID)

# Define color
hmcol <- colorRampPalette(brewer.pal(9, "Blues"))(255)

# Heatmap of the sample-to-sample distances, with batch on x axis
png(filename = "./plot/DESeq2/disPlot_after_batchRemoval.png")
heatmap.2(mat, trace = "none", col = rev(hmcol), margin = c(4,6))
dev.off()

# Save vsd object ---------------------------------------------------------

saveRDS(vsd, file = "./data/paw_RNAseq/DESeq2_preFilter_vst_batchRm_object.rds")
