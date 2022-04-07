# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: Differential expression analysis on macrophages under normoxia
# Data: DESeqDataSet dds

# Load library ------------------------------------------------------------

library(data.table)
library(tibble)
library(DESeq2)
library(dplyr)
library(stringr)
library(ggplot2)

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

# Run Wald test -----------------------------------------------------------

# Full model
design(dds)

# Run differential expression pipeline
dds <- DESeq(dds)

# List the coefficients
resultsNames(dds)

# Save dds
saveRDS(dds, file = "./data/DESeq2_preFilter_Wald_Test_normoxia.rds")

# Save results LPS vs unstim ----------------------------------------------

# Build results table
res <- results(dds, contrast = c("condition", "LPS", "unstim"))
res

# Export results
res_df <- as.data.frame(res) %>%
  rownames_to_column(var = "ensembl_gene_id")
View(res_df)

fwrite(res_df, file = "./results/deseq_wald_test_normoxia.csv")

# Optional: Plot pairwise normalized gene count ---------------------------

# Acod1
geneCounts <- plotCounts(dds,
  gene = "ENSG00000102794",
  intgroup = c("condition", "donor"),
  returnData = TRUE
)

ggplot(geneCounts, aes(x = condition, y = count, color = donor, group = donor)) +
  scale_y_log10() +
  geom_point(size = 3) +
  geom_line()
