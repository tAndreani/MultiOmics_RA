# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: Differential expression analysis on different stages of RA vs. healthy control
# Data: DESeqDataSet dds

# Load library ------------------------------------------------------------

library(data.table)
library(tibble)
library(DESeq2)
library(dplyr)
library(stringr)
library(ggplot2)

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

# Run Wald test -----------------------------------------------------------

# Full model
design(dds)

# Run differential expression pipeline
dds <- DESeq(dds)

# List the coefficients
resultsNames(dds)

# Save dds
saveRDS(dds, file = "./data/DESeq2_preFilter_Wald_Test.rds")

# condition_Arthralgia_vs_Healthy -----------------------------------------

# Build results table
res <- results(dds, contrast = c("condition", "Arthralgia", "Healthy"))
res

# Export results
res_df <- as.data.frame(res) %>%
  rownames_to_column(var = "ensembl_gene_id") %>%
  mutate(condition = "Arthralgia")
View(res_df)

fwrite(res_df, file = "./results/deseq_wald_test_Arthralgia_vs_Healthy.csv")

rm(res)

# condition_UA_vs_Healthy -------------------------------------------------

# Build results table
res <- results(dds, contrast = c("condition", "UA", "Healthy"))
res

# Export results
res_df <- as.data.frame(res) %>%
  rownames_to_column(var = "ensembl_gene_id") %>%
  mutate(condition = "UA")
View(res_df)

fwrite(res_df, file = "./results/deseq_wald_test_UA_vs_Healthy.csv")

rm(res)

# condition_Early_RA_vs_Healthy -------------------------------------------

# Build results table
res <- results(dds, contrast = c("condition", "Early_RA", "Healthy"))
res

# Export results
res_df <- as.data.frame(res) %>%
  rownames_to_column(var = "ensembl_gene_id") %>%
  mutate(condition = "Early RA")
View(res_df)

fwrite(res_df, file = "./results/deseq_wald_test_Early_RA_vs_Healthy.csv")

# Optional: Plot pairwise normalized gene count ---------------------------

# MMP3
geneCounts <- plotCounts(dds,
                         gene = "ENSG00000149968.12",
                         intgroup = "condition",
                         returnData = TRUE)

ggplot(geneCounts, aes(x = condition, y = count)) +
  scale_y_log10() +
  geom_point(size = 3) +
  geom_line()
