# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: Tables of log2FC and BH adjusted p-values
# Data: DESeq2 object dds

# Reference: https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/08b_time_course_analyses.html
# Reference: https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html
# Bioconductor post: https://support.bioconductor.org/p/101002/#101015

# Load library ------------------------------------------------------------

library(data.table)
library(DESeq2)
library(dplyr)
library(stringr)
library(ggplot2)
library(gplots)
library(RColorBrewer)

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

# Run LRT test ------------------------------------------------------------

# Full model
design(dds)

# LRT test with reduced model
dds <- DESeq(dds, test = "LRT", reduced = ~Batch + Treatment + Timepoint)

# List the coefficients
resultsNames(dds)

# Cook’s distance: Indication of an outlier count
par(mar = c(8, 5, 2, 2))
boxplot(log10(assays(dds)[["cooks"]]), range = 0, las = 2)

# Wald tests for the log2 fold changes at individual time points ----------

## 2 days
res_2_days <- results(dds,
  name = "Treatment_CIA_vs_Ctrl",
  test = "Wald", alpha = 0.05
)
res_2_days
summary(res_2_days)

## How many genes with adjusted p value < 0.05 and |log2FC| > 1
sum(abs(res_2_days$log2FoldChange) > 1 & res_2_days$padj < 0.05, na.rm = TRUE)

## 2 weeks
res_2_weeks <- results(dds,
  contrast = list(c("Treatment_CIA_vs_Ctrl", "TreatmentCIA.Timepoint2_weeks")),
  test = "Wald", alpha = 0.05
)
res_2_weeks
summary(res_2_weeks)

## 3 weeks
res_3_weeks <- results(dds,
  contrast = list(c("Treatment_CIA_vs_Ctrl", "TreatmentCIA.Timepoint3_weeks")),
  test = "Wald", alpha = 0.05
)
res_3_weeks
summary(res_3_weeks)

## 4 weeks
res_4_weeks <- results(dds,
  contrast = list(c("Treatment_CIA_vs_Ctrl", "TreatmentCIA.Timepoint4_weeks")),
  test = "Wald", alpha = 0.05
)
res_4_weeks
summary(res_4_weeks)

## 6 weeks
res_6_weeks <- results(dds,
  contrast = list(c("Treatment_CIA_vs_Ctrl", "TreatmentCIA.Timepoint6_weeks")),
  test = "Wald", alpha = 0.05
)
summary(res_6_weeks)

## 7 weeks
res_7_weeks <- results(dds,
  contrast = list(c("Treatment_CIA_vs_Ctrl", "TreatmentCIA.Timepoint7_weeks")),
  test = "Wald", alpha = 0.05
)
summary(res_7_weeks)

## 8 weeks
res_8_weeks <- results(dds,
  contrast = list(c("Treatment_CIA_vs_Ctrl", "TreatmentCIA.Timepoint8_weeks")),
  test = "Wald", alpha = 0.05
)
summary(res_8_weeks)

## 10 weeks
res_10_weeks <- results(dds,
  contrast = list(c("Treatment_CIA_vs_Ctrl", "TreatmentCIA.Timepoint10_weeks")),
  test = "Wald", alpha = 0.05
)
summary(res_10_weeks)

# Loop to save FC and padj for individual time points ---------------------

## No filter for data here!

# List of results
res_list <- list(
  res_10_weeks, res_2_days, res_2_weeks, res_3_weeks,
  res_4_weeks, res_6_weeks, res_7_weeks, res_8_weeks
)

## Note: Manually save res_2_days first
tp <- "2_days"

df <- res_2_days %>%
  as.data.frame() %>%
  mutate(Timepoint = rep(tp, nrow(res_2_days)), gene_id = rownames(res_2_days)) %>%
  select(Timepoint, gene_id, log2FoldChange, pvalue, padj)
View(df)

fwrite(df, file = "./result/DESeq2_Wald_Test_Result/DESeq2_FC_padj_noFilter_2_days.csv")

rm(df)
rm(tp)

## Loop to save other time points
for (i in res_list[-2]) {
  tp <- i@elementMetadata$description[2] %>%
    word(., 5, sep = " ") %>%
    word(., 2, sep = "Timepoint")

  print(tp)

  df <- i %>%
    as.data.frame() %>%
    mutate(Timepoint = rep(tp, nrow(i)), gene_id = rownames(i)) %>%
    select(Timepoint, gene_id, log2FoldChange, pvalue, padj)

  fwrite(df, paste0("./result/DESeq2_Wald_Test_Result/DESeq2_FC_padj_noFilter_", tp, ".csv"))
}

# Save as a merged data frame for IPA -------------------------------------

# Rename column names
colnames(res_2_days) <- paste0(colnames(res_2_days), rep("_2_days", 6))
colnames(res_2_weeks) <- paste0(colnames(res_2_weeks), rep("_2_weeks", 6))
colnames(res_3_weeks) <- paste0(colnames(res_3_weeks), rep("_3_weeks", 6))
colnames(res_4_weeks) <- paste0(colnames(res_4_weeks), rep("_4_weeks", 6))
colnames(res_6_weeks) <- paste0(colnames(res_6_weeks), rep("_6_weeks", 6))
colnames(res_7_weeks) <- paste0(colnames(res_7_weeks), rep("_7_weeks", 6))
colnames(res_8_weeks) <- paste0(colnames(res_8_weeks), rep("_8_weeks", 6))
colnames(res_10_weeks) <- paste0(colnames(res_10_weeks), rep("_10_weeks", 6))

# Loop to bind columns
res_list <- list(
  res_10_weeks, res_2_days, res_2_weeks, res_3_weeks,
  res_4_weeks, res_6_weeks, res_7_weeks, res_8_weeks
)

df <- res_10_weeks

for (i in res_list[-1]) {
  print(colnames(i))

  df <- cbind(df, i)
}

dim(df) # 16803 x 48

colnames(df)

# Only save log2FoldChange, pvalue and padj
res_df <- df %>%
  as.data.frame() %>%
  select(starts_with("log2FoldChange_") | starts_with("pvalue_") | starts_with("padj_"))

colnames(res_df)

# Add a column with gene ID
res_df$gene_id <- rownames(res_df)

# Save for IPA
res_df <- res_df[, c(25, 1:24)]
View(res_df)

fwrite(res_df, "./result/IPA/DESeq2_FC_padj_noFilter_IPA.csv")

# Save dds after LRT ------------------------------------------------------

saveRDS(dds, file = "./data/paw_RNAseq/DESeq2_preFilter_LRT_object.rds")
