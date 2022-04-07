# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: DESeqDataSet
# Data: CIA paw RNAseq countData and colData

# Load library ------------------------------------------------------------

library(DESeq2)
library(stringr)

library(BiocParallel)
register(MulticoreParam(4))

# Load data ---------------------------------------------------------------

# Count matrix
countData <- read.csv("./data/paw_RNAseq/combined_countFile.csv", row.names = 1)

# Table of sample information
colData <- read.csv("./data/metadata/DESeq2_colData.csv", row.names = 1, stringsAsFactors = TRUE)

# Tidy data ---------------------------------------------------------------

# Tidy up factor names in Timepoint
colData$Timepoint <- str_replace_all(colData$Timepoint, " ", "_")

# Set all columns as factors in colData
colData$Timepoint <- as.factor(colData$Timepoint)
colData$Batch <- as.factor(colData$Batch)
colData$Score <- as.factor(colData$Score)

# The columns of the count matrix and the rows of colData must be in the same order
print("Are columns of the count matrix equal to rows of colData?")

print(all(rownames(colData) %in% colnames(countData))) # TRUE

countData <- countData[, rownames(colData)]

print(all(rownames(colData) == colnames(countData))) # TRUE

# Construct a DESeqDataSet ------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = colData,
  design = ~ Batch + Treatment + Timepoint + Treatment:Timepoint
)

# The first level of a factor be the reference level
dds$Treatment <- factor(dds$Treatment, levels = c("Ctrl", "CIA"))
dds$Timepoint <- factor(dds$Timepoint, levels = c(
  "2_days", "2_weeks", "3_weeks", "4_weeks",
  "6_weeks", "7_weeks", "8_weeks", "10_weeks"
))

## Quality check: Average read counts per sample
summary(round(colSums(assay(dds)) / 1e6, 1)) # Average: 8.858 million read counts per sample

# Save dds
saveRDS(dds, file = "./data/paw_RNAseq/DESeq2_dds.rds")
