# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: DESeqDataSet on macrophages under normoxia
# Data: Array Studio count and design table

# Load library ------------------------------------------------------------

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

library(DESeq2)

# Load data ---------------------------------------------------------------

# Design table
design_tb <- read.csv("./data/Corinne_design_table.csv")
View(design_tb)

# Count table
count_tb <- read.csv("./data/Corinne_ensembl_count.csv")
View(count_tb)

# Prepare colData ---------------------------------------------------------

colData <- design_tb %>%
  filter(Cell.type == "macrophage", Conditions == "normoxia", 
         Treatment %in% c("unstim", "LPS")) %>%
  select(observationid, Donor.ID, Treatment) %>%
  dplyr::rename(condition = Treatment, donor = Donor.ID) 

fwrite(colData, file = "./data/colData_macrophage_normoxia.csv")

# Prepare count table -----------------------------------------------------

countData <- count_tb[, colnames(count_tb) %in% c("geneid", colData$observationid)]
View(countData)

fwrite(countData, file = "./data/countData_macrophage_normoxia.csv")

# Tidy data ---------------------------------------------------------------

# Tidy colData
colData <- colData %>%
  column_to_rownames(var = "observationid") %>%
  mutate(donor = as.factor(str_replace_all(donor, " ", "_")),
         condition = as.factor(condition))

# Tidy countData
countData <- countData %>%
  column_to_rownames(var = "geneid")

# The columns of the count matrix and the rows of colData must be in the same order
print("Are columns of the count matrix equal to rows of colData?")

print(all(rownames(colData) %in% colnames(countData))) # TRUE

print(all(rownames(colData) == colnames(countData))) # TRUE

# Construct a DESeqDataSet ------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = round(countData),
  colData = colData,
  design = ~ donor + condition
)

# The first level of a factor be the reference level
dds$condition <- factor(dds$condition, levels = c("unstim", "LPS"))

## Quality check: Average read counts per sample
summary(round(colSums(assay(dds)) / 1e6, 1)) # Average: 36.40 million read counts per sample

# Save dds
saveRDS(dds, file = "./data/DESeq2_dds_macrophage_normoxia.rds")
