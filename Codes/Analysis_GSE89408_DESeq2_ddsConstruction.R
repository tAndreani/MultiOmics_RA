# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: DESeqDataSet on healthy, arthralgia, UA and early RA samples
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
design_tb <- read.csv("./data/GSE89408_meta_data.csv")
View(design_tb)

# Count table
count_tb <- read.csv("./data/GSE89408_count_table.csv")
View(count_tb)

# Prepare colData ---------------------------------------------------------

colData <- design_tb %>%
  mutate(tissue = str_split_fixed(SubjectID, "_", 2)[, 2]) %>%
  mutate(condition = case_when(
    str_detect(tissue, "healthy tissue") ~ "Healthy",
    str_detect(tissue, "osteoarthritis tissue") ~ "Osteoarthritis",
    str_detect(tissue, "arthralgia tissue") ~ "Arthralgia",
    str_detect(tissue, "undifferentiated arthritis tissue") ~ "UA",
    str_detect(DiseaseStage, "early-stage") ~ "Early_RA",
    TRUE ~ "Established RA"
  )) %>%
  select(sampleid, condition) %>%
  filter(condition %in% c("Healthy", "Arthralgia", "UA", "Early_RA"))

View(colData)
dim(colData) # 100 x 2
  
# Statistics of meta data
colData %>%
  group_by(condition) %>%
  summarise(count = n())

## Arthralgia 10, Early_RA 57, Healthy 27, UA 6

fwrite(colData, file = "./data/colData_early_disease.csv")

# Prepare count table -----------------------------------------------------

countData <- count_tb[, colnames(count_tb) %in% c("id", colData$sampleid)]
View(countData)
dim(countData) # 60699 x 101

fwrite(countData, file = "./data/countData_early_disease.csv")

# Tidy data ---------------------------------------------------------------

# Tidy colData
colData <- colData %>%
  column_to_rownames(var = "sampleid") %>%
  mutate(condition = as.factor(condition))

# Tidy countData
countData <- countData %>%
  column_to_rownames(var = "id")

# The columns of the count matrix and the rows of colData must be in the same order
print("Are columns of the count matrix equal to rows of colData?")

print(all(rownames(colData) %in% colnames(countData))) # TRUE

print(all(rownames(colData) == colnames(countData))) # TRUE

# Construct a DESeqDataSet ------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = round(countData),
  colData = colData,
  design = ~condition
)

# The first level of a factor be the reference level
dds$condition <- factor(dds$condition,
  levels = c("Healthy", "Arthralgia", "Early_RA", "UA")
)
                                                  
## Quality check: Average read counts per sample
summary(round(colSums(assay(dds)) / 1e6, 1)) # Average: 34.07 million read counts per sample

# Save dds
saveRDS(dds, file = "./data/DESeq2_dds_early_disease.rds")
