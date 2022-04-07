# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: DESeq2 colData (a table with information about the samples)
# Data: CIA mouse metadata

# Load library ------------------------------------------------------------

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

# Load data ---------------------------------------------------------------

# Count matrix
countData <- read.csv("./data/paw_RNAseq/combined_countFile.csv", row.names = 1)
View(countData)

# Batch data
batchData <- read.csv("./data/metadata/paw_sample_merged_all_batches.csv")
View(batchData)

# Individual paw score at sacrifice
Score_Ind <- read.csv("./data/metadata/CIA_Score_Individual_Sacrifice.csv",
  row.names = 1
)
View(Score_Ind)

# Construct colData -------------------------------------------------------

# The columns of the count matrix and the rows of the colData must be in the same order
df_1 <- data.frame(sampleID = colnames(countData))
View(df_1)

# Extract mouse IDs (The number right after "X")
df_1$mouseID <- str_split_fixed(df_1$sampleID, "_", 2)[, 1] %>%
  str_replace("X", "")

# Merge with batch data
batchData$mouseID <- as.character(batchData$mouseID)

df_2 <- df_1 %>%
  left_join(batchData[, c(2, 5, 6, 8)], by = "mouseID")
View(df_2)

# Merge with score of hind right paw
Score_Ind$mouseID <- as.character(Score_Ind$mouseID)

colData <- Score_Ind %>%
  filter(Location == "HindRight") %>%
  select(Score, mouseID, Treatment) %>%
  right_join(df_2, by = c("mouseID", "Treatment"))
View(colData)

# Define score of Ctrl paw as 0
colData$Score <- replace_na(colData$Score, 0)

# Re-order columns
colData <- colData[, c(4, 2, 3, 5, 6, 1)]
glimpse(colData)

# Save colData
fwrite(colData, "./data/metadata/DESeq2_colData.csv")
