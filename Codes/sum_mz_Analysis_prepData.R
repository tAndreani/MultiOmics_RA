# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: Peak intensity file for construction of mSet in MetaboAnalystR
# Data: tMSI raw data and annotation (Caution: A few metabolite annotation might be duplicated)

# Load library ------------------------------------------------------------

library(data.table)
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)

# Load data ---------------------------------------------------------------

# Load raw data
rawData <- fread("~/cloud-data/snf-mgln-immunometabolism/CIA_FH2000/tMSI/data/DMD1659/DMD1659_v2_raw.csv")
View(rawData)
dim(rawData) # 44 x 11515

# Load annotation file
annoMol <- read.csv("~/cloud-data/snf-mgln-immunometabolism/CIA_FH2000/tMSI/data/DMD1659/DMD1659_annotation.csv")
View(annoMol)
dim(annoMol) # 92112 x 15

# Mouse metadata
metadata <- fread("~/cloud-data/snf-mgln-immunometabolism/CIA_FH2000/CIA_multiomics_figure/data/metadata/paw_sample_merged_all_batches.csv")
View(metadata)

# Tidy annotation table ---------------------------------------------------

# Tidy annotation table
annoMol.tidy <- distinct(annoMol[, 3:15]) %>%
  mutate(Obs_mz_paste = paste0("X", format(round(Obs..m.z., 8), nsmall = 3))) %>%
  mutate(Obs_mz_X = str_remove_all(Obs_mz_paste, " "))
View(annoMol.tidy)

# Create a match table with unique metabolite IDs
matchTb <- data.frame(Candidate = unique(annoMol.tidy$Candidate), 
                      metaID = paste0("meta_", 1:length(unique(annoMol.tidy$Candidate))))
head(matchTb)
dim(matchTb) # 8176 x 2

# Add unique metabolite IDs to annotation table
annoMol.metaID <- annoMol.tidy %>%
  left_join(matchTb, by = "Candidate") %>%
  distinct()
View(annoMol.metaID)

annoMol.metaID <- annoMol.metaID[, c(16, 15, 1:13)]
dim(annoMol.metaID) # 11514 x 15

# Save annotation table
fwrite(annoMol.metaID, "./data/DMD1659/DMD1659_sum_mz_annotation_tidy.csv")

# Tidy raw data table -----------------------------------------------------

# Clean up m/z and mouseID
rawData.v1 <- rawData %>%
  mutate(mouseID = as.numeric(str_split_fixed(Regions, "_M", 2)[, 2])) %>%
  select(!which(sapply(., is.character))) %>%
  # Remove columns with characters
  pivot_longer(!mouseID, names_to = "column_name", values_to = "peak_area") %>%
  mutate(mz = as.numeric(str_split_fixed(column_name, " ", 10)[, 2])) %>%
  # Extract m/z
  mutate(mz_column = str_remove_all(paste0("X", format(round(mz, 8), nsmall = 3)), " "))

View(rawData.v1)
dim(rawData.v1) # 506484 x 5

# Merge with annotation table
rawData.v2 <- rawData.v1 %>%
  select(mouseID, peak_area, mz_column) %>%
  left_join(annoMol.metaID[, c(1, 2)], by = c("mz_column" = "Obs_mz_X")) %>%
  distinct()
View(rawData.v2)
dim(rawData.v2) # 506484 x 4

length(unique(rawData.v2$mz_column)) # 11511 m/z
length(unique(rawData.v2$metaID)) # 8174 metabolite names

# Sum up m/z of the same metabolite together
rawData.v3 <- rawData.v2 %>%
  group_by(mouseID, metaID) %>%
  summarise(sum_peak_area = sum(peak_area))
View(rawData.v3)

# Add mouse metadata
rawData.v4 <- rawData.v3 %>%
  filter(!is.na(sum_peak_area)) %>%
  left_join(metadata[, c(2,5,6)], by = "mouseID") %>%
  pivot_wider(names_from = metaID, values_from = sum_peak_area)
View(rawData.v4)
dim(rawData.v4) # 44 x 8177

# Tidy column names
rawData.tidy <- rawData.v4 %>%
  rename(Time = Timepoint)

head(colnames(rawData.tidy))

# Save prepared raw data
finalCol <- ncol(rawData.tidy)

fwrite(rawData.tidy[, c(1, 3, 2, 4:finalCol)], "./data/DMD1659/DMD1659_sum_mz_v2_raw_tidy.csv")

# Prepare data for statistical analysis -----------------------------------

dim(rawData.tidy) # 44 x 8177

# Loop to save table with sample as column, metabolite as row
for (i in unique(rawData.tidy$Time)) {
  
  print(i)
  
  # Subset with timepoint
  rawData.subset <- rawData.tidy %>%
    filter(Time == i) %>%
    select(!Time) %>%
    rename(Sample = mouseID, Label = Treatment)
  
  # Transpose table
  trans.table <- t(rawData.subset)
  
  # Define file name
  fileName <- paste0("./data/DMD1659/stat_DMD1659_sum_mz_v2_raw_", str_replace_all(i, " ", "_"), ".csv")
  fileName
  
  # Save table
  write.table(trans.table, file = fileName, quote = FALSE, row.names = TRUE, col.names = FALSE, sep = ",")
}
