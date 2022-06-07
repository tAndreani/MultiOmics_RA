# Introduction ------------------------------------------------------------

# Project: CIA Paw RNAseq pre-processed by RASflow
# Code adapted from: https://angus.readthedocs.io/en/2014/analyzing_drosophila_htseq.html
# This R script combines all countFiles from individual tsv files

# Load library ------------------------------------------------------------

library(stringr)
library(data.table)

# Extract sample IDs ------------------------------------------------------

# Define countFile directory path
countFilePath <- "~/cloud-data/snf-mgln-immunometabolism/CIA_FH2000/RNAseq/RASflowOutput/CIA_Paw_Combi_Run_20210607_114851/genome/countFile"

# List count tsv file
tsvDir <- list.files(countFilePath, pattern = "_count.tsv")
tsvDir

# Extract sample IDs
samples <- str_split_fixed(tsvDir, "_count.tsv", 2)[, 1]
samples

# Create a function to read countFile -------------------------------------

read.sample <- function(sample.name) {
  file.name <- paste0(countFilePath, "/", sample.name, "_count.tsv")
  result <- read.delim(file.name, col.names = c("gene_id", "count"), sep = "\t", colClasses = c("character", "numeric"))
}

# Read the first sample ---------------------------------------------------

sample.1 <- read.sample(samples[1])
View(sample.1)
dim(sample.1)

# Read all the rest of the samples in loop --------------------------------

# Use data frame of the first sample as the base
all.data <- sample.1

# Loop to read the rest of the samples
for (c in 2:length(samples)) {
  
  temp.data <- read.sample(samples[c])

  print(paste0("Finish reading sample: ", samples[c]))

  if (temp.data$gene_id != all.data$gene_id) {
    break
  } else {
    print(paste0("Adding sample ", samples[c], " to combined table"))

    all.data <- cbind(all.data, temp.data$count)
  }
}

# Define column names with sample IDs
colnames(all.data)[2:ncol(all.data)] <- samples

# Check combined table
View(all.data)
tail(all.data)
dim(all.data) # 53646 x 172

# Save combined table
fwrite(all.data, 
       file = "~/cloud-data/snf-mgln-immunometabolism/CIA_FH2000/RNAseq/CIA_Paw_RNAseq/data/combined_countFile.csv")
