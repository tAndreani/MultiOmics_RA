# Introduction ------------------------------------------------------------

# This R snippet reads flowSet, compensates and applies biexponential transformation
# Configuration: >=8 Cores, >=32 Gb RAM
# Data: fcs files from CIA study
# Output: Saved GatingSet after compensation and transformation

install.packages("styler")

# Load packages -----------------------------------------------------------

library(flowCore)
library(ncdfFlow)
library(flowWorkspace)
library(openCyto)
library(tidyverse)

# List directories --------------------------------------------------------

fcsDir <- list.files("./fcs", full = TRUE)
fcsDir

# Loop through multiple directories of fcs files --------------------------

## Note: Compensation and transformation must be applied directly to GatingSet instead of flowSet

for (i in c(1:8)) {
  
  print(paste0("Start read fcs files in ", fcsDir[i]))

  # Define fcs file paths
  fcsFiles <- list.files(path = fcsDir[i], pattern = "ms", full = TRUE)

  # Read flowSet
  fs <- read.ncdfFlowSet(fcsFiles,
    truncate_max_range = FALSE,
    channels = c(
      "FSC-H", "FSC-A", "SSC-A", "FL1-A",
      "FL2-A", "FL3-A", "FL4-A", "FL5-A",
      "FL6-A", "FL7-A", "FL8-A", "FL9-A",
      "FL10-A", "FL11-A", "FL12-A", "FL13-A"
    ), mc.cores = 8
  )
  
  print("Start working on GatingSet")
  
  # Create GatingSet object
  gs <- GatingSet(fs)
  
  # Define metadata file path
  mdPath <- list.files(path = fcsDir[i], pattern = "metadata", full = TRUE)
  print(paste0("Read metadata: ", mdPath))
  
  # Load metadata
  md <- read.csv(mdPath)
  
  # Merge pData and metadata
  pdMerge <- inner_join(pData(gs), md, by = "name")
  
  print("Update pData of GatingSet")
  
  # Update pData
  pData(gs)$sample_id <- pdMerge$sample_id
  pData(gs)$condition <- pdMerge$condition
  pData(gs)$patient_id <- pdMerge$patient_id

  # Define compensation file path
  compPath <- list.files(path = fcsDir[i], pattern = "Compensation", full = TRUE)
  
  print(paste0("Read compensation matrix: ", compPath))
  
  # Extract pre-calculated spillover matrix
  compDf <- read.csv(compPath, row.names = 1)
  
  # Rename column names
  colnames(compDf) <- rownames(compDf)
  
  # Create a compensation object
  comp <- compensation(compDf)
  
  # Apply spillover matrix to GatingSet
  gs <- compensate(gs, comp)
  
  print("Tidy up channel names")
  
  # Correspondence between channels and marker names
  df_chl_marker <- gh_pop_get_data(gs)@parameters@data %>%
    select(name, desc)
  
  # Tidy up channel description
  Split <- strsplit(as.character(df_chl_marker$desc), " ", fixed = TRUE)
  marker <- sapply(Split, "[", 1)
  
  marker[11] <- "CD45R"
  
  # Rename
  colnames(gs) <- marker
  print(colnames(gs))
  
  print("Start biexponential transformation")
  
  # Define channels for transformation
  chnls <- colnames(gs)[4:16]
  print(chnls)
  
  # Define biexponential transformation object
  biexpTrans <- flowjo_biexp_trans(
    channelRange = 4096, maxValue = 262144,
    pos = 4.5, neg = 0, widthBasis = -10
  )
  
  # Create a transformerList
  tf <- transformerList(chnls, biexpTrans)
  
  # Add transformerList to GatingSet
  gs <- transform(gs, tf)
  
  # Define GatingSet file path
  gsPath <- paste0("./gatingSet/After_Comp_Trans/", str_split_fixed(fcsDir[i], "/", 3)[, 3])
  
  print(paste0("Save compensated and transformed GatingSet at: ", gsPath))
  
  save_gs(gs, path = gsPath, verbose = TRUE)
  
  rm(fcsFiles)
  rm(fs)
  rm(gs)
  rm(tf)
  
}
