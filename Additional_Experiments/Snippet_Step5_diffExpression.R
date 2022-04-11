# Introduction ------------------------------------------------------------

# This R snippet performs differential analysis of marker expression stratified by cell population
# Data: sce object
# Configuration: >= 8 Cores, >= 64 Gb RAM
# Output: Heatmap with FDR data on cell state markers of B cells, CD4 T cells and CD8 T cells

# Load packages -----------------------------------------------------------

install.packages("styler")

library(CATALYST)
library(tidyverse)
library(data.table)
library(cowplot)
library(diffcyt)

# Restore sce -------------------------------------------------------------

sce <- readRDS("./sce_CIA_wo_6weeks_typeMarkers.rds")

# Manual cluster merging and annotation -----------------------------------

# Import manual annotation
merging_table1 <- read.csv("./clustering/CIA_cluster_merging1_wo_6 weeks_typeMarkers_20210202.csv")
View(merging_table1)

# Convert to factor with merged clusters in desired order
merging_table1$new_cluster <- factor(merging_table1$new_cluster)

# Apply manual merging
sce <- mergeClusters(sce, k = "meta12", table = merging_table1, id = "merging1")

# Heatmap of median marker expression with median expression data  --------

FDR_cutoff <- 0.05

# Loop for timepoints
patient_id_Name <- unique(sce$patient_id)
patient_id_Name

for (i in seq(1:length(patient_id_Name))) {
  
  tp <- patient_id_Name[i]
  
  print(paste0("Heatmap on ", tp))
  
  # Subset
  sceSubset <- filterSCE(sce, patient_id == tp)
  
  # Extract meta data
  ei <- metadata(sceSubset)$experiment_info
  
  # Fit a regular linear model
  ds_formula1 <- createFormula(ei, cols_fixed = "condition")
  
  contrast <- createContrast(c(0, 1))
  
  ds_res1 <- diffcyt(sceSubset, 
                     formula = ds_formula1, contrast = contrast,
                     analysis_type = "DS", method_DS = "diffcyt-DS-LMM",
                     clustering_to_use = "merging1", verbose = TRUE)
  
  # Output data
  diffData <- topTable(ds_res1, all = TRUE, order_by = "cluster_id", 
                       show_meds = TRUE, format_vals = TRUE, digits = 2) %>%
    as.data.frame()
  
  tableName <- paste0("./statistics/diffcyt_", tp, "_20210210.csv")
  
  fwrite(diffData, tableName)
  
  # Heatmap
  # hp <- plotDiffHeatmap(sceSubset, rowData(ds_res1$res), all = TRUE, fdr = FDR_cutoff, sort_by = "none")
  
  # plotName <- paste0("./plot/clusteringPlot_wo_6weeks_typeMarkers/diffHeatmap_", tp, "_20210202.png")
  
  # hp
  # png(plotName, width = 700, height = 400, units = "px")
  # print(hp)
  # dev.off()
}

# Count the number of differential findings 
table(rowData(ds_res1$res)$p_adj < FDR_cutoff)

# Extract p value and FDR
rowData(ds_res1$res)

