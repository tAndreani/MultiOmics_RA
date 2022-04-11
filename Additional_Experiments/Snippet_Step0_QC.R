# This R snippet applies flowAI for automatic quality control of flow cytometry data
# Data: fcs files containing live cells after gating
# Configuration: >= 8 Cores, >= 64 Gb RAM
# Output: Newly created directory resultsQC containing QC statistics in html format and fcs files after QC

library(flowAI)

# Define file paths
fcsFiles <- list.files(path = "./fcs_liveCells", pattern = ".fcs", full = TRUE)
fcsFiles

resQC <- flow_auto_qc(fcsFiles, output = 0)
