# Introduction ------------------------------------------------------------

# This R snippet applies OpenCyto for automated gating to select CD45+ cells
# Configuration: >= 16 Cores, >= 64 Gb RAM
# Data: GatingSet files from CIA study, after compensation and transformation
# Output: Saved csv file with population statistics after autoGating
# Output: Saved fcs files with gated CD45+ cell events (with inverse transformation)
# Output: Saved GatingSet with gated CD45+ cell events

# Samples with abnormal number of dead cells ------------------------------

# week7_20201110
## All fcs files! Therefore save the fcs files in a separate folder

# Load packages -----------------------------------------------------------

install.packages("styler")

library(flowWorkspace)
library(openCyto)
library(data.table)
library(ggcyto)
library(tidyverse)
library(parallel)

rm(list = ls()) # Clear all variables in environment

# List directories --------------------------------------------------------

gsDir <- list.files("./gatingSet/After_Comp_Trans", full = TRUE)
gsDir

# Define which directory to apply autoGating for live cells
i <- 7
gsDir[i]

# Load GatingSet ----------------------------------------------------------

gs <- load_gs(path = gsDir[i], verbose = TRUE)

# Load gating template ----------------------------------------------------

## Only upload gating template into the directory after GatingSet is loaded
# Define file path of gating template
gtFile <- list.files(path = "./gatingSet", pattern = ".csv", full = TRUE)
gtFile

# Get gating template
gt_lc <- gatingTemplate(gtFile)
gt_lc

plot(gt_lc)

# View gating template
View(fread(gtFile))

# Gating for live cells ---------------------------------------------------

# Apply gating template
gt_gating(gt_lc, gs, parallel_type = "multicore", mc.cores = 8)

# Get names of all nodes in GatingSet
gs_get_pop_paths(gs)

# Plot gating scheme
plot(gs[[1]])

# To be saved: Dot plot of selected gate of all fcs files -----------------

# Dot plot with nonDebris gate on all fcs files
pdfName <- paste0("./plot/gatingPlot/nonDebris_", str_split_fixed(gsDir[i], "/", 4)[, 4], ".pdf")
pdfName

p <- ggcyto(gs, aes(x = `FSC-A`, y = `SSC-A`), subset = "root") +
  geom_hex(bins = 6400) +
  geom_gate("nonDebris") +
  geom_stats(size = 6, color = "white", fill = "black", adjust = 0.3) +
  labs_cyto("marker")

pdfPlot <- p +
  facet_wrap(~name, ncol = 4) +
  ggcyto_par_set(
    limits = list(x = c(0, 2e6), y = c(0, 2e6)),
    hex_fill = scale_fill_viridis_c(option = "plasma", alpha = .8)
  )

# All plots saved together in pdf
ggsave(pdfName, plot = pdfPlot, width = 18, height = 18, units = "in")

# Single demo plot saved in png
ggsave("./plot/gatingPlot/nonDebris_ms9_week7_20201110.png", plot = pdfPlot, width = 9, height = 8, units = "cm")

# Dot plot with singlets gate on all fcs files
pdfName <- paste0("./plot/gatingPlot/singlets_", str_split_fixed(gsDir[i], "/", 4)[, 4], ".pdf")
pdfName

p <- ggcyto(gs, aes(x = `FSC-A`, y = `FSC-H`), subset = "nonDebris") +
  geom_hex(bins = 6400) +
  geom_gate("singlets") +
  geom_stats(size = 6, color = "white", fill = "black", adjust = 0.3) +
  labs_cyto("marker")

pdfPlot <- p +
  facet_wrap(~name, ncol = 4) +
  ggcyto_par_set(
    limits = list(x = c(0, 2e6), y = c(0, 2e6)),
    hex_fill = scale_fill_viridis_c(option = "plasma", alpha = .8)
  )

# All plots saved together in pdf
ggsave(pdfName, plot = pdfPlot, width = 18, height = 18, units = "in")

# Single demo plot saved in png
ggsave("./plot/gatingPlot/singlets_ms9_week7_20201110.png", plot = pdfPlot, width = 9, height = 8, units = "cm")

# Dot plot with live cell gate on all fcs files
pdfName <- paste0("./plot/gatingPlot/liveCells_", str_split_fixed(gsDir[i], "/", 4)[, 4], ".pdf")
pdfName

p <- ggcyto(gs, aes(x = DCM, y = `FSC-A`), subset = "singlets") +
  geom_hex(bins = 6400) +
  geom_gate("liveCells") +
  geom_stats(size = 6, color = "white", fill = "black", adjust = 0.3) +
  axis_x_inverse_trans() +
  labs_cyto("marker")

pdfPlot <- p +
  facet_wrap(~name, ncol = 4) +
  ggcyto_par_set(
    limits = list(x = c(-1e3, 5e3), y = c(0, 1e6)),
    hex_fill = scale_fill_viridis_c(option = "plasma", alpha = .8)
  )

# All plots saved together in pdf
ggsave(pdfName, plot = pdfPlot, width = 18, height = 18, units = "in")

# Single demo plot saved in png
ggsave("./plot/gatingPlot/liveCells_ms9_week7_20201110.png", plot = pdfPlot, width = 9, height = 8, units = "cm")

# R session crashed here!!
# Dot plot with CD45 cell gate on all fcs files
pdfName <- paste0("./plot/gatingPlot/CD45Cells_", str_split_fixed(gsDir[i], "/", 4)[, 4], ".pdf")
pdfName

p <- ggcyto(gs[[16]], aes(x = CD45, y = `SSC-A`), subset = "liveCells") +
  geom_hex(bins = 6400) +
  geom_gate("CD45Cells") +
  geom_stats(size = 6, color = "white", fill = "black", adjust = 0.01) +
  axis_x_inverse_trans() +
  labs_cyto("marker")

pdfPlot <- p +
  facet_wrap(~name, ncol = 4) +
  ggcyto_par_set(
    limits = list(x = c(-4e3, 4e3), y = c(0, 1.5e6)),
    hex_fill = scale_fill_viridis_c(option = "plasma", alpha = .8))

# All plots saved together in pdf
ggsave(pdfName, plot = pdfPlot, width = 18, height = 18, units = "in")

# Single demo plot saved in png
ggsave("./plot/gatingPlot/CD45Cells_ms9_week7_20201110.png", plot = pdfPlot, width = 9, height = 8, units = "cm")

# Extra visualization -----------------------------------------------------

# Autoplot all gates for a single sample
autoplot(gs[[2]], strip.text = "gate", 
         bins = 3200, merge = F, 
         axis_inverse_trans = F) +
  ggcyto_par_set(facet = facet_wrap(~name, scales = "free"), limits = list(x = c(0, 5e6), y = c(0, 5e6))) +
  theme_bw()

# Density plot of nonDebris gate on all fcs files
p <- ggcyto(gs, aes(x = `FSC-A`)) +
  geom_density() +
  geom_gate("nonDebris") +
  geom_stats(size = 6, color = "white", fill = "black", adjust = 0.3) +
  labs_cyto("marker")

p +
  facet_wrap(~name, ncol = 4)

# Density plot of concatenated flowFrames, facet by condition
fs <- gs_pop_get_data(gs, "nonDebris")

ggplot(fs, aes(x = `FSC-A`)) +
  geom_density(fill = "blue", alpha = 0.5) +
  facet_wrap(~condition) +
  xlim(0, 5e6)

# Change facetting (Default is facet_wrap(~name))
p + facet_grid(condition ~ day)

# Output gating stats -----------------------------------------------------

# Get all the population statistics
ps <- data.frame(gs_pop_get_count_fast(gs))

# Calculate "percent of parent"
ps$percent_of_parent <- ps$Count / ps$ParentCount

# Add metadata to population stats
psm <- merge(ps, pData(gs), by = "name")

# Save csv file
statName <- paste0("./statistics/CD45Cells_autoGating_", str_split_fixed(gsDir[i], "/", 4)[, 4], ".csv")
statName

write.csv(psm, file = statName)

# Export fcs files with gated CD45Cells events (inverse trans) ------------

# Check how many samples in GatingSet
gs

# Loop to write FCS file from a flowFrame
for (j in c(1:22)) {
  
  # Retrieve flowFrame with inverse transformation of data
  fs_CD45Cells <- gh_pop_get_data(gs[[j]], "CD45Cells", inverse.transform = TRUE)
  
  # Define output file name
  outFile <- paste0("./fcs_CD45Cells/week7_CD45Cells/CD45Cells-", fs_CD45Cells@description$`$FIL`)
  print(paste0("Save fcs file as: ", outFile))
  
  write.FCS(fs_CD45Cells, outFile)

}

# Save GatingSet with gated live cell events ------------------------------

gsPath <- paste0("./gatingSet/CD45Cells/", str_split_fixed(gsDir[i], "/", 4)[, 4])
gsPath

save_gs(gs, path = gsPath, verbose = TRUE)

# Double check on subsetting of CD45Cells ---------------------------------

CD45Cells_count <- read.FCS("./fcs_CD45Cells/CD45Cells-01-healthy_ms164-F8.fcs", transform=FALSE)
CD45Cells_count # 246920 cells

root_count <- read.FCS("./fcs/day2_20200924/01-healthy_ms164-F8.fcs", transform=FALSE)
root_count # 505772 cells

colnames(CD45Cells_count) # "FSC-H" "FSC-A" "SSC-A" "CD62L" "PD-1" "CD45" "CD11b" "CD19" "CD3" "CD44" "CD45R" "CD138" "DCM" "CD8" "CXCR5" "CD4"
