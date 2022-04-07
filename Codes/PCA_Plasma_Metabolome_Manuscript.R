# Introduction ------------------------------------------------------------

# Image: metabolomicsr
# Experiment: Untargeted metabolomics on EDTA plasma from CIA mice
# This R script generates PCA plot

# Load library ------------------------------------------------------------

library(data.table)
library(tidyverse)
library(MetaboAnalystR)
library(pheatmap)
library(RColorBrewer)

setwd("~/cloud-data/snf-mgln-immunometabolism/CIA_FH2000/Metabolomics")

# Load meta data ----------------------------------------------------------

# Load mouse arthritis score at sacrifice
arth_sum <- fread("~/cloud-data/snf-mgln-immunometabolism/CIA_FH2000/CIA_multiomics_figure/data/metadata/CIA_Ctrl_Score_Sum_Sacrifice_20210118.csv")
View(arth_sum)

# Prepare data for time-series/two-factor design module -------------------

# Create objects for storing processed data for time-series analysis
mSet <- InitDataObjects("conc", "ts", paired = FALSE)

mSet <- SetDesignType(mSet, "time") ## The time points group must be labeled as <b>Time</b>

mSet <- Read.TextData(
  mSetObj = mSet,
  "./data/Batch_Norm_Impute_Intensity/Plasma_Batch_Norm_Impute_Intensity_with_MetaData.csv",
  "rowts", "disc"
)
View(mSet)

# Sanity check
mSet <- SanityCheckData(mSet)

# Minimum value replacing
mSet <- ReplaceMin(mSet)

# Check if the sample size is too small
mSet <- IsSmallSmplSize(mSet)

# Normalization to sample median, with Log
mSet <- PreparePrenormData(mSet)
mSet <- Normalization(mSet, "MedianNorm", "LogNorm", "NULL", ratio = FALSE, ratioNum = 20)

# Plot sample normalization summary
mSet <- PlotNormSummary(mSet, "./plot/heatmap/normConc_timeseries_", "png", 72, width = NA)
mSet <- PlotSampleNormSummary(mSet, "./plot/heatmap/normSample_timeseries_", "png", 72, width = NA)

# For the following analysis, change working directory to save the output
setwd("~/cloud-data/snf-mgln-immunometabolism/CIA_FH2000/Metabolomics/result")

# Set plot parameter ------------------------------------------------------

theme_set(theme_bw())
theme_update(
  text = element_text(family = "Helvetica", size = 7),
  axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
  plot.title = element_text(size = 7),
  legend.title = element_text(size = 7),
  legend.key.size = unit(4, "mm")
)

# PCA ---------------------------------------------------------------------

# Calculate PCA
mSet <- PCA.Anal(mSet)

# Create data frame with PCA scores
score <- as.data.frame(mSet$analSet$pca$x)
View(score)

# Add meta data to PCA scores
score <- score %>%
  rownames_to_column(var = "mouseID") %>%
  mutate(mouseID = as.numeric(mouseID), Time = mSet$dataSet$facB) %>%
  left_join(arth_sum, by = "mouseID")

# Rank the order of timepoints
score$Time <- factor(score$Time, levels = c(
  "2 days", "2 weeks", "3 weeks", "4 weeks",
  "6 weeks", "7 weeks", "8 weeks", "10 weeks"
))

# Rank the mouse treatments
score$Treatment <- factor(score$Treatment, levels = c("Ctrl", "CIA"))

# Define label names
xlabel <- paste0("PC1: ", round(100 * mSet$analSet$pca$variance[1], 1), "% variance")
xlabel # PC1 (23.8%)
ylabel <- paste0("PC2: ", round(100 * mSet$analSet$pca$variance[2], 1), "% variance")
ylabel # PC2 (8.7%)

# Create data frame to extract convex hull points
pca_hull <-
  score %>%
  group_by(Treatment) %>%
  slice(chull(PC1, PC2))

# PCA plot with timepoint
pca_tp <- ggplot(score, aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = Time, shape = Treatment), size = 3) +
  scale_shape_manual(values = c(1, 16)) +
  labs(x = xlabel, y = ylabel, title = "PCA on plasma metabolome") +
  geom_polygon(data = pca_hull, aes(fill = Treatment), alpha = 0.2, show.legend = TRUE) +
  scale_fill_manual(values = c("#D95F02", "#1B9E77")) +
  scale_color_brewer(palette = "Paired")

# Save plot
imgName <- "./plot/clustering/PCA_Treatment_Timepoint.tif"
ggsave(imgName,
       plot = pca_tp, device = "tiff",
       units = "mm", width = 100, height = 85, dpi = 300
)
