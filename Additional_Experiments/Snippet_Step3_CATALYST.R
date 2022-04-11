# Introduction ------------------------------------------------------------

# This R snippet applies CATALYST for FlowSOM clustering
# Data: fcs files containing CD45+ cells after gating (without 6 and 7 weeks)
# Configuration: >= 8 Cores, >= 64 Gb RAM
# Output: Saved sce object after clustering
# Output: Automatic clustering based on self-organizing maps
# Output: Heatmaps, t-SNE and UMAP plots
# Output: Cell population abundance and marker expression

# Load packages -----------------------------------------------------------

install.packages("styler")

library(CATALYST)
library(flowCore)
library(ncdfFlow)
library(tidyverse)
library(FlowSOM)
library(data.table)
library(cowplot) # For arranging multiple plots into a grid

# One-time only: Prepare metadata for fcs files ---------------------------

# Extract file names
file_name <- list.files(path = "./fcs_CD45Cells", pattern = ".fcs")

# Create data frame
df <- data.frame(file_name)
View(df)
dim(df) # 136 x 1

# Extract sample IDs
df <- df %>% mutate(name_1 = str_split_fixed(file_name, "-", 4)[, 3]) %>%
  mutate(name_2 = str_split_fixed(name_1, "_", 2)[,2]) %>%
  mutate(sample_id = str_replace_all(name_2, "ms", ""))
# Note: 189 - 194 should be Ctrl

# Load timepoints
CIA_Timepoints <- read.csv("/cloud-data/snf-mgln-immunometabolism/CIA_FH2000/FACS/clustering/CIA_Timepoints.csv")
View(CIA_Timepoints)
dim(CIA_Timepoints) # 195 x 3

CIA_Timepoints$sample_id <- as.character(CIA_Timepoints$sample_id)

# Merge files to create metadata file
md <- left_join(df, CIA_Timepoints, by = "sample_id")
View(md)
dim(md) # 136 x 6
unique(md$patient_id)

# Remove extra columns
md <- md %>% select(!c(name_1, name_2))

# Save metadata
write.csv(md, file = "./clustering/CD45Cells_fcs_wo_6weeks_metadata.csv")

# Load metadata
md <- read.csv("/cloud-data/snf-mgln-immunometabolism/CIA_FH2000/FACS/clustering/CD45Cells_fcs_wo_6weeks_metadata.csv", row.names = 1)
View(md)

# Load fcs files into flowSet ---------------------------------------------

# Define file paths of CD45 fcs files
fcsFiles <- list.files(path = "./fcs_CD45Cells", pattern = ".fcs", full = TRUE)
fcsFiles

# Read flowSet
fs <- read.ncdfFlowSet(fcsFiles,
  truncate_max_range = FALSE,
  channels = c(
    "CD62L", "PD-1", "CD45", "CD11b",
    "CD19", "CD3", "CD44", "CD45R",
    "CD138", "CD8", "CXCR5", "CD4"
  ),
  mc.cores = 16
)

fs

# Descriptive infomation
colnames(fs)
sampleNames(fs)
summary(fs[[1]])

# Import panel ------------------------------------------------------------

panel <- read.csv("/cloud-data/snf-mgln-immunometabolism/CIA_FH2000/FACS/clustering/Marker_Panel.csv")
View(panel)

# Check that all panel columns are in the flowSet object
all(panel$fcs_colname %in% colnames(fs))

# Prepare SingleCellExperiment (SCE) class --------------------------------

# Specify levels for conditions & sample IDs to assure desired ordering
## For the statistical modeling, make the condition variable a factor with Ref being the reference level
glimpse(md)

md$condition <- factor(md$condition, levels = c("Ctrl", "CIA"))

md$sample_id <- factor(md$sample_id,
  levels = md$sample_id[order(md$condition)]
)

md$patient_id <- factor(md$patient_id, 
                        levels = c("2 days", "2 weeks", "3 weeks",
                                   "4 weeks", "8 weeks", "10 weeks"))

# Construct SingleCellExperiment
# arcsinh transforms with a cofactor of 150 for flow cytometry, and 5 for CyTOF
sce <- prepData(fs, panel, md,
  features = panel$fcs_colname,
  transform = TRUE,
  cofactor = 150,
  FACS = TRUE
)

View(sce)
metadata(sce)

# Double-check cofactor
int_metadata(sce)$cofactor

# Marker names
rowData(sce)

# View experimental design table
exprDesign <- ei(sce)
exprDesign

write.csv(exprDesign, file = "./clustering/CD45Cells_fcs_wo_6weeks_cellCount.csv")

# Tidy up the table structure
exprDesign$sample_id <- as.factor(exprDesign$sample_id)
exprDesign$patient_id <- factor(exprDesign$patient_id,
  levels = c(
    "2 days", "2 weeks", "3 weeks",
    "4 weeks", "6 weeks", "8 weeks", "10 weeks"
  )
)

sum(exprDesign$n_cells) # In total: 30,497,842 CD45 cells (wo 6 and 7 weeks)

# To save memory, remove fs
rm(fs)

# Cell count plot ---------------------------------------------------------

# Cell count
n_cells(sce)

# Cell count bar plot
cellCount <- plotCounts(sce, group_by = "sample_id", color_by = "condition") +
  ggtitle("Count of single cells") +
  labs(x = "Sample ID", y = "Number of gated live cells", title = "Count of gated live cells")

ggsave("./plot/clusteringPlot/cellCount_barPlot_20210101.png", plot = cellCount, width = 20, height = 14, units = "cm")

# Cell count bar plot, group by timepoints
cellCount2 <- exprDesign %>% ggplot(aes(x = patient_id, y = n_cells, fill = condition)) +
  geom_boxplot(position = position_dodge()) +
  labs(x = NULL, y = "Cell number", title = "CD45+ cell count") +
  theme(legend.title = element_blank()) +
  scale_fill_brewer(palette = "Dark2")

ggsave("./plot/clusteringPlot/cellCount_boxPlot_20210117.png", plot = cellCount2, width = 20, height = 14, units = "cm")

# MDS plot ----------------------------------------------------------------

# MDS plot on the median marker expressions
mdsPlot <- pbMDS(sce, color_by = "patient_id", shape_by = "condition", label_by = NULL) +
  scale_color_discrete(name = NULL)

ggsave("./plot/clusteringPlot_wo_6weeks_typeMarkers/MDSplot_wo_6weeks_20210202.png", plot = mdsPlot, width = 16, height = 16, units = "cm")

# PCA plot ----------------------------------------------------------------

pcaPlot <- plotDR(sce, "PCA", color_by = "patient_id") +
  scale_color_discrete(name = NULL)

ggsave("./plot/clusteringPlot/PCAplot_20210117.png", plot = pcaPlot, width = 16, height = 16, units = "cm")

# Heatmap on median marker expressions ------------------------------------

medianHp <- plotExprHeatmap(sce,
  by = "sample_id",
  row_anno = c("condition", "patient_id"),
  scale = "last",
  hm_pal = rev(hcl.colors(10, "YlGnBu"))
)

medianHp

png("./plot/clusteringPlot/exprHeatmap_20210117.png", width = 700, height = 1100, units = "px")
print(medianHp)
dev.off()

# Marker ranking with NRS -------------------------------------------------

## Empty black circles indicate the mean NR scores from all the samples
nrsPlot <- plotNRS(sce, features = "type", color_by = "patient_id") +
  theme(legend.title = element_blank())
nrsPlot

ggsave("./plot/clusteringPlot_wo_6weeks_typeMarkers/NRSplot_wo_6weeks_20210202.png", plot = nrsPlot, width = 16, height = 16, units = "cm")

# Cell population identification ------------------------------------------

## With FlowSOM and ConsensusClusterPlus
# cluster() will default to using "type" markers for clustering
# xdim and ydim specify the grid size of SOM
# ~25 minutes of run time for 34,070,837 cells

# Previously, use seed 1234 for 20 clusters
# Now use seed 2021 for 24 clusters
# Now use seed 2021 for 12 cluster with type markers
set.seed(2021)
sce <- cluster(sce,
  features = "type",
  xdim = 10, ydim = 10,
  maxK = 12, seed = 2021,
  verbose = TRUE
)

# View all clusters in SOM
View(cluster_codes(sce))

# Save sce
## This object can be restored into new R sessions using the readRDS function
saveRDS(file = "./sce_CIA_wo_6weeks_typeMarkers.rds", sce)
sce <- readRDS("./sce_CIA_wo_6weeks_typeMarkers.rds") # Restore sce

# Extract the delta area plot
deltaArea <- delta_area(sce)
ggsave("./plot/clusteringPlot_wo_6weeks_typeMarkers/deltaAreaPlot_wo_6weeks_20210202.png", plot = deltaArea, width = 16, height = 16, units = "cm")

# Heatmap of the median marker expression of clusters
clusterHp <- plotExprHeatmap(sce,
  features = "type",
  by = "cluster_id", k = "meta12",
  bars = TRUE, perc = TRUE
)

png("./plot/clusteringPlot_wo_6weeks_typeMarkers/clusterHeatmap_wo_6weeks_20210202.png", width = 600, height = 600, units = "px")
print(clusterHp)
dev.off()

# Output expression data for clusters -------------------------------------

# Define FlowSOM clusters
sce$cluster_id <- cluster_ids(sce, "meta12")

# Construct data frame of expression matrix include cell metadata
cd <- colData(sce)
es <- assay(sce, "exprs")
df <- data.frame(t(es), cd, check.names = FALSE)
glimpse(df)

# Save csv file
setDTthreads(threads = 16)

fwrite(df, file = "./clustering/clusterExpression_wo_6weeks_meta12_20210202.csv", 
       verbose = TRUE)

# Manual cluster merging and annotation -----------------------------------

# Import manual annotation
merging_table1 <- read.csv("./clustering/CIA_cluster_merging1_wo_6 weeks_typeMarkers_20210202.csv")
View(merging_table1)

# Convert to factor with merged clusters in desired order
merging_table1$new_cluster <- factor(merging_table1$new_cluster)

# Apply manual merging
sce <- mergeClusters(sce, k = "meta12", table = merging_table1, id = "merging1")

# Heatmap of median marker expression
mergeCluterHp <- plotExprHeatmap(sce, features = "type", by = "cluster_id", 
                                 k = "merging1", bars = TRUE, perc = TRUE)
mergeCluterHp
png("./plot/clusteringPlot_wo_6weeks_typeMarkers/mergeClusterHeatmap_20210202.png", width = 720, height = 400, units = "px")
print(mergeCluterHp)
dev.off()

# Dimension reduction -----------------------------------------------------

## t-SNE ------------------------------------------------------------------

# t-SNE: Run at most 500 cells per sample
set.seed(1984)
sce <- runDR(sce, "TSNE", cells = 400, features = "type")

# t-SNE: Plot on all clusters together
tsneCluster <- plotDR(sce, "TSNE", color_by = "merging1") +
  theme(legend.title = element_blank())

ggsave("./plot/clusteringPlot_wo_6weeks_typeMarkers/tSNE_Cluster_wo_6weeks_20210202.png", plot = tsneCluster, width = 20, height = 16, units = "cm")

# t-SNE: Loop for each marker
mkName <- rowData(sce)[,2]

for (i in seq(1:length(mkName))) {
  print(paste0("t-SNE plot on ", mkName[i]))
  
  tsnePlot <- plotDR(sce, "TSNE", color_by = mkName[i])
  
  plotName <- paste0("./plot/clusteringPlot_wo_6weeks_typeMarkers/tSNE_", mkName[i], "_20210202.png")
  
  ggsave(plotName, plot = tsnePlot, width = 20, height = 16, units = "cm")
}

# t-SNE Facet: By condition
tsneCond <- plotDR(sce, "TSNE", color_by = "merging1", facet_by = "condition") +
  theme(legend.title = element_blank())

ggsave("./plot/clusteringPlot_wo_6weeks_typeMarkers/tSNE_Condition_20210202.png", plot = tsneCond, width = 20, height = 16, units = "cm")

# t-SNE Facet: By timepoints
## Order by timepoints
sce$patient_id <- factor(sce$patient_id,
                         levels = c(
                           "2 days", "2 weeks", "3 weeks",
                           "4 weeks", "6 weeks", "8 weeks", "10 weeks"
                         )
)

tsneTime <- plotDR(sce, "TSNE", color_by = "merging1", facet_by = "patient_id") +
  theme(legend.title = element_blank())

tsneTime

ggsave("./plot/clusteringPlot_wo_6weeks/tSNE_Timepoint_20210118.png", plot = tsneTime, width = 20, height = 16, units = "cm")

# t-SNE Facet: By condition, and loop for timepoints
patient_id_Name <- c("2 days", "2 weeks", "3 weeks", "4 weeks", "6 weeks", "8 weeks", "10 weeks")

for (i in seq(1:length(patient_id_Name))) {
  print(paste0("t-SNE plot on ", patient_id_Name[i]))
  
  # Subset
  sceSubset <- subset(sce, , patient_id == patient_id_Name[i])
  
  # t-SNE plot
  tsnePlot <- plotDR(sceSubset, "TSNE", color_by = "merging1", facet_by = "condition") +
    labs(title = patient_id_Name[i]) +
    theme(legend.title = element_blank())
  
  plotName <- paste0("./plot/clusteringPlot_wo_6weeks/tSNE_", patient_id_Name[i], "_Condition_20210118.png")
  
  ggsave(plotName, plot = tsnePlot, width = 20, height = 16, units = "cm")
}

# UMAP --------------------------------------------------------------------

# UMAP: Run at most 1000 cells per sample
set.seed(1984)
sce <- runDR(sce, "UMAP", cells = 400, features = "type")

# UMAP: Plot on all clusters together
umapCluster <- plotDR(sce, "UMAP", color_by = "merging1") +
  theme(legend.title = element_blank())

ggsave("./plot/clusteringPlot_wo_6weeks_typeMarkers/UMAP_Cluster_wo_6weeks_20210202.png", plot = umapCluster, width = 20, height = 16, units = "cm")

# UMAP: Loop for each marker
mkName <- rowData(sce)[,2]

for (i in seq(1:length(mkName))) {
  print(paste0("UMAP plot on ", mkName[i]))
  
  umapPlot <- plotDR(sce, "UMAP", color_by = mkName[i])
  
  plotName <- paste0("./plot/clusteringPlot_wo_6weeks_typeMarkers/UMAP_", mkName[i], "_20210202.png")
  
  ggsave(plotName, plot = umapPlot, width = 20, height = 16, units = "cm")
}

# UMAP Facet: By condition
umapCond <- plotDR(sce, "UMAP", color_by = "merging1", facet_by = "condition") +
  theme(legend.title = element_blank())

ggsave("./plot/clusteringPlot_wo_6weeks_typeMarkers/UMAP_Condition_20210202.png", plot = umapCond, width = 20, height = 16, units = "cm")

# UMAP Facet: By timepoints
umapTime <- plotDR(sce, "UMAP", color_by = "merging1", facet_by = "patient_id") +
  theme(legend.title = element_blank())

umapTime

ggsave("./plot/clusteringPlot/UMAP_Timepoint_20210102.png", plot = umapTime, width = 20, height = 16, units = "cm")

# Differential cell population abundance ----------------------------------

# Get cell counts by cluster & sample
ns <- table(cluster_id = cluster_ids(sce, "merging1"), sample_id = sample_ids(sce))
df_count <- as.data.frame(ns)
View(df_count)
colnames(df_count)[3] <- "cellCount"

write.csv(df_count, "./statistics/CD45Cells_FlowSOM_12x_wo_6weeks_typeMarkers_cellCount.csv")

# Get frequencies by cluster & sample
fq <- prop.table(ns, 2) * 100
df <- as.data.frame(fq)
View(df)
write.csv(df, "./statistics/CD45Cells_FlowSOM_12x_wo_6weeks_typeMarkers_Freq.csv")

# Relative abundance of the cell populations in each sample
## Barplot
freqBarplot <- plotAbundances(sce, k = "merging1", by = "sample_id") +
  theme(legend.title = element_blank())

ggsave("./plot/clusteringPlot_wo_6weeks/Freq_barPlot_Condition_20210131.png", plot = freqBarplot, width = 60, height = 15, units = "cm")

## Boxplot loop for each cell population
setDTthreads(threads = 8)
df <- fread("./statistics/CD45Cells_FlowSOM_20x_Freq.csv", header = T, nThread = 8, verbose = TRUE)
glimpse(df)
df <- df[,-1]

# Load timepoints
CIA_Timepoints <- read.csv("/cloud-data/snf-mgln-immunometabolism/CIA_FH2000/FACS/clustering/CIA_Timepoints.csv")
glimpse(CIA_Timepoints)

# Add timepoint information
df <- df %>% left_join(CIA_Timepoints, by = "sample_id")
glimpse(df)

# Tidy data
df$condition <- factor(df$condition, levels = c("Ctrl", "CIA"))
df$patient_id <- factor(df$patient_id, 
                        levels = c("2 days", "2 weeks", "3 weeks", "4 weeks", 
                                   "6 weeks", "8 weeks", "10 weeks"))

# Loop by cell population
clusterName <- unique(df$cluster_id)
clusterName

for (i in seq(1:length(clusterName))) {
  print(paste0("Population plot on ", clusterName[i]))
  
  # Subset
  dfSubset <- df %>% filter(cluster_id == clusterName[i])

  # Boxplot
  freqBoxplot <- dfSubset %>% ggplot(aes(x = patient_id, y = Freq, fill = condition)) +
    geom_boxplot(aes(color = condition),
                 fill = "white",
                 position = position_dodge(0.9), outlier.shape = NA) +
    geom_point(position = position_jitterdodge(dodge.width = 0.9, seed = 2019), size = 0.5) +
    labs(x = "Day after immunization", y = "Fraction of total population", title = clusterName[i]) +
    theme_bw() +
    theme(legend.title = element_blank()) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2")
  
  # Save
  plotName <- paste0("./plot/clusteringPlot/Freq_boxPlot_cluster", i, "_20210117.png")
  
  ggsave(plotName, plot = freqBoxplot, width = 16, height = 10, units = "cm")
  
  print(paste0(plotName, " is done"))
}

# Plot median state-marker expressions stratified by cell population
## Cell type markers
p <- plotPbExprs(sce, k = "merging1", features = "type", facet_by = "cluster_id", shape_by = NULL) +
  scale_color_brewer(palette = "Dark2")

p$facet$params$ncol <- 2
p

ggsave("./plot/clusteringPlot_wo_6weeks_typeMarkers/Freq_boxPlot_typeMarker_20210202.png", plot = p, width = 36, height = 20, units = "cm")

## Cell state markers
p <- plotPbExprs(sce, k = "merging1", features = "state", facet_by = "cluster_id", shape_by = NULL) +
  scale_color_brewer(palette = "Dark2")

p$facet$params$ncol <- 2
p

ggsave("./plot/clusteringPlot_wo_6weeks_typeMarkers/Freq_boxPlot_stateMarker_20210202.png", plot = p, width = 36, height = 20, units = "cm")

# Distributions of state marker intensities -------------------------------

# Distributions of marker intensities (arcsinh-transformed) of state markers
plotClusterExprs(sce, k = "merging1", features = "state")

# Output expression data for cell population ------------------------------

# List different timepoints
timePoint <- unique(sce$patient_id)
timePoint

for (i in seq(1:length(timePoint))) {

  # Subset
  sceSubset <- subset(sce, , patient_id == timePoint[i])
  print(unique(sceSubset$patient_id))
  
  # Define cell population
  sceSubset$cluster_id <- cluster_ids(sceSubset, "merging1")
  
  # Construct data frame of expression matrix include cell metadata
  cd <- colData(sceSubset)
  es <- assay(sceSubset, "exprs")
  df <- data.frame(t(es), cd, check.names = FALSE)

  # Save csv file
  csvName <- paste0("./clustering/markerExpression_wo_6weeks_", timePoint[i], "_20210202.csv")
  
  setDTthreads(threads = 16)
  fwrite(df, file = csvName, verbose = TRUE)
  
  print(csvName)
}
