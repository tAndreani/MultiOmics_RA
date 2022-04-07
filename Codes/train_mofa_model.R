# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: MOFA model on RNAseq and metabolomics data
# Data: DESeq2 object after variance stabilizing transformation (VST)
# Data: Metabolon peak area data

# Reference: https://huber-group-embl.github.io/mofaCLL/analysisProcedure.html
# Updated tutorial: https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/CLL.html

# Set environment variables -----------------------------------------------

Sys.setenv(PATH = paste("/cloud-home/I0442220/.cache/basilisk/1.2.1/MOFA2-1.0.1/mofa_env",
                        Sys.getenv()["PATH"],
                        sep = ";"
))

# Load library ------------------------------------------------------------

library(data.table)
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)

library(DESeq2)
library(MOFA2)
library(ggplot2)
library(RColorBrewer)
library(corrplot)

# Prepare RNAseq dataset --------------------------------------------------

# Load pre-filtered and VST transformed data
rna.vst <- readRDS("~/cloud-data/snf-mgln-immunometabolism/CIA_FH2000/CIA_multiomics_figure/data/paw_RNAseq/DESeq2_preFilter_vst_batchRm_object.rds")
nrow(rna.vst) # 16803 genes

# Select top 5000 most variable genes
exprMat <- assay(rna.vst)
nTop <- 5000
sds <- genefilter::rowSds(exprMat)
exprMat <- exprMat[order(sds, decreasing = T)[1:nTop], ]

View(exprMat) # Rows: ENSEMBL gene names, columns: mouse IDs

# Data distribution of RNAseq data
boxplot(exprMat, outline = FALSE, col = "cornflowerblue", main = "Transformed RNAseq data")

# Tidy up column names
colnames(exprMat) <- str_split_fixed(colnames(exprMat), "_", 2)[, 1]

# Prepare metabolomics dataset --------------------------------------------

# Load peak area data
peakData <- fread("~/cloud-data/snf-mgln-immunometabolism/CIA_FH2000/Metabolomics/data/Batch_Norm_Impute_Intensity/Plasma_Batch_Norm_Impute_Intensity_with_ChemAnno_Score.csv")
View(peakData)

# Tidy data
metaMat <- peakData %>%
  select(mouseID, Intensity, CHEMICAL_NAME) %>%
  mutate(mouseID = paste0("X", mouseID), Intensity = log(Intensity)) %>%
  # Log transformation
  spread(key = mouseID, value = Intensity) %>%
  column_to_rownames("CHEMICAL_NAME") %>%
  as.matrix()

View(metaMat) # Rows: metabolites, columns: mouse IDs

# Data distribution of metabolomics data
boxplot(metaMat, outline = FALSE, col = "green", main = "Transformed metabolomics data")

# Create the MOFA obejct --------------------------------------------------

# List of data matrix
mofaData <- list(Metabolites = metaMat, mRNA = exprMat)
lapply(mofaData, dim)

# Extract samples that appears in both assays
sampleList <- lapply(mofaData, colnames)
useSamples <- intersect(sampleList$Metabolites, sampleList$mRNA)

# Only keep samples that appears in both assays
f2 <- function(x) x[, useSamples]
mofaData <- lapply(mofaData, f2)

# Build MOFA object -------------------------------------------------------

# Create the MOFA object
MOFAobject <- create_mofa(mofaData)
MOFAobject

# Plot sample size of each modality
fig_data <- plot_data_overview(MOFAobject) +
  scale_fill_brewer(palette = "Set2")

ggsave("./plot/plot_data_overview.png",
  plot = fig_data,
  width = 15, height = 10, units = "cm"
)

# Setup MOFA training parameters ------------------------------------------

# List data options
DataOptions <- get_default_data_options(MOFAobject)
DataOptions

# Define model options
ModelOptions <- get_default_model_options(MOFAobject)
ModelOptions$num_factors <- 25 # Number of factors
ModelOptions

# Define training options
TrainOptions <- get_default_training_options(MOFAobject)
TrainOptions$drop_factor_threshold <- 0.02
TrainOptions$convergence_mode <- "slow"
TrainOptions$seed <- 1234
TrainOptions$verbose <- TRUE
TrainOptions

# Prepare the MOFA object -------------------------------------------------

MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = DataOptions,
                           model_options = ModelOptions,
                           training_options = TrainOptions
)

# Train the MOFA model ----------------------------------------------------

# Define file path for saving the trained model
outfile <- file.path("./result/model_CIA.hdf5")

# Train the MOFA model
MOFAobject <- run_mofa(MOFAobject, outfile, use_basilisk = TRUE)

# Variance explained by MOFA for each omic data ---------------------------

# Load model
model <- load_model("./result/model_CIA.hdf5")
model

# Calculate the variance explained (R2) per factor in each view
calculate_variance_explained(model)

# Plot variance explained
theme_set(theme_classic())
theme_update(
  text = element_text(family = "Helvetica", size = 7),
  plot.title = element_text(size = 7),
  legend.title = element_text(size = 7),
  legend.text = element_text(size = 7),
  legend.key.size = unit(4, "mm")
)

fig_v1 <- plot_variance_explained(model) +
  scale_fill_viridis_c(breaks = c(5, 10, 15, 20, 25))

ggsave("./plot/plot_var_per_factor_viridis.png", plot = fig_v1, width = 15, height = 10, units = "cm")

ggsave("./plot/plot_var_per_factor_viridis.tif",
  plot = fig_v1, device = "tiff",
  units = "mm", width = 90, height = 45, dpi = 300
)

# Total variance
theme_set(theme_bw())
theme_update(
  text = element_text(family = "Helvetica", size = 7),
  plot.title = element_text(size = 7),
  legend.title = element_text(size = 7),
  legend.text = element_text(size = 7),
  legend.key.size = unit(4, "mm")
)

fig_v2 <- plot_variance_explained(model, plot_total = T)[[2]] +
  scale_y_continuous(breaks = seq(0, 70, 10), limits = c(0, 70))

ggsave("./plot/plot_var_per_view.png", plot = fig_v2, width = 15, height = 10, units = "cm")

ggsave("./plot/plot_var_per_view.tif",
       plot = fig_v2, device = "tiff",
       units = "mm", width = 90, height = 45, dpi = 300
)

# Sanity check: Factors should be uncorrelated ----------------------------

# Get factor data
Z <- get_factors(model)

# Compute correlation
r <- abs(cor(
  x = do.call(rbind, Z), y = do.call(rbind, Z),
  method = "pearson", use = "complete.obs"
))

# Corrplot
png(
  filename = "./plot/factor_corrPlot.png",
  width = 600, height = 600, units = "px"
)
corrplot(r,
         method = "color", type = "upper", order = "hclust",
         col = rev(brewer.pal(n = 8, name = "RdYlBu")),
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, tl.cex = 1, # Text label color, rotation and size
         diag = FALSE
)
dev.off()
