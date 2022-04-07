# Introduction ------------------------------------------------------------

# Image: metabolomicsr
# Output: mSet object
# Data: Tidy tMSI raw data combined with metadata

# Load library ------------------------------------------------------------

library(data.table)
library(MetaboAnalystR)

# Double check input raw data ---------------------------------------------

textData <- fread("./data/DMD1659/DMD1659_sum_mz_v2_raw_tidy.csv")
View(textData)
dim(textData) # 44 x 8177

head(colnames(textData)) #"mouseID", "Treatment", "Time", "meta_1", "meta_10", "meta_100"

# Prepare data for time-series/two-factor design module -------------------

# Create objects for storing processed data for time-series analysis
mSet <- InitDataObjects("conc", "ts", paired = FALSE)

mSet <- SetDesignType(mSet, "time") ## The time points group must be labeled as <b>Time</b>

mSet <- Read.TextData(
  mSetObj = mSet, "./data/DMD1659/DMD1659_sum_mz_v2_raw_tidy.csv", "rowts", "disc"
)

View(mSet)

# Sanity check
mSet <- SanityCheckData(mSet) # A total of 0 (0%) missing values were detected

# Minimum value replacing
mSet <- ReplaceMin(mSet)

# Check if the sample size is too small
mSet <- IsSmallSmplSize(mSet)

# Normalization and scaling -----------------------------------------------

mSet <- PreparePrenormData(mSet)

# Log Normalization and Mean Centering
mSet <- Normalization(mSet, "NULL", "LogNorm", "MeanCenter", ratio = FALSE, ratioNum = 20)

# Plot sample normalization summary
mSet <- PlotSampleNormSummary(mSet, "./plot/DMD1659/normSample_timeseries_", 
                              "png", 300, width = NA)

# Save mSet ---------------------------------------------------------------

saveRDS(mSet, file = "./data/DMD1659/sum_mz_mSet_MetaboAnalyst_norm.rds")

