# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: Violin plot correlating factor values with arthritis score
# Data: MOFA model on RNAseq and metabolomics data

# Load library ------------------------------------------------------------

library(data.table)
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)
library(ggpubr)

library(MOFA2)
library(ggplot2)
library(psych)

# Load model --------------------------------------------------------------

model <- load_model("./result/model_CIA.hdf5")
model

# Set plot theme ----------------------------------------------------------

theme_update(
  text = element_text(family = "Helvetica", size = 7),
  plot.title = element_text(size = 7),
  legend.title = element_text(size = 7),
  legend.text = element_text(size = 7),
  legend.key.size = unit(4, "mm"),
  legend.position = "bottom"
)

# Add sample metadata to the model ----------------------------------------

# Load DESeq2 colData with individual paw scores at sacrifice
DESeq2_colData <- read.csv("~/cloud-data/snf-mgln-immunometabolism/CIA_FH2000/CIA_multiomics_figure/data/metadata/DESeq2_colData.csv")
View(DESeq2_colData)

colnames(DESeq2_colData)[6] <- "RNAseq_Paw_Score"

# Load total arthritis score at sacrifice
score <- fread("~/cloud-data/snf-mgln-immunometabolism/CIA_FH2000/CIA_multiomics_figure/data/metadata/CIA_Ctrl_Score_Sum_Sacrifice_20210118.csv")
View(score)

# Merge data and subset for samples used in MOFA
sample_metadata <- score %>%
  dplyr::rename(Total_Score = Score) %>%
  left_join(DESeq2_colData, by = c("mouseID", "Treatment")) %>%
  mutate(sample = paste0("X", mouseID)) %>%
  filter(sample %in% samples_names(model)$group1) %>%
  select(sample, Treatment, Timepoint, RNAseq_Paw_Score, Total_Score, Batch)
View(sample_metadata)

dim(sample_metadata) # 171 x 6

rowSums(is.na(sample_metadata)) # 0 NA value

# Check if the sample names fit
all(sort(sample_metadata$sample) == sort(unlist(samples_names(model))))

# Add metadata to the model
samples_metadata(model) <- sample_metadata

head(model@samples_metadata, n = 3)

# Get factor values -------------------------------------------------------

# Get factors
allFactors <- get_factors(model, factors = "all", as.data.frame = TRUE) %>%
  as_tibble()
head(allFactors)

fwrite(allFactors, "./result/sample_factor_values.csv")

# Merge factor values with metadata
plotTab <- allFactors %>%
  left_join(as_tibble(model@samples_metadata), by = c("sample", "group"))
fwrite(plotTab, "./result/sample_factor_values_metadata.csv")

# Violin plot to correlate factor values with RNAseq paw scores -----------

# Facet by all factors
fig_F8 <- plot_factor(model, factors = 1:3, color_by = "RNAseq_Paw_Score", dot_size = 1, dodge = T, add_violin = T, violin_alpha = 0.25) + 
  scale_fill_discrete(name = "Paw Score") +
  ylim(-6, 15) +
  facet_wrap(vars(factor))

ggsave("./plot/dotPlot_factors_vs_RNAseq_paw_score.tif",
       plot = fig_F8, device = "tiff",
       units = "mm", width = 110, height = 60, dpi = 300
)

# ANOVA to compare groups with RNAseq paw scores --------------------------

plotTab$RNAseq_Paw_Score <- factor(plotTab$RNAseq_Paw_Score)

for (i in unique(plotTab$factor)) {
  
  df <- plotTab %>%
    filter(factor == i)
  
  print(i)
  
  res.aov <- aov(value ~ RNAseq_Paw_Score, data = df)
  print(summary(res.aov))
  
}

# Factor1: <2e-16 ***, Factor2: 0.282, Factor3: 2.7e-08 ***
