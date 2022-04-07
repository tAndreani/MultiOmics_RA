# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: Dumbbell chart on genes with largest weights
# Data: MOFA model on RNAseq and metabolomics data

# Load library ------------------------------------------------------------

library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(ggrepel)

library(AnnotationDbi)
library(org.Mm.eg.db)

# Set plot parameter ------------------------------------------------------

theme_set(theme_bw())
theme_update(
  text = element_text(family = "Helvetica", size = 9),
  plot.title = element_text(size = 9),
  legend.title = element_text(size = 9),
  legend.text = element_text(size = 9),
  legend.key.size = unit(4, "mm")
)

# Get weight --------------------------------------------------------------

allWeights <- get_weights(model, views = "all", factors = "all", as.data.frame = TRUE) %>%
  as_tibble()
head(allWeights)

fwrite(allWeights, "./result/feature_weights.csv")

# Map ENSEMBL IDs to Entrez IDs -------------------------------------------

# Tidy up Ensembl IDs
result_mRNA <- allWeights %>%
  filter(view == "mRNA") %>%
  separate(feature, c("ensembl_gene_id", "version"))
View(result_mRNA)

# Mapping to Entrez ID
result_mRNA$entrez <- mapIds(org.Mm.eg.db,
                             keys = result_mRNA$ensembl_gene_id,
                             keytype = "ENSEMBL", column = "ENTREZID",
                             multiVals = "first"
)

print(paste0(sum(!is.na(result_mRNA$entrez)), " genes mapped")) # 27966 genes mapped

# Mapping to gene symbol
result_mRNA$symbol <- mapIds(org.Mm.eg.db,
                             keys = result_mRNA$ensembl_gene_id,
                             keytype = "ENSEMBL", column = "SYMBOL",
                             multiVals = "first"
)

print(paste0(sum(!is.na(result_mRNA$symbol)), " genes mapped")) # 27966 genes mapped

# Save result table
fwrite(result_mRNA, file = "./result/mRNA_feature_weights_annotated.csv")

# Prepare table -----------------------------------------------------------

# Load result table
result_mRNA <- fread("./result/mRNA_feature_weights_annotated.csv")
View(result_mRNA)

# Define factor
fc <- "Factor1"

# Subset genes from selected factor
mat <- result_mRNA %>%
  dplyr::filter(view == "mRNA", factor == fc, symbol != "") %>%
  arrange(value) %>%
  mutate(symbol = str_to_upper(symbol))
View(mat)

# Extract features with largest weight
W <- mat$value
names(W) <- mat$symbol

featureName <- c(names(W)[tail(order(W), n = 10)], # Top 10 negative features
                 names(W)[head(order(W), n = 10)]) # Top 10 positive features
featureName

# Dumbbell chart ----------------------------------------------------------

# Rank features by factor values
df <- mat %>%
  filter(symbol %in% featureName) %>%
  arrange(value)

# Define the order of features
sorted_level <- df$symbol
sorted_level

# Dumbbell chart
p <- df %>%
  ggplot(aes(x = value, y = symbol)) +
  geom_segment(aes(x = 0, xend = value, y = symbol, yend = symbol), size = 1, color = "#7570B3") +
  labs(x = "Feature weight", y = NULL, title = fc) +
  scale_y_discrete(limits = sorted_level)

ggsave("./plot/dumbbellPlot_F1_20_Gene_v2.tif", plot = p, device = "tiff",
       units = "mm", width = 90, height = 90, dpi = 300)
