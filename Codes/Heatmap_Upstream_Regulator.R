# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: Heatmap on IPA upstream regulators
# Data: z score and p-value from IPA upstream analysis

# Load library ------------------------------------------------------------

library(readxl)
library(data.table)
library(tidyr)
library(dplyr)
library(stringr)

library(RColorBrewer)
library(pheatmap)

# Load data ---------------------------------------------------------------

reg_table <- read_excel("./data/Human_Macrophage_Normoxia_LPS_Upstream_Regulator.xls")
View(reg_table)

colnames(reg_table)[5:6] <- c("zScore", "pval")

# Extract top regulators --------------------------------------------------

# Regulators with biggest positive z score
up_regulator <- reg_table %>% 
  filter(pval < 0.05, `Molecule Type`!= "chemical toxicant") %>% 
  arrange(desc(zScore)) %>% 
  pull(`Upstream Regulators`) %>% 
  unique()

head(up_regulator, 20)

# Regulators with biggest negative z score
down_regulator <- reg_table %>% 
  filter(pval < 0.05, !`Molecule Type`%in% c("chemical drug", "chemical reagent", "fusion gene/product", "biologic drug")) %>% 
  arrange(zScore) %>% 
  pull(`Upstream Regulators`) %>% 
  unique()

head(down_regulator, 20)

# Combine top regulators
combi_regulator <- c(head(up_regulator, 20), head(down_regulator, 20))
combi_regulator

## Export a table for upstream regulators (top 20 up or down)
up_regulator <- reg_table %>% 
  filter(pval < 0.05, `Molecule Type`!= "chemical toxicant") %>% 
  arrange(desc(zScore)) %>%
  slice_head(n = 20) %>%
  select(!`Expr Log Ratio`)
View(up_regulator)

down_regulator <- reg_table %>% 
  filter(pval < 0.05, !`Molecule Type`%in% c("chemical drug", "chemical reagent", "fusion gene/product", "biologic drug")) %>% 
  arrange(zScore) %>%
  slice_head(n = 20) %>%
  select(!`Expr Log Ratio`)
View(down_regulator)

combi_table <- rbind(up_regulator, down_regulator)
combi_table

fwrite(combi_table, "./results/supplement_table_top40_upstream_regulators.csv")

# Prepare matrix ----------------------------------------------------------

mat_input <- reg_table %>%
  filter(pval < 0.05, `Upstream Regulators` %in% combi_regulator) %>%
  select(`Upstream Regulators`, zScore) 

mat <- mat_input %>%
  select(-`Upstream Regulators`) %>%
  as.matrix()

rownames(mat) <- mat_input$`Upstream Regulators`
colnames(mat) <- "LPS"
mat

min(mat, na.rm = TRUE) # -8.166
max(mat, na.rm = TRUE) # 15.202

mat_row_ordered <- as.matrix(mat[combi_regulator, ])
mat_row_ordered

colnames(mat_row_ordered) <- "LPS"

# Heatmap with NA value ---------------------------------------------------

# Break for gradient
breakList <- seq(-10, 10, length.out = 100)

# Heatmap
pheatmap(mat_row_ordered,
  na_col = "grey45", scale = "none", cluster_cols = FALSE, cluster_rows = FALSE, angle_col = 45,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "PiYG")))(length(breakList))
)

# Export TIFF 200 x 1360 heatmap_normoxia_upstream_regulator.tiff

