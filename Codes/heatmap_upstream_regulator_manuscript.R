# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: Heatmap on IPA upstream regulators
# Data: z score and BH p-value (-log10) from IPA upstream analysis

# Load library ------------------------------------------------------------

library(readxl)
library(data.table)
library(tidyr)
library(dplyr)
library(stringr)

library(RColorBrewer)
library(pheatmap)

# Load data ---------------------------------------------------------------

# Table of -log10 BH p-values
BH_pValue <- read_excel("data/paw_RNAseq/CIA_Paw_RNAseq_Upstream_Regulator_BH_pValue.xls")
View(BH_pValue)

# Table of z score
z_score <- read_excel("data/paw_RNAseq/CIA_Paw_RNAseq_Upstream_Regulator_z_Score.xls",
                      na = "NA"
)
View(z_score)

# Molecule type
mol_type <- read_excel("data/paw_RNAseq/CIA_Paw_RNAseq_Upstream_Regulator_3_week_Molecule_Type.xls")
View(mol_type)

# Merge table -------------------------------------------------------------

# Change to long format
bh_long <- BH_pValue %>%
  pivot_longer(cols = 2:9, names_to = "Timepoint", values_to = "bh_pval")
View(bh_long)

zScore_long <- z_score %>%
  pivot_longer(cols = 2:9, names_to = "Timepoint", values_to = "zScore")
View(zScore_long)

# Merge table
reg_table <- full_join(bh_long, zScore_long, by = c("Upstream Regulators", "Timepoint")) %>%
  mutate(Timepoint = str_replace_all(Timepoint, "_", " ")) %>%
  rename(`Upstream Regulator` = `Upstream Regulators`) %>%
  left_join(mol_type[, c(1, 3)], by = "Upstream Regulator")
View(reg_table)

# Extract top regulators --------------------------------------------------

# Regulators with biggest positive z score
up_regulator <- reg_table %>% 
  filter(bh_pval > -log10(0.05)) %>% 
  arrange(desc(zScore)) %>% 
  pull(`Upstream Regulator`) %>% 
  unique()

head(up_regulator, 10)
head(up_regulator, 23)

# Regulators with biggest negative z score
down_regulator <- reg_table %>% 
  filter(bh_pval > -log10(0.05), 
         !`Molecule Type`%in% c("fusion gene/product", "chemical drug", "chemical - kinase inhibitor")) %>% 
  arrange(zScore) %>% 
  pull(`Upstream Regulator`) %>% 
  unique()

head(down_regulator, 10)
head(down_regulator, 20)

# Combine top regulators
combi_regulator <- c(head(up_regulator, 10), head(down_regulator, 10))
combi_regulator

## Export a table for upstream regulators (top 20 up or down)
up_regulator_table <- reg_table %>%
  filter(`Upstream Regulator` %in% head(up_regulator, 23)[-c(12, 13, 16)]) %>%
  arrange(desc(zScore)) %>%
  pivot_wider(
    names_from = Timepoint,
    names_sep = ".",
    values_from = c(zScore, bh_pval)
  )
View(up_regulator_table)

down_regulator_table <- reg_table %>% 
  filter(`Upstream Regulator` %in% c(head(down_regulator, 20))) %>% 
  arrange(zScore) %>%
  pivot_wider(
    names_from = Timepoint,
    names_sep = ".",
    values_from = c(zScore, bh_pval)
  )
View(down_regulator_table)

combi_table <- rbind(up_regulator_table, down_regulator_table)
combi_table

fwrite(combi_table, "./result/supplement_table_top40_upstream_regulators_cia_paw.csv")
  
# Prepare matrix ----------------------------------------------------------

# Spread table
mat_input <- reg_table %>%
  filter(bh_pval > -log10(0.05), `Upstream Regulator` %in% combi_regulator) %>%
  select(-bh_pval) %>%
  spread(Timepoint, zScore) %>%
  mutate(`2 days` = NA)

mat_input$`2 days` <- as.numeric(mat_input$`2 days`)

# Create matrix 
mat <- mat_input[, c(1, 10, 4:9, 3)] %>%
  select(-`Upstream Regulator`) %>%
  as.matrix()

rownames(mat) <- mat_input$`Upstream Regulator`
mat

mat_row_ordered <- mat[combi_regulator,]
mat_row_ordered

min(mat_row_ordered, na.rm = TRUE) # -7.821
max(mat_row_ordered, na.rm = TRUE) # 13.224

# Heatmap with NA value ---------------------------------------------------

# Break for gradient
breakList <- seq(-10, 10, length.out = 100)

# Heatmap
pheatmap(mat_row_ordered,
  scale = "none", cluster_cols = FALSE, cluster_rows = FALSE,
  angle_col = 45, na_col = "grey45",
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "PiYG")))(length(breakList))
)

# Export TIFF 490 x 680: heatmap_CIA_upstream_regulator_v1.tiff
