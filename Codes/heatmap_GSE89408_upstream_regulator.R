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
BH_pValue <- read_excel("data/GSE89408_upstream_regulator_BH_pValue.xls")
View(BH_pValue)

# Table of z score
z_score <- read_excel("data/GSE89408_upstream_regulator_zScore.xls", na = "N/A")
View(z_score)

# Molecule type
mol_type <- read_excel("data/GSE89408_upstream_regulator_early_RA_molecule_type.xls")
View(mol_type)

# Merge table -------------------------------------------------------------

# Change to long format
bh_long <- BH_pValue %>%
  pivot_longer(cols = 2:4, names_to = "Timepoint", values_to = "bh_pval")
View(bh_long)

zScore_long <- z_score %>%
  pivot_longer(cols = 2:4, names_to = "Timepoint", values_to = "zScore")
View(zScore_long)

# Merge table
reg_table <- full_join(bh_long, zScore_long, by = c("Upstream Regulators", "Timepoint")) %>%
  rename(`Upstream Regulator` = `Upstream Regulators`) %>%
  left_join(mol_type[, c(1, 3)], by = "Upstream Regulator")
View(reg_table)

# Extract top regulators --------------------------------------------------

# Regulators with biggest positive z score
up_regulator <- reg_table %>% 
  filter(bh_pval > -log10(0.05), !`Molecule Type` %in% c("chemical - endogenous non-mammalian", "chemical reagent", "chemical toxicant")) %>% 
  arrange(desc(zScore)) %>% 
  pull(`Upstream Regulator`) %>% 
  unique()

head(up_regulator, 10)
head(up_regulator, 21)

# Regulators with biggest negative z score
down_regulator <- reg_table %>%
  filter(
    bh_pval > -log10(0.05),
    !`Molecule Type` %in% c(
      "chemical - protease inhibitor",
      "chemical drug", "chemical reagent", "biologic drug",
      "chemical - endogenous mammalian", "chemical - kinase inhibitor"
    )
  ) %>%
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
  filter(`Upstream Regulator` %in% head(up_regulator, 21)[-19]) %>%
  mutate(BH_pval = 10^-bh_pval) %>%
  arrange(desc(zScore)) %>%
  select(-bh_pval) %>%
  pivot_wider(
    names_from = Timepoint,
    names_sep = ".",
    values_from = c(zScore, BH_pval)
  )
View(up_regulator_table)

down_regulator_table <- reg_table %>%
  filter(`Upstream Regulator` %in% c(head(down_regulator, 20), "SIRT1")) %>%
  mutate(BH_pval = 10^-bh_pval) %>%
  arrange(zScore) %>%
  select(-bh_pval) %>%
  pivot_wider(
    names_from = Timepoint,
    names_sep = ".",
    values_from = c(zScore, BH_pval)
  )
View(down_regulator_table)

combi_table <- rbind(up_regulator_table, down_regulator_table)
combi_table

fwrite(combi_table, "./results/supplement_table_top40_upstream_regulators_hs_synovium.csv")

# Prepare matrix ----------------------------------------------------------

mat_input <- reg_table %>%
  filter(bh_pval > -log10(0.05), `Upstream Regulator` %in% combi_regulator) %>%
  dplyr::select(-bh_pval, -`Molecule Type`) %>%
  spread(Timepoint, zScore)

mat <- mat_input[, c(1, 2, 4, 3)] %>%
  dplyr::select(-`Upstream Regulator`) %>%
  as.matrix()

rownames(mat) <- mat_input$`Upstream Regulator`
mat

mat_row_ordered <- mat[combi_regulator,]
mat_row_ordered

min(mat, na.rm = TRUE) # -5.331
max(mat, na.rm = TRUE) # 9.43

# Heatmap with NA value ---------------------------------------------------

# Break for gradient
breakList <- seq(-10, 10, length.out = 100)

# Heatmap
pheatmap(mat_row_ordered,
  scale = "none", cluster_cols = FALSE, cluster_rows = FALSE,
  angle_col = 90, na_col = "grey45",
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "PiYG")))(length(breakList))
)

# Export TIFF 260 x 680: heatmap_upstream_regulator_v1.tiff

