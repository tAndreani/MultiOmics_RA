# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: Venn diagram on DEG genes
# Data: Annotated DESeq2 result

# Load library ------------------------------------------------------------

library(data.table)
library(purrr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggVennDiagram)

# Set plot parameter ------------------------------------------------------

theme_set(theme_bw())
theme_update(
  text = element_text(family = "Helvetica", size = 7),
  axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
  plot.title = element_text(size = 7),
  legend.title = element_text(size = 7),
  legend.key.size = unit(4, "mm")
)

# Load data ---------------------------------------------------------------

# Load DESeq2 results on macrophages under normoxia
mac_norm <- fread("./results/Annotated_deseq_wald_test_normoxia.csv")
View(mac_norm)

# Load DESeq2 results on CIA paws
csvPath <- list.files("~/cloud-data/snf-mgln-immunometabolism/CIA_FH2000/CIA_multiomics_figure/result/Annotated_DESeq2_Wald_Test_Result",
  full.names = TRUE
)

paw <- map_df(csvPath, read.csv)
View(paw)

# Tidy table --------------------------------------------------------------

# Subset for 3 weeks
paw_3week <- paw %>%
  filter(Timepoint == "3_weeks") %>%
  mutate(symbol = str_to_upper(symbol))
View(paw_3week)

# List of up-regulated genes at 3 weeks of CIA
CIA_3weeks_up <- paw_3week %>%
  filter(is.na(entrez) == FALSE) %>%
  filter(padj < 0.05 & log2FoldChange > log2(1.5)) %>%
  pull(symbol) %>%
  unique()

# List of up-regulated genes on macrophages under normoxia
mac_norm_up <- mac_norm %>%
  filter(is.na(entrez) == FALSE) %>%
  filter(padj < 0.05 & log2FoldChange > log2(1.5)) %>%
  pull(symbol) %>%
  unique()

# Create a list for up-regulated genes
up <- list(CIA_3weeks_up = CIA_3weeks_up, mac_norm_up = mac_norm_up)
names(up) <- c("CIA paw (3 weeks)", "LPS stimulated macrophage")
View(up)

# List of down-regulated genes at 3 weeks of CIA
CIA_3weeks_down <- paw_3week %>%
  filter(is.na(entrez) == FALSE) %>%
  filter(padj < 0.05 & log2FoldChange < log2(1 / 1.5)) %>%
  pull(symbol) %>%
  unique()

# List of down-regulated genes on macrophages under normoxia
mac_norm_down <- mac_norm %>%
  filter(is.na(entrez) == FALSE) %>%
  filter(padj < 0.05 & log2FoldChange < log2(1 / 1.5)) %>%
  pull(symbol) %>%
  unique()

# Create a list for down-regulated genes
down <- list(CIA_3weeks_down = CIA_3weeks_down, mac_norm_down = mac_norm_down)
names(down) <- c("CIA paw (3 weeks)", "LPS stimulated macrophage")
View(down)

# Venn diagram ------------------------------------------------------------

# Venn diagram for significantly up-regulated genes
p1 <- ggVennDiagram(up, label_alpha = 0, label_percent_digit = 1) +
  scale_fill_gradient(low = "white", high = "red", guide = "none")

ggsave("./plot/venn_up_normoxia_v1.tif", p1,
  device = "tiff", units = "mm",
  width = 200, height = 100, dpi = 300
)

# Venn diagram for significantly down-regulated genes
p2 <- ggVennDiagram(down, label_alpha = 0, label_percent_digit = 1) +
  scale_fill_gradient(low = "white", high = "cyan", guide = "none")

ggsave("./plot/venn_down_normoxia_v1.tif", p2,
  device = "tiff", units = "mm",
  width = 200, height = 100, dpi = 300
)

# Save shared significant up-regulated genes ------------------------------

up_symbol <- intersect(CIA_3weeks_up, mac_norm_up)

# Save ensembl IDs
up_ensembl <- mac_norm %>% 
  filter(symbol %in% up_symbol) %>% 
  pull(ensembl_gene_id) %>% 
  unique()

saveRDS(up_ensembl, file = "./results/shared_upregulated_genes_ensembl.rds")

# Save Entrez IDs
up_entrez <- mac_norm %>% 
  filter(symbol %in% up_symbol) %>% 
  pull(entrez) %>% 
  unique()

saveRDS(up_entrez, file = "./results/shared_upregulated_genes_entrez.rds")
