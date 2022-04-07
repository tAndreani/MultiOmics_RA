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

# Load DESeq2 results on CIA paws
csvPath <- list.files("~/cloud-data/snf-mgln-immunometabolism/CIA_FH2000/CIA_multiomics_figure/result/Annotated_DESeq2_Wald_Test_Result",
                      full.names = TRUE
)

paw <- map_df(csvPath, read.csv)
View(paw)

# Load DESeq2 results on human synovium
csvPath <- list.files("./results", pattern = "deseq_wald_test_annotated_", full.names = TRUE)
csvPath

synovium <- map_df(csvPath, read.csv)
View(synovium)

unique(synovium$condition)

# Subset table for conditions ---------------------------------------------

# Subset for 3 weeks of CIA samples
paw_3week <- paw %>%
  filter(Timepoint == "3_weeks") %>%
  mutate(symbol = str_to_upper(symbol))
View(paw_3week)

# List up-regulated genes -------------------------------------------------

# List of up-regulated genes at 3 weeks of CIA
CIA_3weeks_up <- paw_3week %>%
  filter(is.na(entrez) == FALSE) %>%
  filter(padj < 0.05 & log2FoldChange > log2(1.5)) %>%
  pull(symbol) %>%
  unique()

head(CIA_3weeks_up)

# List of up-regulated genes of human synoviums
synovium_up <- synovium %>%
  filter(is.na(entrez) == FALSE) %>%
  filter(padj < 0.05 & log2FoldChange > log2(1.5)) %>%
  unique()
head(synovium_up)

# Create a combined list for up-regulated genes
list_up <- list(
  CIA_3weeks_up = CIA_3weeks_up,
  Arthralgia = synovium_up[synovium_up$condition == "Arthralgia", ]$symbol,
  UA = synovium_up[synovium_up$condition == "UA", ]$symbol,
  Early_RA = synovium_up[synovium_up$condition == "Early RA", ]$symbol
)

# Rename list 
names(list_up)[1] <- "CIA (3 wks)"
names(list_up)[4] <- "Early RA"

# List down-regulated genes -----------------------------------------------

# List of down-regulated genes at 3 weeks of CIA
CIA_3weeks_down <- paw_3week %>%
  filter(is.na(entrez) == FALSE) %>%
  filter(padj < 0.05 & log2FoldChange < log2(1 / 1.5)) %>%
  pull(symbol) %>%
  unique()

head(CIA_3weeks_down)

# List of down-regulated genes of human synoviums
synovium_down <- synovium %>%
  filter(is.na(entrez) == FALSE) %>%
  filter(padj < 0.05 & log2FoldChange < log2(1 / 1.5)) %>%
  unique()
head(synovium_down)

# Create a combined list for down-regulated genes
list_down <- list(
  CIA_3weeks_down = CIA_3weeks_down,
  Arthralgia = synovium_down[synovium_down$condition == "Arthralgia", ]$symbol,
  UA = synovium_down[synovium_down$condition == "UA", ]$symbol,
  Early_RA = synovium_down[synovium_down$condition == "Early RA", ]$symbol
)

# Rename list 
names(list_down)[1] <- "CIA (3 wks)"
names(list_down)[4] <- "Early RA"

# Venn diagram ------------------------------------------------------------

# Venn diagram for significantly up-regulated genes
p1 <- ggVennDiagram(list_up, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "red", guide = "none")

ggsave("./plot/venn_up_CIA.tif", p1,
  device = "tiff", units = "mm",
  width = 320, height = 160, dpi = 300
)

# Venn diagram for significantly down-regulated genes
p2 <- ggVennDiagram(list_down, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "dodgerblue", guide = "none")

ggsave("./plot/venn_down_CIA.tif", p2, 
       device = "tiff", units = "mm",
       width = 300, height = 120, dpi = 300)

# Save shared significant up-regulated genes ------------------------------

# Intersect CIA_3weeks_up and early RA
up_symbol <- intersect(CIA_3weeks_up, synovium_up[synovium_up$condition == "Early RA", ]$symbol)
head(up_symbol)

# Save ensembl IDs
up_ensembl <- synovium %>% 
  filter(symbol %in% up_symbol) %>% 
  pull(gene_id) %>% 
  unique()

saveRDS(up_ensembl, file = "./results/shared_upregulated_genes_ensembl.rds")

# Save Entrez IDs
up_entrez <- synovium %>% 
  filter(symbol %in% up_symbol) %>% 
  pull(entrez) %>% 
  unique()

saveRDS(up_entrez, file = "./results/shared_upregulated_genes_entrez.rds")
