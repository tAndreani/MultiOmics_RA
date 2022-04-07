# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: Dumbbell chart on significant GO pathways
# Data: mRNA feature weight from MOFA model

# Load library ------------------------------------------------------------

library(gage)
library(fgsea)

library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)

# Set plot parameter ------------------------------------------------------

theme_set(theme_bw())
theme_update(
  text = element_text(family = "Helvetica", size = 8),
  axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
  plot.title = element_text(size = 8),
  legend.title = element_text(size = 7),
  legend.text = element_text(size = 7),
  legend.key.size = unit(2, "mm"),
  legend.spacing = unit(0, "mm"),
  legend.position = "bottom",
  legend.direction = "vertical"
)

# Generate up-to-date GO (Gene Ontology) gene sets ------------------------

## https://rdrr.io/bioc/gage/man/go.gsets.html
go.mmu <- go.gsets(species = "Mouse")

# Subset only biological process
go.bp <- go.mmu$go.sets[go.mmu$go.subs$BP]
View(go.bp)

length(go.bp) # 16025 GO biological process pathways

# Prepare gene list for GSEA on selected Factor 

# Define factor
fc <- "Factor1"

# Load weight data
input <- fread("./result/mRNA_feature_weights_annotated.csv")
View(input)

# Remove rows with unmapped Entrez IDs
input_dropNA <- input[!is.na(input$entrez), ]

# Select factor and remove Entrez IDs that are mapped to multiple genes
res_factor <- input_dropNA %>%
  filter(factor == fc) %>%
  group_by(entrez) %>%
  filter(n() == 1)

print(paste0(nrow(res_factor), " unique genes as GSEA input")) # 4661 unique genes

# Prepare geneList by sorting fold change in decreasing order
deseq2.fc <- res_factor$value
names(deseq2.fc) <- res_factor$entrez
deseq2.fc <- sort(deseq2.fc, decreasing = TRUE)
print(tail(deseq2.fc))

# Save as geneList
save(deseq2.fc, file = "./result/geneList_F3.csv")

# fgsea on metabolism pathways
fgsea_met <- fgseaMultilevel(pathways = go.bp, stats = deseq2.fc, minSize = 12, maxSize = 600) # GO minSize = 20

print(paste0("Number of significant metabolic pathways (padj < 0.05): ", sum(fgsea_met$padj < 0.05, na.rm = TRUE)))

View(fgsea_met)

# Save result on metabolism pathways
df_met <- fgsea_met %>%
  as.data.frame()

fwrite(df_met, file = "./result/GSEA_go_bp_F1_result.csv")

# Dumbbell charts on GSEA pathways ----------------------------------------

# Load GSEA results
df_met <- fread("./result/GSEA_go_bp_F1_result.csv")
View(df_met)

# For manuscript
pathwayList <- c(
  "GO:0006096 glycolytic process",
  "GO:0006099 tricarboxylic acid cycle",
  "GO:0006119 oxidative phosphorylation",
  "GO:0022904 respiratory electron transport chain",
  "GO:0006809 nitric oxide biosynthetic process",
  "GO:0072593 reactive oxygen species metabolic process",
  "GO:0002224 toll-like receptor signaling pathway",
  "GO:0042119 neutrophil activation",
  "GO:0030225 macrophage differentiation",
  "GO:0050900 leukocyte migration",
  "GO:0072676 lymphocyte migration"
)

# Tidy table
gsea_table <- df_met %>%
  filter(pathway %in% pathwayList, padj < 0.05) %>%
  mutate(
    pathwayName = str_to_sentence(substring(pathway, 12)),
    minusLog10padj = -log10(padj),
    Direction = case_when(NES < 0 ~ "Depleted", TRUE ~ "Enriched")
  ) %>%
  arrange(padj)
View(gsea_table)

# Define the order of the pathway names
sorted_level <- gsea_table %>%
  dplyr::pull(pathwayName)

sorted_level <- sorted_level[c(1:7, 12:13)]

# Dumbbell chart
p <- gsea_table %>%
  filter(pathwayName %in% sorted_level) %>%
  ggplot(aes(x = minusLog10padj, y = pathwayName)) +
  geom_segment(aes(x = 0, xend = minusLog10padj, y = pathwayName, yend = pathwayName), size = 1) +
  geom_point(aes(color = Direction, size = size)) +
  labs(x = "-log10(FDR)", y = NULL, size = "Gene module size",
       title = "Factor 1: Gene Ontology biological process", subtitle = "FDR < 0.05") +
  scale_color_manual(values = c("#74ADD1", "#F46D43")) +
  scale_y_discrete(limits = rev(sorted_level)) +
  scale_size_continuous(breaks = c(25, 50, 100, 150)) 

ggsave("./plot/dumbbellPlot_GSEA_go_bp_F1_v2.tif",
       plot = p, device = "tiff",
       units = "mm", width = 110, height = 90, dpi = 300
)
