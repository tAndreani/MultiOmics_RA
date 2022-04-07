# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: GSEA result on GO Biological Process pathways
# Data: Annotated DESeq2 results under normoxia and hypoxia

# Load library ------------------------------------------------------------

library(gage)
library(fgsea)

library(data.table)
library(dplyr)
library(stringr)

# Generate up-to-date GO (Gene Ontology) gene sets ------------------------

## https://rdrr.io/bioc/gage/man/go.gsets.html
go.hs <- go.gsets(species = "Human")

# Subset only biological process
go.bp <- go.hs$go.sets[go.hs$go.subs$BP]
View(go.bp)

length(go.bp) # 16104 GO biological process pathways

# GSEA on GO biological process pathway -----------------------------------

# List annotated DESeq2 result csv file path
csvPath <- list.files("./results", full.names = TRUE, pattern = "Annotated")
csvPath

# Loop to perform GSEA in different time points
for (singleFile in csvPath) {
  
  # Load csv file
  input <- fread(singleFile)
  
  # Define oxidative condition
  tp <- str_split_fixed(singleFile, "_", 5)[,5]
  print(tp)
  
  # Remove rows with unmapped Entrez IDs
  input_dropNA <- input[!is.na(input$entrez), ]
  
  # Remove Entrez IDs that are mapped to multiple genes
  res_DESeq2 <- input_dropNA %>%
    group_by(entrez) %>%
    filter(n() == 1)
  
  print(paste0(nrow(res_DESeq2), " unique genes as GSEA input"))
  
  # Prepare geneList by sorting fold change in decreasing order
  deseq2.fc <- res_DESeq2$log2FoldChange
  names(deseq2.fc) <- res_DESeq2$entrez
  deseq2.fc <- sort(deseq2.fc, decreasing = TRUE)
  print(head(deseq2.fc))
  
  # fgsea on GO biological process pathway
  fgsea_go <- fgseaMultilevel(pathways = go.bp, stats = deseq2.fc, minSize = 15, maxSize = 600, eps = 0)
  
  print(paste0("Number of significant GO pathways (padj < 0.05): ", sum(fgsea_go$padj < 0.05, na.rm = TRUE), " at ", tp))
  
  # Save result
  df_go <- fgsea_go %>%
    as.data.frame() %>%
    mutate(Oxidative_Condition = rep(tp, nrow(fgsea_go)))
  
  fwrite(df_go, file = paste0("./results/fgsea_go_bp_", tp))
  
  rm(deseq2.fc)
}

#"Number of significant GO pathways (padj < 0.05): 777 at hypoxia.csv"
#"Number of significant GO pathways (padj < 0.05): 888 at normoxia.csv"
