# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: GSEA result on GO Biological Process pathways
# Data: Annotated DESeq2 result

# Load library ------------------------------------------------------------

library(gage)
library(fgsea)

library(data.table)
library(dplyr)
library(stringr)

# Generate up-to-date GO (Gene Ontology) gene sets ------------------------

## https://rdrr.io/bioc/gage/man/go.gsets.html
go.mmu <- go.gsets(species = "Mouse")

# Subset only biological process
go.bp <- go.mmu$go.sets[go.mmu$go.subs$BP]
View(go.bp)

length(go.bp) # 16025 GO biological process pathways

# GSEA on GO biological process pathway -----------------------------------

# List annotated DESeq2 result csv file path
csvPath <- list.files("./result/Annotated_DESeq2_Wald_Test_Result", full.names = TRUE)
csvPath

# Loop to perform GSEA in different time points
for (singleFile in csvPath) {
  
  # Load csv file
  input <- fread(singleFile)
  
  # Define time point
  tp <- unique(input$Timepoint)
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
  fgsea_go <- fgseaMultilevel(pathways = go.bp, stats = deseq2.fc, minSize = 15, maxSize = 600)
  
  print(paste0("Number of significant GO pathways (padj < 0.05): ", sum(fgsea_go$padj < 0.05, na.rm = TRUE), " at ", tp))
  
  # Save result
  df_go <- fgsea_go %>%
    as.data.frame() %>%
    mutate(Timepoint = rep(tp, nrow(fgsea_go)))
  
  fwrite(df_go, file = paste0("./result/fgseaResults/fgsea_go_bp_", tp, ".csv"))
  
  rm(deseq2.fc)
}

# [1] "Number of significant GO pathways (padj < 0.05): 345 at 2_days"
# [1] "Number of significant GO pathways (padj < 0.05): 292 at 2_weeks"
# [1] "Number of significant GO pathways (padj < 0.05): 837 at 3_weeks"
# [1] "Number of significant GO pathways (padj < 0.05): 916 at 4_weeks"
# [1] "Number of significant GO pathways (padj < 0.05): 578 at 6_weeks"
# [1] "Number of significant GO pathways (padj < 0.05): 680 at 7_weeks"
# [1] "Number of significant GO pathways (padj < 0.05): 526 at 8_weeks"
# [1] "Number of significant GO pathways (padj < 0.05): 951 at 10_weeks"
