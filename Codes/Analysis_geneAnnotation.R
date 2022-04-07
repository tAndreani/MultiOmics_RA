# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: Annotation of ENSEMBL gene IDs in DESeq2 Wald test results
# Data: DESeq2 Wald test table

# Load library ------------------------------------------------------------

library(AnnotationDbi)
library(org.Mm.eg.db)

library(data.table)
library(tidyr)
library(purrr)

# Print keytypes in mouse database
keytypes(org.Mm.eg.db)

# Use AnnotationDbi to map to symbol and Entrez ID ------------------------

# 2 methods available: mapIds can only return 1 column, whereas AnnotationDbi::select can return multiple columns

## Trial with a few genes
mapIds(org.Mm.eg.db, keys = resultFC$ensembl_gene_id[1:5], keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
AnnotationDbi::select(org.Mm.eg.db, keys = resultFC$ensembl_gene_id[1:5], keytype = "ENSEMBL", columns = c("ENTREZID", "SYMBOL"))

# Loop to annotate DESeq2 results -----------------------------------------

csvPath <- list.files("./result/DESeq2_Wald_Test_Result/", full.names = TRUE)
csvPath

for (i in csvPath) {

  # Load data
  resultFC <- fread(i)

  # Extract time point
  tp <- unique(resultFC$Timepoint)
  print(tp)

  # Tidy up Ensembl IDs
  resultFC <- resultFC %>%
    separate(gene_id, c("ensembl_gene_id", "version"))

  # Mapping to Entrez ID
  resultFC$entrez <- mapIds(org.Mm.eg.db,
    keys = resultFC$ensembl_gene_id,
    keytype = "ENSEMBL", column = "ENTREZID",
    multiVals = "first"
  )

  print(paste0(sum(!is.na(resultFC$entrez)), " genes mapped")) # 15290 genes mapped

  # Mapping to gene symbol
  resultFC$symbol <- mapIds(org.Mm.eg.db,
    keys = resultFC$ensembl_gene_id,
    keytype = "ENSEMBL", column = "SYMBOL",
    multiVals = "first"
  )

  print(paste0(sum(!is.na(resultFC$symbol)), " genes mapped")) # 15290 genes mapped

  # Count how many genes that are not mapped but have significant FDR
  padj_unmap <- resultFC[is.na(resultFC$entrez), ]$padj %>%
    discard(is.na)

  print(paste0(length(padj_unmap[padj_unmap < 0.05]), " unmapped genes are significant at ", tp))

  # Save result table
  fwrite(resultFC, file = paste0("./result/Annotated_DESeq2_Wald_Test_Result/DESeq2_annotated_noFilter_", tp, ".csv"))

  rm(resultFC)
}

# 1 unmapped genes are significant at 2_days
# 4 unmapped genes are significant at 2_weeks
# 64 unmapped genes are significant at 3_weeks
# 17 unmapped genes are significant at 4_weeks
# 12 unmapped genes are significant at 6_weeks
# 10 unmapped genes are significant at 7_weeks
# 3 unmapped genes are significant at 8_weeks
# 11 unmapped genes are significant at 10_weeks

# Optional: Using biomaRt to annotate genes -------------------------------

library(biomaRt)

# Step1: Identifying the database you need
listEnsembl()

ensembl <- useEnsembl(biomart = "genes")

# Step 2: Choosing a dataset
datasets <- listDatasets(ensembl)
View(datasets)

ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

# Mapping genes
ensemblID <- resultFC$ensembl_gene_id

trial <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  filters = "ensembl_gene_id",
  values = ensemblID,
  mart = ensembl
)
View(trial)

sum(!is.na(trial$entrezgene_id)) # 21945 genes mapped, which is less than AnnotationDbi
