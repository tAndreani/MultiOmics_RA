# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: Annotation of ENSEMBL gene IDs in DESeq2 Wald test results
# Data: DESeq2 Wald test table

# Load library ------------------------------------------------------------

library(AnnotationDbi)
library(org.Hs.eg.db)

library(data.table)
library(tidyr)
library(purrr)

# Print keytypes in human database
keytypes(org.Hs.eg.db)

# Loop to annotate DESeq2 results -----------------------------------------

csvPath <- list.files("./results", full.names = TRUE, pattern = "deseq_wald_test")
csvPath

for (i in csvPath) {

  # Load data
  resultFC <- fread(i)

  # Extract condition
  tp <- unique(resultFC$condition)
  print(tp)

  # Tidy up Ensembl IDs
  resultFC <- resultFC %>%
    separate(ensembl_gene_id, c("gene_id", "version"))

  # Mapping to Entrez ID
  resultFC$entrez <- mapIds(org.Hs.eg.db,
    keys = resultFC$gene_id,
    keytype = "ENSEMBL", column = "ENTREZID",
    multiVals = "first"
  )

  print(paste0(sum(!is.na(resultFC$entrez)), " genes mapped")) # 18075 genes mapped

  # Mapping to gene symbol
  resultFC$symbol <- mapIds(org.Hs.eg.db,
    keys = resultFC$gene_id,
    keytype = "ENSEMBL", column = "SYMBOL",
    multiVals = "first"
  )

  print(paste0(sum(!is.na(resultFC$symbol)), " genes mapped")) # 18075 genes mapped

  # Count how many genes that are not mapped but have significant FDR
  padj_unmap <- resultFC[is.na(resultFC$entrez), ]$padj %>%
    discard(is.na)

  print(paste0(length(padj_unmap[padj_unmap < 0.05]), " unmapped genes are significant at ", tp))

  # Save result table
  fwrite(resultFC, file = paste0("./results/deseq_wald_test_annotated_", tp, ".csv"))

  rm(resultFC)
}
