# Introduction ------------------------------------------------------------

# Image: R_Lib Common
# Output: Annotation of ENSEMBL gene IDs in DESeq2 Wald test results
# Data: DESeq2 Wald test tables under normoxia

# Load library ------------------------------------------------------------

library(AnnotationDbi)
library(org.Hs.eg.db)

library(data.table)
library(tidyr)
library(purrr)

# Print keytypes in human database
keytypes(org.Hs.eg.db)

# Load DESeq2 results -----------------------------------------------------

resultFC <- fread("./results/deseq_wald_test_normoxia.csv")
View(resultFC)

head(resultFC)

# ID mapping --------------------------------------------------------------

# Mapping to Entrez ID
resultFC$entrez <- mapIds(org.Hs.eg.db,
                          keys = resultFC$ensembl_gene_id,
                          keytype = "ENSEMBL", column = "ENTREZID",
                          multiVals = "first"
)

print(paste0(sum(!is.na(resultFC$entrez)), " genes mapped")) # 12579 genes mapped

# Mapping to gene symbol
resultFC$symbol <- mapIds(org.Hs.eg.db,
                          keys = resultFC$ensembl_gene_id,
                          keytype = "ENSEMBL", column = "SYMBOL",
                          multiVals = "first"
)

print(paste0(sum(!is.na(resultFC$symbol)), " genes mapped")) # 12579 genes mapped

# Count how many genes that are not mapped but have significant FDR
padj_unmap <- resultFC[is.na(resultFC$entrez), ]$padj %>%
  discard(is.na)

print(paste0(length(padj_unmap[padj_unmap < 0.05]), " unmapped genes are significant")) # 377 for normoxia, 405 for hypoxia

# Save result table
fwrite(resultFC, file = "./results/Annotated_deseq_wald_test_normoxia.csv")

