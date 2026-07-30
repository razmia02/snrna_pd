library(Seurat)

# Quick smoke test for snRNA-seq workflow
counts <- matrix(rpois(5000, lambda = 2), nrow = 100, ncol = 50)
rownames(counts) <- paste0("Gene", 1:100)
colnames(counts) <- paste0("Cell", 1:50)

seurat_obj <- CreateSeuratObject(counts = counts, project = "snRNA_test")
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)

print("CI Check Passed: Basic snRNA-seq workflow executed successfully!")
