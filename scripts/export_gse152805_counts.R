## Export GSE152805 counts and metadata for pySCENIC
suppressMessages({library(Seurat); library(Matrix)})

rds <- "C:/Users/liuha/Desktop/\u5b5f\u5fb7\u5c14/revision_20261002/outputs/GSE152805_seurat.rds"
out <- "C:/Users/liuha/OpenClawVDR/tools/scenic"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

s <- readRDS(rds)
s <- JoinLayers(s)
cat("cells:", ncol(s), "genes:", nrow(s), "\n")

cnt <- LayerData(s, assay = "RNA", layer = "counts")
cnt <- as(cnt, "CsparseMatrix")
writeMM(cnt, file.path(out, "counts.mtx"))
writeLines(rownames(cnt), file.path(out, "genes.txt"))
writeLines(colnames(cnt), file.path(out, "cells.txt"))

meta <- data.frame(cell = colnames(s),
                   cluster = as.character(s$seurat_clusters),
                   donor = as.character(s$donor),
                   stringsAsFactors = FALSE)
write.csv(meta, file.path(out, "cellmeta.csv"), row.names = FALSE)
cat("written counts.mtx, genes.txt, cells.txt, cellmeta.csv\n")
cat("DONE_EXPORT\n")
