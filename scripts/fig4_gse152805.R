## Recompute Figure 4 panel data from the original GSE152805 Seurat object
suppressMessages({library(Seurat); library(data.table)})

out <- "C:/Users/liuha/OpenClawVDR/repo/tables"
rds <- "C:/Users/liuha/Desktop/\u5b5f\u5fb7\u5c14/revision_20261002/outputs/GSE152805_seurat.rds"
merged <- readRDS(rds)
merged <- JoinLayers(merged)
cat("cells:", ncol(merged), " genes:", nrow(merged), "\n")
cat("meta cols:", paste(head(colnames(merged@meta.data), 12), collapse = ", "), "\n")

data <- LayerData(merged, assay = "RNA", layer = "data")
cl <- as.character(merged$seurat_clusters)
donor <- as.character(merged$donor)
cat("clusters:", length(unique(cl)), " donors:", paste(sort(unique(donor)), collapse = ","), "\n")

panel <- list(
  HomC  = c("ACAN","COL2A1","SOX9","COMP","PRG4","TIMP3","COL11A1"),
  HTC   = c("COL10A1","IBSP","MMP13","RUNX2","IHH","ALPL","PHOSPHO1"),
  preHTC= c("COL10A1","COL9A1","SPP1"),
  FC    = c("COL1A1","COL3A1","DCN","FN1","LUM"),
  ProC  = c("MKI67","TOP2A","UBE2C","CENPF","CCNB1"),
  EC    = c("VEGFA","HIF1A","ANGPTL4","SLC2A1"),
  RegC  = c("CCL3","CCL4","CCL4L2","CXCL8","CXCL1","CXCL2","CXCL3","CXCL6","IL6","NFKBIA",
            "TNFAIP3","SOCS3","CEBPB","FOS","FOSB","JUN","JUNB","ATF3","NFKB1","RELA","ICAM1","CXCL5")
)
SIGS <- list(
  regulatory = panel$RegC,
  catabolic  = c("MMP13","COL10A1","ADAMTS5","IBSP","RUNX2","IL6","CXCL8"),
  NLRP3      = c("NLRP3","PYCARD","CASP1","IL1B","IL18","GSDMD","TXNIP","AIM2","NLRC4","CASP8")
)

zc <- function(g) {
  g <- intersect(g, rownames(data))
  if (!length(g)) return(NULL)
  m <- colMeans(data[g, , drop = FALSE])
  # z-score the score itself across cells is not needed; use mean of per-gene z
  Z <- t(scale(t(data[g, , drop = FALSE])))
  Z[is.na(Z)] <- 0
  colMeans(Z)
}

cell <- data.table(cell = colnames(merged), cluster = cl, donor = donor,
                   vdr = as.integer(data["VDR", ] > 0),
                   total = merged$nCount_RNA)
for (k in names(SIGS)) {
  s <- zc(SIGS[[k]])
  cell[[k]] <- if (is.null(s)) NA_real_ else as.numeric(s)
}
cat("signatures available:", paste(names(SIGS)[!is.na(cell[1, names(SIGS), with = FALSE])], collapse = ", "), "\n")

## state assignment by best-scoring panel
sc <- sapply(names(panel), function(k) {
  g <- intersect(panel[[k]], rownames(data))
  if (!length(g)) return(rep(NA_real_, length(unique(cl))))
  Zc <- t(scale(t(data[g, , drop = FALSE]))); Zc[is.na(Zc)] <- 0
  sapply(sort(unique(cl)), function(c) mean(Zc[, cl == c, drop = FALSE]))
})
rownames(sc) <- sort(unique(cl))
state_of <- colnames(sc)[apply(sc, 1, which.max)]
names(state_of) <- rownames(sc)
cell[, cluster_state := state_of[cluster]]
print(table(cell$cluster_state))

fwrite(cell[, .(cell, cluster, cluster_state, donor, vdr, total,
                regulatory, catabolic, NLRP3)],
       file.path(out, "fig4_gse152805_cells.csv"))

## panel A/B: per donor, RegC/HomC vs other
cell[, grp := cluster_state %in% c("RegC", "HomC")]
ab <- cell[, .(n = .N, det = sum(vdr)), by = .(donor, grp)]
ab[, frac := det / n]
print(dcast(ab, donor ~ grp, value.var = c("frac", "n")))
fwrite(ab, file.path(out, "fig4_gse152805_donor_detection.csv"))

## panel C/D: within-state VDR+ vs VDR-
res <- list()
for (st in unique(cell$cluster_state)) {
  if (is.na(st)) next
  sub <- cell[cluster_state == st]
  if (sum(sub$vdr == 1) < 8 || sum(sub$vdr == 0) < 8) next
  for (k in names(SIGS)) {
    a <- sub[vdr == 1][[k]]; b <- sub[vdr == 0][[k]]
    if (all(is.na(a)) || all(is.na(b))) next
    tt <- t.test(a, b)
    res[[length(res) + 1]] <- data.table(state = st, signature = k,
                                         diff = mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE),
                                         p = tt$p.value, n_pos = length(a), n_neg = length(b))
  }
}
eff <- rbindlist(res)
print(dcast(eff, state ~ signature, value.var = "diff"))
fwrite(eff, file.path(out, "fig4_gse152805_within_state.csv"))
cat("\nDONE_FIG4\n")
