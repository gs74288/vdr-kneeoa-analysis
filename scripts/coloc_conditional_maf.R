## Recompute conditional colocalisation with the same specification as the windowed analysis
suppressMessages(library(coloc))

H  <- "C:/Users/liuha/OpenClawVDR/repo/tables/coloc_harmonised_input_GCST90566800.csv"
C  <- "C:/Users/liuha/OpenClawVDR/tools/cojo/coloc_conditional_input.csv"
OUT<- "C:/Users/liuha/OpenClawVDR/repo/tables"

h <- read.csv(H, stringsAsFactors = FALSE)
c <- read.csv(C, stringsAsFactors = FALSE)
cat("harmonised:", nrow(h), "conditional:", nrow(c), "\n")
cat("harmonised cols:", paste(colnames(h), collapse=","), "\n")

m <- merge(c[, c("snp","b2_orig","se2_orig","b2_cond","se2_cond")],
           h[, c("snp","b1","se1","maf","ncase","ncontrol")], by="snp")
cat("merged:", nrow(m), "| maf available:", sum(!is.na(m$maf)), "\n")

N_EQTL <- 31684
N_CASE <- median(m$ncase, na.rm = TRUE)
N_CTRL <- median(m$ncontrol, na.rm = TRUE)
N_OA   <- N_CASE + N_CTRL
S_OA   <- N_CASE / N_OA
cat(sprintf("N_CASE=%.0f N_CTRL=%.0f N_OA=%.0f s=%.4f\n", N_CASE, N_CTRL, N_OA, S_OA))

d1 <- list(snp=m$snp, beta=m$b1, varbeta=m$se1^2, N=N_EQTL, type="quant", sdY=1, MAF=m$maf)

run <- function(beta, varbeta, tag, use_maf=TRUE) {
  d2 <- list(snp=m$snp, beta=beta, varbeta=varbeta, N=N_OA, type="cc", s=S_OA)
  if (use_maf) d2$MAF <- m$maf
  r <- suppressMessages(coloc.abf(d1, d2, p12=1e-5))$summary
  cat(sprintf("%-28s n=%d  H0=%.2g H1=%.2g H2=%.2g H3=%.4f H4=%.3e\n",
              tag, r[["nsnps"]], r[["PP.H0.abf"]], r[["PP.H1.abf"]], r[["PP.H2.abf"]],
              r[["PP.H3.abf"]], r[["PP.H4.abf"]]))
  data.frame(window=tag, nsnps=r[["nsnps"]], H0=r[["PP.H0.abf"]], H1=r[["PP.H1.abf"]],
             H2=r[["PP.H2.abf"]], PP.H3=r[["PP.H3.abf"]], PP.H4=r[["PP.H4.abf"]])
}

rows <- list()
rows[[1]] <- run(m$b2_orig, m$se2_orig^2, "disease-peak 100kb, uncond", TRUE)
rows[[2]] <- run(m$b2_cond, m$se2_cond^2, "disease-peak 100kb, cond", TRUE)
rows[[3]] <- run(m$b2_orig, m$se2_orig^2, "uncond, no MAF", FALSE)
rows[[4]] <- run(m$b2_cond, m$se2_cond^2, "cond, no MAF", FALSE)

out <- do.call(rbind, rows)
write.csv(out, file.path(OUT, "coloc_conditional_recomputed.csv"), row.names=FALSE)
cat("\nwritten coloc_conditional_recomputed.csv\n")
cat(sprintf("sum check: %s\n", paste(sprintf("%.6f", rowSums(out[,c("H0","H1","H2","PP.H3","PP.H4")])), collapse=" ")))
