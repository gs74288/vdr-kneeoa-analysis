## Colocalisation before and after conditioning the OA signal on rs11168351
suppressMessages(library(coloc))

IN  <- "C:/Users/liuha/OpenClawVDR/tools/cojo/coloc_conditional_input.csv"
OUT <- "C:/Users/liuha/OpenClawVDR/repo/tables"
d <- read.csv(IN, stringsAsFactors = FALSE)
cat("variants:", nrow(d), "\n")

N_EQTL <- 31684
N_OA <- 1182675
S_OA <- 0.1382

d1 <- list(snp = d$snp, beta = d$b1, varbeta = d$se1^2, N = N_EQTL, type = "quant", sdY = 1)

run <- function(beta, varbeta, tag) {
  d2 <- list(snp = d$snp, beta = beta, varbeta = varbeta, N = N_OA, type = "cc", s = S_OA)
  s <- suppressMessages(coloc.abf(d1, d2, p12 = 1e-5))$summary
  cat(sprintf("%-12s n=%d  H1=%.3g  H3=%.4f  H4=%.3e\n",
              tag, s[["nsnps"]], s[["PP.H1.abf"]], s[["PP.H3.abf"]], s[["PP.H4.abf"]]))
  s
}

a <- run(d$b2_orig, d$se2_orig^2, "unconditional")
b <- run(d$b2_cond, d$se2_cond^2, "conditional")

out <- rbind(as.data.frame(t(a)), as.data.frame(t(b)))
out$model <- c("unconditional", "conditional")
write.csv(out, file.path(OUT, "coloc_conditional_summary.csv"), row.names = FALSE)
cat("\nwritten coloc_conditional_summary.csv\n")
