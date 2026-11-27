## Colocalisation: VDR cis-eQTL vs a FinnGen endpoint in chr12:46.8-48.9 Mb
suppressMessages(library(coloc))
args <- commandArgs(trailingOnly = TRUE)
IN   <- args[1]
NCASE<- as.numeric(args[2])
NCTRL<- as.numeric(args[3])
TAG  <- args[4]

d <- read.csv(IN, stringsAsFactors = FALSE)
cat(sprintf("%s: %d variants\n", TAG, nrow(d)))

N_EQTL <- 31684
N_OA   <- NCASE + NCTRL
S_OA   <- NCASE / N_OA
cat(sprintf("N_case=%.0f N_ctrl=%.0f s=%.4f\n", NCASE, NCTRL, S_OA))

d1 <- list(snp = d$snp, beta = d$b1, varbeta = d$se1^2, N = N_EQTL, type = "quant",
           sdY = 1, MAF = d$maf)
d2 <- list(snp = d$snp, beta = d$b2, varbeta = d$se2^2, N = N_OA, type = "cc",
           s = S_OA, MAF = d$maf)

r <- suppressMessages(coloc.abf(d1, d2, p12 = 1e-5))$summary
cat(sprintf("n=%d  H0=%.2g H1=%.2g H2=%.2g H3=%.4f H4=%.3e\n",
            r[["nsnps"]], r[["PP.H0.abf"]], r[["PP.H1.abf"]], r[["PP.H2.abf"]],
            r[["PP.H3.abf"]], r[["PP.H4.abf"]]))

out <- data.frame(endpoint = TAG, nsnps = r[["nsnps"]], H0 = r[["PP.H0.abf"]],
                  H1 = r[["PP.H1.abf"]], H2 = r[["PP.H2.abf"]],
                  H3 = r[["PP.H3.abf"]], H4 = r[["PP.H4.abf"]])
write.csv(out, file.path("C:/Users/liuha/OpenClawVDR/repo/tables",
                         sprintf("coloc_finngen_%s.csv", TAG)), row.names = FALSE)
cat("written\n")
