## NHANES III: sample flow and weighted baseline by radiographic OA status
suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu = "adjust")

csv <- "C:/Users/liuha/Desktop/\u5b5f\u5fb7\u5c14/revision_20261002/clinical_nhanes3_20261011/nhanes3_knee_vitD_design.csv"
df <- read.csv(csv, stringsAsFactors = FALSE)
num <- function(x) suppressWarnings(as.numeric(x))
n0 <- nrow(df)

df$vitD <- num(df$vitD_nmol); df$kll <- num(df$KL_left); df$klr <- num(df$KL_right)
df$age <- num(df$age); df$bmi <- num(df$bmi); df$month <- num(df$exam_month)
df$wt <- num(df$wt); df$psu <- as.integer(df$psu); df$stratum <- as.integer(df$stratum)

has_rad <- !is.na(df$kll) & df$kll %in% 0:4 & !is.na(df$klr) & df$klr %in% 0:4
has_vit <- !is.na(df$vitD)
df1 <- df[has_rad & has_vit, ]
df1$bmi[df1$bmi < 12 | df1$bmi > 70] <- NA
df1$klmax <- pmax(df1$kll, df1$klr)
df1$roa <- as.integer(df1$klmax >= 2)
df1$female <- as.integer(df1$sex == 2)
df1$smoke <- ifelse(df1$smoke100 == 1, 1L, ifelse(df1$smoke100 == 2, 0L, NA_integer_))
df1$season <- ifelse(df1$month %in% c(12,1,2), "winter",
              ifelse(df1$month %in% c(3,4,5), "spring",
              ifelse(df1$month %in% c(6,7,8), "summer", "autumn")))
df1$race <- factor(df1$race); df1$season <- factor(df1$season)
cc <- complete.cases(df1[, c("roa","age","female","race","bmi","season","smoke","wt","psu","stratum")])

cat(sprintf("records in file            : %d\n", n0))
cat(sprintf("with gradable radiograph    : %d\n", sum(has_rad)))
cat(sprintf("and measured 25(OH)D       : %d\n", nrow(df1)))
cat(sprintf("excluded: missing covariates: %d\n", sum(!cc)))
cat(sprintf("analysed                   : %d\n", sum(cc)))

d <- df1[cc, ]
des <- svydesign(ids = ~psu, strata = ~stratum, weights = ~wt, nest = TRUE, data = d)

cat("\n=== weighted baseline by radiographic OA ===\n")
mk <- function(var, label) {
  f <- as.formula(paste("~", var))
  a <- svyby(f, ~roa, des, svymean, na.rm = TRUE)
  cat(sprintf("%-22s OA: %s | no OA: %s\n", label,
              paste(sprintf("%.1f", a[2, -1]), collapse = " "),
              paste(sprintf("%.1f", a[1, -1]), collapse = " ")))
}
mk("age", "age (years)")
mk("bmi", "BMI (kg/m2)")
a <- svyby(~vitD, ~roa, des, svymean, na.rm = TRUE)
cat(sprintf("%-22s OA: %.1f | no OA: %.1f\n", "25(OH)D (nmol/L)", a[2,2], a[1,2]))
p <- svyby(~female, ~roa, des, svymean, na.rm = TRUE)
cat(sprintf("%-22s OA: %.1f%% | no OA: %.1f%%\n", "female (%)", 100*p[2,2], 100*p[1,2]))
