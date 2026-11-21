# Extract the two panels of Figure 6 from the NHANES design file
suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu = "adjust")

csv <- "C:/Users/liuha/Desktop/\u5b5f\u5fb7\u5c14/revision_20261002/clinical_nhanes3_20261011/nhanes3_knee_vitD_design.csv"
out <- "C:/Users/liuha/OpenClawVDR/repo/tables/fig6_data"

dir.create(out, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(csv, stringsAsFactors = FALSE)
num <- function(x) suppressWarnings(as.numeric(x))
df$vitD <- num(df$vitD_nmol); df$kll <- num(df$KL_left); df$klr <- num(df$KL_right)
df$age <- num(df$age); df$bmi <- num(df$bmi); df$month <- num(df$exam_month)
df$wt <- num(df$wt); df$psu <- as.integer(df$psu); df$stratum <- as.integer(df$stratum)
df <- df[!is.na(df$vitD) & df$kll %in% 0:4 & df$klr %in% 0:4, ]
df$bmi[df$bmi < 12 | df$bmi > 70] <- NA
df$klmax <- pmax(df$kll, df$klr)
df$roa <- as.integer(df$klmax >= 2)
df$female <- as.integer(df$sex == 2)
df$smoke <- ifelse(df$smoke100 == 1, 1L, ifelse(df$smoke100 == 2, 0L, NA_integer_))
df$season <- ifelse(df$month %in% c(12,1,2), "winter",
             ifelse(df$month %in% c(3,4,5), "spring",
             ifelse(df$month %in% c(6,7,8), "summer", "autumn")))
df$vitD25 <- df$vitD/25
df$race <- factor(df$race); df$season <- factor(df$season)

d <- df[complete.cases(df[, c("roa","vitD25","age","female","race","bmi","season","smoke","wt","psu","stratum")]), ]
cat("N:", nrow(d), "\n")

des <- svydesign(ids = ~psu, strata = ~stratum, weights = ~wt, nest = TRUE, data = d)

## panel A: weighted mean 25(OH)D by KL grade
d$klf <- factor(d$klmax, levels = 0:4)
desA <- svydesign(ids = ~psu, strata = ~stratum, weights = ~wt, nest = TRUE, data = d)
m <- svyby(~vitD, ~klf, desA, svymean, na.rm = TRUE)
n_by <- as.data.frame(table(d$klmax))
names(n_by) <- c("grade", "n")
panA <- data.frame(grade = as.integer(as.character(m$klf)),
                   mean = as.numeric(m$vitD),
                   se = as.numeric(SE(m)),
                   n = n_by$n[match(as.integer(as.character(m$klf)), n_by$grade)])
write.csv(panA, file.path(out, "panelA_mean_by_grade.csv"), row.names = FALSE)
print(panA)

## panel B: progressive adjustment, OR per 25 nmol/L
models <- list(
  list("Crude",                 roa ~ vitD25),
  list("+ age, sex",            roa ~ vitD25 + age + female),
  list("+ race",                roa ~ vitD25 + age + female + race),
  list("+ BMI",                 roa ~ vitD25 + age + female + race + bmi),
  list("+ season",              roa ~ vitD25 + age + female + race + bmi + season),
  list("+ smoking",             roa ~ vitD25 + age + female + race + bmi + season + smoke)
)
rows <- list()
for (mm in models) {
  fit <- svyglm(mm[[2]], design = des, family = quasibinomial())
  b <- unname(coef(fit)["vitD25"]); se <- unname(SE(fit)["vitD25"])
  dfr <- max(1, degf(des)); tq <- qt(0.975, dfr)
  rows[[length(rows) + 1]] <- data.frame(
    model = mm[[1]], OR = exp(b),
    lo = exp(b - tq * se), hi = exp(b + tq * se),
    p = 2 * pt(-abs(b/se), dfr))
}
panB <- do.call(rbind, rows)
panB$model <- factor(panB$model, levels = panB$model)
write.csv(panB, file.path(out, "panelB_progressive_OR.csv"), row.names = FALSE)
print(panB)
cat("\nwritten to", out, "\n")
