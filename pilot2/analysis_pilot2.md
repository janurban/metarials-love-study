Pilot Study 2 (Polish replication of Tang et al., 2025): Reproducible
Analysis
================

This document reproduces every analysis reported in *Pilot Study 2:
Materials and Scale Validation*, which is summarised in the manuscript
(Study 2, Method: “Pilot Study 2, N = 129”). It re-runs the original
SPSS syntax (`Analizy.spv`) and the JASP confirmatory factor analyses in
R. Knit this file in RStudio (**Knit** button) or with
`rmarkdown::render("analysis_pilot2.Rmd")` to regenerate
`analysis_pilot2.md`.

**Data source.** By default, the document downloads
`pilot2/data_pilot2.sav` from the GitHub repository
[janurban/metarials-love-study](https://github.com/janurban/metarials-love-study).
It can also read a local `.sav` file by passing a path with
`rmarkdown::render("analysis_pilot2.Rmd", params = list(data = "path/to/file.sav"))`.

**Required packages.** `haven`, `dplyr`, `knitr`, `lavaan`. Install them
with
`install.packages(c("haven", "dplyr", "knitr", "lavaan", "rmarkdown"))`.

# 0. Setup

``` r
library(haven)   # read the SPSS .sav file
library(dplyr)   # data handling
library(knitr)   # tables
library(lavaan)  # CFA and McDonald's omega
```

Helper functions. Standardised mean differences are Cohen’s *d*
(one-sample: (*M* – 4) / *SD*; paired: *d<sub>z</sub>* =
*M*<sub>diff</sub> / *SD*<sub>diff</sub>), with 95% CIs from the
noncentral *t* distribution. These are the same values SPSS prints with
`/ES DISPLAY(TRUE)`.

``` r
# number formatting
fmt  <- function(x, d = 2) formatC(x, format = "f", digits = d)
no0  <- function(x, d = 3) sub("^(-?)0\\.", "\\1.", fmt(x, d))   # APA style: drop leading zero
fmtp <- function(p) ifelse(p < .001, "< .001", paste("=", no0(p, 3)))
pct  <- function(n, N) fmt(100 * n / N, 1)

# 95% CI for Cohen's d via the noncentral t distribution (as in SPSS)
d_ci <- function(t, df, n, level = .95) {
  a  <- (1 - level) / 2
  lo <- uniroot(function(ncp) pt(t, df, ncp, lower.tail = FALSE) - a, c(t - 10, t + 10))$root
  hi <- uniroot(function(ncp) pt(t, df, ncp) - a, c(t - 10, t + 10))$root
  c(lo, hi) / sqrt(n)
}

# One-sample t test against the scale midpoint
one_sample <- function(x, mu = 4) {
  x  <- x[!is.na(x)]
  n  <- length(x)
  tt <- t.test(x, mu = mu)
  d  <- (mean(x) - mu) / sd(x)
  ci <- d_ci(unname(tt$statistic), n - 1, n)
  data.frame(M = mean(x), SD = sd(x), t = unname(tt$statistic), df = n - 1,
             p = tt$p.value, d = d, d_low = ci[1], d_high = ci[2])
}

# Cronbach's alpha
cronbach_alpha <- function(X) {
  X <- as.matrix(X); k <- ncol(X)
  k / (k - 1) * (1 - sum(apply(X, 2, var)) / var(rowSums(X)))
}

# McDonald's omega from a one-factor ML model:
# omega = (sum of loadings)^2 / [(sum of loadings)^2 + sum of residual variances].
# This is the estimator used by SPSS RELIABILITY /MODEL=OMEGA.
mcdonald_omega <- function(data, items) {
  fit <- cfa(paste("f =~", paste(items, collapse = " + ")), data = data, std.lv = TRUE)
  pe  <- parameterEstimates(fit)
  lam <- pe$est[pe$op == "=~"]
  th  <- pe$est[pe$op == "~~" & pe$lhs == pe$rhs & pe$lhs %in% items]
  sum(lam)^2 / (sum(lam)^2 + sum(th))
}

txt <- list()  # formatted numbers used in the text
```

# 1. Data

``` r
data_source <- params$data
if (grepl("^https?://", data_source, ignore.case = TRUE)) {
  # A github.com/.../blob/... link points to the web page, not the file:
  # rewrite it to the raw-file address
  data_source <- sub("^https?://github\\.com/([^/]+)/([^/]+)/blob/",
                     "https://raw.githubusercontent.com/\\1/\\2/", data_source)
  data_file <- tempfile(fileext = ".sav")
  download_status <- download.file(data_source, data_file, mode = "wb", quiet = TRUE)
  if (download_status != 0 || !file.exists(data_file) || file.size(data_file) == 0) {
    stop("Could not download the data file from ", data_source)
  }
  # A valid SPSS system file starts with "$FL2" (or "$FL3" if compressed with zlib)
  if (!readChar(data_file, 4, useBytes = TRUE) %in% c("$FL2", "$FL3")) {
    stop("The file downloaded from ", data_source, " is not an SPSS .sav file. ",
         "Check that the URL points to the raw file and that the repository is public.")
  }
} else {
  data_file <- data_source
}

dat_raw <- read_sav(data_file)
dim(dat_raw)
```

    ## [1] 129  50

``` r
# Labelled categorical variables -> factors; all other variables -> plain numeric/character
cat_vars <- c("Eco_products1", "Eco_products2", "Gender", "Language", "Education",
              paste0("ComprehensionCheck", 1:8))
dat <- dat_raw %>%
  mutate(
    across(all_of(cat_vars), ~ as_factor(.x, levels = "labels"), .names = "{.col}_f"),
    Time_overall = as.numeric(Time_overall),
    Time_vignette = as.numeric(Time_vignette)
  ) %>%
  zap_labels()
```

The shared data file already contains only the 129 native Polish
speakers (the four participants with another native language were
removed before the file was saved), so the exclusion step below changes
nothing. It is kept to document the criterion.

``` r
dat <- dat %>% filter(Language == 1)
N <- nrow(dat)
N
```

    ## [1] 129

Derived variables in the data file are checked against their
definitions:

``` r
key <- c(2, 1, 2, 2, 2, 1, 1, 1)  # correct answer for checks 1-8 (1 = True, 2 = False)
score_ok <- sapply(1:8, function(i)
  all(as.numeric(dat[[paste0("ComprehensionCheck", i)]] == key[i]) ==
        dat[[paste0("ComprehensionCheck", i, "_score")]]))
c(item_scores       = all(score_ok),
  sum_score         = all(rowSums(dat[paste0("ComprehensionCheck", 1:8, "_score")]) ==
                            dat$Overall_comprehension_score),
  virtue_index      = isTRUE(all.equal(rowMeans(dat[c("Virtue1", "Virtue2", "Virtue3")]),
                                       dat$Virtue_index, check.attributes = FALSE)),
  vice_index        = isTRUE(all.equal(rowMeans(dat[c("Vice1", "Vice2", "Vice3")]),
                                       dat$Vice_index, check.attributes = FALSE)),
  intention3_recode = all(8 - dat$Intention3 == dat$Intention3_recoded),
  no_missing        = !anyNA(select(dat, -Open_ended_question, -Study_discipline)))
```

    ##       item_scores         sum_score      virtue_index        vice_index intention3_recode 
    ##              TRUE              TRUE              TRUE              TRUE              TRUE 
    ##        no_missing 
    ##              TRUE

# 2. Participants

``` r
txt$age_m  <- fmt(mean(dat$Age)); txt$age_sd <- fmt(sd(dat$Age))
txt$age_min <- min(dat$Age);      txt$age_max <- max(dat$Age)

tab_freq <- function(f) {
  n <- table(f)
  data.frame(Category = names(n), n = as.vector(n), Percent = pct(as.vector(n), N))
}
```

**Gender**

``` r
kable(tab_freq(dat$Gender_f))
```

| Category                     |   n | Percent |
|:-----------------------------|----:|:--------|
| Female                       |  77 | 59.7    |
| Male                         |  48 | 37.2    |
| Non-binary or gender-diverse |   3 | 2.3     |
| Prefer not to say            |   1 | 0.8     |

**Highest completed education**

``` r
kable(tab_freq(dat$Education_f))
```

| Category                               |   n | Percent |
|:---------------------------------------|----:|:--------|
| Primary                                |   5 | 3.9     |
| Secondary without high school diploma  |   2 | 1.6     |
| Secondary with high school diploma     |  44 | 34.1    |
| Incomplete higher (currently studying) |  55 | 42.6    |
| Higher                                 |  23 | 17.8    |

Age: *M* = 23.98, *SD* = 8.70, range 18–56.

Field of study among participants currently studying (incomplete higher
education):

``` r
students <- dat %>% filter(Education == 4)
disc <- trimws(students$Study_discipline)
disc[disc == ""] <- "(not reported)"
disc_tab <- as.data.frame(table(disc), responseName = "n", stringsAsFactors = FALSE)
disc_tab <- disc_tab[order(-disc_tab$n, disc_tab$disc, method = "radix"), ]
rownames(disc_tab) <- NULL
kable(disc_tab, col.names = c("Field of study", "n"))
```

| Field of study                                               |   n |
|:-------------------------------------------------------------|----:|
| Psychologia                                                  |  39 |
| Ekonomia                                                     |   2 |
| Informatyka                                                  |   2 |
| Pedagogika                                                   |   2 |
| (not reported)                                               |   1 |
| Administracja publiczna w systemie bezpieczeństwa narodowego |   1 |
| Antropologia społeczna                                       |   1 |
| Bezpieczeństwo wewnętrzne                                    |   1 |
| Finanse i rachunkowość                                       |   1 |
| Lekarski                                                     |   1 |
| Marketing i komunikacja rynkowa                              |   1 |
| Politologia                                                  |   1 |
| Prawo                                                        |   1 |
| Produkcja filmowa                                            |   1 |

``` r
txt$n_students <- nrow(students)
txt$n_psych    <- sum(disc == "Psychologia")
txt$pct_psych  <- pct(sum(disc == "Psychologia"), nrow(students))
```

39 of 55 students (70.9%) studied psychology.

# 3. Completion times

`Time_overall` and `Time_vignette` are stored in seconds.

``` r
mmss <- function(s) sprintf("%d:%02d", floor(s / 60), round(s %% 60))
times <- data.frame(
  Variable = c("Whole survey", "Vignette viewing"),
  M_s   = c(mean(dat$Time_overall),   mean(dat$Time_vignette)),
  Mdn_s = c(median(dat$Time_overall), median(dat$Time_vignette)),
  SD_s  = c(sd(dat$Time_overall),     sd(dat$Time_vignette)))
times_out <- times %>% mutate(
  `M (s)` = fmt(M_s), `Mdn (s)` = fmt(Mdn_s), `SD (s)` = fmt(SD_s),
  `M (min)` = fmt(M_s / 60), `Mdn (min)` = fmt(Mdn_s / 60), `SD (min)` = fmt(SD_s / 60),
  `M (min:s)` = mmss(M_s), `Mdn (min:s)` = mmss(Mdn_s), `SD (min:s)` = mmss(SD_s)) %>%
  select(-M_s, -Mdn_s, -SD_s)
kable(times_out)
```

| Variable | M (s) | Mdn (s) | SD (s) | M (min) | Mdn (min) | SD (min) | M (min:s) | Mdn (min:s) | SD (min:s) |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| Whole survey | 416.35 | 239.00 | 679.86 | 6.94 | 3.98 | 11.33 | 6:56 | 3:59 | 11:20 |
| Vignette viewing | 37.22 | 11.22 | 160.56 | 0.62 | 0.19 | 2.68 | 0:37 | 0:11 | 2:41 |

# 4. Familiarity with ecological products

**“Have you ever seen ecological products marked with the previously
shown green-leaf logo in a store?”**

``` r
kable(tab_freq(dat$Eco_products1_f))
```

| Category       |   n | Percent |
|:---------------|----:|:--------|
| Yes            | 105 | 81.4    |
| No             |   6 | 4.7     |
| Don’t remember |  18 | 14.0    |

**“Have you ever purchased ecological products?”**

``` r
kable(tab_freq(dat$Eco_products2_f))
```

| Category       |   n | Percent |
|:---------------|----:|:--------|
| Yes            |  89 | 69.0    |
| No             |   5 | 3.9     |
| Don’t remember |  35 | 27.1    |

# 5. Comprehension check (Table 1)

``` r
cc <- dat$Overall_comprehension_score
txt$cc_m <- fmt(mean(cc)); txt$cc_sd <- fmt(sd(cc))
txt$cc_pct <- pct(mean(cc), 8)
txt$cc_all <- sum(cc == 8); txt$cc_all_pct <- pct(sum(cc == 8), N)
```

Distribution of the comprehension sum score (0–8):

``` r
kable(as.data.frame(table(Score = cc), responseName = "n"))
```

| Score |   n |
|:------|----:|
| 2     |   1 |
| 3     |   1 |
| 5     |   3 |
| 6     |  12 |
| 7     |  27 |
| 8     |  85 |

**Table 1.** Comprehension (true/false) check: response accuracy (*N* =
129).

``` r
items_cc <- c(
  "The production of ecological products involves the intensive use of pesticides and artificial fertilisers in order to ensure high yields. (False)",
  "Ecological products must meet defined technical standards and undergo comprehensive quality control. (True)",
  "Ecological products differ from conventional foods only in price and packaging appearance. (False)",
  "The use of the ecological product logo is not regulated by any regulations. (False)",
  "Ecological products may compromise safety standards in favour of environmental benefits. (False)",
  "The ecological product logo is also meant to remind people of the importance of protecting the environment and preventing pollution. (True)",
  "The logo of ecological products is intended to strengthen harmony between people and the natural world. (True)",
  "Ecological products are produced in excellent conditions with minimal use of chemical substances. (True)")

table1 <- data.frame(Item = paste0(1:8, ". ", items_cc))
correct <- sapply(1:8, function(i) sum(dat[[paste0("ComprehensionCheck", i, "_score")]]))
table1$`Correct, n (%)`   <- paste0(correct, " (", pct(correct, N), "%)")
table1$`Incorrect, n (%)` <- paste0(N - correct, " (", pct(N - correct, N), "%)")
kable(table1)
```

| Item | Correct, n (%) | Incorrect, n (%) |
|:---|:---|:---|
| 1\. The production of ecological products involves the intensive use of pesticides and artificial fertilisers in order to ensure high yields. (False) | 127 (98.4%) | 2 (1.6%) |
| 2\. Ecological products must meet defined technical standards and undergo comprehensive quality control. (True) | 126 (97.7%) | 3 (2.3%) |
| 3\. Ecological products differ from conventional foods only in price and packaging appearance. (False) | 123 (95.3%) | 6 (4.7%) |
| 4\. The use of the ecological product logo is not regulated by any regulations. (False) | 121 (93.8%) | 8 (6.2%) |
| 5\. Ecological products may compromise safety standards in favour of environmental benefits. (False) | 119 (92.2%) | 10 (7.8%) |
| 6\. The ecological product logo is also meant to remind people of the importance of protecting the environment and preventing pollution. (True) | 119 (92.2%) | 10 (7.8%) |
| 7\. The logo of ecological products is intended to strengthen harmony between people and the natural world. (True) | 106 (82.2%) | 23 (17.8%) |
| 8\. Ecological products are produced in excellent conditions with minimal use of chemical substances. (True) | 120 (93.0%) | 9 (7.0%) |

Comprehension score: *M* = 7.45 of 8, *SD* = 0.98 (93.1% correct); 85
participants (65.9%) answered all eight items correctly.

# 6. Perceived comprehensibility and credibility of the vignette

One-sample *t* tests against the scale midpoint (4). SPSS:
`T-TEST /TESTVAL=4 /VARIABLES=Vignette_Comprehension Vignette_Credibility`.

``` r
ost <- rbind(one_sample(dat$Vignette_Comprehension),
             one_sample(dat$Vignette_Credibility))
ost_out <- data.frame(
  Rating = c("Comprehensibility", "Credibility"),
  M = fmt(ost$M), SD = fmt(ost$SD), t = fmt(ost$t), df = ost$df, p = fmtp(ost$p),
  d = fmt(ost$d), `95% CI d` = paste0("[", fmt(ost$d_low), ", ", fmt(ost$d_high), "]"),
  check.names = FALSE)
kable(ost_out)
```

| Rating            | M    | SD   | t     |  df | p       | d    | 95% CI d       |
|:------------------|:-----|:-----|:------|----:|:--------|:-----|:---------------|
| Comprehensibility | 5.81 | 1.43 | 14.41 | 128 | \< .001 | 1.27 | \[1.04, 1.50\] |
| Credibility       | 5.30 | 1.60 | 9.22  | 128 | \< .001 | 0.81 | \[0.61, 1.01\] |

# 7. Ecological products as virtue vs. vice foods

## 7.1 Reliability of the virtue and vice indices

``` r
virtue_items <- c("Virtue1", "Virtue2", "Virtue3")
vice_items   <- c("Vice1", "Vice2", "Vice3")
rel_vv <- data.frame(
  Index = c("Virtue (health)", "Vice (hedonic)"),
  alpha = no0(c(cronbach_alpha(dat[virtue_items]), cronbach_alpha(dat[vice_items]))),
  omega = no0(c(mcdonald_omega(dat, virtue_items), mcdonald_omega(dat, vice_items))))
kable(rel_vv, col.names = c("Index", "Cronbach's alpha", "McDonald's omega"))
```

| Index           | Cronbach’s alpha | McDonald’s omega |
|:----------------|:-----------------|:-----------------|
| Virtue (health) | .869             | .875             |
| Vice (hedonic)  | .873             | .879             |

## 7.2 Virtue vs. vice: paired *t* test

SPSS: `T-TEST PAIRS=Virtue_index WITH Vice_index (PAIRED)`.

``` r
pt_vv <- t.test(dat$Virtue_index, dat$Vice_index, paired = TRUE)
diff  <- dat$Virtue_index - dat$Vice_index
dz    <- mean(diff) / sd(diff)
dz_ci <- d_ci(unname(pt_vv$statistic), N - 1, N)
paired_out <- data.frame(
  `Virtue M (SD)` = paste0(fmt(mean(dat$Virtue_index)), " (", fmt(sd(dat$Virtue_index)), ")"),
  `Vice M (SD)`   = paste0(fmt(mean(dat$Vice_index)), " (", fmt(sd(dat$Vice_index)), ")"),
  r = no0(cor(dat$Virtue_index, dat$Vice_index)),
  `M diff [95% CI]` = paste0(fmt(mean(diff)), " [", fmt(pt_vv$conf.int[1]), ", ",
                             fmt(pt_vv$conf.int[2]), "]"),
  t = fmt(pt_vv$statistic), df = N - 1, p = fmtp(pt_vv$p.value),
  d_z = fmt(dz), `95% CI d_z` = paste0("[", fmt(dz_ci[1]), ", ", fmt(dz_ci[2]), "]"),
  check.names = FALSE)
kable(paired_out)
```

|  | Virtue M (SD) | Vice M (SD) | r | M diff \[95% CI\] | t | df | p | d_z | 95% CI d_z |
|:---|:---|:---|:---|:---|:---|---:|:---|:---|:---|
| t | 5.46 (1.33) | 4.38 (1.39) | .698 | 1.09 \[0.90, 1.27\] | 11.64 | 128 | \< .001 | 1.02 | \[0.81, 1.24\] |

## 7.3 Healthier than conventional products

One-sample *t* test against the midpoint (4).

``` r
hc <- one_sample(dat$Healthiness_comparative)
hc_out <- data.frame(M = fmt(hc$M), SD = fmt(hc$SD), t = fmt(hc$t), df = hc$df, p = fmtp(hc$p),
                     d = fmt(hc$d), `95% CI d` = paste0("[", fmt(hc$d_low), ", ", fmt(hc$d_high), "]"),
                     check.names = FALSE)
kable(hc_out)
```

| M    | SD   | t     |  df | p       | d    | 95% CI d       |
|:-----|:-----|:------|----:|:--------|:-----|:---------------|
| 5.43 | 1.60 | 10.13 | 128 | \< .001 | 0.89 | \[0.69, 1.09\] |

# 8. Purchase intention scale

Item 3 (“I am not willing to buy ecological products”) is reverse-scored
(`Intention3_recoded` = 8 – `Intention3`).

``` r
pi_items <- c("Intention1", "Intention2", "Intention3_recoded", paste0("Intention", 4:10))
pi_text <- c(
  "1. I am willing to spend more on ecological products because they are worth the extra cost.",
  "2. During my next shopping trip I plan to buy food marked with the ecological product logo.",
  "3. I am not willing to buy ecological products. (reversed)",
  "4. I intend to buy ecological products in the future.",
  "5. I always willingly buy more ecological products for my family's needs.",
  "6. I always try to seek out ecological products, even if it requires travelling out of town.",
  "7. If I had to do it again, I would buy ecological products.",
  "8. I make an effort to buy ecological products.",
  "9. I plan to buy ecological products.",
  "10. I intend to buy ecological products.")
pi_source <- c(rep("Tang et al. (2025, Exp. 1a); White & Peloza (2009)", 3),
               rep("Kumar et al. (2023)", 3),
               rep("Tang et al. (2025, Exp. 2 & 3); Lee & Yun (2015); Loebnitz et al. (2022)", 4))
```

## 8.1 Three-factor CFA (items grouped by source)

Maximum-likelihood estimation, as in JASP 0.19.1 (which uses `lavaan`).

``` r
model3 <- '
  Tang1a =~ Intention1 + Intention2 + Intention3_recoded
  Kumar  =~ Intention4 + Intention5 + Intention6
  Tang23 =~ Intention7 + Intention8 + Intention9 + Intention10
'
fit3 <- cfa(model3, data = dat, estimator = "ML", std.lv = TRUE)

# Is the solution admissible?
lavInspect(fit3, "post.check")
```

    ## [1] FALSE

``` r
# Latent correlations: values > 1 make the latent covariance matrix non-positive definite
cor_lv <- unclass(lavInspect(fit3, "cor.lv"))
kable(fmt(cor_lv, 3))
```

|        | Tang1a | Kumar | Tang23 |
|:-------|:-------|:------|:-------|
| Tang1a | 1.000  | 1.027 | 0.913  |
| Kumar  | 1.027  | 1.000 | 1.029  |
| Tang23 | 0.913  | 1.029 | 1.000  |

``` r
round(eigen(cor_lv)$values, 3)
```

    ## [1]  2.980  0.087 -0.067

Two latent correlations exceed 1 and one eigenvalue of the latent
covariance matrix is negative. The three-factor solution is therefore
inadmissible: the three groups of items cannot be separated empirically.
`lavaan` gives the corresponding warning (suppressed here), and JASP
reports that the latent covariance matrix is not positive definite.

## 8.2 One-factor CFA (Table 2)

Fit of the one-factor model (lavaan’s default χ² is based on *N*, as in
JASP):

``` r
model1 <- paste("PI =~", paste(pi_items, collapse = " + "))
fit1 <- cfa(model1, data = dat, estimator = "ML", std.lv = TRUE)

fm <- fitMeasures(fit1, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea",
                          "rmsea.ci.lower", "rmsea.ci.upper"))
fm <- setNames(as.numeric(fm), names(fm))  # drop lavaan's print class
fit_out <- data.frame(
  `chi2` = fmt(fm["chisq"]), df = fm["df"], p = fmtp(fm["pvalue"]),
  CFI = no0(fm["cfi"]), TLI = no0(fm["tli"]), RMSEA = no0(fm["rmsea"]),
  `RMSEA 90% CI` = paste0("[", no0(fm["rmsea.ci.lower"]), ", ", no0(fm["rmsea.ci.upper"]), "]"),
  check.names = FALSE, row.names = NULL)
kable(fit_out)
```

| chi2   |  df | p       | CFI  | TLI  | RMSEA | RMSEA 90% CI   |
|:-------|----:|:--------|:-----|:-----|:------|:---------------|
| 106.81 |  35 | \< .001 | .933 | .914 | .126  | \[.099, .154\] |

**Table 2.** Standardised factor loadings (λ), explained variance (*R*²)
and item descriptives (*N* = 129).

``` r
std <- standardizedSolution(fit1) %>% filter(op == "=~")
r2  <- lavInspect(fit1, "r2")[pi_items]

table2 <- data.frame(
  Source = pi_source, Item = pi_text,
  lambda = no0(std$est.std), `R2` = no0(r2), p = fmtp(std$pvalue),
  M  = fmt(sapply(dat[pi_items], mean), 2),
  SD = fmt(sapply(dat[pi_items], sd), 3),
  check.names = FALSE, row.names = NULL)
kable(table2)
```

| Source | Item | lambda | R2 | p | M | SD |
|:---|:---|:---|:---|:---|:---|:---|
| Tang et al. (2025, Exp. 1a); White & Peloza (2009) | 1\. I am willing to spend more on ecological products because they are worth the extra cost. | .737 | .544 | \< .001 | 3.79 | 1.694 |
| Tang et al. (2025, Exp. 1a); White & Peloza (2009) | 2\. During my next shopping trip I plan to buy food marked with the ecological product logo. | .757 | .573 | \< .001 | 3.18 | 1.548 |
| Tang et al. (2025, Exp. 1a); White & Peloza (2009) | 3\. I am not willing to buy ecological products. (reversed) | .686 | .470 | \< .001 | 4.72 | 1.776 |
| Kumar et al. (2023) | 4\. I intend to buy ecological products in the future. | .861 | .741 | \< .001 | 4.49 | 1.682 |
| Kumar et al. (2023) | 5\. I always willingly buy more ecological products for my family’s needs. | .739 | .546 | \< .001 | 3.12 | 1.723 |
| Kumar et al. (2023) | 6\. I always try to seek out ecological products, even if it requires travelling out of town. | .573 | .328 | \< .001 | 2.01 | 1.482 |
| Tang et al. (2025, Exp. 2 & 3); Lee & Yun (2015); Loebnitz et al. (2022) | 7\. If I had to do it again, I would buy ecological products. | .843 | .710 | \< .001 | 4.15 | 1.625 |
| Tang et al. (2025, Exp. 2 & 3); Lee & Yun (2015); Loebnitz et al. (2022) | 8\. I make an effort to buy ecological products. | .848 | .720 | \< .001 | 3.43 | 1.749 |
| Tang et al. (2025, Exp. 2 & 3); Lee & Yun (2015); Loebnitz et al. (2022) | 9\. I plan to buy ecological products. | .893 | .798 | \< .001 | 4.20 | 1.669 |
| Tang et al. (2025, Exp. 2 & 3); Lee & Yun (2015); Loebnitz et al. (2022) | 10\. I intend to buy ecological products. | .937 | .878 | \< .001 | 4.23 | 1.623 |

## 8.3 Reliability of the 10-item scale

``` r
txt$pi_alpha <- no0(cronbach_alpha(dat[pi_items]))
txt$pi_omega <- no0(mcdonald_omega(dat, pi_items))
```

Cronbach’s α = .943; McDonald’s ω = .944.

# 9. Comparison with the report

All statistics in the Pilot Study 2 report and in the SPSS output
(`Analizy.spv`) are reproduced, apart from these points:

1.  **Completion times (Section 2 of the report).** The report gives *M*
    = 6.56 min, *Mdn* = 3.59 min, *SD* = 11.19 min for the whole survey.
    These are min:s values (6:56, 3:59, 11:20) written as decimal
    minutes. In decimal minutes the values are *M* = 6.94, *Mdn* = 3.98,
    *SD* = 11.33 (see Section 3 above).
2.  **Education (Section 1.1).** The percentages for “higher education
    degree” and “incomplete higher education” are swapped. The data show
    higher education *n* = 23 (17.8%) and incomplete higher education
    *n* = 55 (42.6%). The latter matches the report’s own *n* = 55 for
    current students.
3.  **Psychology students.** The data contain 39 of 55 students (70.9%)
    in psychology, plus one student with no field reported. The report
    gives 40 (72.7%). Please check which figure is correct.
4.  **Variable labels.** The labels of `Virtue_index` and `Vice_index`
    in the .sav file list α = .871 / ω = .877 and α = .871 / ω = .876.
    Those values do not match the current data. The data, the SPSS
    output and the report all give virtue α = .869, ω = .875 and vice α
    = .873, ω = .879.
5.  **Item 7 wording in Table 2** reads “organic products” in the report
    and “ecological products” in the data file and questionnaire.

# Session information

``` r
sessionInfo()
```

    ## R version 4.5.1 (2025-06-13)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.5 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=cs_CZ.UTF-8       
    ##  [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=cs_CZ.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=cs_CZ.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
    ## [10] LC_TELEPHONE=C             LC_MEASUREMENT=cs_CZ.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Europe/Prague
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] lavaan_0.6-17 knitr_1.45    dplyr_1.1.4   haven_2.5.4  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] crayon_1.5.3     vctrs_0.6.5      cli_3.6.5        rlang_1.1.6      xfun_0.41       
    ##  [6] forcats_1.0.0    renv_1.0.3       generics_0.1.4   glue_1.8.0       pbivnorm_0.6.0  
    ## [11] htmltools_0.5.7  stats4_4.5.1     hms_1.1.3        rmarkdown_2.25   quadprog_1.5-8  
    ## [16] evaluate_0.23    tibble_3.3.0     MASS_7.3-65      tzdb_0.5.0       fastmap_1.1.1   
    ## [21] yaml_2.3.12      lifecycle_1.0.4  compiler_4.5.1   pkgconfig_2.0.3  digest_0.6.37   
    ## [26] R6_2.6.1         readr_2.1.5      tidyselect_1.2.1 parallel_4.5.1   mnormt_2.1.1    
    ## [31] pillar_1.10.2    magrittr_2.0.4   withr_3.0.2      tools_4.5.1
