Pilot Study 3 (Polish replication of Tang et al., 2025): Reproducible
Analysis
================

This document reproduces every analysis reported in *Pilot Study 3:
Validation of the Affect Induction*, which is summarised in the
manuscript (Study 2, Method: “Pilot Study 3, N = 194; d = 0.48, 95% CI
\[0.19, 0.78\]”). It re-runs the original SPSS syntax (`Analizy.spv`) in
R. Knit this file in RStudio (**Knit** button) or with
`rmarkdown::render("analysis_pilot3.Rmd")` to regenerate
`analysis_pilot3.md`.

**Data source.** By default, the document downloads
`pilot3/data_pilot3_anonymized.sav` from the GitHub repository
[janurban/metarials-love-study](https://github.com/janurban/metarials-love-study).
It can also read a local `.sav` file by passing a path with
`rmarkdown::render("analysis_pilot3.Rmd", params = list(data = "path/to/file.sav"))`.

**Required packages.** `haven`, `dplyr`, `knitr`. Install them with
`install.packages(c("haven", "dplyr", "knitr", "rmarkdown"))`.

# 0. Setup

``` r
library(haven)   # read the SPSS .sav file
library(dplyr)   # data handling
library(knitr)   # tables
```

Helper functions. Group comparisons follow the SPSS command
`T-TEST GROUPS=Condition(2 1) /ES DISPLAY(TRUE) /HOMOGENEITY DISPLAY(TRUE)`:

- Differences are love minus neutral.
- Levene’s test is based on the mean, as in SPSS.
- The report uses the Student *t* test when Levene’s *p* ≥ .05 and the
  Welch *t* test when Levene’s *p* \< .05.
- Cohen’s *d* always uses the pooled *SD*. Its 95% CI comes from the
  noncentral *t* distribution of the pooled-variance *t*, which is what
  SPSS prints in the “Independent Samples Effect Sizes” table.

``` r
# number formatting
fmt  <- function(x, d = 2) formatC(x, format = "f", digits = d)
no0  <- function(x, d = 3) sub("^(-?)0\\.", "\\1.", fmt(x, d))   # APA style: drop leading zero
fmtp <- function(p) ifelse(p < .001, "< .001", paste("=", no0(p, 3)))
pct  <- function(n, N) fmt(100 * n / N, 1)

# Levene's test based on the mean (SPSS default)
levene_mean <- function(y, g) {
  z <- abs(y - ave(y, g))
  a <- anova(lm(z ~ factor(g)))
  c(F = a$`F value`[1], df1 = a$Df[1], df2 = a$Df[2], p = a$`Pr(>F)`[1])
}

# 95% CI for Cohen's d (two independent groups) via the noncentral t distribution
d_ci_2 <- function(t, df, n1, n2, level = .95) {
  a  <- (1 - level) / 2
  lo <- uniroot(function(ncp) pt(t, df, ncp, lower.tail = FALSE) - a, c(t - 10, t + 10))$root
  hi <- uniroot(function(ncp) pt(t, df, ncp) - a, c(t - 10, t + 10))$root
  c(lo, hi) / sqrt(n1 * n2 / (n1 + n2))
}

# Love vs. neutral comparison for one variable
compare_groups <- function(y, cond) {
  x1 <- y[cond == "Love"];    x1 <- x1[!is.na(x1)]
  x2 <- y[cond == "Neutral"]; x2 <- x2[!is.na(x2)]
  n1 <- length(x1); n2 <- length(x2)
  lev     <- levene_mean(c(x1, x2), rep(1:2, c(n1, n2)))
  student <- t.test(x1, x2, var.equal = TRUE)
  welch   <- t.test(x1, x2, var.equal = FALSE)
  use     <- if (lev["p"] < .05) welch else student
  sd_pooled <- sqrt(((n1 - 1) * var(x1) + (n2 - 1) * var(x2)) / (n1 + n2 - 2))
  d  <- (mean(x1) - mean(x2)) / sd_pooled
  ci <- d_ci_2(unname(student$statistic), n1 + n2 - 2, n1, n2)
  data.frame(
    M_love = mean(x1), SD_love = sd(x1), n_love = n1,
    M_neutral = mean(x2), SD_neutral = sd(x2), n_neutral = n2,
    levene_F = unname(lev["F"]), levene_p = unname(lev["p"]),
    test = if (lev["p"] < .05) "Welch" else "Student",
    t = unname(use$statistic), df = unname(use$parameter), p = use$p.value,
    d = d, d_low = ci[1], d_high = ci[2])
}
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

    ## [1] 204  22

Variables are recoded to English labels. `Education` is stored as
LimeSurvey answer codes: `AO01` = primary, `AO02` = secondary without
high school diploma (*matura*), `AO03` = secondary with diploma, `AO05`
= higher, and `-oth` = “Other”. Participants who chose “Other” typed in
their field of study (`Education_open`), so the report counts them as
current students.

``` r
dat_all <- dat_raw %>%
  zap_labels() %>%
  mutate(
    Condition_f = factor(Condition, levels = c(2, 1), labels = c("Love", "Neutral")),
    Gender_f    = factor(Sex, levels = 1:4,
                         labels = c("Female", "Male", "Non-binary or gender-diverse",
                                    "Prefer not to say")),
    Education_f = factor(trimws(Education),
                         levels = c("AO01", "AO02", "AO03", "AO05", "-oth"),
                         labels = c("Primary", "Secondary without high school diploma",
                                    "Secondary with high school diploma", "Higher",
                                    "Other (currently studying)")),
    Text = NA_character_,
    Text_chars = NA_integer_)
```

## 1.1 Exclusions

Participants whose native language was not Polish were excluded
(`Excluded = 1`; SPSS: `FILTER BY filter_$` with `Excluded = 2`). The
`Excluded` flag is identical to `Language1 = 2`.

``` r
all(xor(dat_all$Excluded == 1, dat_all$Language1 == 1))   # flag matches language
```

    ## [1] TRUE

``` r
dat_all %>%
  filter(Excluded == 1) %>%
  count(`Native language` = trimws(Language_other), name = "n") %>%
  kable()
```

| Native language      |   n |
|:---------------------|----:|
| Niemiecki            |   1 |
| Rosyjski             |   1 |
| Ukraiński            |   4 |
| Ukraiński i rosyjski |   1 |
| francuski            |   1 |
| ukraiński            |   2 |

``` r
dat <- dat_all %>% filter(Excluded == 2)
N <- nrow(dat)
c(before = nrow(dat_all), excluded = sum(dat_all$Excluded == 1), analysed = N)
```

    ##   before excluded analysed 
    ##      204       10      194

# 2. Participants

``` r
tab_freq <- function(f) {
  n <- table(f)
  data.frame(Category = names(n), n = as.vector(n), Percent = pct(as.vector(n), N))
}
age <- c(M = mean(dat$Age), Mdn = median(dat$Age), SD = sd(dat$Age),
         Min = min(dat$Age), Max = max(dat$Age))
fmt(age)
```

    ##       M     Mdn      SD     Min     Max 
    ## "25.80" "22.00"  "9.68" "18.00" "64.00"

**Gender**

``` r
kable(tab_freq(dat$Gender_f))
```

| Category                     |   n | Percent |
|:-----------------------------|----:|:--------|
| Female                       | 122 | 62.9    |
| Male                         |  67 | 34.5    |
| Non-binary or gender-diverse |   2 | 1.0     |
| Prefer not to say            |   3 | 1.5     |

**Highest completed education**

``` r
kable(tab_freq(dat$Education_f))
```

| Category                              |   n | Percent |
|:--------------------------------------|----:|:--------|
| Primary                               |   3 | 1.5     |
| Secondary without high school diploma |   7 | 3.6     |
| Secondary with high school diploma    |  24 | 12.4    |
| Higher                                |  57 | 29.4    |
| Other (currently studying)            | 103 | 53.1    |

**Condition** (random assignment)

``` r
kable(tab_freq(dat$Condition_f))
```

| Category |   n | Percent |
|:---------|----:|:--------|
| Love     |  77 | 39.7    |
| Neutral  | 117 | 60.3    |

# 3. Completion and writing times

`Survey_time` (whole study) and `Writing_time` (time on the writing-task
page, including the 30-s minimum) are stored in seconds.

``` r
time_row <- function(x, label) data.frame(
  Variable = label,
  `M (s)` = fmt(mean(x)), `Mdn (s)` = fmt(median(x)), `SD (s)` = fmt(sd(x)),
  `M (min)` = fmt(mean(x) / 60), `Mdn (min)` = fmt(median(x) / 60),
  `SD (min)` = fmt(sd(x) / 60),
  check.names = FALSE)
kable(rbind(
  time_row(dat$Survey_time, "Whole study"),
  time_row(dat$Writing_time, "Writing page, all"),
  time_row(dat$Writing_time[dat$Condition_f == "Love"], "Writing page, love"),
  time_row(dat$Writing_time[dat$Condition_f == "Neutral"], "Writing page, neutral"),
  time_row(dat$Writing_time - 30, "Writing page beyond the 30-s minimum, all")))
```

| Variable | M (s) | Mdn (s) | SD (s) | M (min) | Mdn (min) | SD (min) |
|:---|:---|:---|:---|:---|:---|:---|
| Whole study | 613.23 | 290.52 | 1986.58 | 10.22 | 4.84 | 33.11 |
| Writing page, all | 454.82 | 176.99 | 1967.26 | 7.58 | 2.95 | 32.79 |
| Writing page, love | 333.02 | 192.33 | 418.10 | 5.55 | 3.21 | 6.97 |
| Writing page, neutral | 534.97 | 163.84 | 2511.61 | 8.92 | 2.73 | 41.86 |
| Writing page beyond the 30-s minimum, all | 424.82 | 146.99 | 1967.26 | 7.08 | 2.45 | 32.79 |

Length of the written texts was not reported in the anonymized data
file. The original text fields are omitted, so character counts are not
available.

``` r
data.frame(
  Condition = c("Love", "Neutral"),
  `Character counts` = c("Not available", "Not available"),
  check.names = FALSE
) %>% kable()
```

| Condition | Character counts |
|:----------|:-----------------|
| Love      | Not available    |
| Neutral   | Not available    |

# 4. Manipulation check

“Right now, I feel loved” (1 = strongly disagree, 7 = strongly agree).

``` r
mc <- compare_groups(dat$ManipCheck, dat$Condition_f)
```

``` r
res_table <- function(res, labels) data.frame(
  Measure = labels,
  `Love M (SD)`    = paste0(fmt(res$M_love), " (", fmt(res$SD_love), ")"),
  `Neutral M (SD)` = paste0(fmt(res$M_neutral), " (", fmt(res$SD_neutral), ")"),
  `Levene F (p)`   = paste0(fmt(res$levene_F), " (", no0(res$levene_p), ")"),
  Test = res$test,
  t = fmt(res$t), df = ifelse(res$test == "Welch", fmt(res$df, 1), fmt(res$df, 0)),
  p = fmtp(res$p), d = fmt(res$d),
  `95% CI d` = paste0("[", fmt(res$d_low), ", ", fmt(res$d_high), "]"),
  check.names = FALSE)
kable(res_table(mc, "Feeling loved"))
```

| Measure | Love M (SD) | Neutral M (SD) | Levene F (p) | Test | t | df | p | d | 95% CI d |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| Feeling loved | 5.90 (1.47) | 5.13 (1.66) | 2.27 (.134) | Student | 3.30 | 192 | = .001 | 0.48 | \[0.19, 0.78\] |

``` r
c(n_love = mc$n_love, n_neutral = mc$n_neutral)
```

    ##    n_love n_neutral 
    ##        77       117

# 5. Evaluation of the writing task

``` r
task_vars <- c(Task_comprehension = "Instructions comprehensible",
               Task_easy          = "Task easy",
               Task_effort        = "Task effortful",
               Task_engagement    = "Task engaging",
               Task_time          = "30-s minimum appropriate")
task_res <- do.call(rbind, lapply(names(task_vars),
                                  function(v) compare_groups(dat[[v]], dat$Condition_f)))
kable(res_table(task_res, unname(task_vars)))
```

| Measure | Love M (SD) | Neutral M (SD) | Levene F (p) | Test | t | df | p | d | 95% CI d |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| Instructions comprehensible | 6.65 (0.89) | 6.30 (1.30) | 7.31 (.007) | Welch | 2.23 | 191.8 | = .027 | 0.30 | \[0.01, 0.59\] |
| Task easy | 4.81 (1.75) | 5.51 (1.74) | 0.33 (.567) | Student | -2.77 | 192 | = .006 | -0.41 | \[-0.70, -0.11\] |
| Task effortful | 3.82 (1.85) | 2.92 (1.82) | 0.01 (.921) | Student | 3.33 | 192 | = .001 | 0.49 | \[0.20, 0.78\] |
| Task engaging | 5.88 (1.17) | 5.14 (1.66) | 12.81 (.000) | Welch | 3.68 | 191.0 | \< .001 | 0.50 | \[0.21, 0.79\] |
| 30-s minimum appropriate | 4.77 (1.95) | 4.91 (2.20) | 2.49 (.116) | Student | -0.48 | 192 | = .631 | -0.07 | \[-0.36, 0.22\] |

The report also compares the condition means with the scale midpoint (4)
descriptively. For ease, both means are above 4; for effort, both are
below 4.

# 6. Comparison with the report

All statistics in the Pilot Study 3 report match the SPSS output
(`Analizy.spv`) and the analyses above, apart from these points:

1.  **Manipulation check, *p* value.** The report gives *p* \< .001. The
    two-tailed *p* is .001 (*p* = .00115). The one-tailed *p* is .0006.
2.  **Manipulation check, *SD* in the love condition.** The report gives
    1.46. The value is 1.47 (1.4653).
3.  **Engagement, 95% CI for *d*.** The report gives \[0.21, 0.80\]. The
    upper limit is 0.79 (0.7949).
4.  **Education.** The report gives higher education *n* = 58 (29.9%)
    and secondary education without the diploma *n* = 6 (3.1%). The data
    give *n* = 57 (29.4%) and *n* = 7 (3.6%). The other categories
    match.
5.  **Writing time.** The report says participants wrote for a median of
    3.0 min “after the forced-pacing minimum had elapsed”. The median of
    2.95 min (177 s) is the total time on the writing page, *including*
    the 30-s minimum. The median time beyond the minimum is about 2.45
    min.
6.  **Typo in the summary.** It begins “Pilot Study 2 indicated…”; this
    should be Pilot Study 3.

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
