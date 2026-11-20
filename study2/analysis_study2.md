Study 2 (Polish replication of Tang et al., 2025, Experiment 1a):
Reproducible Analysis
================

This document reproduces every analysis of Study 2 reported in the
manuscript (Study 2: Method and Results; General Discussion) and in
Supplementary Materials S2.1–S2.4. It re-runs the original SPSS syntax
(`Analizy.spv`) in R. Knit this file in RStudio (**Knit** button) or
with `rmarkdown::render("analysis_study2.Rmd")` to regenerate
`analysis_study2.md` in the same folder.

**Data source.** By default, the document downloads
`study2/data_study2_anonymized.sav` from the GitHub repository
[janurban/metarials-love-study](https://github.com/janurban/metarials-love-study).
It can also read a local `.sav` file by passing a path with
`rmarkdown::render("analysis_study2.Rmd", params = list(data = "path/to/file.sav"))`.
The analyses do not use the free-text responses (`control_condition`,
`experimental_love_condition`), so they run on the anonymized file
without them.

**Required packages.** `haven`, `dplyr`, `knitr`, `lavaan`. Install them
with
`install.packages(c("haven", "dplyr", "knitr", "lavaan", "rmarkdown"))`.

# 0. Setup

``` r
library(haven)   # read the SPSS .sav file
library(dplyr)   # data handling
library(knitr)   # tables
library(lavaan)  # McDonald's omega
```

Helper functions. Group comparisons follow the SPSS command
`T-TEST GROUPS=condition(2 1) /ES DISPLAY(TRUE) /HOMOGENEITY DISPLAY(TRUE)`:

- Differences are feeling loved minus control.
- Levene’s test is based on the mean, as in SPSS.
- The manuscript uses the Student *t* test when Levene’s *p* ≥ .05 and
  the Welch *t* test when Levene’s *p* \< .05.
- Cohen’s *d* uses the pooled *SD*. Its 95% CI comes from the noncentral
  *t* distribution of the pooled-variance *t*, as in SPSS’s “Independent
  Samples Effect Sizes” table.

Fleiss’s kappa follows the SPSS command `FLEISS MULTIRATER KAPPA`. Cases
with a missing rating are excluded (`CLASSMISSING=EXCLUDE`). The 95% CIs
use the asymptotic standard errors of Fleiss, Nee and Landis (1979),
which SPSS prints.

``` r
# number formatting
fmt  <- function(x, d = 2) formatC(x, format = "f", digits = d)
no0  <- function(x, d = 3) sub("^(-?)0\\.", "\\1.", fmt(x, d))   # APA style: drop leading zero
fmtp <- function(p) ifelse(p < .001, "< .001", paste("=", no0(p, 3)))
pct  <- function(n, N) fmt(100 * n / N, 1)
ci_txt <- function(lo, hi, f = fmt) paste0("[", f(lo), ", ", f(hi), "]")

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

# Feeling loved vs. control comparison for one variable
compare_groups <- function(y, cond) {
  x1 <- y[cond == "Feeling loved"]; x1 <- x1[!is.na(x1)]
  x2 <- y[cond == "Control"];       x2 <- x2[!is.na(x2)]
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
    M_control = mean(x2), SD_control = sd(x2), n_control = n2,
    levene_F = unname(lev["F"]), levene_p = unname(lev["p"]),
    test = if (lev["p"] < .05) "Welch" else "Student",
    t = unname(use$statistic), df = unname(use$parameter), p = use$p.value,
    d = d, d_low = ci[1], d_high = ci[2])
}

# Cronbach's alpha
cronbach_alpha <- function(X) {
  X <- as.matrix(X); k <- ncol(X)
  k / (k - 1) * (1 - sum(apply(X, 2, var)) / var(rowSums(X)))
}

# McDonald's omega from a one-factor ML model (the estimator of SPSS RELIABILITY /MODEL=OMEGA)
mcdonald_omega <- function(data, items) {
  fit <- cfa(paste("f =~", paste(items, collapse = " + ")), data = data, std.lv = TRUE)
  pe  <- parameterEstimates(fit)
  lam <- pe$est[pe$op == "=~"]
  th  <- pe$est[pe$op == "~~" & pe$lhs == pe$rhs & pe$lhs %in% items]
  sum(lam)^2 / (sum(lam)^2 + sum(th))
}

# Fleiss's kappa: overall and category-specific, with SPSS's asymptotic CIs
fleiss_kappa <- function(ratings, categories, level = .95) {
  R <- as.matrix(ratings)
  R <- R[complete.cases(R), , drop = FALSE]
  n <- nrow(R); m <- ncol(R)
  nij <- sapply(categories, function(k) rowSums(R == k))   # cases x categories
  p  <- colSums(nij) / (n * m); q <- 1 - p; pq <- p * q
  P_i   <- rowSums(nij * (nij - 1)) / (m * (m - 1))
  P_e   <- sum(p^2)
  kappa <- (mean(P_i) - P_e) / (1 - P_e)
  se    <- sqrt(2) / (sum(pq) * sqrt(n * m * (m - 1))) *
           sqrt(sum(pq)^2 - sum(pq * (q - p)))
  k_j   <- 1 - colSums(nij * (m - nij)) / (n * m * (m - 1) * pq)
  se_j  <- sqrt(2 / (n * m * (m - 1)))
  cond  <- colSums(nij * (nij - 1)) / (colSums(nij) * (m - 1))   # conditional probability
  z <- qnorm(1 - (1 - level) / 2)
  list(n = n, raters = m, kappa = kappa, se = se, ci = kappa + c(-1, 1) * z * se,
       by_category = data.frame(code = categories, cond_prob = cond, kappa = k_j,
                                se = se_j, ci_low = k_j - z * se_j, ci_high = k_j + z * se_j))
}

# Pearson correlation with Fisher-z 95% CI
cor_ci <- function(x, y) {
  ct <- cor.test(x, y)
  data.frame(r = unname(ct$estimate), df = unname(ct$parameter),
             lo = ct$conf.int[1], hi = ct$conf.int[2], p = ct$p.value)
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

    ## [1] 1020   46

Variables are recoded to English labels. Written responses were coded
1–6 (1–2 valid, 3–6 invalid; Table S2.2). `EXCLUSIONS` marks valid (1)
and invalid (2) responses. In the love-content codes, `0` means “does
not describe feeling loved”; this code has no value label in the data
file.

``` r
excl_labels <- c("1. Valid (specific event)", "2. Valid (general statement)",
                 "3. No event recalled", "4. Refusal or comment only",
                 "5. Off-topic content", "6. Nonsensical or copied text")
source_labels <- c("Does not describe feeling loved", "Romantic partner", "Family members",
                   "People outside the family", "Animal", "God or a spiritual being",
                   "Oneself", "Several sources or unspecified")
form_labels <- c("Does not describe feeling loved", "Trust and acceptance",
                 "Support in needs and goals", "Symbolic and physical expressions",
                 "Shared time and presence", "Other")

dat <- dat_raw %>%
  zap_labels() %>%
  mutate(
    Condition_f = factor(condition, levels = c(2, 1), labels = c("Feeling loved", "Control")),
    Sex_f       = factor(sex, levels = 1:4,
                         labels = c("Woman", "Man", "Non-binary or gender-diverse",
                                    "Prefer not to say")),
    Age_group_f = factor(age, levels = 1:5,
                         labels = c("15–17 years", "18–35 years", "36–50 years",
                                    "51–70 years", "71 years or older")),
    Education_f = factor(education, levels = 1:4,
                         labels = c("Primary or lower secondary", "Basic vocational",
                                    "Secondary or post-secondary", "Higher")),
    Final_code_f = factor(Final_decision_exclusion, levels = 1:6, labels = excl_labels),
    Valid = EXCLUSIONS == 1)

N <- nrow(dat)
pi_items <- c("Purchase_intention1", "Purchase_intention2", "Purchase_intention3_RECODED",
              paste0("Purchase_intention", 4:10))
```

Derived variables in the data file are checked against their
definitions:

``` r
c(intention3_recode = all(8 - dat$Purchase_intention3 == dat$Purchase_intention3_RECODED),
  index_3_items     = isTRUE(all.equal(rowMeans(dat[pi_items[1:3]]),
                                       dat$Purchase_intention_index_3item,
                                       check.attributes = FALSE)),
  index_10_items    = isTRUE(all.equal(rowMeans(dat[pi_items]),
                                       dat$Purchase_intention_index_10items,
                                       check.attributes = FALSE)),
  exclusions_match_codes = all((dat$Final_decision_exclusion <= 2) == dat$Valid),
  unique_ids        = !anyDuplicated(dat$ID),
  no_missing_outcomes = !anyNA(dat[c("manipulation_check", pi_items)]))
```

    ##      intention3_recode          index_3_items         index_10_items exclusions_match_codes 
    ##                   TRUE                   TRUE                   TRUE                   TRUE 
    ##             unique_ids    no_missing_outcomes 
    ##                   TRUE                   TRUE

# 2. Sample characteristics (Table S2.1)

``` r
tab_freq <- function(f, NN = length(f)) {
  n <- table(f)
  data.frame(Category = names(n), n = as.vector(n), Percent = pct(as.vector(n), NN))
}
kable(rbind(
  cbind(Characteristic = "Sex", tab_freq(dat$Sex_f)),
  cbind(Characteristic = "Age group", tab_freq(droplevels(dat$Age_group_f))),
  cbind(Characteristic = "Education", tab_freq(dat$Education_f)),
  cbind(Characteristic = "Condition", tab_freq(dat$Condition_f))))
```

| Characteristic | Category                     |   n | Percent |
|:---------------|:-----------------------------|----:|:--------|
| Sex            | Woman                        | 510 | 50.0    |
| Sex            | Man                          | 510 | 50.0    |
| Sex            | Non-binary or gender-diverse |   0 | 0.0     |
| Sex            | Prefer not to say            |   0 | 0.0     |
| Age group      | 18–35 years                  | 316 | 31.0    |
| Age group      | 36–50 years                  | 337 | 33.0    |
| Age group      | 51–70 years                  | 367 | 36.0    |
| Education      | Primary or lower secondary   |   9 | 0.9     |
| Education      | Basic vocational             |  83 | 8.1     |
| Education      | Secondary or post-secondary  | 356 | 34.9    |
| Education      | Higher                       | 572 | 56.1    |
| Condition      | Feeling loved                | 510 | 50.0    |
| Condition      | Control                      | 510 | 50.0    |

``` r
age <- c(M = mean(dat$year), SD = sd(dat$year), Min = min(dat$year), Max = max(dat$year))
fmt(age)
```

    ##       M      SD     Min     Max 
    ## "45.01" "13.49" "18.00" "70.00"

Completion and reading times. `survey_finish_time` and
`time_read_vignette` are in seconds. `time_writing_task` is the time
spent writing after the 30-s lock was released, also in seconds.

``` r
time_row <- function(x, label, scale = 1, unit = "s") data.frame(
  Variable = label, Unit = unit,
  M = fmt(mean(x) / scale), SD = fmt(sd(x) / scale), Mdn = fmt(median(x) / scale),
  Min = fmt(min(x) / scale), Max = fmt(max(x) / scale))
kable(rbind(
  time_row(dat$survey_finish_time, "Whole study", 60, "min"),
  time_row(dat$time_read_vignette, "Reading the vignette"),
  time_row(dat$time_writing_task, "Writing after the 30-s lock")))
```

| Variable                    | Unit | M     | SD    | Mdn   | Min   | Max    |
|:----------------------------|:-----|:------|:------|:------|:------|:-------|
| Whole study                 | min  | 6.84  | 1.83  | 6.66  | 4.02  | 11.57  |
| Reading the vignette        | s    | 45.43 | 14.47 | 43.00 | 23.00 | 89.00  |
| Writing after the 30-s lock | s    | 53.44 | 43.38 | 46.00 | 0.00  | 188.00 |

# 3. Reliability of the purchase intention measures

``` r
rel <- data.frame(
  Measure = c("Primary: 3 items (Tang et al., Exp. 1a)", "Extended: 10 items"),
  alpha = no0(c(cronbach_alpha(dat[pi_items[1:3]]), cronbach_alpha(dat[pi_items]))),
  omega = no0(c(mcdonald_omega(dat, pi_items[1:3]), mcdonald_omega(dat, pi_items))))
kable(rel, col.names = c("Measure", "Cronbach's alpha", "McDonald's omega"))
```

| Measure                                 | Cronbach’s alpha | McDonald’s omega |
|:----------------------------------------|:-----------------|:-----------------|
| Primary: 3 items (Tang et al., Exp. 1a) | .825             | .832             |
| Extended: 10 items                      | .965             | .965             |

# 4. Validity of the written responses (S2.2)

## 4.1 Agreement among the five coders (Table S2.3)

``` r
judges <- paste0("Judge", 1:5, "_exclusions")
fk <- fleiss_kappa(dat[judges], 1:6)
c(responses = fk$n, raters = fk$raters)
```

    ## responses    raters 
    ##      1018         5

``` r
kable(data.frame(kappa = no0(fk$kappa), SE = fmt(fk$se, 4),
                 `95% CI` = ci_txt(fk$ci[1], fk$ci[2], no0), check.names = FALSE))
```

| kappa | SE     | 95% CI         |
|:------|:-------|:---------------|
| .799  | 0.0066 | \[.786, .812\] |

``` r
kable(data.frame(
  Category = excl_labels,
  `Conditional probability` = no0(fk$by_category$cond_prob, 2),
  kappa = no0(fk$by_category$kappa, 2),
  `95% CI` = ci_txt(fk$by_category$ci_low, fk$by_category$ci_high, function(x) no0(x, 2)),
  check.names = FALSE))
```

| Category                       | Conditional probability | kappa | 95% CI       |
|:-------------------------------|:------------------------|:------|:-------------|
| 1\. Valid (specific event)     | .94                     | .80   | \[.79, .82\] |
| 2\. Valid (general statement)  | .75                     | .71   | \[.69, .73\] |
| 3\. No event recalled          | .96                     | .95   | \[.93, .97\] |
| 4\. Refusal or comment only    | .90                     | .90   | \[.88, .92\] |
| 5\. Off-topic content          | .54                     | .52   | \[.51, .54\] |
| 6\. Nonsensical or copied text | .53                     | .53   | \[.51, .55\] |

## 4.2 Final codes by condition (Table S2.4 and note to Table 1)

``` r
code_tab <- function(x) {
  n <- table(factor(x, levels = 1:6))
  c(as.vector(n), sum(n[1:2]), sum(n[3:6]))
}
cnt <- sapply(list(`Feeling loved` = dat$Final_decision_exclusion[dat$Condition_f == "Feeling loved"],
                   Control         = dat$Final_decision_exclusion[dat$Condition_f == "Control"],
                   Total           = dat$Final_decision_exclusion), code_tab)
tot <- c(sum(dat$Condition_f == "Feeling loved"), sum(dat$Condition_f == "Control"), N)
final_tab <- data.frame(
  `Final code` = c(excl_labels, "Valid (Codes 1–2)", "Invalid (Codes 3–6)"),
  sapply(1:3, function(j) paste0(cnt[, j], " (", pct(cnt[, j], tot[j]), ")")),
  check.names = FALSE)
names(final_tab)[2:4] <- paste0(c("Feeling loved", "Control", "Total"), ", n (%)")
kable(final_tab)
```

| Final code | Feeling loved, n (%) | Control, n (%) | Total, n (%) |
|:---|:---|:---|:---|
| 1\. Valid (specific event) | 270 (52.9) | 427 (83.7) | 697 (68.3) |
| 2\. Valid (general statement) | 140 (27.5) | 12 (2.4) | 152 (14.9) |
| 3\. No event recalled | 69 (13.5) | 59 (11.6) | 128 (12.5) |
| 4\. Refusal or comment only | 13 (2.5) | 3 (0.6) | 16 (1.6) |
| 5\. Off-topic content | 14 (2.7) | 8 (1.6) | 22 (2.2) |
| 6\. Nonsensical or copied text | 4 (0.8) | 1 (0.2) | 5 (0.5) |
| Valid (Codes 1–2) | 410 (80.4) | 439 (86.1) | 849 (83.2) |
| Invalid (Codes 3–6) | 100 (19.6) | 71 (13.9) | 171 (16.8) |

``` r
valid_love <- sum(dat$Final_decision_exclusion[dat$Condition_f == "Feeling loved"] %in% 1:2)
c(general_statements_love = sum(dat$Final_decision_exclusion == 2 & dat$Condition_f == "Feeling loved"),
  valid_love = valid_love)
```

    ## general_statements_love              valid_love 
    ##                     140                     410

## 4.3 Check: final code vs. modal code

The manuscript states that the final code is the modal code among the
five coders. The table below lists the responses where the five codes
have no unique mode (a tie), or where the final code differs from the
unique mode.

``` r
mode_info <- function(x) {
  x <- x[!is.na(x)]
  tab <- table(x)
  top <- as.numeric(names(tab)[tab == max(tab)])
  paste(top, collapse = "/")
}
dat$Modal_code <- apply(dat[judges], 1, mode_info)
dat %>%
  filter(grepl("/", Modal_code) | Modal_code != as.character(Final_decision_exclusion)) %>%
  transmute(ID, Condition = Condition_f,
            `Coder codes` = apply(across(all_of(judges)), 1, paste, collapse = ", "),
            `Modal code(s)` = Modal_code, `Final code` = Final_decision_exclusion) %>%
  kable()
```

|  ID | Condition     | Coder codes   | Modal code(s) | Final code |
|----:|:--------------|:--------------|:--------------|-----------:|
| 145 | Feeling loved | 4, 5, 1, 5, 1 | 1/5           |          4 |
| 229 | Feeling loved | 2, 3, 2, 5, 5 | 2/5           |          3 |
| 258 | Feeling loved | 2, 2, 5, 1, 1 | 1/2           |          2 |
| 587 | Feeling loved | 1, 2, 5, 1, 6 | 1             |          2 |
| 733 | Feeling loved | 6, 5, 5, 6, 3 | 5/6           |          5 |

# 5. Preregistered tests (Table 1)

``` r
outcomes <- c(Purchase_intention_index_3item   = "Purchase intention, 3 items (H1)",
              Purchase_intention_index_10items = "Purchase intention, 10 items (H1)",
              manipulation_check               = "Felt love (H2)")
run_tests <- function(d, sample_label) {
  res <- do.call(rbind, lapply(names(outcomes), function(v) compare_groups(d[[v]], d$Condition_f)))
  data.frame(
    Sample = sample_label, Measure = unname(outcomes),
    `Feeling loved M (SD)` = paste0(fmt(res$M_love), " (", fmt(res$SD_love), ")"),
    `Control M (SD)`       = paste0(fmt(res$M_control), " (", fmt(res$SD_control), ")"),
    `n (FL / C)` = paste0(res$n_love, " / ", res$n_control),
    `Levene p` = no0(res$levene_p), Test = res$test,
    `t (df)` = paste0(fmt(res$t), " (", ifelse(res$test == "Welch", fmt(res$df, 2), fmt(res$df, 0)), ")"),
    p = fmtp(res$p),
    `d [95% CI]` = paste0(fmt(res$d), " ", ci_txt(res$d_low, res$d_high)),
    check.names = FALSE)
}
table1 <- rbind(run_tests(dat, "Full (primary)"),
                run_tests(filter(dat, Valid), "Restricted (sensitivity)"))
kable(table1)
```

| Sample | Measure | Feeling loved M (SD) | Control M (SD) | n (FL / C) | Levene p | Test | t (df) | p | d \[95% CI\] |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| Full (primary) | Purchase intention, 3 items (H1) | 4.26 (1.48) | 4.19 (1.51) | 510 / 510 | .776 | Student | 0.78 (1018) | = .435 | 0.05 \[-0.07, 0.17\] |
| Full (primary) | Purchase intention, 10 items (H1) | 4.21 (1.47) | 4.15 (1.46) | 510 / 510 | .738 | Student | 0.62 (1018) | = .538 | 0.04 \[-0.08, 0.16\] |
| Full (primary) | Felt love (H2) | 5.25 (1.70) | 5.05 (1.70) | 510 / 510 | .685 | Student | 1.88 (1018) | = .061 | 0.12 \[-0.01, 0.24\] |
| Restricted (sensitivity) | Purchase intention, 3 items (H1) | 4.36 (1.45) | 4.26 (1.52) | 410 / 439 | .367 | Student | 0.96 (847) | = .338 | 0.07 \[-0.07, 0.20\] |
| Restricted (sensitivity) | Purchase intention, 10 items (H1) | 4.31 (1.43) | 4.21 (1.46) | 410 / 439 | .834 | Student | 0.98 (847) | = .327 | 0.07 \[-0.07, 0.20\] |
| Restricted (sensitivity) | Felt love (H2) | 5.56 (1.47) | 5.10 (1.69) | 410 / 439 | .013 | Welch | 4.20 (842.80) | \< .001 | 0.29 \[0.15, 0.42\] |

# 6. Descriptive statistics and correlations (Table S2.7, exploratory correlations)

``` r
measures <- c(manipulation_check = "1. Felt love",
              Purchase_intention_index_3item = "2. Intention, 3 items (Tang et al., Exp. 1a)",
              Purchase_intention_index_10items = "3. Intention, 10 items")
desc_cor <- function(d, sample_label) {
  r12 <- cor_ci(d$manipulation_check, d$Purchase_intention_index_3item)
  r13 <- cor_ci(d$manipulation_check, d$Purchase_intention_index_10items)
  r23 <- cor_ci(d$Purchase_intention_index_3item, d$Purchase_intention_index_10items)
  data.frame(
    Sample = sample_label, Measure = unname(measures),
    M  = fmt(sapply(d[names(measures)], mean)),
    SD = fmt(sapply(d[names(measures)], sd)),
    `r with 1 [95% CI]` = c("—", paste(no0(r12$r, 2), ci_txt(r12$lo, r12$hi, function(x) no0(x, 2))),
                            paste(no0(r13$r, 2), ci_txt(r13$lo, r13$hi, function(x) no0(x, 2)))),
    `r with 2` = c("", "—", no0(r23$r, 2)),
    `df` = c("", r12$df, r13$df),
    `p` = c("", fmtp(r12$p), fmtp(r13$p)),
    check.names = FALSE, row.names = NULL)
}
kable(rbind(desc_cor(dat, "Full (N = 1,020)"),
            desc_cor(filter(dat, Valid), "Restricted (N = 849)")))
```

| Sample | Measure | M | SD | r with 1 \[95% CI\] | r with 2 | df | p |
|:---|:---|:---|:---|:---|:---|:---|:---|
| Full (N = 1,020) | 1\. Felt love | 5.15 | 1.70 | — |  |  |  |
| Full (N = 1,020) | 2\. Intention, 3 items (Tang et al., Exp. 1a) | 4.22 | 1.50 | .24 \[.18, .29\] | — | 1018 | \< .001 |
| Full (N = 1,020) | 3\. Intention, 10 items | 4.18 | 1.46 | .26 \[.20, .32\] | .95 | 1018 | \< .001 |
| Restricted (N = 849) | 1\. Felt love | 5.33 | 1.61 | — |  |  |  |
| Restricted (N = 849) | 2\. Intention, 3 items (Tang et al., Exp. 1a) | 4.31 | 1.49 | .23 \[.16, .29\] | — | 847 | \< .001 |
| Restricted (N = 849) | 3\. Intention, 10 items | 4.26 | 1.45 | .25 \[.19, .31\] | .95 | 847 | \< .001 |

# 7. Item-level statistics (Table S2.8)

``` r
item_text <- c(
  "1. I am willing to spend more on ecological products because they are worth the extra cost.",
  "2. During my next shopping trip, I plan to buy food marked with the ecological product logo.",
  "3. I am NOT willing to buy ecological products. (reverse-coded)",
  "4. I intend to buy ecological products in the future.",
  "5. I always willingly buy more ecological products for my family's needs.",
  "6. I always try to seek out ecological products, even if it requires travelling out of town.",
  "7. If I had to do it again, I would buy ecological products.",
  "8. I make an effort to buy ecological products.",
  "9. I plan to buy ecological products.",
  "10. I intend to buy ecological products.")
kable(data.frame(Item = item_text,
                 M  = fmt(sapply(dat[pi_items], mean)),
                 SD = fmt(sapply(dat[pi_items], sd)),
                 row.names = NULL))
```

| Item | M | SD |
|:---|:---|:---|
| 1\. I am willing to spend more on ecological products because they are worth the extra cost. | 3.99 | 1.74 |
| 2\. During my next shopping trip, I plan to buy food marked with the ecological product logo. | 4.04 | 1.63 |
| 3\. I am NOT willing to buy ecological products. (reverse-coded) | 4.65 | 1.84 |
| 4\. I intend to buy ecological products in the future. | 4.47 | 1.59 |
| 5\. I always willingly buy more ecological products for my family’s needs. | 4.06 | 1.70 |
| 6\. I always try to seek out ecological products, even if it requires travelling out of town. | 3.45 | 1.78 |
| 7\. If I had to do it again, I would buy ecological products. | 4.33 | 1.59 |
| 8\. I make an effort to buy ecological products. | 4.13 | 1.68 |
| 9\. I plan to buy ecological products. | 4.33 | 1.63 |
| 10\. I intend to buy ecological products. | 4.38 | 1.59 |

# 8. Content of the written responses in the feeling-loved condition (S2.3, Table S2.5)

Six coders coded the source and the form of love in the valid responses
of the feeling-loved condition (*n* = 410). Agreement includes the code
“does not describe feeling loved” (`0`).

``` r
love_valid <- filter(dat, Condition_f == "Feeling loved", Valid)
nrow(love_valid)
```

    ## [1] 410

``` r
fk_source <- fleiss_kappa(love_valid[paste0("Judge_source_love", 1:6)], 0:7)
fk_form   <- fleiss_kappa(love_valid[paste0("Judge_type_love", 1:6)], 0:5)
kable(data.frame(
  Dimension = c("Source of love", "Form of love"),
  n = c(fk_source$n, fk_form$n),
  kappa = no0(c(fk_source$kappa, fk_form$kappa), 2),
  `95% CI` = c(ci_txt(fk_source$ci[1], fk_source$ci[2], function(x) no0(x, 2)),
               ci_txt(fk_form$ci[1], fk_form$ci[2], function(x) no0(x, 2))),
  check.names = FALSE))
```

| Dimension      |   n | kappa | 95% CI       |
|:---------------|----:|:------|:-------------|
| Source of love | 410 | .81   | \[.80, .83\] |
| Form of love   | 410 | .60   | \[.58, .61\] |

``` r
content_tab <- function(final, fk, labels, dim_label) {
  n <- table(factor(final, levels = fk$by_category$code))
  data.frame(Dimension = dim_label, Category = labels, n = as.vector(n),
             Percent = pct(as.vector(n), length(final)),
             kappa = no0(fk$by_category$kappa, 2))
}
kable(rbind(
  content_tab(love_valid$Final_decision_source_love, fk_source, source_labels, "Source of love"),
  content_tab(love_valid$Final_decision_type_love,   fk_form,   form_labels,   "Form of love")))
```

| Dimension      | Category                          |   n | Percent | kappa |
|:---------------|:----------------------------------|----:|:--------|:------|
| Source of love | Does not describe feeling loved   |   0 | 0.0     | .28   |
| Source of love | Romantic partner                  | 130 | 31.7    | .91   |
| Source of love | Family members                    | 144 | 35.1    | .84   |
| Source of love | People outside the family         |  17 | 4.1     | .86   |
| Source of love | Animal                            |   6 | 1.5     | 1.00  |
| Source of love | God or a spiritual being          |   4 | 1.0     | 1.00  |
| Source of love | Oneself                           |   2 | 0.5     | .64   |
| Source of love | Several sources or unspecified    | 107 | 26.1    | .69   |
| Form of love   | Does not describe feeling loved   |   0 | 0.0     | .28   |
| Form of love   | Trust and acceptance              |  24 | 5.9     | .61   |
| Form of love   | Support in needs and goals        |  81 | 19.8    | .71   |
| Form of love   | Symbolic and physical expressions | 120 | 29.3    | .68   |
| Form of love   | Shared time and presence          | 110 | 26.8    | .50   |
| Form of love   | Other                             |  75 | 18.3    | .52   |

# 9. Comparison with the manuscript

All statistics reported for Study 2 in the manuscript and the supplement
are reproduced, apart from these points:

1.  **Table 1, restricted sample, *d* values.** The *d* values are
    truncated: “0.0” for the 10-item measure should be 0.07, and “0.2”
    for felt love should be 0.29. The text gives 0.29 correctly.
2.  **Table S2.7 note, CI for the correlation of felt love with the
    10-item measure in the restricted sample.** The note gives \[.19,
    .32\]. The upper limit is .31 (.315).
3.  **Table S2.8, item 5.** The note gives *SD* = 1.69. The value is
    1.70 (1.695).
4.  **Final codes vs. modal codes (Section 4.3).** The manuscript says
    the final code is the modal code among the five coders. For three
    responses (IDs 145, 229 and 587), the final code is not a modal
    code. Two of these (IDs 145 and
    229) were tied between a valid and an invalid code, and the final
         code is a third, invalid code, presumably set by adjudication.
         Adding one sentence on how ties were resolved would make the
         procedure fully reproducible.
5.  **Table numbers in S2.2–S2.3.** The text refers to Table S2.3 for
    the coding scheme, Table S2.4 for category agreement, Table S2.5 for
    final codes and Table S2.6 for source and form. The tables are
    numbered S2.2, S2.3, S2.4 and S2.5.
6.  **ω for the 10-item measure.** SPSS prints .9649 and the `lavaan` ML
    estimate is .9652. Both round to the reported .965.

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
