# Feeling Loved and Purchase Intention for Ecological Products: Replication of Tang et al. (2025, Study 1a) — Study 1 Analysis

Invalid Date

- [Overview](#overview)
  - [Preregistered hypotheses](#preregistered-hypotheses)
  - [Notes on correspondence between the preregistration, the script,
    and the
    draft](#notes-on-correspondence-between-the-preregistration-the-script-and-the-draft)
- [Setup](#setup)
- [Data](#data)
  - [Anonymization (not run)](#anonymization-not-run)
  - [Coding of open-ended texts (not
    run)](#coding-of-open-ended-texts-not-run)
  - [Loading the anonymized data](#loading-the-anonymized-data)
  - [Loading the anonymized text
    codes](#loading-the-anonymized-text-codes)
  - [Participant flow](#participant-flow)
    - [Drop-out by condition
      (exploratory)](#drop-out-by-condition-exploratory)
  - [Completion times (appendix; not reported in the
    paper)](#completion-times-appendix-not-reported-in-the-paper)
  - [Analysis sample and
    sociodemographics](#analysis-sample-and-sociodemographics)
    - [Repeated participation](#repeated-participation)
- [Measures](#measures)
  - [Purchase intention](#purchase-intention)
  - [Self-perception (self-esteem;
    mediator)](#self-perception-self-esteem-mediator)
  - [Reliability](#reliability)
  - [Validity of the manipulation
    texts](#validity-of-the-manipulation-texts)
- [Descriptive statistics](#descriptive-statistics)
- [Confirmatory analyses](#confirmatory-analyses)
  - [H2: Manipulation check](#h2-manipulation-check)
  - [H1: Effect of the manipulation on purchase
    intention](#h1-effect-of-the-manipulation-on-purchase-intention)
  - [Replication verdict](#replication-verdict)
- [Extension: Mediation by self-perception
  (H3)](#extension-mediation-by-self-perception-h3)
  - [Effect of the manipulation on
    self-perception](#effect-of-the-manipulation-on-self-perception)
  - [Mediation model](#mediation-model)
- [Sensitivity analyses
  (preregistered)](#sensitivity-analyses-preregistered)
  - [Six-item intention scale](#six-item-intention-scale)
  - [Excluding participants with invalid
    texts](#excluding-participants-with-invalid-texts)
- [Exploratory analyses (not
  preregistered)](#exploratory-analyses-not-preregistered)
  - [Content of the manipulation
    texts](#content-of-the-manipulation-texts)
  - [Inter-rater reliability (coder 1 vs. coder
    2)](#inter-rater-reliability-coder-1-vs-coder-2)
  - [Small-telescopes test (Simonsohn,
    2015)](#small-telescopes-test-simonsohn-2015)
- [Supplementary tables](#supplementary-tables)
  - [Table A1: Means, standard deviations, and
    correlations](#table-a1-means-standard-deviations-and-correlations)
  - [Table A2: Descriptive statistics by
    condition](#table-a2-descriptive-statistics-by-condition)
- [References](#references)
- [Session information](#session-information)

# Overview

This document reproduces all analyses of Study 1, a preregistered
replication of Tang et al. (2025, Study 1a), and links each analysis to
the preregistration (Replication Recipe) and to the draft of the paper.
All numbers reported in the text below are computed in-line from the
data, so the document can be re-knitted whenever the data or the
analysis change.

**Original study.** Tang, D., Chen, J., Jiang, Q., Huang, H., & Zheng,
M. (2025). How feeling loved inspires the purchase of green foods.
*Journal of Environmental Psychology, 102*, 102523.
https://doi.org/10.1016/j.jenvp.2025.102523. Original effect: *d* =
0.53, 95% CI \[0.19, 0.87\], *N* = 141.

**Materials and preregistration:**
https://osf.io/vfy5g/overview?view_only=5e910d70bd6b4f669ff74c3a16c4d0bf

## Preregistered hypotheses

| Label | Hypothesis | Status |
|----|----|----|
| H1 | Participants in the feeling-loved condition (vs. control) will have a higher intention to purchase ecological (organic) food. | Confirmatory |
| H2 | Participants in the feeling-loved condition (vs. control) will feel more loved (manipulation check). | Confirmatory |
| H3 | The effect of the manipulation on intention is positively mediated by self-esteem (labelled “self-perception” in the draft and in the data). | Confirmatory extension |

**Criterion for a successful replication (preregistered):**
statistically significant (two-tailed, α = .05) differences in both
intention (H1) and the feeling of being loved (H2), both in the
predicted direction. **Criterion for H3:** self-esteem positively
mediates the effect of the manipulation on intention.

## Notes on correspondence between the preregistration, the script, and the draft

1.  **Primary intention measure.** As preregistered, all confirmatory
    tests use the *three-item* intention scale from Tang et al. (2025).
    The *six-item* scale (three items from Tang et al. plus three from
    Kumar et al., 2023) is used in a sensitivity analysis (see
    *Sensitivity analyses*). Note that the draft currently reports the
    six-item scale as the main measure; the Results section of the draft
    needs to be updated accordingly.
2.  **Hypothesis labels.** The preregistration labels the mediation
    hypothesis H4 in the hypothesis section and H3 in the analysis plan.
    This document uses H3.
3.  **Sensitivity analysis excluding invalid texts.** Preregistered and
    implemented under *Sensitivity analyses*, based on the coding of the
    manipulation texts (see *Validity of the manipulation texts*). A
    text is treated as *valid* when its final code matches the
    instructions of the participant’s condition: a love-related code in
    the feeling-loved condition, everyday activities in the control
    condition. All other codes are treated as invalid; this includes
    otherwise valid texts that describe the content of the other
    condition (flagged as `wrong_condition` when processing the codes).
    Every completed questionnaire must have a final code; this is
    checked before the analyses.
4.  **Sign convention.** All differences are computed as **experimental
    − control**, so positive values indicate an effect in the predicted
    direction. The original script used the formula interface of
    `t.test()` and `cohen.d()`, which computes control − experimental
    (hence the negative *t* and *d* in the draft). Absolute values are
    identical.
5.  **Minimum text length.** The informed consent states a minimum of
    130 characters, whereas the confirmation page and the writing task
    require 300 characters. The Methods section should report 300
    characters (the enforced minimum) and may note the discrepancy.

# Setup

``` r
if (!require("pacman")) install.packages("pacman")
# mediation (via MASS) and car mask some dplyr functions, so dplyr is loaded last
# and dplyr::select / car::recode are called with explicit namespaces.
pacman::p_load(car, psych, GPArotation, mediation, pwr, effsize, knitr, tidyr, dplyr)
```

``` r
# Formatting helpers for in-text reporting
fmt   <- function(x, digits = 2) formatC(x, format = "f", digits = digits)
fmt0  <- function(x, digits = 2) sub("^(-?)0\\.", "\\1.", fmt(x, digits))   # no leading zero (APA)
fmt_p <- function(p) ifelse(p < .001, "< .001", paste("=", fmt0(p, 3)))
pct   <- function(x, digits = 0) paste0(fmt(100 * x, digits), "%")
supported <- function(ok) ifelse(ok, "supported", "not supported")

# Welch two-sample t-test and Cohen's d, computed as experimental (1) minus control (0)
compare_groups <- function(var, data, conf.level = .95) {
  x1 <- data[[var]][data$treatment == 1]
  x0 <- data[[var]][data$treatment == 0]
  tt <- t.test(x1, x0)                       # two-tailed Welch test
  dd <- effsize::cohen.d(x1, x0, conf.level = conf.level)
  list(
    t = tt, d = dd,
    ok = tt$p.value < .05 & tt$statistic > 0,   # significant and in the predicted direction
    n1 = length(x1), m1 = mean(x1), sd1 = sd(x1),
    n0 = length(x0), m0 = mean(x0), sd0 = sd(x0),
    text = sprintf("*t*(%s) = %s, *p* %s, *d* = %s, 95%% CI [%s, %s]",
                   fmt(tt$parameter), fmt(tt$statistic), fmt_p(tt$p.value),
                   fmt(dd$estimate), fmt(dd$conf.int[1]), fmt(dd$conf.int[2]))
  )
}

# One row of a results table from compare_groups() output
result_row <- function(label, res) {
  data.frame(Test = label,
             `Experimental M (SD)` = sprintf("%.2f (%.2f)", res$m1, res$sd1),
             `Control M (SD)`      = sprintf("%.2f (%.2f)", res$m0, res$sd0),
             Result = res$text,
             check.names = FALSE)
}

# Reliability: Cronbach's alpha, McDonald's omega hierarchical and omega total
reliability <- function(items) {
  om <- tryCatch(suppressWarnings(psych::omega(items, plot = FALSE)),
                 error = function(e) NULL)
  if (is.null(om)) {
    return(c(alpha = psych::alpha(items)$total$raw_alpha, omega_h = NA, omega_t = NA))
  }
  c(alpha = unname(om$alpha), omega_h = unname(om$omega_h), omega_t = unname(om$omega.tot))
}

# Summary table of a mediation::mediate() object
med_summary <- function(med) {
  data.frame(
    Effect   = c("Indirect effect (ACME)", "Direct effect (ADE)", "Total effect", "Proportion mediated"),
    Estimate = c(med$d.avg, med$z.avg, med$tau.coef, med$n.avg),
    Lower    = c(med$d.avg.ci[1], med$z.avg.ci[1], med$tau.ci[1], med$n.avg.ci[1]),
    Upper    = c(med$d.avg.ci[2], med$z.avg.ci[2], med$tau.ci[2], med$n.avg.ci[2]),
    p        = c(med$d.avg.p, med$z.avg.p, med$tau.p, med$n.avg.p),
    row.names = NULL
  )
}
med_ok <- function(med) med$d.avg > 0 & med$d.avg.ci[1] > 0   # positive indirect effect, CI excludes 0
```

# Data

## Anonymization (not run)

The raw data exported from LimeSurvey were anonymized with the code
below before publication. The code is shown for transparency and is not
evaluated, because the raw data are not public. The written texts
(manipulation task) were deleted from the public data set because they
may contain identifying information.

``` r
setwd("~honza/Google Drive/CUK/uceni/FSV/Kvant praktika/Love study 2025-2026 kvant praktika/s1 data and analyses")
data <- readxl::read_excel("results-survey997596.xlsx")

# renaming the two id variables in the dataset
data <- data %>%
  rename_with(~ c("syst.id", "subject"), c("id...1", "id...10"))

# anonymizing / deleting potentially identifying information
data$subject   <- as.numeric(as.factor(data$subject)) # anonymized participant identifier
data$ipaddr    <- NULL                                # IP address
data$refurl    <- NULL                                # redirecting URL
data$treatment <- NULL                                # written text, experimental condition
data$control   <- NULL                                # written text, control condition

saveRDS(data, "data_s1_anonymized.RDS")
writexl::write_xlsx(data, "data_s1_anonymized.xlsx")
```

## Coding of open-ended texts (not run)

The texts written in the manipulation task were coded independently by
two coders (coder 1 and coder 2) into content categories; disagreements
were resolved by a third coder (coder 3). The code below creates the
final code and translates the Czech category labels into English. The
texts themselves are removed before the codes are published in
`text_codes.RDS`, which is loaded below (*Loading the anonymized text
codes*).

``` r
library(readxl)
setwd("~honza/Google Drive/CUK/uceni/FSV/Kvant praktika/Love study 2025-2026 kvant praktika/s1 data and analyses")
dat <- readxl::read_excel("text_coding2.xlsx", na = c("", "NA"))

dat <- dat[order(dat$id),]

# Create final text_code from three coders
# - coders 1 and 2 agree (incl. both NA)  -> coder1
# - coders 1 and 2 disagree               -> coder3

# Treat two NAs as agreement; one NA vs. a code as disagreement
agree <- (is.na(dat$coder1) & is.na(dat$coder2)) |
         (!is.na(dat$coder1) & !is.na(dat$coder2) & dat$coder1 == dat$coder2)

dat$final_code <- ifelse(agree, dat$coder1, dat$coder3)

table(dat$final_code)

# Check: disagreements that coder3 did not resolve
# Add row-level flags so unresolved cases can be filtered or reviewed later.
dat$unresolved <- as.integer(!agree & is.na(dat$coder3))

# Recode the Czech category labels (as entered by coders) into English labels.
text_recode_en <- function(x) {
  dplyr::recode(
    x,
    "Irelevantní, nesmyslné" = "irrelevant_or_nonsense_(not_valid)",
    "Každodenní aktivity (legitimní výrok)" = "everyday_activities_(valid)",
    "Neví, cítí se špatně, necítí se milován(a)" = "uncertain_or_unloved_(not_valid)",
    "Obecné výroky (ne o sobe)" = "general_statements_(not_valid)",
    "Prátelská/ od bližního" = "friend's_or_other_person's_love_(valid)",
    "Rodinná" = "family_love_(valid)",
    "Romantická" = "romantic_love_(valid)",
    "Více typů lásky" = "multiple_kinds_of_love_(valid)",
    .default = NA_character_
  )
}

dat <- dat %>%
  dplyr::mutate(across(c(coder1, coder2, coder3, final_code), text_recode_en, .names = "{.col}_en"))

dat$text <- NULL  # remove the original texts from the anonymized data

saveRDS(dat, "text_codes.RDS")
writexl::write_xlsx(dat, "text_codes.xlsx")
```

## Loading the anonymized data

``` r
read_rds_from_url <- function(url) {
  tmp <- tempfile(fileext = ".rds")
  download.file(url, destfile = tmp, mode = "wb", quiet = TRUE)
  readRDS(tmp)
}

data <- read_rds_from_url(params$data_url)

# LimeSurvey stores answer codes as "A1", "A2", ...; convert them to numbers
data <- data %>%
  mutate(across(c(gender:educ, random, manipCheck, `intention[SQ001]`:repeated),
                ~ as.numeric(gsub("A", "", .x))))
```

## Loading the anonymized text codes

The codes are merged with the survey data under *Validity of the
manipulation texts*, after the analysis sample has been defined.

``` r
codes <- read_rds_from_url(params$codes_url)
dplyr::glimpse(codes)
```

    #> Rows: 4,824
    #> Columns: 11
    #> $ id            <dbl> 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 1…
    #> $ condition     <chr> NA, "treatment", "treatment", NA, "control", "control", …
    #> $ coder1        <chr> NA, "Rodinná", "Rodinná", NA, "Každodenní aktivity (legi…
    #> $ coder2        <chr> NA, "Rodinná", "Rodinná", NA, "Každodenní aktivity (legi…
    #> $ coder3        <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
    #> $ final_code    <chr> NA, "Rodinná", "Rodinná", NA, "Každodenní aktivity (legi…
    #> $ unresolved    <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
    #> $ coder1_en     <chr> NA, "family_love_(valid)", "family_love_(valid)", NA, "e…
    #> $ coder2_en     <chr> NA, "family_love_(valid)", "family_love_(valid)", NA, "e…
    #> $ coder3_en     <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
    #> $ final_code_en <chr> NA, "family_love_(valid)", "family_love_(valid)", NA, "e…

## Participant flow

Participants who were screened out before randomization (variable
`random` missing) are excluded from the flow. Participants who were
randomized but did not reach the last page (page 17) are classified as
drop-outs. As preregistered, the main analyses include all participants
who completed the study.

``` r
data$status <- "screenout"
data$status <- ifelse(!is.na(data$random), "dropout", data$status)
data$status <- ifelse(data$lastpage == 17, "complete", data$status)

n_screenout <- sum(data$status == "screenout")
data <- data %>% filter(status != "screenout")   # from here on: randomized participants only

flow <- data %>%
  count(status) %>%
  mutate(prop = n / sum(n))
flow
```

    #> # A tibble: 2 × 3
    #>   status       n  prop
    #>   <chr>    <int> <dbl>
    #> 1 complete  1000 0.596
    #> 2 dropout    679 0.404

``` r
n_entered  <- nrow(data)
n_complete <- sum(data$status == "complete")
n_dropout  <- sum(data$status == "dropout")
```

In total, 1679 participants entered the study (after excluding NA
participants who were screened out before randomization) and 1000 (60%)
completed it; 679 participants dropped out after randomization.

### Drop-out by condition (exploratory)

The draft notes that most drop-outs happened in the writing task.
Because condition was assigned before the writing task, differential
drop-out across conditions could bias the comparison of completers. This
analysis was not preregistered. It includes only randomized
participants, so screen-outs (who were never assigned to a condition)
are not part of it.

``` r
data <- data %>%
  mutate(condition = factor(random, levels = c(1, 2), labels = c("control", "experimental")))

dropout_tab <- data %>%
  group_by(Condition = condition) %>%
  summarise(Randomized = n(),
            Completed  = sum(status == "complete"),
            `Dropped out` = sum(status == "dropout"),
            `Drop-out rate` = pct(mean(status == "dropout"), 1),
            .groups = "drop")
kable(dropout_tab, caption = "Drop-out by condition (randomized participants only; screen-outs excluded)")
```

| Condition    | Randomized | Completed | Dropped out | Drop-out rate |
|:-------------|-----------:|----------:|------------:|:--------------|
| control      |        845 |       516 |         329 | 38.9%         |
| experimental |        834 |       484 |         350 | 42.0%         |

Drop-out by condition (randomized participants only; screen-outs
excluded)

``` r
attrition_tab  <- table(Condition = data$condition, Status = data$status)
attrition_test <- chisq.test(attrition_tab)
attrition_test
```

    #> 
    #>  Pearson's Chi-squared test with Yates' continuity correction
    #> 
    #> data:  attrition_tab
    #> X-squared = 1.4781, df = 1, p-value = 0.2241

``` r
# Last page reached by participants who dropped out
data %>%
  filter(status == "dropout") %>%
  count(condition, lastpage) %>%
  pivot_wider(names_from = condition, values_from = n, values_fill = 0) %>%
  arrange(lastpage) %>%
  kable(caption = "Last page reached by drop-outs, by condition")
```

| lastpage | control | experimental |
|---------:|--------:|-------------:|
|       10 |       0 |          347 |
|       11 |     326 |            0 |
|       12 |       0 |            1 |
|       14 |       1 |            1 |
|       16 |       2 |            1 |

Last page reached by drop-outs, by condition

Of the participants randomized to the control condition, 329 (38.9%)
dropped out, compared with 350 (42.0%) in the experimental condition,
χ²(1) = 1.48, *p* = .224.

## Completion times (appendix; not reported in the paper)

``` r
time_vars <- c("interviewtime", "consentTime", "confirmTextTime", "productInfoTime",
               "treatmentTime", "controlTime", "manipCheckTime", "intentionTime",
               "selfperceptionTime", "repeatedTime", "debriefingTime")

time_table <- data %>%
  filter(lastpage == 17) %>%
  summarise(across(all_of(time_vars),
                   list(M      = ~ mean(.x, na.rm = TRUE),
                        SD     = ~ sd(.x, na.rm = TRUE),
                        Min    = ~ fivenum(.x, na.rm = TRUE)[1],
                        Q1     = ~ fivenum(.x, na.rm = TRUE)[2],
                        Median = ~ fivenum(.x, na.rm = TRUE)[3],
                        Q3     = ~ fivenum(.x, na.rm = TRUE)[4],
                        Max    = ~ fivenum(.x, na.rm = TRUE)[5]),
                   .names = "{.col}___{.fn}")) %>%
  pivot_longer(everything(), names_to = c("Variable", ".value"),
               names_pattern = "(.*)___(.*)")

kable(time_table, digits = 2,
      caption = "Completion times in seconds (completers only). interviewtime = total completion time.")
```

| Variable           |      M |      SD |   Min |     Q1 | Median |     Q3 |      Max |
|:-------------------|-------:|--------:|------:|-------:|-------:|-------:|---------:|
| interviewtime      | 748.70 | 2140.60 | 66.61 | 390.19 | 536.42 | 747.62 | 65612.56 |
| consentTime        |  31.41 |   99.38 |  1.66 |   4.70 |   9.76 |  28.83 |  1967.11 |
| confirmTextTime    |  10.27 |   31.01 |  1.80 |   5.21 |   7.76 |  10.68 |   853.69 |
| productInfoTime    |  71.58 |  181.38 |  2.03 |  26.68 |  39.33 |  57.64 |  2849.18 |
| treatmentTime      | 345.55 |  269.31 |  4.61 | 194.88 | 282.35 | 418.06 |  2636.81 |
| controlTime        | 507.05 | 2898.10 |  3.76 | 193.40 | 283.88 | 406.52 | 65346.78 |
| manipCheckTime     |  20.43 |  229.04 |  2.15 |   7.42 |   9.70 |  12.91 |  7185.28 |
| intentionTime      |  47.89 |  140.22 |  4.20 |  22.26 |  32.19 |  43.14 |  2776.99 |
| selfperceptionTime |  66.08 |  137.06 |  5.30 |  36.58 |  50.95 |  69.88 |  3604.74 |
| repeatedTime       |   6.90 |   13.21 |  1.99 |   4.47 |   5.60 |   7.07 |   395.83 |
| debriefingTime     |  10.94 |   14.19 |  1.35 |   4.32 |   6.48 |  12.68 |   202.39 |

Completion times in seconds (completers only). interviewtime = total
completion time.

## Analysis sample and sociodemographics

``` r
data.complete <- data %>% filter(lastpage == 17)

# Age was recorded in a drop-down menu starting at 15 (code 1 = 15 years)
data.complete$ageNum2 <- data.complete$ageNum + 14

data.complete$gender <- car::recode(data.complete$gender,
                                    "1 = 'male'; 2 = 'female'; 3 = 'non-binary/gender-diverse'")
data.complete$educ   <- car::recode(data.complete$educ,
                                    "1 = 'elementary/secondary without diploma'; 3 = 'secondary with diploma'; 4 = 'university'")

gender_tab <- prop.table(table(data.complete$gender))
educ_tab   <- prop.table(table(data.complete$educ))

data.frame(Gender = names(gender_tab), Proportion = round(as.numeric(gender_tab), 3))
```

    #>                      Gender Proportion
    #> 1                    female      0.516
    #> 2                      male      0.481
    #> 3 non-binary/gender-diverse      0.003

``` r
data.frame(Education = names(educ_tab), Proportion = round(as.numeric(educ_tab), 3))
```

    #>                              Education Proportion
    #> 1 elementary/secondary without diploma      0.303
    #> 2               secondary with diploma      0.429
    #> 3                           university      0.268

The analysis sample comprised 1000 participants (age *M* = 46.32, *SD* =
14.94; 52% women). In terms of education, 30% had elementary education
or secondary education without a high school diploma, 43% had a high
school diploma, and 27% had university education.

### Repeated participation

Participants were asked whether they had participated in the study
before (1 = yes, 2 = no). No exclusion based on this item was
preregistered; the answers are reported for information only.

``` r
table(Repeated = data.complete$repeated, useNA = "ifany")
```

    #> Repeated
    #>   1   2 
    #>  32 968

# Measures

``` r
data.complete <- data.complete %>%
  rename_with(~ paste0("intention_", 1:6), `intention[SQ001]`:`intention[SQ006]`) %>%
  rename_with(~ paste0("selfperc_", 1:10), `selfperception[SQ001]`:`selfperception[SQ010]`)

# Treatment indicator: 1 = feeling-loved (experimental) condition, 0 = control condition
data.complete$treatment <- data.complete$random - 1
```

## Purchase intention

Items 1–3 were adapted from Tang et al. (2025), items 4–6 from Kumar et
al. (2023). Item 3 is reverse-coded. Scale scores are item means (1–7).
As preregistered, the **three-item scale (`intention`) is the primary
measure**; the six-item scale (`intention6`) is used in a sensitivity
analysis.

``` r
intention_items <- data.complete %>% dplyr::select(intention_1:intention_6)
intention_items[, 3] <- 8 - intention_items[, 3]               # reverse-code item 3

data.complete$intention  <- rowMeans(intention_items[, 1:3])   # three-item scale (primary; Tang et al., 2025)
data.complete$intention6 <- rowMeans(intention_items)          # six-item scale (sensitivity analysis)

rel_int3 <- reliability(intention_items[, 1:3])
rel_int6 <- reliability(intention_items)
```

## Self-perception (self-esteem; mediator)

Rosenberg’s self-esteem scale (Czech adaptation by Blatný & Osecká,
1994), answered on a seven-point scale as in Tang et al. (2025). Items
2, 5, 6, 8, and 9 are reverse-coded.

``` r
selfperc <- data.complete %>% dplyr::select(selfperc_1:selfperc_10)
rev_items <- c(2, 5, 6, 8, 9)
selfperc[, rev_items] <- 8 - selfperc[, rev_items]

data.complete$selfperc <- rowMeans(selfperc)

rel_self <- reliability(selfperc)
```

## Reliability

``` r
rel_table <- rbind(
  "Intention (3 items, Tang et al.; primary)" = rel_int3,
  "Intention (6 items; sensitivity analysis)" = rel_int6,
  "Self-perception (10 items)"                = rel_self
)
kable(rel_table, digits = 2,
      col.names = c("Cronbach's alpha", "Omega hierarchical", "Omega total"))
```

|  | Cronbach’s alpha | Omega hierarchical | Omega total |
|:---|---:|---:|---:|
| Intention (3 items, Tang et al.; primary) | 0.65 | 0.06 | 0.73 |
| Intention (6 items; sensitivity analysis) | 0.89 | 0.90 | 0.92 |
| Self-perception (10 items) | 0.86 | 0.63 | 0.90 |

The three-item intention scale (primary measure) had α = .65,
ω<sub>t</sub> = .73. Extending the scale with the three items from Kumar
et al. (2023) increased reliability, α = .89, ω<sub>h</sub> = .90,
ω<sub>t</sub> = .92. The self-perception scale had α = .86,
ω<sub>h</sub> = .63, ω<sub>t</sub> = .90.

Note: omega hierarchical is not well defined for a three-item scale
(psych::omega extracts three group factors by default), so only alpha
and omega total are interpreted for the three-item version.

## Validity of the manipulation texts

The text codes are merged with the analysis sample. Each code is
classified as valid or invalid relative to the participant’s condition:
in the feeling-loved condition, the four love-related codes are valid;
in the control condition, only “everyday activities” is valid. Validity
is derived separately for coder 1, coder 2, and the final code.

``` r
love_codes     <- c("family_love_(valid)", "romantic_love_(valid)",
                    "friend's_or_other_person's_love_(valid)", "multiple_kinds_of_love_(valid)")
everyday_codes <- "everyday_activities_(valid)"

# TRUE = valid for the participant's condition, FALSE = invalid, NA = no code
is_valid <- function(code, treatment) {
  dplyr::case_when(is.na(code)    ~ NA,
                   treatment == 1 ~ code %in% love_codes,
                   treatment == 0 ~ code %in% everyday_codes)
}

codes_small <- codes %>%
  mutate(join_id = as.character(.data[[params$codes_id]])) %>%
  dplyr::select(join_id, coder1_en, coder2_en, coder3_en, final_code_en,
                unresolved, any_of(c("wrong_condition")))
stopifnot(!anyDuplicated(codes_small$join_id))

data.complete <- data.complete %>%
  mutate(join_id = as.character(.data[[params$data_id]]),
         coded   = join_id %in% codes_small$join_id) %>%
  left_join(codes_small, by = "join_id") %>%
  mutate(valid_c1    = is_valid(coder1_en, treatment),
         valid_c2    = is_valid(coder2_en, treatment),
         valid_final = is_valid(final_code_en, treatment))

n_coded   <- sum(data.complete$coded)
n_uncoded <- sum(!data.complete$coded)
```

1000 of 1000 completers were matched to a text code (0 without a code).

**Checks of the codes.** The chunk below runs two checks on all
randomized participants (completers and drop-outs) and stops the knit if
either fails. Matching is by identifier only; condition labels are not
compared across files.

1.  *Completeness.* Every completed questionnaire must have a final
    code; missing codes are allowed only for unfinished questionnaires.
2.  *Condition flag* (only if `text_codes.RDS` contains
    `wrong_condition`). The flag created when processing the codes must
    agree with the condition-specific validity derived here.

``` r
code_status <- data %>%
  mutate(join_id = as.character(.data[[params$data_id]])) %>%
  left_join(codes_small, by = "join_id") %>%
  mutate(Code = case_when(!join_id %in% codes_small$join_id ~ "not in codes file",
                          is.na(final_code_en)              ~ "in codes file, no final code",
                          TRUE                              ~ "final code present"))

# Codes-file identifiers without a randomized participant (e.g., screen-outs)
n_orphan <- length(setdiff(codes_small$join_id, code_status$join_id))

# 1. Completeness check
# Only matched IDs should be checked for missing final codes; unmatched IDs are reviewed manually.
completeness_tab <- table(Status = code_status$status, Code = code_status$Code)
completeness_tab
```

    #>           Code
    #> Status     final code present in codes file, no final code
    #>   complete               1000                            0
    #>   dropout                   0                          679

``` r
missing_complete <- code_status %>%
  filter(status == "complete",
         join_id %in% codes_small$join_id,
         Code != "final code present")
if (nrow(missing_complete) > 0) {
  stop(sprintf("Completeness check failed: %d completed questionnaire(s) without a final code (%s: %s).",
               nrow(missing_complete), params$data_id,
               paste(head(missing_complete$join_id, 20), collapse = ", ")))
}

# 2. Condition flag check
if ("wrong_condition" %in% names(data.complete)) {
  derived_flag <- grepl("_\\(valid\\)$", data.complete$final_code_en) & data.complete$valid_final %in% FALSE
  flag_mismatch <- data.complete %>%
    filter(coded, !is.na(wrong_condition), (wrong_condition == 1) != derived_flag)
  if (nrow(flag_mismatch) > 0) {
    stop(sprintf("Condition flag check failed: %d case(s) where wrong_condition disagrees with the derived validity (%s: %s).",
                 nrow(flag_mismatch), params$data_id, paste(head(flag_mismatch$join_id, 20), collapse = ", ")))
  }
}
```

All checks passed. All 1000 completed questionnaires have a final code;
of the 679 drop-outs, 0 have a final code. 3145 identifier(s) in the
codes file do not belong to a randomized participant.

# Descriptive statistics

Means and standard deviations by condition are reported in Table A2 and
means, standard deviations, and correlations in Table A1 (see
*Supplementary tables*).

# Confirmatory analyses

All tests are two-tailed Welch *t*-tests for independent samples, as
preregistered. Cohen’s *d* uses the pooled *SD*; positive values
indicate higher scores in the feeling-loved condition.

## H2: Manipulation check

``` r
h2 <- compare_groups("manipCheck", data.complete)
h2$t
```

    #> 
    #>  Welch Two Sample t-test
    #> 
    #> data:  x1 and x0
    #> t = 5.0612, df = 997.73, p-value = 4.958e-07
    #> alternative hypothesis: true difference in means is not equal to 0
    #> 95 percent confidence interval:
    #>  0.3427470 0.7768315
    #> sample estimates:
    #> mean of x mean of y 
    #>  5.557851  4.998062

``` r
h2$d
```

    #> 
    #> Cohen's d
    #> 
    #> d estimate: 0.3197748 (small)
    #> 95 percent confidence interval:
    #>     lower     upper 
    #> 0.1948114 0.4447382

Participants in the feeling-loved condition reported feeling more loved
(*M* = 5.56, *SD* = 1.71) than participants in the control condition
(*M* = 5.00, *SD* = 1.79), *t*(997.73) = 5.06, *p* \< .001, *d* = 0.32,
95% CI \[0.19, 0.44\]. **H2 was supported.**

## H1: Effect of the manipulation on purchase intention

``` r
h1 <- compare_groups("intention", data.complete)   # three-item scale (preregistered primary measure)
h1$t
```

    #> 
    #>  Welch Two Sample t-test
    #> 
    #> data:  x1 and x0
    #> t = 0.25177, df = 990.92, p-value = 0.8013
    #> alternative hypothesis: true difference in means is not equal to 0
    #> 95 percent confidence interval:
    #>  -0.1456719  0.1885531
    #> sample estimates:
    #> mean of x mean of y 
    #>  4.095730  4.074289

``` r
h1$d
```

    #> 
    #> Cohen's d
    #> 
    #> d estimate: 0.01594207 (negligible)
    #> 95 percent confidence interval:
    #>      lower      upper 
    #> -0.1082330  0.1401172

Purchase intention did not differ substantially between the
feeling-loved condition (*M* = 4.10, *SD* = 1.36) and the control
condition (*M* = 4.07, *SD* = 1.33), *t*(990.92) = 0.25, *p* = .801, *d*
= 0.02, 95% CI \[-0.11, 0.14\]. **H1 was not supported.**

## Replication verdict

``` r
data.frame(
  Criterion = c("H1 (intention, 3 items)", "H2 (manipulation check)"),
  d   = round(c(h1$d$estimate, h2$d$estimate), 2),
  p   = signif(c(h1$t$p.value, h2$t$p.value), 3),
  Met = c(h1$ok, h2$ok)
)
```

    #>                 Criterion    d        p   Met
    #> 1 H1 (intention, 3 items) 0.02 8.01e-01 FALSE
    #> 2 H2 (manipulation check) 0.32 4.96e-07  TRUE

Following the preregistered criterion (significant effects in the
predicted direction for both H1 and H2), the replication was **not
successful**. The replication effect on intention (*d* = 0.02) was
considerably smaller than the original effect (*d* = 0.53, 95% CI
\[0.19, 0.87\]).

# Extension: Mediation by self-perception (H3)

## Effect of the manipulation on self-perception

``` r
h3a <- compare_groups("selfperc", data.complete)
h3a$t
```

    #> 
    #>  Welch Two Sample t-test
    #> 
    #> data:  x1 and x0
    #> t = -1.2386, df = 992.34, p-value = 0.2158
    #> alternative hypothesis: true difference in means is not equal to 0
    #> 95 percent confidence interval:
    #>  -0.21851559  0.04941123
    #> sample estimates:
    #> mean of x mean of y 
    #>  4.708471  4.793023

``` r
h3a$d
```

    #> 
    #> Cohen's d
    #> 
    #> d estimate: -0.07840219 (negligible)
    #> 95 percent confidence interval:
    #>       lower       upper 
    #> -0.20262295  0.04581857

``` r
r_sm <- cor.test(data.complete$selfperc, data.complete$intention)
r_sm
```

    #> 
    #>  Pearson's product-moment correlation
    #> 
    #> data:  data.complete$selfperc and data.complete$intention
    #> t = 2.2926, df = 998, p-value = 0.02208
    #> alternative hypothesis: true correlation is not equal to 0
    #> 95 percent confidence interval:
    #>  0.01043462 0.13377371
    #> sample estimates:
    #>        cor 
    #> 0.07238089

Self-perception did not differ substantially between the feeling-loved
condition (*M* = 4.71, *SD* = 1.08) and the control condition (*M* =
4.79, *SD* = 1.07), *t*(992.34) = -1.24, *p* = .216, *d* = -0.08, 95% CI
\[-0.20, 0.05\]. Self-perception and intention were correlated at *r* =
.07, *p* = .022.

## Mediation model

The mediation model regresses self-perception on treatment (mediator
model) and intention (three-item scale) on treatment and self-perception
(outcome model). The indirect effect (ACME) is tested with a
nonparametric percentile bootstrap (5000 resamples) using the
*mediation* package.

``` r
set.seed(params$seed)
model_m <- lm(selfperc ~ treatment, data = data.complete)
model_y <- lm(intention ~ treatment + selfperc, data = data.complete)

med <- mediation::mediate(model_m, model_y,
                          treat = "treatment", mediator = "selfperc",
                          boot = TRUE, boot.ci.type = "perc",
                          sims = params$n_boot)
summary(model_m)
```

    #> 
    #> Call:
    #> lm(formula = selfperc ~ treatment, data = data.complete)
    #> 
    #> Residuals:
    #>     Min      1Q  Median      3Q     Max 
    #> -3.4930 -0.7085 -0.0930  0.7915  2.2915 
    #> 
    #> Coefficients:
    #>             Estimate Std. Error t value Pr(>|t|)    
    #> (Intercept)  4.79302    0.04748 100.957   <2e-16 ***
    #> treatment   -0.08455    0.06824  -1.239    0.216    
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> 
    #> Residual standard error: 1.078 on 998 degrees of freedom
    #> Multiple R-squared:  0.001536,   Adjusted R-squared:  0.0005354 
    #> F-statistic: 1.535 on 1 and 998 DF,  p-value: 0.2156

``` r
summary(model_y)
```

    #> 
    #> Call:
    #> lm(formula = intention ~ treatment + selfperc, data = data.complete)
    #> 
    #> Residuals:
    #>     Min      1Q  Median      3Q     Max 
    #> -3.2946 -0.7571 -0.0023  0.9201  3.0139 
    #> 
    #> Coefficients:
    #>             Estimate Std. Error t value Pr(>|t|)    
    #> (Intercept)  3.63943    0.19783  18.397   <2e-16 ***
    #> treatment    0.02911    0.08499   0.343   0.7320    
    #> selfperc     0.09073    0.03939   2.303   0.0215 *  
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> 
    #> Residual standard error: 1.342 on 997 degrees of freedom
    #> Multiple R-squared:  0.005356,   Adjusted R-squared:  0.003361 
    #> F-statistic: 2.684 on 2 and 997 DF,  p-value: 0.06876

``` r
summary(med)
```

    #> 
    #> Causal Mediation Analysis 
    #> 
    #> Nonparametric Bootstrap Confidence Intervals with the Percentile Method
    #> 
    #>                  Estimate 95% CI Lower 95% CI Upper p-value
    #> ACME           -0.0076712   -0.0264833    0.0045530  0.2424
    #> ADE             0.0291118   -0.1368822    0.1919641  0.7128
    #> Total Effect    0.0214406   -0.1443095    0.1843410  0.7764
    #> Prop. Mediated -0.3577883   -1.4762658    1.3164411  0.8380
    #> 
    #> Sample Size Used: 1000 
    #> 
    #> 
    #> Simulations: 5000

``` r
kable(med_summary(med), digits = 3,
      col.names = c("Effect", "Estimate", "95% CI lower", "95% CI upper", "p"),
      caption = sprintf("Mediation by self-perception, intention (3 items) as outcome (N = %d, %d bootstrap resamples)",
                        med$nobs, params$n_boot))
```

| Effect                 | Estimate | 95% CI lower | 95% CI upper |     p |
|:-----------------------|---------:|-------------:|-------------:|------:|
| Indirect effect (ACME) |   -0.008 |       -0.026 |        0.005 | 0.242 |
| Direct effect (ADE)    |    0.029 |       -0.137 |        0.192 | 0.713 |
| Total effect           |    0.021 |       -0.144 |        0.184 | 0.776 |
| Proportion mediated    |   -0.358 |       -1.476 |        1.316 | 0.838 |

Mediation by self-perception, intention (3 items) as outcome (N = 1000,
5000 bootstrap resamples)

``` r
h3_ok  <- med_ok(med)
coef_y <- summary(model_y)$coefficients
```

The indirect effect of the manipulation on intention through
self-perception was *b* = -0.008, 95% CI \[-0.026, 0.005\].
Self-perception was associated with intention when controlling for
treatment, *b* = 0.091, *SE* = 0.039, *p* = .021. **H3 was not
supported.** Note that the proportion mediated is not interpretable when
the total effect is close to zero.

# Sensitivity analyses (preregistered)

## Six-item intention scale

As preregistered, H1 and the mediation model (H3) are re-estimated with
the six-item intention scale.

``` r
h1_6 <- compare_groups("intention6", data.complete)
h1_6$t
```

    #> 
    #>  Welch Two Sample t-test
    #> 
    #> data:  x1 and x0
    #> t = 0.11721, df = 992.29, p-value = 0.9067
    #> alternative hypothesis: true difference in means is not equal to 0
    #> 95 percent confidence interval:
    #>  -0.1609896  0.1814425
    #> sample estimates:
    #> mean of x mean of y 
    #>  4.045110  4.034884

``` r
h1_6$d
```

    #> 
    #> Cohen's d
    #> 
    #> d estimate: 0.007419511 (negligible)
    #> 95 percent confidence interval:
    #>      lower      upper 
    #> -0.1167540  0.1315931

``` r
set.seed(params$seed)
model_y6 <- lm(intention6 ~ treatment + selfperc, data = data.complete)
med6 <- mediation::mediate(model_m, model_y6,
                           treat = "treatment", mediator = "selfperc",
                           boot = TRUE, boot.ci.type = "perc",
                           sims = params$n_boot)
summary(med6)
```

    #> 
    #> Causal Mediation Analysis 
    #> 
    #> Nonparametric Bootstrap Confidence Intervals with the Percentile Method
    #> 
    #>                  Estimate 95% CI Lower 95% CI Upper p-value
    #> ACME           -0.0064121   -0.0242598    0.0041204  0.2780
    #> ADE             0.0166386   -0.1557150    0.1880005  0.8320
    #> Total Effect    0.0102265   -0.1631021    0.1810370  0.8912
    #> Prop. Mediated -0.6270124   -1.3045225    1.3477228  0.9292
    #> 
    #> Sample Size Used: 1000 
    #> 
    #> 
    #> Simulations: 5000

``` r
kable(med_summary(med6), digits = 3,
      col.names = c("Effect", "Estimate", "95% CI lower", "95% CI upper", "p"),
      caption = "Mediation by self-perception, intention (6 items) as outcome")
```

| Effect                 | Estimate | 95% CI lower | 95% CI upper |     p |
|:-----------------------|---------:|-------------:|-------------:|------:|
| Indirect effect (ACME) |   -0.006 |       -0.024 |        0.004 | 0.278 |
| Direct effect (ADE)    |    0.017 |       -0.156 |        0.188 | 0.832 |
| Total effect           |    0.010 |       -0.163 |        0.181 | 0.891 |
| Proportion mediated    |   -0.627 |       -1.305 |        1.348 | 0.929 |

Mediation by self-perception, intention (6 items) as outcome

With the six-item scale, intention again did not differ substantially
between the feeling-loved condition (*M* = 4.05, *SD* = 1.39) and the
control condition (*M* = 4.03, *SD* = 1.37), *t*(992.29) = 0.12, *p* =
.907, *d* = 0.01, 95% CI \[-0.12, 0.13\], so H1 was not supported. The
indirect effect through self-perception was *b* = -0.006, 95% CI
\[-0.024, 0.004\], so H3 was not supported.

## Excluding participants with invalid texts

The preregistration specifies a sensitivity analysis that replicates all
analyses after excluding participants whose written text is unrelated to
a situation of being loved (feeling-loved condition) or to everyday
activities (control condition). Here, participants are retained only if
the final code of their text is valid for their condition (see *Validity
of the manipulation texts*); participants with an invalid final code, or
without a final code (uncoded or unresolved disagreement), are excluded.

``` r
data.valid <- data.complete %>% filter(valid_final %in% TRUE)

exclusion_tab <- data.complete %>%
  mutate(Status = case_when(valid_final %in% TRUE  ~ "valid (retained)",
                            valid_final %in% FALSE ~ "invalid (excluded)",
                            TRUE                   ~ "no final code (excluded)")) %>%
  count(Condition = condition, Status) %>%
  group_by(Condition) %>%
  mutate(Percent = pct(n / sum(n), 1)) %>%
  ungroup()
kable(exclusion_tab, caption = "Participants retained and excluded, by condition")
```

| Condition    | Status             |   n | Percent |
|:-------------|:-------------------|----:|:--------|
| control      | invalid (excluded) |  43 | 8.3%    |
| control      | valid (retained)   | 473 | 91.7%   |
| experimental | invalid (excluded) |  67 | 13.8%   |
| experimental | valid (retained)   | 417 | 86.2%   |

Participants retained and excluded, by condition

``` r
n_valid    <- nrow(data.valid)
n_excluded <- nrow(data.complete) - n_valid
```

``` r
v_h2  <- compare_groups("manipCheck", data.valid)
v_h1  <- compare_groups("intention",  data.valid)
v_h16 <- compare_groups("intention6", data.valid)
v_h3a <- compare_groups("selfperc",   data.valid)

kable(rbind(result_row("H2: manipulation check",  v_h2),
            result_row("H1: intention (3 items)", v_h1),
            result_row("H1: intention (6 items)", v_h16),
            result_row("Self-perception",          v_h3a)),
      caption = sprintf("Group comparisons after excluding invalid texts (N = %d)", n_valid))
```

| Test | Experimental M (SD) | Control M (SD) | Result |
|:---|:---|:---|:---|
| H2: manipulation check | 5.79 (1.53) | 5.05 (1.77) | *t*(887.71) = 6.62, *p* \< .001, *d* = 0.44, 95% CI \[0.31, 0.57\] |
| H1: intention (3 items) | 4.14 (1.35) | 4.09 (1.34) | *t*(871.84) = 0.61, *p* = .542, *d* = 0.04, 95% CI \[-0.09, 0.17\] |
| H1: intention (6 items) | 4.09 (1.37) | 4.04 (1.36) | *t*(871.98) = 0.55, *p* = .581, *d* = 0.04, 95% CI \[-0.09, 0.17\] |
| Self-perception | 4.78 (1.10) | 4.83 (1.07) | *t*(868.58) = -0.63, *p* = .526, *d* = -0.04, 95% CI \[-0.17, 0.09\] |

Group comparisons after excluding invalid texts (N = 890)

``` r
set.seed(params$seed)
v_model_m <- lm(selfperc ~ treatment, data = data.valid)
v_model_y <- lm(intention ~ treatment + selfperc, data = data.valid)
v_med <- mediation::mediate(v_model_m, v_model_y,
                            treat = "treatment", mediator = "selfperc",
                            boot = TRUE, boot.ci.type = "perc",
                            sims = params$n_boot)
summary(v_med)
```

    #> 
    #> Causal Mediation Analysis 
    #> 
    #> Nonparametric Bootstrap Confidence Intervals with the Percentile Method
    #> 
    #>                  Estimate 95% CI Lower 95% CI Upper p-value
    #> ACME           -0.0044438   -0.0231824    0.0105250  0.5364
    #> ADE             0.0595338   -0.1168263    0.2349551  0.5168
    #> Total Effect    0.0550900   -0.1207006    0.2312927  0.5496
    #> Prop. Mediated -0.0806641   -1.3041046    0.9876924  0.8220
    #> 
    #> Sample Size Used: 890 
    #> 
    #> 
    #> Simulations: 5000

``` r
kable(med_summary(v_med), digits = 3,
      col.names = c("Effect", "Estimate", "95% CI lower", "95% CI upper", "p"),
      caption = "Mediation by self-perception after excluding invalid texts, intention (3 items) as outcome")
```

| Effect                 | Estimate | 95% CI lower | 95% CI upper |     p |
|:-----------------------|---------:|-------------:|-------------:|------:|
| Indirect effect (ACME) |   -0.004 |       -0.023 |        0.011 | 0.536 |
| Direct effect (ADE)    |    0.060 |       -0.117 |        0.235 | 0.517 |
| Total effect           |    0.055 |       -0.121 |        0.231 | 0.550 |
| Proportion mediated    |   -0.081 |       -1.304 |        0.988 | 0.822 |

Mediation by self-perception after excluding invalid texts, intention (3
items) as outcome

After excluding 110 participants with invalid or missing text codes (*N*
= 890), the manipulation check showed a significant difference in the
predicted direction, *t*(887.71) = 6.62, *p* \< .001, *d* = 0.44, 95% CI
\[0.31, 0.57\], and intention (three items) did not differ significantly
in the predicted direction, *t*(871.84) = 0.61, *p* = .542, *d* = 0.04,
95% CI \[-0.09, 0.17\]. In this subsample, the replication would thus be
**not successful**. The indirect effect through self-perception was *b*
= -0.004, 95% CI \[-0.023, 0.011\] (H3 not supported).

# Exploratory analyses (not preregistered)

## Content of the manipulation texts

Frequencies of the final codes, separately for each condition.
Percentages are computed within condition (coded participants only).

``` r
data.coded <- data.complete %>% filter(coded)

freq_tab <- data.coded %>%
  mutate(Code = ifelse(is.na(final_code_en), "(no final code)", final_code_en)) %>%
  count(condition, Code) %>%
  group_by(condition) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup() %>%
  pivot_wider(names_from = condition, values_from = c(n, pct), values_fill = 0) %>%
  dplyr::select(Code, n_control, pct_control, n_experimental, pct_experimental) %>%
  arrange(Code)

freq_tab <- bind_rows(freq_tab,
                      summarise(freq_tab, Code = "Total", across(where(is.numeric), sum)))

kable(freq_tab, digits = 1,
      col.names = c("Final code", "Control n", "Control %", "Experimental n", "Experimental %"),
      caption = "Frequencies of final text codes by condition")
```

| Final code | Control n | Control % | Experimental n | Experimental % |
|:---|---:|---:|---:|---:|
| everyday_activities\_(valid) | 473 | 91.7 | 0 | 0.0 |
| family_love\_(valid) | 0 | 0.0 | 205 | 42.4 |
| friend’s_or_other_person’s_love\_(valid) | 0 | 0.0 | 29 | 6.0 |
| general_statements\_(not_valid) | 0 | 0.0 | 3 | 0.6 |
| irrelevant_or_nonsense\_(not_valid) | 43 | 8.3 | 51 | 10.5 |
| multiple_kinds_of_love\_(valid) | 0 | 0.0 | 80 | 16.5 |
| romantic_love\_(valid) | 0 | 0.0 | 103 | 21.3 |
| uncertain_or_unloved\_(not_valid) | 0 | 0.0 | 13 | 2.7 |
| Total | 516 | 100.0 | 484 | 100.0 |

Frequencies of final text codes by condition

``` r
# Codes without a final code: both coders left the text uncoded, or the disagreement was not resolved
table(Condition = data.coded$condition,
      Unresolved = data.coded$unresolved, useNA = "ifany")
```

    #>               Unresolved
    #> Condition        0
    #>   control      516
    #>   experimental 484

``` r
# Flags created when processing the codes (shown if present in text_codes.RDS)
if ("single_coder" %in% names(data.coded)) {
  print(table(Condition = data.coded$condition, `Coded by one coder only` = data.coded$single_coder))
}
if ("wrong_condition" %in% names(data.coded)) {
  print(table(Condition = data.coded$condition, `Content of the other condition` = data.coded$wrong_condition))
}
```

## Inter-rater reliability (coder 1 vs. coder 2)

Agreement between the two independent coders is reported for (a) the
refined content categories and (b) the binary classification of texts as
valid vs. invalid for the participant’s condition, separately for each
condition and overall. Agreement is computed on texts that both coders
coded; texts coded by only one coder are counted in a separate column
(under the coding protocol, these count as disagreements resolved by
coder 3). Confidence intervals for Cohen’s κ are percentile bootstrap
intervals (2000 resamples). Note that κ can be low despite high
percentage agreement when one category is very frequent (e.g., almost
all texts valid), so both indices should be interpreted together.

``` r
kappa_cohen <- function(a, b) {
  lev <- sort(union(a, b))
  tab <- table(factor(a, levels = lev), factor(b, levels = lev))
  n   <- sum(tab)
  po  <- sum(diag(tab)) / n
  pe  <- sum(rowSums(tab) * colSums(tab)) / n^2
  (po - pe) / (1 - pe)
}

irr_summary <- function(a, b, n_boot = params$n_boot_kappa) {
  a <- as.character(a); b <- as.character(b)
  one_missing <- sum(xor(is.na(a), is.na(b)))
  both <- !is.na(a) & !is.na(b)
  a <- a[both]; b <- b[both]
  boot <- replicate(n_boot, {
    i <- sample.int(length(a), replace = TRUE)
    kappa_cohen(a[i], b[i])
  })
  ci <- unname(quantile(boot, c(.025, .975), na.rm = TRUE))
  data.frame(`n (both coded)` = length(a), `n (one coder only)` = one_missing,
             `Agreement %` = 100 * mean(a == b),
             kappa = kappa_cohen(a, b), `95% CI lower` = ci[1], `95% CI upper` = ci[2],
             check.names = FALSE)
}

set.seed(params$seed)
irr_table <- bind_rows(lapply(c("control", "experimental", "overall"), function(cond) {
  d <- if (cond == "overall") data.coded else filter(data.coded, condition == cond)
  bind_rows(
    data.frame(Condition = cond, Coding = "Refined categories",
               irr_summary(d$coder1_en, d$coder2_en), check.names = FALSE),
    data.frame(Condition = cond, Coding = "Valid vs. invalid",
               irr_summary(d$valid_c1, d$valid_c2), check.names = FALSE)
  )
}))

kable(irr_table, digits = c(0, 0, 0, 0, 1, 2, 2, 2),
      caption = "Inter-rater reliability of the text coding (coder 1 vs. coder 2)")
```

| Condition | Coding | n (both coded) | n (one coder only) | Agreement % | kappa | 95% CI lower | 95% CI upper |
|:---|:---|---:|---:|---:|---:|---:|---:|
| control | Refined categories | 516 | 0 | 95.9 | 0.74 | 0.63 | 0.84 |
| control | Valid vs. invalid | 516 | 0 | 95.9 | 0.74 | 0.63 | 0.84 |
| experimental | Refined categories | 484 | 0 | 70.5 | 0.60 | 0.54 | 0.65 |
| experimental | Valid vs. invalid | 484 | 0 | 94.4 | 0.79 | 0.71 | 0.86 |
| overall | Refined categories | 1000 | 0 | 83.6 | 0.77 | 0.74 | 0.80 |
| overall | Valid vs. invalid | 1000 | 0 | 95.2 | 0.77 | 0.71 | 0.83 |

Inter-rater reliability of the text coding (coder 1 vs. coder 2)

## Small-telescopes test (Simonsohn, 2015)

The preregistration justifies the sample size partly by Simonsohn’s
(2015) small-telescopes heuristic. Here we apply the corresponding test:
*d*<sub>33%</sub> is the effect size the original study (*N* = 141) had
33% power to detect. If the upper bound of the one-sided 95% confidence
interval (i.e., the two-sided 90% CI) of the replication effect lies
below *d*<sub>33%</sub>, the replication suggests the original study
could not have meaningfully detected the effect. The split of the
original sample into groups is not reported by Tang et al. (2025), so
equal group sizes (70 and 71) are assumed.

``` r
d33 <- pwr::pwr.t2n.test(n1 = 70, n2 = 71, sig.level = .05, power = 1/3)$d

d90_3 <- compare_groups("intention",  data.complete, conf.level = .90)$d
d90_6 <- compare_groups("intention6", data.complete, conf.level = .90)$d

data.frame(
  Scale     = c("Intention (3 items; primary)", "Intention (6 items)"),
  d         = round(c(d90_3$estimate, d90_6$estimate), 3),
  CI90_low  = round(c(d90_3$conf.int[1], d90_6$conf.int[1]), 3),
  CI90_high = round(c(d90_3$conf.int[2], d90_6$conf.int[2]), 3),
  d33       = round(d33, 3),
  Smaller_than_d33 = c(d90_3$conf.int[2] < d33, d90_6$conf.int[2] < d33)
)
```

    #>                          Scale     d CI90_low CI90_high   d33 Smaller_than_d33
    #> 1 Intention (3 items; primary) 0.016   -0.088     0.120 0.259             TRUE
    #> 2          Intention (6 items) 0.007   -0.097     0.112 0.259             TRUE

For the three-item scale, the upper bound of the 90% CI of the
replication effect (0.12) was below *d*<sub>33%</sub> = 0.26.

# Supplementary tables

## Table A1: Means, standard deviations, and correlations

``` r
var_names <- c("treatment"  = "Treatment",
               "manipCheck" = "Manipulation check",
               "intention"  = "Intention (3 items)",
               "intention6" = "Intention (6 items)",
               "selfperc"   = "Self-perception")

df_summary <- data.complete[, names(var_names), drop = FALSE] %>%
  mutate(across(everything(), as.numeric))

corr_obj <- psych::corr.test(df_summary, use = "complete", method = "pearson", adjust = "none")
r_mat <- corr_obj$r
p_mat <- corr_obj$p

star_label <- function(p) {
  if (is.na(p)) return("")
  if (p < .001) "***" else if (p < .01) "**" else if (p < .05) "*" else "ns"
}

format_r <- function(r) {
  r_rounded <- round(r, 2)
  if (abs(r_rounded) < .01) r_rounded <- 0
  sub("^(-?)0", "\\1", sprintf("%.2f", r_rounded))
}

k <- length(var_names)
corr_table <- matrix("", k, k, dimnames = list(var_names, var_names))
for (i in seq_len(k)) {
  for (j in seq_len(k)) {
    if (i == j) {
      corr_table[i, j] <- "1.00"
    } else {
      corr_table[i, j] <- paste0(format_r(r_mat[i, j]), star_label(p_mat[i, j]))
    }
  }
}

table_a1 <- data.frame(
  Variable = var_names,
  M  = fmt(colMeans(df_summary)),
  SD = fmt(apply(df_summary, 2, sd)),
  corr_table,
  check.names = FALSE,
  row.names = NULL
)

names(table_a1) <- c("Variable", "M", "SD", names(var_names))

kable(table_a1, align = paste0("lrr", strrep("c", k)),
      caption = sprintf("Table A1. Means, standard deviations, and correlations (N = %d)", nrow(df_summary)))
```

| Variable | M | SD | treatment | manipCheck | intention | intention6 | selfperc |
|:---|---:|---:|:--:|:--:|:--:|:--:|:--:|
| Treatment | 0.48 | 0.50 | 1.00 | .16\*\*\* | .01ns | .00ns | -.04ns |
| Manipulation check | 5.27 | 1.77 | .16\*\*\* | 1.00 | .16\*\*\* | .20\*\*\* | .36\*\*\* |
| Intention (3 items) | 4.08 | 1.34 | .01ns | .16\*\*\* | 1.00 | .94\*\*\* | .07\* |
| Intention (6 items) | 4.04 | 1.38 | .00ns | .20\*\*\* | .94\*\*\* | 1.00 | .06ns |
| Self-perception | 4.75 | 1.08 | -.04ns | .36\*\*\* | .07\* | .06ns | 1.00 |

Table A1. Means, standard deviations, and correlations (N = 1000)

*Note.* Treatment is a dummy indicator (1 = experimental condition, 0 =
control condition). ns = not significant, \* *p* \< .05, \*\* *p* \<
.01, \*\*\* *p* \< .001.

## Table A2: Descriptive statistics by condition

``` r
desc_vars <- c(manipCheck = "Manipulation check",
               intention  = "Intention (3 items)",
               intention6 = "Intention (6 items)",
               selfperc   = "Self-perception")

n_ctrl <- sum(data.complete$treatment == 0)
n_exp  <- sum(data.complete$treatment == 1)

table_a2 <- data.complete %>%
  dplyr::select(treatment, all_of(names(desc_vars))) %>%
  pivot_longer(-treatment, names_to = "var") %>%
  group_by(var, treatment) %>%
  summarise(cell = sprintf("%.2f (%.2f)", mean(value), sd(value)), .groups = "drop") %>%
  pivot_wider(names_from = treatment, values_from = cell) %>%
  mutate(Variable = desc_vars[var]) %>%
  arrange(match(var, names(desc_vars))) %>%
  dplyr::select(Variable, `0`, `1`)

kable(table_a2,
      col.names = c("Variable",
                    sprintf("Control group (n = %d)", n_ctrl),
                    sprintf("Experimental group (n = %d)", n_exp)),
      caption = "Table A2. Descriptive statistics, M (SD), for the control and experimental groups")
```

| Variable            | Control group (n = 516) | Experimental group (n = 484) |
|:--------------------|:------------------------|:-----------------------------|
| Manipulation check  | 5.00 (1.79)             | 5.56 (1.71)                  |
| Intention (3 items) | 4.07 (1.33)             | 4.10 (1.36)                  |
| Intention (6 items) | 4.03 (1.37)             | 4.05 (1.39)                  |
| Self-perception     | 4.79 (1.07)             | 4.71 (1.08)                  |

Table A2. Descriptive statistics, M (SD), for the control and
experimental groups

# References

Blatný, M., & Osecká, L. (1994). Rosenbergova škála sebehodnocení:
Struktura globálního vztahu k sobě. *Československá psychologie, 38*(6),
481–488.

Kumar, S., Gupta, K., Kumar, A., Singh, A., & Singh, R. K. (2023).
Applying the theory of reasoned action to examine consumers’ attitude
and willingness to purchase organic foods. *International Journal of
Consumer Studies, 47*(1), 118–135.

Simonsohn, U. (2015). Small telescopes: Detectability and the evaluation
of replication results. *Psychological Science, 26*(5), 559–569.

Tang, D., Chen, J., Jiang, Q., Huang, H., & Zheng, M. (2025). How
feeling loved inspires the purchase of green foods. *Journal of
Environmental Psychology, 102*, 102523.

# Session information

``` r
sessionInfo()
```

    #> R version 4.5.1 (2025-06-13)
    #> Platform: x86_64-pc-linux-gnu
    #> Running under: Ubuntu 24.04.5 LTS
    #> 
    #> Matrix products: default
    #> BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0 
    #> LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0  LAPACK version 3.12.0
    #> 
    #> locale:
    #>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    #>  [3] LC_TIME=cs_CZ.UTF-8        LC_COLLATE=en_US.UTF-8    
    #>  [5] LC_MONETARY=cs_CZ.UTF-8    LC_MESSAGES=en_US.UTF-8   
    #>  [7] LC_PAPER=cs_CZ.UTF-8       LC_NAME=C                 
    #>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    #> [11] LC_MEASUREMENT=cs_CZ.UTF-8 LC_IDENTIFICATION=C       
    #> 
    #> time zone: Europe/Prague
    #> tzcode source: system (glibc)
    #> 
    #> attached base packages:
    #> [1] stats     graphics  grDevices utils     datasets  methods   base     
    #> 
    #> other attached packages:
    #>  [1] dplyr_1.1.4          tidyr_1.3.1          knitr_1.45          
    #>  [4] effsize_0.8.1        pwr_1.3-0            mediation_4.5.1     
    #>  [7] sandwich_3.1-0       mvtnorm_1.3-3        Matrix_1.7-4        
    #> [10] MASS_7.3-65          GPArotation_2024.2-1 psych_2.4.1         
    #> [13] car_3.1-2            carData_3.0-5        pacman_0.5.1        
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] gtable_0.3.6        xfun_0.41           ggplot2_4.0.3      
    #>  [4] htmlwidgets_1.6.4   lattice_0.22-5      vctrs_0.6.5        
    #>  [7] tools_4.5.1         Rdpack_2.6.4        generics_0.1.4     
    #> [10] parallel_4.5.1      tibble_3.3.0        cluster_2.1.8.1    
    #> [13] pkgconfig_2.0.3     data.table_1.18.6.1 checkmate_2.3.1    
    #> [16] RColorBrewer_1.1-3  S7_0.2.2            lifecycle_1.0.4    
    #> [19] compiler_4.5.1      farver_2.1.2        stringr_1.5.2      
    #> [22] mnormt_2.1.1        codetools_0.2-20    htmltools_0.5.7    
    #> [25] yaml_2.3.12         htmlTable_2.4.2     Formula_1.2-5      
    #> [28] pillar_1.10.2       nloptr_2.0.3        Hmisc_5.1-1        
    #> [31] reformulas_0.4.1    rpart_4.1.24        boot_1.3-31        
    #> [34] abind_1.4-8         nlme_3.1-168        tidyselect_1.2.1   
    #> [37] digest_0.6.37       stringi_1.8.7       purrr_1.0.4        
    #> [40] splines_4.5.1       fastmap_1.1.1       grid_4.5.1         
    #> [43] colorspace_2.1-1    cli_3.6.5           magrittr_2.0.4     
    #> [46] base64enc_0.1-3     utf8_1.2.6          dichromat_2.0-0.1  
    #> [49] withr_3.0.2         foreign_0.8-90      scales_1.4.0       
    #> [52] backports_1.4.1     rmarkdown_2.25      nnet_7.3-20        
    #> [55] lme4_1.1-37         gridExtra_2.3       zoo_1.8-12         
    #> [58] lpSolve_5.6.20      evaluate_0.23       rbibutils_2.3      
    #> [61] rlang_1.1.6         Rcpp_1.1.1          glue_1.8.0         
    #> [64] rstudioapi_0.15.0   minqa_1.2.8         jsonlite_1.8.8     
    #> [67] R6_2.6.1
