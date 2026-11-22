This document reproduces the analyses in the *Report From Pilot Study
PS1 (Love Replication)*. It follows the original analysis script
(`Analysis_v3.R`). Knit this file in RStudio (**Knit** button) or with
`rmarkdown::render("analysis_pilot1.Rmd")` to regenerate
`analysis_pilot1.md` in the same folder.

**Data source.** By default, the document downloads
`pilot1/data_pilot1.csv` from the GitHub repository
[janurban/metarials-love-study](https://github.com/janurban/metarials-love-study).
It can also read a local `.csv` file by passing a path with
`rmarkdown::render("analysis_pilot1.Rmd", params = list(data = "path/to/file.csv"))`.

**Required packages.** `dplyr`, `readr`, `knitr`, `psych`, `GPArotation`
(`psych::omega()` needs `GPArotation`). Install them with
`install.packages(c("dplyr", "readr", "knitr", "psych", "GPArotation", "rmarkdown"))`.

# 0. Setup

``` r
library(dplyr)   # data handling
library(readr)   # parse_number()
library(knitr)   # tables
library(psych)   # alpha and omega
```

``` r
fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)
pct <- function(n, N, d = 1) formatC(100 * n / N, format = "f", digits = d)

# Frequency table with percentages
tab_freq <- function(x, labels) {
  f <- factor(x, levels = seq_along(labels), labels = labels)
  n <- table(f, useNA = "ifany")
  data.frame(Category = names(n), n = as.vector(n), Percent = pct(as.vector(n), length(x)))
}
```

# 1. Data

The file is the LimeSurvey export. Columns are selected by name, not by
position. The checks below fail with an error if the export has a
different structure.

``` r
data_source <- params$data
if (grepl("^https?://", data_source, ignore.case = TRUE)) {
  # A github.com/.../blob/... link points to the web page, not the file:
  # rewrite it to the raw-file address
  data_source <- sub("^https?://github\\.com/([^/]+)/([^/]+)/blob/",
                     "https://raw.githubusercontent.com/\\1/\\2/", data_source)
  data_file <- tempfile(fileext = ".csv")
  download_status <- download.file(data_source, data_file, mode = "wb", quiet = TRUE)
  if (download_status != 0 || !file.exists(data_file) || file.size(data_file) == 0) {
    stop("Could not download the data file from ", data_source)
  }
} else {
  data_file <- data_source
}

data0 <- read.csv(data_file, header = TRUE, check.names = FALSE,
                  stringsAsFactors = FALSE, encoding = "UTF-8")

expected <- c("id", "submitdate", "lastpage", "startlanguage", "seed", "text",
              paste0("experience[SQ00", 1:2, "]"),
              paste0("quiz[SQ00", 1:9, "]"),
              paste0("intention[SQ00", 1:6, "]"),
              "gender", "age", "educ", "comment")
missing_cols <- setdiff(expected, names(data0))
if (length(missing_cols) > 0) stop("Columns not found: ", paste(missing_cols, collapse = ", "))
dim(data0)
```

    ## [1] 166  27

## 1.1 Sample

As in the original script, the first two rows are test cases and are
removed. Only participants who reached the last page (`lastpage == 8`)
are analysed.

``` r
data_entered <- data0[-c(1, 2), ]
data <- data_entered %>% filter(lastpage == 8)
N <- nrow(data)
c(entered = nrow(data_entered), completed = N,
  completion_rate_pct = round(100 * N / nrow(data_entered)))
```

    ##             entered           completed completion_rate_pct 
    ##                 164                 102                  62

## 1.2 Variables

In the LimeSurvey export:

- `experience[SQ001]` is familiarity and `experience[SQ002]` is past
  purchase. Codes: 1 = yes, 2 = no, 3 = I do not know.
- `quiz[SQ001]` is the item without text. Its instruction told
  participants to choose “True”, and it is not analysed (see the note in
  the questionnaire).
- `quiz[SQ002]` to `quiz[SQ009]` are quiz items 1–8. Codes: 1 = true, 2
  = false.
- `intention[SQ001]` to `intention[SQ006]` are the six intention items
  (1–7).

Answer codes such as `A1` are turned into numbers with
`readr::parse_number()`, as in the original script.

``` r
to_num <- function(x) as.numeric(parse_number(as.character(x)))

data <- data %>%
  transmute(
    id,
    familiarity  = to_num(`experience[SQ001]`),
    pastPurchase = to_num(`experience[SQ002]`),
    emptyItem    = to_num(`quiz[SQ001]`),
    across(all_of(paste0("quiz[SQ00", 2:9, "]")), to_num),
    across(all_of(paste0("intention[SQ00", 1:6, "]")), to_num),
    gender = to_num(gender),
    age    = to_num(age),
    educ   = to_num(educ))
names(data) <- c("id", "familiarity", "pastPurchase", "emptyItem",
                 paste0("quiz", 1:8), paste0("intention", 1:6), "gender", "age", "educ")

# Missing values among the analysed variables
colSums(is.na(select(data, -id)))
```

    ##  familiarity pastPurchase    emptyItem        quiz1        quiz2        quiz3 
    ##            0            0            0            0            0            0 
    ##        quiz4        quiz5        quiz6        quiz7        quiz8   intention1 
    ##            0            0            0            0            0            0 
    ##   intention2   intention3   intention4   intention5   intention6       gender 
    ##            0            0            0            0            0            0 
    ##          age         educ 
    ##            0            0

# 2. Participants

``` r
c(M = mean(data$age, na.rm = TRUE), SD = sd(data$age, na.rm = TRUE),
  Min = min(data$age, na.rm = TRUE), Max = max(data$age, na.rm = TRUE))
```

    ##        M       SD      Min      Max 
    ## 25.97059 10.42793 18.00000 66.00000

**Gender**

``` r
kable(tab_freq(data$gender, c("Man", "Woman", "Non-binary or gender diverse",
                              "Would rather not disclose")))
```

| Category                     |   n | Percent |
|:-----------------------------|----:|:--------|
| Man                          |  29 | 28.4    |
| Woman                        |  69 | 67.6    |
| Non-binary or gender diverse |   2 | 2.0     |
| Would rather not disclose    |   2 | 2.0     |

**Highest completed education**

``` r
kable(tab_freq(data$educ, c("Primary", "Secondary without high school diploma",
                            "Secondary with high school diploma", "University")))
```

| Category                              |   n | Percent |
|:--------------------------------------|----:|:--------|
| Primary                               |   4 | 3.9     |
| Secondary without high school diploma |   1 | 1.0     |
| Secondary with high school diploma    |  63 | 61.8    |
| University                            |  34 | 33.3    |

# 3. Familiarity and past purchase

**“Have you ever seen Ecological Products with the green-leaf label in a
store?”**

``` r
kable(tab_freq(data$familiarity, c("Yes", "No", "I do not know")))
```

| Category      |   n | Percent |
|:--------------|----:|:--------|
| Yes           |  74 | 72.5    |
| No            |   9 | 8.8     |
| I do not know |  19 | 18.6    |

**“Have you ever bought ecological products?”**

``` r
kable(tab_freq(data$pastPurchase, c("Yes", "No", "I do not know")))
```

| Category      |   n | Percent |
|:--------------|----:|:--------|
| Yes           |  71 | 69.6    |
| No            |   0 | 0.0     |
| I do not know |  31 | 30.4    |

# 4. Knowledge quiz

Items 1, 3, 4 and 5 are false statements, and items 2, 6, 7 and 8 are
true statements. Each correct answer scores 1. This is the same scoring
as the original script, which reverses items 1, 3, 4 and 5 (`3 - x`) and
then computes `2 - x`.

``` r
quiz_items <- paste0("quiz", 1:8)
key <- c(2, 1, 2, 2, 2, 1, 1, 1)   # correct answer: 1 = true, 2 = false
quiz_text <- c(
  "1. The production of Ecological Products involves the heavy use of pesticides and artificial fertilizers to ensure high yields. (False)",
  "2. Ecological Products must comply with specific technical standards and comprehensive quality control. (True)",
  "3. Ecological Products are distinguished from conventional foods only by their price and packaging design. (False)",
  "4. The use of the Ecological Product logo is not controlled by any specific regulation. (False)",
  "5. Ecological Products may compromise on safety to achieve environmental benefits. (False)",
  "6. The Ecological Product Logo is also intended to remind people of the importance of environmental protection and preventing pollution. (True)",
  "7. The Ecological Product logo aims to strengthen harmony between humans and the natural world. (True)",
  "8. Ecological products are produced in an excellent ecological environment with minimal use of chemicals. (True)")

correct <- sweep(as.matrix(data[quiz_items]), 2, key, `==`) * 1
data$quiz_score <- rowSums(correct)

# Same result as the original scoring
quiz_orig <- as.matrix(data[quiz_items])
quiz_orig[, c(1, 3, 4, 5)] <- 3 - quiz_orig[, c(1, 3, 4, 5)]
stopifnot(isTRUE(all.equal(unname(rowSums(2 - quiz_orig)), unname(data$quiz_score))))

c(M = mean(data$quiz_score, na.rm = TRUE), SD = sd(data$quiz_score, na.rm = TRUE))
```

    ##         M        SD 
    ## 7.5686275 0.6676366

Distribution of the number of correct answers (0–8):

``` r
kable(as.data.frame(table(`Correct answers` = data$quiz_score), responseName = "n"))
```

| Correct.answers |   n |
|:----------------|----:|
| 5               |   1 |
| 6               |   7 |
| 7               |  27 |
| 8               |  67 |

Accuracy for each item:

``` r
n_correct <- colSums(correct, na.rm = TRUE)
kable(data.frame(Item = quiz_text,
                 `Correct, n (%)` = paste0(n_correct, " (", pct(n_correct, N), "%)"),
                 check.names = FALSE, row.names = NULL))
```

| Item                                                                                                                                             | Correct, n (%) |
|:-------------------------------------------------------------------------------------------------------------------------------------------------|:---------------|
| 1\. The production of Ecological Products involves the heavy use of pesticides and artificial fertilizers to ensure high yields. (False)         | 100 (98.0%)    |
| 2\. Ecological Products must comply with specific technical standards and comprehensive quality control. (True)                                  | 99 (97.1%)     |
| 3\. Ecological Products are distinguished from conventional foods only by their price and packaging design. (False)                              | 100 (98.0%)    |
| 4\. The use of the Ecological Product logo is not controlled by any specific regulation. (False)                                                 | 100 (98.0%)    |
| 5\. Ecological Products may compromise on safety to achieve environmental benefits. (False)                                                      | 93 (91.2%)     |
| 6\. The Ecological Product Logo is also intended to remind people of the importance of environmental protection and preventing pollution. (True) | 99 (97.1%)     |
| 7\. The Ecological Product logo aims to strengthen harmony between humans and the natural world. (True)                                          | 83 (81.4%)     |
| 8\. Ecological products are produced in an excellent ecological environment with minimal use of chemicals. (True)                                | 98 (96.1%)     |

Check of the item without text (participants were told to choose “True”
= 1):

``` r
table(emptyItem = data$emptyItem, useNA = "ifany")
```

    ## emptyItem
    ##  1  2 
    ## 76 26

# 5. Purchase intention

Item 3 (“I am unwilling to buy Ecological Products”) is reverse-coded
(`8 - x`). Items 1–3 are the original measure (Tang et al., 2025, Study
1a); items 4–6 come from Kumar et al. (2023).

``` r
int_items <- paste0("intention", 1:6)
Intention <- data[int_items]
Intention$intention3 <- 8 - Intention$intention3

data$intention_3item <- rowMeans(Intention[, 1:3])
data$intention_6item <- rowMeans(Intention)

int_text <- c("1. I am willing to spend more on Ecological Products, as it is worth the extra charge.",
              "2. I plan to purchase food with the Ecological Products logo next time I shop.",
              "3. I am unwilling to buy Ecological Products. (reverse-coded)",
              "4. I intend to purchase Ecological products in the future.",
              "5. I am always interested in buying more Ecological Products for the family's needs.",
              "6. I always intend to look for Ecological Products, although outside the city.")
kable(data.frame(Item = int_text,
                 M  = fmt(sapply(Intention, mean, na.rm = TRUE)),
                 SD = fmt(sapply(Intention, sd, na.rm = TRUE)),
                 row.names = NULL))
```

| Item                                                                                    | M    | SD   |
|:----------------------------------------------------------------------------------------|:-----|:-----|
| 1\. I am willing to spend more on Ecological Products, as it is worth the extra charge. | 4.36 | 1.40 |
| 2\. I plan to purchase food with the Ecological Products logo next time I shop.         | 3.88 | 1.26 |
| 3\. I am unwilling to buy Ecological Products. (reverse-coded)                          | 5.56 | 1.64 |
| 4\. I intend to purchase Ecological products in the future.                             | 5.22 | 1.48 |
| 5\. I am always interested in buying more Ecological Products for the family’s needs.   | 3.89 | 1.47 |
| 6\. I always intend to look for Ecological Products, although outside the city.         | 3.37 | 1.51 |

## 5.1 Descriptive statistics and floor or ceiling effects

A floor or ceiling effect would show as many participants at the lowest
(1) or highest (7) possible scale mean.

``` r
scale_desc <- function(x, label) data.frame(
  Scale = label, M = fmt(mean(x)), SD = fmt(sd(x)), Min = fmt(min(x)), Max = fmt(max(x)),
  `% at floor (1)` = pct(sum(x == 1), length(x)),
  `% at ceiling (7)` = pct(sum(x == 7), length(x)),
  check.names = FALSE)
kable(rbind(scale_desc(data$intention_3item, "3 items (Tang et al., 2025, S1a)"),
            scale_desc(data$intention_6item, "6 items (extended)")))
```

| Scale                            | M    | SD   | Min  | Max  | % at floor (1) | % at ceiling (7) |
|:---------------------------------|:-----|:-----|:-----|:-----|:---------------|:-----------------|
| 3 items (Tang et al., 2025, S1a) | 4.60 | 1.14 | 1.00 | 6.67 | 2.0            | 0.0              |
| 6 items (extended)               | 4.38 | 1.06 | 1.17 | 6.33 | 0.0            | 0.0              |

## 5.2 Internal consistency

The original script uses `psych::omega()` with its defaults (three group
factors, oblimin rotation). For the 3-item scale, this solution is
underidentified, and `psych` gives a warning. The values below are the
ones `psych::omega()` returns.

``` r
om3 <- psych::omega(Intention[, 1:3], plot = FALSE)
om6 <- psych::omega(Intention, plot = FALSE)
kable(data.frame(
  Scale = c("3 items (Tang et al., 2025, S1a)", "6 items (extended)"),
  alpha   = fmt(c(om3$alpha, om6$alpha)),
  omega_h = fmt(c(om3$omega_h, om6$omega_h)),
  omega_t = fmt(c(om3$omega.tot, om6$omega.tot))),
  col.names = c("Scale", "Cronbach's alpha", "Omega hierarchical", "Omega total"))
```

| Scale                            | Cronbach’s alpha | Omega hierarchical | Omega total |
|:---------------------------------|:-----------------|:-------------------|:------------|
| 3 items (Tang et al., 2025, S1a) | 0.71             | 0.06               | 0.72        |
| 6 items (extended)               | 0.82             | 0.69               | 0.88        |

# 6. Comparison with the report

The report gives:

- N = 163 entered and N = 104 completed (63%).
- Age M = 25.90, SD = 10.34; 62% with a high school diploma, 33% with a
  university degree; 68% women.
- 71% had seen Ecological Products, 9% had not and 20% did not know.
- 68% had bought them and 32% did not know.
- Quiz: M = 7.57 correct answers, SD = 0.67.
- 3-item scale: M = 4.61, SD = 1.13, α = 0.63, ω<sub>t</sub> = 0.68.
- 6-item scale: M = 4.41, SD = 1.07, α = 0.80, ω<sub>h</sub> = 0.85.

Two points to check against the output above:

1.  **Quiz score.** The comments in the original script give *M* = 7.5,
    *SD* = 0.82. The report gives *M* = 7.57, *SD* = 0.67. Only one pair
    can match the data.
2.  **Omega for the 6-item scale.** The report gives ω<sub>h</sub> =
    0.85, which is higher than α = 0.80. Hierarchical omega is usually
    lower than α. The value may be ω<sub>t</sub> (omega total), so the
    label “ω<sub>h</sub>” may need to be corrected.

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
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=cs_CZ.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=cs_CZ.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=cs_CZ.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=cs_CZ.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Europe/Prague
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] psych_2.4.1 knitr_1.45  readr_2.1.5 dplyr_1.1.4
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] vctrs_0.6.5          nlme_3.1-168         cli_3.6.5           
    ##  [4] rlang_1.1.6          xfun_0.41            GPArotation_2024.2-1
    ##  [7] generics_0.1.4       glue_1.8.0           htmltools_0.5.7     
    ## [10] hms_1.1.3            rmarkdown_2.25       grid_4.5.1          
    ## [13] evaluate_0.23        tibble_3.3.0         tzdb_0.5.0          
    ## [16] fastmap_1.1.1        yaml_2.3.12          lifecycle_1.0.4     
    ## [19] compiler_4.5.1       pkgconfig_2.0.3      lattice_0.22-5      
    ## [22] digest_0.6.37        R6_2.6.1             tidyselect_1.2.1    
    ## [25] pillar_1.10.2        mnormt_2.1.1         parallel_4.5.1      
    ## [28] magrittr_2.0.4       withr_3.0.2          tools_4.5.1
