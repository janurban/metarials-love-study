# This is an analytical R script for Replication study 1 of Tang et al.'s (2025) Study 1a

# # Data anonymization
# if (!require("pacman")) install.packages("pacman")
# pacman::p_load(car, lme4, stringr, reshape, car, TAM, psych, plyr, dplyr, stringr,
#                tidyr, readxl, sjPlot, correlation, effsize, mediation, effsize, writexl) 

# # Reading data from local file
# setwd("/home/honza/Google Drive/CUK/uceni/FSV/Kvant praktika/Love study 2025-2026 kvant praktika/s1 data and analyses")
# data <- read_excel("results-survey997596.xlsx")

# # renaming the two id variables in the dataset
# data <- data %>%
#   rename_with(~ c("syst.id", "subject") , c("id...1", "id...10"))

# # Anonymizing/deleting potentially identifying information
# data$subject <- as.numeric(as.factor(data$subject)) # anonymizing person's identifier
# data$ipaddr <- NULL # deleting IP address
# data$refurl <- NULL # deleting redirecting URL
# data$treatment <- NULL # deleting the written text in the treatment condition
# data$control <- NULL # deleting the written text in the control condition

# # Saving anonymized data
# setwd("/home/honza/github/metarials-love-study")
# saveRDS(data, "data_s1_anonymized.RDS")
# write_xlsx(data, "data_s1_anonymized.xlsx") 


# Analytical script that works with anonymized data starts here
# This is an analytical R script for Replication study 1 of Tang et al.'s (2025) Study 1a

if (!require("pacman")) install.packages("pacman")
pacman::p_load(car, lme4, stringr, reshape, car, TAM, psych, plyr, dplyr, stringr,
               tidyr, readxl, sjPlot, correlation, effsize, mediation, effsize, writexl) 

# Function that reads data from Github
read_rds_from_url <- function(url) {
  tmp <- tempfile(fileext = ".rds")
  download.file(url, destfile = tmp, mode = "wb", quiet = TRUE)
  readRDS(tmp)
}

data <- read_rds_from_url(
  "https://raw.githubusercontent.com/janurban/metarials-love-study/main/data_s1_anonymized.RDS"
)

data <-   data %>% mutate(across(c(gender:educ, random, manipCheck, `intention[SQ001]`:repeated), ~ as.numeric(gsub("A", "", .x))))

# Time statistics
summary_table <- data %>% filter(lastpage == 17) %>%
  summarise(
    across(
      .cols = c(interviewtime, consentTime, confirmTextTime, productInfoTime, treatmentTime, controlTime,
        manipCheckTime, intentionTime, selfperceptionTime, repeatedTime, debriefingTime),
      .fns = list(
        M      = ~ mean(.x, na.rm = TRUE),
        SD     = ~ sd(.x, na.rm = TRUE),
        Min    = ~ fivenum(.x, na.rm  = T)[1],
        Q1     = ~ fivenum(.x, na.rm  = T)[2],
        Median = ~ fivenum(.x, na.rm  = T)[3],
        Q3     = ~ fivenum(.x, na.rm  = T)[4],
        Max    = ~ fivenum(.x, na.rm  = T)[5]
      ),
      .names = "{.col}___{.fn}"
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("Variable", ".value"),
    names_pattern = "(.*)___(.*)"
  ) %>% as.data.frame()
print(summary_table, digits = 4, row.names = FALSE)


data$status <- "screenout"
data$status <- ifelse (!is.na(data$random), "dropout", data$status)
data$status <- ifelse (data$lastpage == 17, "complete", data$status)

data <- data %>% filter(status != "screenout") 

table(data$status)
# complete  dropout 
#     1000      679
 
prop.table(table(data$status))
#  complete   dropout 
# 0.5955926 0.4044074 

data.complete <- data %>% filter(lastpage == 17)

# Sociodemographics
data.complete$ageNum2 <- data.complete$ageNum + 14
mean(data.complete$ageNum2) # 46.318
sd(data.complete$ageNum2) # 14.94

data.complete$gender<- car::recode(data.complete$gender, "1 = 'M'; 2 = 'F'; 3 = 'other.'")
prop.table(table(data.complete$gender)) 
# F      M      other. 
# 0.516  0.481  0.003 

prop.table(table(data.complete$educ))
# 1 = elementary + high school without leaving exam
# 3 = high school with leaving exam
# 4 = university education
#     1     3     4 
# 0.303 0.429 0.268 

# renaming variables
colnames(data.complete)
data.complete <- data.complete %>%
  rename_with(~ paste0("intention_", c(1:6)), `intention[SQ001]`:`intention[SQ006]`)
data.complete <- data.complete %>%
  rename_with(~ paste0("selfperc_", c(1:10)), `selfperception[SQ001]`:`selfperception[SQ010]`)

# # Recoding values to numeric
# data.complete <-   data.complete %>% mutate(across(c(random, manipCheck, intention_1:repeated), ~ as.numeric(gsub("A", "", .x))))
# intention

intention0 <- data.complete %>% dplyr::select(intention_1:intention_6)
intention <- intention0
intention[,3] <- 8 - intention0[,3] 
data.complete$intention <- rowMeans(intention) # intention score for the use in the main analysis
data.complete$intention3 <- rowMeans(intention[, c(1:3)]) # intention score for the exploratory analysis with 3 items only
psych::omega(intention) # reliability for the 6-item scale
psych::omega(intention[, c(1:3)]) # reliability for the 3-item scale

selfperc0 <- data.complete %>% dplyr::select(selfperc_1:selfperc_10)
selfperc <- selfperc0
selfperc[,c(2, 5, 6, 8, 9 )] <- 8 - selfperc0[,c(2, 5, 6, 8, 9 )] 
data.complete$selfperc <- rowMeans(selfperc)
psych::omega(selfperc) # reliability for the 10-item scale

data.complete$treatment <- as.factor(data.complete$random-1) # treatment(1 = experimental, 0 = control)

t.test(manipCheck ~ treatment, data = data.complete)
cohen.d(manipCheck ~ treatment, data = data.complete) #d = -0.3197748

t.test(intention ~ treatment, data = data.complete)
cohen.d(intention ~ treatment, data = data.complete) #d = -0.007419511

# exploratory analysis with 3 intention items
t.test(intention3 ~ treatment, data = data.complete)
cohen.d(intention3 ~ treatment, data = data.complete) #d = -0.007419511

t.test(selfperc ~ treatment, data = data.complete)
cohen.d(selfperc ~ treatment, data = data.complete) #d = 0.07840219 

cor.test(data.complete$selfperc, data.complete$intention)

# Run mediation model
set.seed(1234)
model_m <- lm(selfperc ~ treatment, data = data.complete)
model_y <- lm(intention ~ treatment + selfperc, data = data.complete)

results <- mediation::mediate(model_m, 
                              model_y, 
                              treat = "treatment", 
                              mediator = "selfperc", 
                              boot = TRUE,            
                              sims = 5000)
summary(model_m)
summary(model_y)  
summary(results)

#  Descriptive statistics (M, SD, correlations)
# 1. Variables for the final descriptive and correlation tables
vars <- c("treatment", "manipCheck", "intention", "selfperc")

df_summary <- data.complete[, vars, drop = FALSE] |>
  mutate(across(everything(), ~ as.numeric(as.character(.x))))

# 2. Means and SDs
summary_table <- psych::describe(df_summary, fast = TRUE) |>
  as.data.frame() |>
  tibble::rownames_to_column("Variable") |>
  transmute(
    Variable,
    M = round(mean, 2),
    SD = round(sd, 2)
  )

# 3. Correlation matrix with stars instead of p-values
corr_obj <- psych::corr.test(df_summary, use = "complete", method = "pearson", adjust = "none")

r_matrix <- round(corr_obj$r, 2)
p_matrix <- corr_obj$p

format_r <- function(r) {
  r_rounded <- round(r, 2)
  if (is.na(r_rounded)) return(NA_character_)
  if (abs(r_rounded) < 0.01) r_rounded <- 0
  sub("^0", "", sprintf("%.2f", r_rounded))
}

star_label <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "ns"
}

corr_table <- matrix(NA_character_, nrow = nrow(r_matrix), ncol = ncol(r_matrix))
for (i in seq_len(nrow(r_matrix))) {
  for (j in seq_len(ncol(r_matrix))) {
    if (i == j) {
      corr_table[i, j] <- "1.00"
    } else {
      corr_table[i, j] <- paste0(format_r(r_matrix[i, j]), star_label(p_matrix[i, j]))
    }
  }
}
rownames(corr_table) <- rownames(r_matrix)
colnames(corr_table) <- colnames(r_matrix)

# 4. Final table layout: variable, M, SD, and each pairwise correlation with stars
combined_table <- summary_table |>
  mutate(
    treatment = corr_table[, "treatment"],
    manipCheck = corr_table[, "manipCheck"],
    intention = corr_table[, "intention"],
    selfperc = corr_table[, "selfperc"]
  )

colnames(combined_table) <- c("Variable", "M", "SD", "treatment", "manipCheck", "intention", "selfperc")
print(combined_table, row.names = FALSE)





