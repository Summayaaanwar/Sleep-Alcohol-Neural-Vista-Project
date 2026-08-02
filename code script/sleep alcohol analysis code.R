############################################################
# NHANES 2017-2018: Alcohol Use Patterns and Sleep Disorders
# CORRECTED, REPRODUCIBLE ANALYSIS
#
# Main decisions:
# 1) Start from ALL adults >=18 years
# 2) Keep ONE clean master analysis sheet
# 3) Preserve raw NHANES variable codes
# 4) Create derived variables once
# 5) Remove functional_impairment from master sheet and all models
# 6) Remove multiple_sleep_problems from master sheet and all models
# 7) Use model-specific complete-case analysis
# 8) Use NHANES survey design: WTMEC2YR, SDMVPSU, SDMVSTRA
# 9) Include p-values in all model outputs
# 10) Export tables + figures
############################################################

rm(list = ls())
gc()
options(repos = c(CRAN = "https://cloud.r-project.org"))

############################################################
# 00 — PACKAGES
############################################################

required_pkgs <- c(
  "dplyr",
  "tidyr",
  "stringr",
  "readr",
  "haven",
  "survey",
  "ggplot2",
  "tibble"
)

optional_pkgs <- c(
  "DiagrammeR",
  "DiagrammeRsvg",
  "rsvg",
  "htmlwidgets"
)

install_if_missing <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p, dependencies = TRUE)
    }
  }
}

install_if_missing(required_pkgs)

library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(haven)
library(survey)
library(ggplot2)
library(tibble)

options(survey.lonely.psu = "adjust")
set.seed(1234)

############################################################
# 01 — FOLDER STRUCTURE
############################################################

base_dir    <- "E:/Neural Vista/projct alcohol and sleep"
dir_raw     <- file.path(base_dir, "raw_data")
dir_outputs <- file.path(base_dir, "outputs")
dir_step    <- file.path(dir_outputs, "data")
dir_tables  <- file.path(dir_outputs, "tables")
dir_figures <- file.path(dir_outputs, "figures")
dir_logs    <- file.path(dir_outputs, "logs")

for (p in c(base_dir, dir_raw, dir_step, dir_outputs, dir_tables, dir_figures, dir_logs)) {
  dir.create(p, showWarnings = FALSE, recursive = TRUE)
}

setwd(base_dir)
cat("Project folder:", base_dir, "\n")

############################################################
# 02 — GLOBAL SETTINGS
############################################################

CYCLE <- "2017-2018"
AGE_MIN <- 18

############################################################
# 03 — HELPER FUNCTIONS
############################################################

# Check whether downloaded file is a valid SAS XPT file
is_valid_xpt <- function(path) {
  hdr <- tryCatch(readChar(path, nchars = 40, useBytes = TRUE), error = function(e) "")
  grepl("HEADER RECORD", hdr, fixed = TRUE)
}

# Read a required local XPT file. Raw files are never downloaded, deleted, or changed.
load_local_xpt <- function(file_name, raw_dir = dir_raw) {
  local_path <- file.path(raw_dir, file_name)
  if (!file.exists(local_path)) stop("Missing raw file: ", local_path)
  if (!is_valid_xpt(local_path)) stop("Invalid XPT file: ", local_path)
  message("Reading: ", local_path)
  read_xpt(local_path)
}

# Safe numeric conversion
to_num <- function(x) suppressWarnings(as.numeric(x))

# Keep only the documented valid values for a categorical variable.
keep_codes <- function(x, valid) {
  x <- to_num(x)
  x[!is.na(x) & !(x %in% valid)] <- NA_real_
  x
}

# Keep only a documented continuous range.
keep_range <- function(x, lower, upper) {
  x <- to_num(x)
  x[!is.na(x) & (x < lower | x > upper)] <- NA_real_
  x
}

# Clean coefficient names for exported summary tables
clean_level_name <- function(x) {
  x %>%
    str_replace("^.*?([A-Za-z].*)$", "\\1") %>%
    str_replace("^alcohol_category", "") %>%
    str_replace("^alcohol_user", "") %>%
    str_replace("^alcohol_risk", "") %>%
    str_replace("^gender", "") %>%
    str_replace("^sleep_trouble", "") %>%
    str_replace("^daytime_sleepy", "") %>%
    str_replace("^insufficient_sleep", "") %>%
    str_replace("^snoring_frequent", "") %>%
    str_replace("^breathing_pauses", "") %>%
    str_trim()
}

# Extract p-values from svyglm
extract_pvals <- function(model) {
  sm <- summary(model)
  cm <- sm$coefficients
  
  if ("Pr(>|t|)" %in% colnames(cm)) {
    return(as.numeric(cm[, "Pr(>|t|)"]))
  }
  
  if ("Pr(>|z|)" %in% colnames(cm)) {
    return(as.numeric(cm[, "Pr(>|z|)"]))
  }
  
  if ("t value" %in% colnames(cm)) {
    stat <- as.numeric(cm[, "t value"])
    df <- tryCatch(survey::degf(model$survey.design), error = function(e) NA_real_)
    if (!is.na(df)) {
      return(2 * pt(abs(stat), df = df, lower.tail = FALSE))
    }
  }
  
  if ("z value" %in% colnames(cm)) {
    stat <- as.numeric(cm[, "z value"])
    return(2 * pnorm(abs(stat), lower.tail = FALSE))
  }
  
  rep(NA_real_, nrow(cm))
}

# Significance stars
sig_stars <- function(p) {
  dplyr::case_when(
    !is.na(p) & p < 0.001 ~ "***",
    !is.na(p) & p < 0.01  ~ "**",
    !is.na(p) & p < 0.05  ~ "*",
    !is.na(p) ~ "NS",
    TRUE ~ "Not applicable"
  )
}

# Extract OR table from logistic survey model
extract_or <- function(model, model_name, design_obj, outcome_name) {
  b <- coef(model)
  se <- sqrt(diag(vcov(model)))
  df_test <- survey::degf(design_obj)
  statistic <- as.numeric(b / se)
  pvals <- 2 * pt(abs(statistic), df = df_test, lower.tail = FALSE)
  critical <- qt(0.975, df = df_test)
  
  tibble(
    outcome = outcome_name,
    model = model_name,
    n = nrow(design_obj$variables),
    term = names(b),
    estimate = as.numeric(b),
    se = as.numeric(se),
    test_statistic = statistic,
    df_test = df_test,
    p_method = "Design-based Wald t test (df = PSU - strata)",
    OR = exp(b),
    CI_low = exp(b - critical * se),
    CI_high = exp(b + critical * se),
    p = as.numeric(pvals),
    significance = sig_stars(as.numeric(pvals))
  )
}

# Extract beta table from linear survey model
extract_beta <- function(model, model_name, design_obj, outcome_name) {
  b <- coef(model)
  se <- sqrt(diag(vcov(model)))
  df_test <- survey::degf(design_obj)
  statistic <- as.numeric(b / se)
  pvals <- 2 * pt(abs(statistic), df = df_test, lower.tail = FALSE)
  critical <- qt(0.975, df = df_test)
  
  tibble(
    outcome = outcome_name,
    model = model_name,
    n = nrow(design_obj$variables),
    term = names(b),
    beta = as.numeric(b),
    se = as.numeric(se),
    test_statistic = statistic,
    df_test = df_test,
    p_method = "Design-based Wald t test (df = PSU - strata)",
    CI_low = as.numeric(b - critical * se),
    CI_high = as.numeric(b + critical * se),
    p = as.numeric(pvals),
    significance = sig_stars(as.numeric(pvals))
  )
}

# Factor conversion for common variables
factorize_common <- function(df) {
  df %>%
    mutate(
      sleep_trouble      = factor(sleep_trouble, levels = c("No", "Yes")),
      alcohol_category   = factor(alcohol_category, levels = c("Non-drinker", "Drinkers (non-binge)", "Binge drinker")),
      alcohol_user       = factor(alcohol_user, levels = c("No", "Yes")),
      alcohol_risk       = factor(alcohol_risk, levels = c("None", "Low/Moderate", "High")),
      age_group          = factor(age_group, levels = c("18-29", "30-44", "45-59", "60+")),
      gender             = factor(gender, levels = c("Female", "Male")),
      race_ethnicity     = factor(
        race_ethnicity,
        levels = c("Non-Hispanic White", "Non-Hispanic Black", "Non-Hispanic Asian",
                   "Mexican American", "Other Hispanic", "Other / Multi-racial")
      ),
      education          = factor(
        education,
        levels = c("College graduate+", "Some college/AA", "High school/GED",
                   "9-11th grade", "<9th grade")
      ),
      poverty_category   = factor(poverty_category),
      partnered          = factor(partnered, levels = c("No", "Yes")),
      smoking_status     = factor(smoking_status, levels = c("Never", "Former", "Current")),
      bmi_category       = factor(bmi_category, levels = c("Normal", "Underweight", "Overweight", "Obese")),
      diabetes           = factor(diabetes, levels = c("No", "Yes")),
      hypertension       = factor(hypertension, levels = c("No", "Yes")),
      heart_disease      = factor(heart_disease, levels = c("No", "Yes")),
      arthritis          = factor(arthritis, levels = c("No", "Yes")),
      daytime_sleepy     = factor(daytime_sleepy, levels = c("No", "Yes")),
      insufficient_sleep = factor(insufficient_sleep, levels = c("No", "Yes")),
      snoring_frequent   = factor(snoring_frequent, levels = c("No", "Yes")),
      breathing_pauses   = factor(breathing_pauses, levels = c("No", "Yes"))
    )
}



# Model-specific complete-case dataset
make_analysis_df <- function(df, vars_needed) {
  df %>%
    filter(
      if_all(all_of(vars_needed), ~ !is.na(.x)),
      is.finite(WTMEC2YR),
      WTMEC2YR > 0
    ) %>%
    factorize_common()
}

# NHANES survey design
make_design <- function(df) {
  svydesign(
    ids = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~WTMEC2YR,
    nest = TRUE,
    data = df
  )
}

# Export weighted proportions
export_weighted_prop <- function(design_obj, var_name, file_name) {
  est <- svymean(as.formula(paste0("~", var_name)), design_obj, na.rm = TRUE)
  df_test <- degf(design_obj)
  critical <- qt(0.975, df = df_test)
  estimate <- as.numeric(coef(est))
  se_value <- as.numeric(SE(est))
  
  out <- tibble(
    variable = var_name,
    category = clean_level_name(names(coef(est))),
    n_unweighted = nrow(design_obj$variables),
    estimate = estimate,
    se = se_value,
    CI_low = pmax(0, estimate - critical * se_value),
    CI_high = pmin(1, estimate + critical * se_value),
    design_df = df_test,
    p_applicability = "Descriptive proportion; hypothesis-test p-value not applicable"
  )
  
  write_csv(out, file.path(dir_tables, file_name))
  out
}

# Model 4 helper for secondary binary outcomes
run_binary_model4 <- function(outcome_var, df) {
  needed <- c(
    outcome_var,
    "alcohol_category",
    "WTMEC2YR", "SDMVPSU", "SDMVSTRA",
    "age_group", "gender", "race_ethnicity", "education",
    "poverty_ratio", "partnered", "DMDHHSIZ",
    "smoking_status", "bmi_category", "waist_circumference",
    "depression_score_no_sleep",
    "diabetes", "hypertension", "heart_disease", "arthritis"
  )
  
  needed <- needed[needed %in% names(df)]
  
  ana <- make_analysis_df(df, needed)
  des <- make_design(ana)
  
  f <- as.formula(
    paste0(
      "I(", outcome_var, " == 'Yes') ~ alcohol_category + ",
      "age_group + gender + race_ethnicity + education + poverty_ratio + partnered + DMDHHSIZ + ",
      "smoking_status + bmi_category + waist_circumference + ",
      "depression_score_no_sleep + diabetes + hypertension + heart_disease + arthritis"
    )
  )
  
  mod <- svyglm(f, design = des, family = quasibinomial())
  extract_or(mod, "Model 4", des, outcome_var)
}

############################################################
# 04 — IMPORT RAW NHANES FILES
############################################################

DEMO_J <- load_local_xpt("DEMO_J.XPT")
ALQ_J  <- load_local_xpt("ALQ_J.XPT")
SLQ_J  <- load_local_xpt("SLQ_J.XPT")
BMX_J  <- load_local_xpt("BMX_J.XPT")
DPQ_J  <- load_local_xpt("DPQ_J.XPT")

SMQ_J  <- load_local_xpt("SMQ_J.XPT")
MCQ_J  <- load_local_xpt("MCQ_J.XPT")
DIQ_J  <- load_local_xpt("DIQ_J.XPT")
BPQ_J  <- load_local_xpt("BPQ_J.XPT")

############################################################
# 05 — SELECT NEEDED VARIABLES AND MERGE
############################################################

DEMO_clean <- DEMO_J %>%
  select(
    SEQN, SDMVPSU, SDMVSTRA, WTMEC2YR,
    RIAGENDR, RIDAGEYR, RIDRETH3,
    DMDEDUC2, INDFMPIR, DMDMARTL, DMDHHSIZ
  )

ALQ_clean <- ALQ_J %>%
  select(
    SEQN, ALQ111, ALQ121, ALQ130, ALQ142, ALQ270, ALQ280, ALQ151
  )

SLQ_clean <- SLQ_J %>%
  select(
    SEQN, SLQ050, SLD012, SLD013, SLQ120, SLQ030, SLQ040
  )

BMX_clean <- BMX_J %>%
  select(
    SEQN, BMXBMI, BMXWAIST
  )

DPQ_clean <- DPQ_J %>%
  select(
    SEQN, DPQ010:DPQ090
  )

SMQ_clean <- SMQ_J %>%
  select(
    SEQN, SMQ020, SMQ040, SMD641, SMD650, SMQ078
  )

MCQ_clean <- MCQ_J %>%
  select(
    SEQN, MCQ160A, MCQ160B, MCQ160C, MCQ160D, MCQ160E
  )

DIQ_clean <- DIQ_J %>%
  select(
    SEQN, DIQ010
  )

BPQ_clean <- BPQ_J %>%
  select(
    SEQN, BPQ020
  )

step1_merged <- DEMO_clean %>%
  left_join(ALQ_clean, by = "SEQN") %>%
  left_join(SLQ_clean, by = "SEQN") %>%
  left_join(BMX_clean, by = "SEQN") %>%
  left_join(DPQ_clean, by = "SEQN") %>%
  left_join(SMQ_clean, by = "SEQN") %>%
  left_join(MCQ_clean, by = "SEQN") %>%
  left_join(DIQ_clean, by = "SEQN") %>%
  left_join(BPQ_clean, by = "SEQN") %>%
  filter(RIDAGEYR >= AGE_MIN)

saveRDS(step1_merged, file.path(dir_step, "step1_merged_adults.rds"))

############################################################
# 06 — VARIABLE-SPECIFIC CLEANING
############################################################

# IMPORTANT: NHANES missing codes are variable-specific. Never recode every
# 7, 9, 77, or 99 across the full data set; many are valid measurements/codes.
step2_clean <- step1_merged %>%
  mutate(across(-SEQN, to_num)) %>%
  mutate(
    # Demographics
    RIAGENDR = keep_codes(RIAGENDR, 1:2),
    RIDAGEYR = keep_range(RIDAGEYR, 18, 80),
    RIDRETH3 = keep_codes(RIDRETH3, c(1, 2, 3, 4, 6, 7)),
    DMDEDUC2 = keep_codes(DMDEDUC2, 1:5),
    INDFMPIR = keep_range(INDFMPIR, 0, 5),
    DMDMARTL = keep_codes(DMDMARTL, 1:6),
    DMDHHSIZ = keep_codes(DMDHHSIZ, 1:7),

    # Alcohol: 7 and 9 are valid for the frequency variables below.
    ALQ111 = keep_codes(ALQ111, 1:2),
    ALQ121 = keep_codes(ALQ121, 0:10),
    ALQ130 = keep_codes(ALQ130, 1:15),
    ALQ142 = keep_codes(ALQ142, 0:10),
    ALQ270 = keep_codes(ALQ270, 0:10),
    ALQ280 = keep_codes(ALQ280, 0:10),
    ALQ151 = keep_codes(ALQ151, 1:2),

    # Sleep
    SLQ050 = keep_codes(SLQ050, 1:2),
    SLD012 = keep_range(SLD012, 2, 14),
    SLD013 = keep_range(SLD013, 2, 14),
    SLQ030 = keep_codes(SLQ030, 0:3),
    SLQ040 = keep_codes(SLQ040, 0:3),
    SLQ120 = keep_codes(SLQ120, 0:4),

    # Smoking and chronic conditions
    SMQ020 = keep_codes(SMQ020, 1:2),
    SMQ040 = keep_codes(SMQ040, 1:3),
    MCQ160A = keep_codes(MCQ160A, 1:2),
    MCQ160B = keep_codes(MCQ160B, 1:2),
    MCQ160C = keep_codes(MCQ160C, 1:2),
    MCQ160D = keep_codes(MCQ160D, 1:2),
    MCQ160E = keep_codes(MCQ160E, 1:2),
    DIQ010 = keep_codes(DIQ010, 1:3),
    BPQ020 = keep_codes(BPQ020, 1:2)
  )

phq_items <- paste0("DPQ", sprintf("%03d", seq(10, 90, 10)))
phq_items <- phq_items[phq_items %in% names(step2_clean)]

step2_clean <- step2_clean %>%
  mutate(across(all_of(phq_items), ~ {
    v <- to_num(.x)
    v[!(v %in% 0:3)] <- NA
    v
  }))

# Fix NHANES pseudo-zero values in alcohol variables
fix_tiny_zero <- function(x, tol = 1e-12) {
  x <- to_num(x)
  x[!is.na(x) & abs(x) < tol] <- 0
  x
}

alq_zero_vars <- c("ALQ121", "ALQ142", "ALQ270", "ALQ280")
alq_zero_vars <- alq_zero_vars[alq_zero_vars %in% names(step2_clean)]

step2_clean <- step2_clean %>%
  mutate(across(all_of(alq_zero_vars), fix_tiny_zero))

saveRDS(step2_clean, file.path(dir_step, "step2_cleaned_na.rds"))

############################################################
# 07 — DERIVE ANALYSIS VARIABLES
############################################################

step3_derived <- step2_clean %>%
  mutate(
    age_years = RIDAGEYR,
    age_group = case_when(
      age_years >= 18 & age_years <= 29 ~ "18-29",
      age_years >= 30 & age_years <= 44 ~ "30-44",
      age_years >= 45 & age_years <= 59 ~ "45-59",
      age_years >= 60 ~ "60+",
      TRUE ~ NA_character_
    ),
    gender = case_when(
      RIAGENDR == 1 ~ "Male",
      RIAGENDR == 2 ~ "Female",
      TRUE ~ NA_character_
    ),
    race_ethnicity = case_when(
      RIDRETH3 == 1 ~ "Mexican American",
      RIDRETH3 == 2 ~ "Other Hispanic",
      RIDRETH3 == 3 ~ "Non-Hispanic White",
      RIDRETH3 == 4 ~ "Non-Hispanic Black",
      RIDRETH3 == 6 ~ "Non-Hispanic Asian",
      RIDRETH3 == 7 ~ "Other / Multi-racial",
      TRUE ~ NA_character_
    ),
    education = case_when(
      DMDEDUC2 == 1 ~ "<9th grade",
      DMDEDUC2 == 2 ~ "9-11th grade",
      DMDEDUC2 == 3 ~ "High school/GED",
      DMDEDUC2 == 4 ~ "Some college/AA",
      DMDEDUC2 == 5 ~ "College graduate+",
      TRUE ~ NA_character_
    ),
    poverty_ratio = INDFMPIR,
    poverty_category = case_when(
      !is.na(poverty_ratio) & poverty_ratio < 1  ~ "<1.0 (Below poverty)",
      !is.na(poverty_ratio) & poverty_ratio < 2  ~ "1.0-1.99",
      !is.na(poverty_ratio) & poverty_ratio < 4  ~ "2.0-3.99",
      !is.na(poverty_ratio) & poverty_ratio >= 4 ~ ">=4.0",
      TRUE ~ NA_character_
    ),
    marital_status = case_when(
      DMDMARTL == 1 ~ "Married",
      DMDMARTL == 2 ~ "Widowed",
      DMDMARTL == 3 ~ "Divorced",
      DMDMARTL == 4 ~ "Separated",
      DMDMARTL == 5 ~ "Never married",
      DMDMARTL == 6 ~ "Living with partner",
      TRUE ~ NA_character_
    ),
    partnered = case_when(
      DMDMARTL %in% c(1, 6) ~ "Yes",
      DMDMARTL %in% c(2, 3, 4, 5) ~ "No",
      TRUE ~ NA_character_
    ),
    bmi = BMXBMI,
    waist_circumference = BMXWAIST,
    bmi_category = case_when(
      !is.na(bmi) & bmi < 18.5 ~ "Underweight",
      !is.na(bmi) & bmi < 25   ~ "Normal",
      !is.na(bmi) & bmi < 30   ~ "Overweight",
      !is.na(bmi) & bmi >= 30  ~ "Obese",
      TRUE ~ NA_character_
    ),
    smoking_status = case_when(
      SMQ020 == 2 ~ "Never",
      SMQ020 == 1 & SMQ040 %in% c(1, 2) ~ "Current",
      SMQ020 == 1 & SMQ040 == 3 ~ "Former",
      TRUE ~ NA_character_
    ),
    depression_score = rowSums(across(all_of(phq_items)), na.rm = FALSE),
    depression_score_no_sleep = rowSums(across(all_of(phq_items)), na.rm = FALSE) - DPQ030,
    depressed = case_when(
      !is.na(depression_score) & depression_score >= 10 ~ "Yes",
      !is.na(depression_score) ~ "No",
      TRUE ~ NA_character_
    ),
    diabetes = case_when(
      DIQ010 == 1 ~ "Yes",
      DIQ010 == 2 ~ "No",
      TRUE ~ NA_character_
    ),
    hypertension = case_when(
      BPQ020 == 1 ~ "Yes",
      BPQ020 == 2 ~ "No",
      TRUE ~ NA_character_
    ),
    heart_disease = case_when(
      MCQ160B == 1 | MCQ160C == 1 | MCQ160D == 1 | MCQ160E == 1 ~ "Yes",
      MCQ160B == 2 & MCQ160C == 2 & MCQ160D == 2 & MCQ160E == 2 ~ "No",
      TRUE ~ NA_character_
    ),
    arthritis = case_when(
      MCQ160A == 1 ~ "Yes",
      MCQ160A == 2 ~ "No",
      TRUE ~ NA_character_
    ),
    sleep_trouble = case_when(
      SLQ050 == 1 ~ "Yes",
      SLQ050 == 2 ~ "No",
      TRUE ~ NA_character_
    ),
    sleep_hours = case_when(
      !is.na(SLD012) & !is.na(SLD013) ~ ((5 * SLD012) + (2 * SLD013)) / 7,
      !is.na(SLD012) & is.na(SLD013)  ~ SLD012,
      is.na(SLD012) & !is.na(SLD013)  ~ SLD013,
      TRUE ~ NA_real_
    ),
    insufficient_sleep = case_when(
      !is.na(sleep_hours) & sleep_hours < 7 ~ "Yes",
      !is.na(sleep_hours) ~ "No",
      TRUE ~ NA_character_
    ),
    snoring_frequent = case_when(
      SLQ030 == 3 ~ "Yes",
      SLQ030 %in% 0:2 ~ "No",
      TRUE ~ NA_character_
    ),
    breathing_pauses = case_when(
      SLQ040 %in% 1:3 ~ "Yes",
      SLQ040 == 0 ~ "No",
      TRUE ~ NA_character_
    ),
    daytime_sleepy = case_when(
      SLQ120 %in% c(3, 4) ~ "Yes",
      SLQ120 %in% 0:2 ~ "No",
      TRUE ~ NA_character_
    ),
    alcohol_user = case_when(
      ALQ111 == 1 ~ "Yes",
      ALQ111 == 2 ~ "No",
      TRUE ~ NA_character_
    ),
    drinks_per_day = ALQ130,
    drinking_frequency = ALQ121,
    # Ordered score used only for Spearman correlation: larger = more frequent.
    # ALQ121 code 1 is most frequent and code 10 is least frequent.
    drinking_frequency_rank = case_when(
      ALQ121 == 0 ~ 0,
      ALQ121 %in% 1:10 ~ 11 - ALQ121,
      TRUE ~ NA_real_
    ),
    binge_past_year = case_when(
      !is.na(ALQ142) & ALQ142 > 0 ~ "Yes",
      !is.na(ALQ142) & ALQ142 == 0 ~ "No",
      TRUE ~ NA_character_
    ),
    binge_niaaa = case_when(
      !is.na(ALQ270) & ALQ270 > 0 ~ "Yes",
      !is.na(ALQ270) & ALQ270 == 0 ~ "No",
      TRUE ~ NA_character_
    ),
    heavy_episodic = case_when(
      !is.na(ALQ280) & ALQ280 > 0 ~ "Yes",
      !is.na(ALQ280) & ALQ280 == 0 ~ "No",
      TRUE ~ NA_character_
    ),
    binge_any = case_when(
      binge_past_year == "Yes" | binge_niaaa == "Yes" | heavy_episodic == "Yes" ~ "Yes",
      rowSums(cbind(binge_past_year == "Yes", binge_niaaa == "Yes", heavy_episodic == "Yes"), na.rm = TRUE) == 0 &
        rowSums(cbind(!is.na(binge_past_year), !is.na(binge_niaaa), !is.na(heavy_episodic))) > 0 ~ "No",
      TRUE ~ NA_character_
    ),
    alcohol_category = case_when(
      alcohol_user == "No" ~ "Non-drinker",
      alcohol_user == "Yes" & binge_any == "Yes" ~ "Binge drinker",
      alcohol_user == "Yes" & binge_any == "No" ~ "Drinkers (non-binge)",
      TRUE ~ NA_character_
    ),
    alcohol_risk = case_when(
      alcohol_user == "No" ~ "None",
      alcohol_user == "Yes" & (binge_any == "Yes" | (!is.na(drinks_per_day) & drinks_per_day >= 4)) ~ "High",
      alcohol_user == "Yes" & binge_any == "No" ~ "Low/Moderate",
      TRUE ~ NA_character_
    )
  )

############################################################
# 08 — EXPORT FINAL MASTER ANALYSIS SHEET
############################################################

analysis_vars <- c(
  "SEQN", "SDMVPSU", "SDMVSTRA", "WTMEC2YR",
  "ALQ111", "ALQ121", "ALQ130", "ALQ142", "ALQ270", "ALQ280", "ALQ151",
  "alcohol_user", "alcohol_category", "alcohol_risk",
  "drinks_per_day", "drinking_frequency", "drinking_frequency_rank",
  "binge_past_year", "binge_niaaa", "heavy_episodic", "binge_any",
  "SLQ050", "SLD012", "SLD013", "SLQ120", "SLQ030", "SLQ040",
  "sleep_trouble", "sleep_hours", "insufficient_sleep", "daytime_sleepy",
  "snoring_frequent", "breathing_pauses",
  "RIAGENDR", "RIDAGEYR", "RIDRETH3", "DMDEDUC2", "INDFMPIR", "DMDMARTL", "DMDHHSIZ",
  "age_years", "age_group", "gender", "race_ethnicity", "education",
  "poverty_ratio", "poverty_category", "marital_status", "partnered",
  "BMXBMI", "BMXWAIST", "bmi", "bmi_category", "waist_circumference",
  "SMQ020", "SMQ040", "SMD641", "SMD650", "SMQ078", "smoking_status",
  phq_items, "depression_score", "depression_score_no_sleep", "depressed",
  "DIQ010", "BPQ020", "MCQ160A", "MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E",
  "diabetes", "hypertension", "heart_disease", "arthritis"
)

analysis_vars <- analysis_vars[analysis_vars %in% names(step3_derived)]

final_analysis_sheet <- step3_derived %>%
  select(all_of(analysis_vars))

saveRDS(final_analysis_sheet, file.path(dir_step, "final_analysis_sheet.rds"))
write_csv(final_analysis_sheet, file.path(dir_step, "final_analysis_sheet.csv"))

# Confirm that legitimate values previously lost by blanket recoding were kept.
valid_code_checks <- tribble(
  ~variable, ~value,
  "RIDAGEYR", 77,
  "RIDRETH3", 7,
  "DMDHHSIZ", 7,
  "ALQ121", 7,
  "ALQ121", 9,
  "ALQ130", 7,
  "ALQ130", 9,
  "ALQ142", 7,
  "ALQ142", 9,
  "ALQ270", 7,
  "ALQ270", 9,
  "ALQ280", 7,
  "ALQ280", 9,
  "SLD012", 7,
  "SLD012", 9,
  "SLD013", 7,
  "SLD013", 9,
  "SLQ030", 0,
  "SLQ040", 0,
  "SLQ120", 0
)

cleaning_qc <- bind_rows(lapply(seq_len(nrow(valid_code_checks)), function(i) {
  v <- valid_code_checks$variable[i]
  z <- valid_code_checks$value[i]
  tibble(
    variable = v,
    value = z,
    raw_count = sum(to_num(step1_merged[[v]]) == z, na.rm = TRUE),
    retained_count = sum(to_num(step2_clean[[v]]) == z, na.rm = TRUE)
  )
})) %>%
  mutate(status = if_else(raw_count == retained_count, "PASS", "FAIL"))

write_csv(cleaning_qc, file.path(dir_tables, "cleaning_qc_valid_values.csv"))

if (any(cleaning_qc$status == "FAIL")) {
  stop("Cleaning QC failed: at least one legitimate code was lost. See cleaning_qc_valid_values.csv")
}
if (anyDuplicated(final_analysis_sheet$SEQN) > 0) {
  stop("Duplicate SEQN values found in final_analysis_sheet.")
}

############################################################
# 09 — QC SUMMARY
############################################################

qc_file <- file.path(dir_logs, "qc_summary.txt")
sink(qc_file)

cat("=== QC SUMMARY ===\n")
cat("Cycle:", CYCLE, "\n")
cat("Adults N:", nrow(final_analysis_sheet), "\n\n")
cat("PHQ-9 summary:\n")
print(summary(final_analysis_sheet$depression_score))
cat("\nSleep hours summary:\n")
print(summary(final_analysis_sheet$sleep_hours))
cat("\nAlcohol user counts:\n")
print(table(final_analysis_sheet$alcohol_user, useNA = "ifany"))
cat("\nAlcohol category counts:\n")
print(table(final_analysis_sheet$alcohol_category, useNA = "ifany"))
cat("\nSleep trouble counts:\n")
print(table(final_analysis_sheet$sleep_trouble, useNA = "ifany"))
cat("\nDaytime sleepy counts:\n")
print(table(final_analysis_sheet$daytime_sleepy, useNA = "ifany"))

sink()
#####Run this once to confirm where the remaining NA comes from:

with(final_analysis_sheet, table(alcohol_user, is.na(alcohol_category), useNA = "ifany"))

###

final_analysis_sheet %>%
  filter(alcohol_user == "Yes", is.na(alcohol_category)) %>%
  summarise(
    n = n(),
    miss_ALQ142 = sum(is.na(ALQ142)),
    miss_ALQ270 = sum(is.na(ALQ270)),
    miss_ALQ280 = sum(is.na(ALQ280)),
    miss_all_binge = sum(is.na(ALQ142) & is.na(ALQ270) & is.na(ALQ280))
  )
############################################################
# 10 — MISSINGNESS TABLE
############################################################

missing_table <- function(df) {
  tibble(
    variable = names(df),
    n_missing = sapply(df, function(x) sum(is.na(x))),
    pct_missing = round(100 * n_missing / nrow(df), 2)
  ) %>%
    arrange(desc(pct_missing))
}

miss_tbl <- missing_table(final_analysis_sheet)
write_csv(miss_tbl, file.path(dir_tables, "missing_table.csv"))

############################################################
# 11 — FLOW COUNTS
############################################################

cov_model0 <- c()
cov_model1 <- c("age_group", "gender", "race_ethnicity")
cov_model2 <- c(cov_model1, "education", "poverty_ratio", "partnered", "DMDHHSIZ")
cov_model3 <- c(cov_model2, "smoking_status", "bmi_category", "waist_circumference")
cov_model4 <- c(cov_model3, "depression_score_no_sleep", "diabetes", "hypertension", "heart_disease", "arthritis")

needed_m0 <- c("sleep_trouble", "alcohol_category", "WTMEC2YR", "SDMVPSU", "SDMVSTRA")
needed_m1 <- c("sleep_trouble", "alcohol_category", cov_model1, "WTMEC2YR", "SDMVPSU", "SDMVSTRA")
needed_m2 <- c("sleep_trouble", "alcohol_category", cov_model2, "WTMEC2YR", "SDMVPSU", "SDMVSTRA")
needed_m3 <- c("sleep_trouble", "alcohol_category", cov_model3, "WTMEC2YR", "SDMVPSU", "SDMVSTRA")
needed_m4 <- c("sleep_trouble", "alcohol_category", cov_model4, "WTMEC2YR", "SDMVPSU", "SDMVSTRA")

flow_counts <- tibble(
  step = c(
    "N_adults",
    "N_positive_MEC_weight",
    "N_non_missing_outcome",
    "N_non_missing_exposure",
    "N_complete_model0",
    "N_complete_model1",
    "N_complete_model2",
    "N_complete_model3",
    "N_complete_model4"
  ),
  n = c(
    nrow(final_analysis_sheet),
    final_analysis_sheet %>% filter(!is.na(WTMEC2YR), WTMEC2YR > 0) %>% nrow(),
    final_analysis_sheet %>% filter(!is.na(WTMEC2YR), WTMEC2YR > 0, !is.na(sleep_trouble)) %>% nrow(),
    final_analysis_sheet %>% filter(!is.na(WTMEC2YR), WTMEC2YR > 0, !is.na(sleep_trouble), !is.na(alcohol_category)) %>% nrow(),
    nrow(make_analysis_df(final_analysis_sheet, needed_m0)),
    nrow(make_analysis_df(final_analysis_sheet, needed_m1)),
    nrow(make_analysis_df(final_analysis_sheet, needed_m2)),
    nrow(make_analysis_df(final_analysis_sheet, needed_m3)),
    nrow(make_analysis_df(final_analysis_sheet, needed_m4))
  )
)

write_csv(flow_counts, file.path(dir_tables, "flow_counts.csv"))

############################################################
# 12 — BUILD MODEL-SPECIFIC SURVEY DESIGNS
############################################################

analysis_m0 <- make_analysis_df(final_analysis_sheet, needed_m0)
analysis_m1 <- make_analysis_df(final_analysis_sheet, needed_m1)
analysis_m2 <- make_analysis_df(final_analysis_sheet, needed_m2)
analysis_m3 <- make_analysis_df(final_analysis_sheet, needed_m3)
analysis_m4 <- make_analysis_df(final_analysis_sheet, needed_m4)

des_m0 <- make_design(analysis_m0)
des_m1 <- make_design(analysis_m1)
des_m2 <- make_design(analysis_m2)
des_m3 <- make_design(analysis_m3)
des_m4 <- make_design(analysis_m4)

saveRDS(des_m0, file.path(dir_step, "survey_design_model0.rds"))
saveRDS(des_m1, file.path(dir_step, "survey_design_model1.rds"))
saveRDS(des_m2, file.path(dir_step, "survey_design_model2.rds"))
saveRDS(des_m3, file.path(dir_step, "survey_design_model3.rds"))
saveRDS(des_m4, file.path(dir_step, "survey_design_model4.rds"))

############################################################
# 13 — WEIGHTED DESCRIPTIVES
############################################################

des_alcohol <- make_design(
  make_analysis_df(final_analysis_sheet, c("alcohol_category", "WTMEC2YR", "SDMVPSU", "SDMVSTRA"))
)

des_sleep <- make_design(
  make_analysis_df(final_analysis_sheet, c("sleep_trouble", "WTMEC2YR", "SDMVPSU", "SDMVSTRA"))
)

des_prev <- make_design(
  make_analysis_df(final_analysis_sheet, c("alcohol_category", "sleep_trouble", "WTMEC2YR", "SDMVPSU", "SDMVSTRA"))
)

prop_alcohol <- export_weighted_prop(des_alcohol, "alcohol_category", "weighted_alcohol_props.csv")
prop_sleep   <- export_weighted_prop(des_sleep, "sleep_trouble", "weighted_sleep_props.csv")

############################################################
# 14 — TABLE 1 BY ALCOHOL CATEGORY
############################################################

# Use variable-specific available-case samples for Table 1. Requiring every
# Table 1 field simultaneously would unnecessarily discard many participants.
table1_continuous_one <- function(v) {
  needed <- c("alcohol_category", v, "WTMEC2YR", "SDMVPSU", "SDMVSTRA")
  ana <- make_analysis_df(final_analysis_sheet, needed)
  des <- make_design(ana)
  est <- as.data.frame(
    svyby(as.formula(paste0("~", v)), ~alcohol_category, des, svymean,
          vartype = c("se", "ci"), na.rm = TRUE)
  )
  se_col <- grep("^se", names(est), value = TRUE)[1]
  low_col <- grep("^ci_l", names(est), value = TRUE)[1]
  high_col <- grep("^ci_u", names(est), value = TRUE)[1]
  counts <- ana %>% count(alcohol_category, name = "n_unweighted")
  p_overall <- tryCatch(
    regTermTest(svyglm(as.formula(paste0(v, " ~ alcohol_category")), design = des),
                ~alcohol_category)$p,
    error = function(e) NA_real_
  )
  tibble(
    variable = v,
    alcohol_category = est$alcohol_category,
    estimate = est[[v]],
    se = est[[se_col]],
    CI_low = est[[low_col]],
    CI_high = est[[high_col]],
    p_overall = as.numeric(p_overall)
  ) %>%
    left_join(counts, by = "alcohol_category")
}

table1_cont_vars <- c(
  "age_years", "bmi", "waist_circumference", "poverty_ratio", "depression_score"
)
table1_cont <- bind_rows(lapply(table1_cont_vars, table1_continuous_one))
write_csv(table1_cont, file.path(dir_tables, "table1_continuous.csv"))

table1_cat_vars <- c(
  "age_group", "gender", "race_ethnicity", "education", "smoking_status",
  "bmi_category", "diabetes", "hypertension", "heart_disease", "arthritis",
  "sleep_trouble", "daytime_sleepy", "insufficient_sleep", "snoring_frequent", "breathing_pauses"
)

table1_cat_list <- lapply(table1_cat_vars, function(v) {
  needed <- c("alcohol_category", v, "WTMEC2YR", "SDMVPSU", "SDMVSTRA")
  ana <- make_analysis_df(final_analysis_sheet, needed)
  table1_design <- make_design(ana)
  est <- svyby(
    as.formula(paste0("~", v)),
    ~alcohol_category,
    table1_design,
    svymean,
    vartype = "se",
    na.rm = TRUE
  )
  
  out <- as.data.frame(est)
  cols_est <- names(out)[grepl(paste0("^", v), names(out))]
  cols_se  <- names(out)[grepl("^se\\.", names(out))]
  counts <- ana %>% count(alcohol_category, name = "n_unweighted")
  p_overall <- tryCatch(
    svychisq(as.formula(paste0("~", v, " + alcohol_category")),
             table1_design, statistic = "F")$p.value,
    error = function(e) NA_real_
  )
  
  res <- tibble()
  
  for (i in seq_along(cols_est)) {
    tmp <- tibble(
      variable = v,
      alcohol_category = out$alcohol_category,
      category = str_remove(cols_est[i], paste0("^", v)),
      estimate = out[[cols_est[i]]],
      se = out[[cols_se[i]]],
      p_overall = as.numeric(p_overall)
    )
    res <- bind_rows(res, tmp)
  }
  
  res %>% left_join(counts, by = "alcohol_category")
})

table1_cat <- bind_rows(table1_cat_list)
write_csv(table1_cat, file.path(dir_tables, "table1_categorical.csv"))

############################################################
# 15 — PRIMARY MODELS 0–4
############################################################

m0 <- svyglm(
  I(sleep_trouble == "Yes") ~ alcohol_category,
  design = des_m0,
  family = quasibinomial()
)

m1 <- svyglm(
  I(sleep_trouble == "Yes") ~ alcohol_category + age_group + gender + race_ethnicity,
  design = des_m1,
  family = quasibinomial()
)

m2 <- svyglm(
  I(sleep_trouble == "Yes") ~ alcohol_category + age_group + gender + race_ethnicity +
    education + poverty_ratio + partnered + DMDHHSIZ,
  design = des_m2,
  family = quasibinomial()
)

m3 <- svyglm(
  I(sleep_trouble == "Yes") ~ alcohol_category + age_group + gender + race_ethnicity +
    education + poverty_ratio + partnered + DMDHHSIZ +
    smoking_status + bmi_category + waist_circumference,
  design = des_m3,
  family = quasibinomial()
)

m4 <- svyglm(
  I(sleep_trouble == "Yes") ~ alcohol_category + age_group + gender + race_ethnicity +
    education + poverty_ratio + partnered + DMDHHSIZ +
    smoking_status + bmi_category + waist_circumference +
    depression_score_no_sleep + diabetes + hypertension + heart_disease + arthritis,
  design = des_m4,
  family = quasibinomial()
)

primary_results <- bind_rows(
  extract_or(m0, "Model 0", des_m0, "sleep_trouble"),
  extract_or(m1, "Model 1", des_m1, "sleep_trouble"),
  extract_or(m2, "Model 2", des_m2, "sleep_trouble"),
  extract_or(m3, "Model 3", des_m3, "sleep_trouble"),
  extract_or(m4, "Model 4", des_m4, "sleep_trouble")
)

write_csv(primary_results, file.path(dir_tables, "primary_sleep_trouble_models.csv"))

primary_overall_tests <- bind_rows(
  tibble(model = "Model 0", n = nrow(des_m0$variables), design_df = degf(des_m0),
         p_overall_alcohol = as.numeric(regTermTest(m0, ~alcohol_category, df = degf(des_m0))$p)),
  tibble(model = "Model 1", n = nrow(des_m1$variables), design_df = degf(des_m1),
         p_overall_alcohol = as.numeric(regTermTest(m1, ~alcohol_category, df = degf(des_m1))$p)),
  tibble(model = "Model 2", n = nrow(des_m2$variables), design_df = degf(des_m2),
         p_overall_alcohol = as.numeric(regTermTest(m2, ~alcohol_category, df = degf(des_m2))$p)),
  tibble(model = "Model 3", n = nrow(des_m3$variables), design_df = degf(des_m3),
         p_overall_alcohol = as.numeric(regTermTest(m3, ~alcohol_category, df = degf(des_m3))$p)),
  tibble(model = "Model 4", n = nrow(des_m4$variables), design_df = degf(des_m4),
         p_overall_alcohol = as.numeric(regTermTest(m4, ~alcohol_category, df = degf(des_m4))$p))
) %>%
  mutate(
    significance = sig_stars(p_overall_alcohol),
    p_method = "Design-based joint Wald F test (denominator df = PSU - strata)"
  )

write_csv(
  primary_overall_tests,
  file.path(dir_tables, "primary_alcohol_overall_pvalues.csv")
)

############################################################
# 16 — SECONDARY CONTINUOUS OUTCOME: SLEEP HOURS
############################################################

needed_sleep_hours <- c(
  "sleep_hours", "alcohol_category",
  "WTMEC2YR", "SDMVPSU", "SDMVSTRA",
  "age_group", "gender", "race_ethnicity", "education",
  "poverty_ratio", "partnered", "DMDHHSIZ",
  "smoking_status", "bmi_category", "waist_circumference",
  "depression_score_no_sleep", "diabetes", "hypertension", "heart_disease", "arthritis"
)

analysis_sleep_hours <- make_analysis_df(final_analysis_sheet, needed_sleep_hours)
design_sleep_hours   <- make_design(analysis_sleep_hours)

sleep_hours_model <- svyglm(
  sleep_hours ~ alcohol_category + age_group + gender + race_ethnicity +
    education + poverty_ratio + partnered + DMDHHSIZ +
    smoking_status + bmi_category + waist_circumference +
    depression_score_no_sleep + diabetes + hypertension + heart_disease + arthritis,
  design = design_sleep_hours
)

sleep_hours_results <- extract_beta(sleep_hours_model, "Model 4", design_sleep_hours, "sleep_hours")
write_csv(sleep_hours_results, file.path(dir_tables, "sleep_hours_model.csv"))

############################################################
# 17 — SECONDARY BINARY OUTCOMES
############################################################

secondary_binary_results <- bind_rows(
  run_binary_model4("insufficient_sleep", final_analysis_sheet),
  run_binary_model4("daytime_sleepy", final_analysis_sheet),
  run_binary_model4("snoring_frequent", final_analysis_sheet),
  run_binary_model4("breathing_pauses", final_analysis_sheet)
)

write_csv(secondary_binary_results, file.path(dir_tables, "secondary_binary_models.csv"))

############################################################
# 18 — DOSE-RESPONSE ANALYSES
############################################################

drinkers_df <- final_analysis_sheet %>%
  filter(alcohol_user == "Yes") %>%
  factorize_common()

needed_dose1 <- c(
  "sleep_hours", "drinks_per_day",
  "WTMEC2YR", "SDMVPSU", "SDMVSTRA",
  "age_group", "gender", "race_ethnicity", "education",
  "poverty_ratio", "partnered", "DMDHHSIZ",
  "smoking_status", "bmi_category", "waist_circumference",
  "depression_score_no_sleep", "diabetes", "hypertension", "heart_disease", "arthritis"
)

ana_dose1 <- drinkers_df %>%
  filter(if_all(all_of(needed_dose1), ~ !is.na(.x)), WTMEC2YR > 0)
des_dose1 <- make_design(ana_dose1)

mod_dose1 <- svyglm(
  sleep_hours ~ drinks_per_day + age_group + gender + race_ethnicity +
    education + poverty_ratio + partnered + DMDHHSIZ +
    smoking_status + bmi_category + waist_circumference +
    depression_score_no_sleep + diabetes + hypertension + heart_disease + arthritis,
  design = des_dose1
)

needed_dose2 <- c(
  "sleep_trouble", "drinks_per_day",
  "WTMEC2YR", "SDMVPSU", "SDMVSTRA",
  "age_group", "gender", "race_ethnicity", "education",
  "poverty_ratio", "partnered", "DMDHHSIZ",
  "smoking_status", "bmi_category", "waist_circumference",
  "depression_score_no_sleep", "diabetes", "hypertension", "heart_disease", "arthritis"
)

ana_dose2 <- drinkers_df %>%
  filter(if_all(all_of(needed_dose2), ~ !is.na(.x)), WTMEC2YR > 0)
des_dose2 <- make_design(ana_dose2)

mod_dose2 <- svyglm(
  I(sleep_trouble == "Yes") ~ drinks_per_day + age_group + gender + race_ethnicity +
    education + poverty_ratio + partnered + DMDHHSIZ +
    smoking_status + bmi_category + waist_circumference +
    depression_score_no_sleep + diabetes + hypertension + heart_disease + arthritis,
  design = des_dose2,
  family = quasibinomial()
)

dose_sleep_hours <- extract_beta(
  mod_dose1, "Dose-response", des_dose1, "sleep_hours_by_drinks_per_day"
) %>%
  transmute(
    outcome, model, n, term,
    coefficient = beta,
    effect = beta,
    effect_measure = "Beta (hours per additional drink)",
    se, test_statistic, df_test, p_method, CI_low, CI_high, p, significance
  )

dose_sleep_trouble <- extract_or(
  mod_dose2, "Dose-response", des_dose2, "sleep_trouble_by_drinks_per_day"
) %>%
  transmute(
    outcome, model, n, term,
    coefficient = estimate,
    effect = OR,
    effect_measure = "Odds ratio per additional drink",
    se, test_statistic, df_test, p_method, CI_low, CI_high, p, significance
  )

dose_response_results <- bind_rows(dose_sleep_hours, dose_sleep_trouble)

write_csv(dose_response_results, file.path(dir_tables, "dose_response_models.csv"))

############################################################
# 19 — SENSITIVITY ANALYSES
############################################################

needed_sens_bin <- c(
  "sleep_trouble", "alcohol_user",
  "WTMEC2YR", "SDMVPSU", "SDMVSTRA",
  "age_group", "gender", "race_ethnicity", "education",
  "poverty_ratio", "partnered", "DMDHHSIZ",
  "smoking_status", "bmi_category", "waist_circumference",
  "depression_score_no_sleep", "diabetes", "hypertension", "heart_disease", "arthritis"
)

ana_sens_bin <- make_analysis_df(final_analysis_sheet, needed_sens_bin)
des_sens_bin <- make_design(ana_sens_bin)

mod_sens_bin <- svyglm(
  I(sleep_trouble == "Yes") ~ alcohol_user + age_group + gender + race_ethnicity +
    education + poverty_ratio + partnered + DMDHHSIZ +
    smoking_status + bmi_category + waist_circumference +
    depression_score_no_sleep + diabetes + hypertension + heart_disease + arthritis,
  design = des_sens_bin,
  family = quasibinomial()
)

needed_sens_risk <- c(
  "sleep_trouble", "alcohol_risk",
  "WTMEC2YR", "SDMVPSU", "SDMVSTRA",
  "age_group", "gender", "race_ethnicity", "education",
  "poverty_ratio", "partnered", "DMDHHSIZ",
  "smoking_status", "bmi_category", "waist_circumference",
  "depression_score_no_sleep", "diabetes", "hypertension", "heart_disease", "arthritis"
)

ana_sens_risk <- make_analysis_df(final_analysis_sheet, needed_sens_risk)
des_sens_risk <- make_design(ana_sens_risk)

mod_sens_risk <- svyglm(
  I(sleep_trouble == "Yes") ~ alcohol_risk + age_group + gender + race_ethnicity +
    education + poverty_ratio + partnered + DMDHHSIZ +
    smoking_status + bmi_category + waist_circumference +
    depression_score_no_sleep + diabetes + hypertension + heart_disease + arthritis,
  design = des_sens_risk,
  family = quasibinomial()
)

needed_sens_phq <- c(
  "sleep_trouble", "alcohol_category", "depression_score",
  "WTMEC2YR", "SDMVPSU", "SDMVSTRA",
  "age_group", "gender", "race_ethnicity", "education",
  "poverty_ratio", "partnered", "DMDHHSIZ",
  "smoking_status", "bmi_category", "waist_circumference",
  "diabetes", "hypertension", "heart_disease", "arthritis"
)

ana_sens_phq <- make_analysis_df(final_analysis_sheet, needed_sens_phq)
des_sens_phq <- make_design(ana_sens_phq)

mod_sens_phq <- svyglm(
  I(sleep_trouble == "Yes") ~ alcohol_category + age_group + gender + race_ethnicity +
    education + poverty_ratio + partnered + DMDHHSIZ +
    smoking_status + bmi_category + waist_circumference +
    depression_score + diabetes + hypertension + heart_disease + arthritis,
  design = des_sens_phq,
  family = quasibinomial()
)

sensitivity_results <- bind_rows(
  extract_or(mod_sens_bin, "Sensitivity", des_sens_bin, "sleep_trouble_alcohol_user"),
  extract_or(mod_sens_risk, "Sensitivity", des_sens_risk, "sleep_trouble_alcohol_risk"),
  extract_or(mod_sens_phq, "Sensitivity", des_sens_phq, "sleep_trouble_PHQ_full")
)

write_csv(sensitivity_results, file.path(dir_tables, "sensitivity_models.csv"))

############################################################
# 20 — INTERACTION ANALYSIS
############################################################

needed_interaction <- c(
  "sleep_trouble", "alcohol_category", "gender",
  "WTMEC2YR", "SDMVPSU", "SDMVSTRA",
  "age_group", "race_ethnicity", "education",
  "poverty_ratio", "partnered", "DMDHHSIZ",
  "smoking_status", "bmi_category", "waist_circumference",
  "depression_score_no_sleep", "diabetes", "hypertension", "heart_disease", "arthritis"
)

ana_int <- make_analysis_df(final_analysis_sheet, needed_interaction)
des_int <- make_design(ana_int)

mod_int_gender <- svyglm(
  I(sleep_trouble == "Yes") ~ alcohol_category * gender + age_group +
    race_ethnicity + education + poverty_ratio + partnered + DMDHHSIZ +
    smoking_status + bmi_category + waist_circumference +
    depression_score_no_sleep + diabetes + hypertension + heart_disease + arthritis,
  design = des_int,
  family = quasibinomial()
)

interaction_results <- extract_or(mod_int_gender, "Interaction", des_int, "sleep_trouble_alcohol_x_gender")
write_csv(interaction_results, file.path(dir_tables, "interaction_models.csv"))

interaction_joint <- regTermTest(
  mod_int_gender,
  ~alcohol_category:gender,
  df = degf(des_int),
  method = "Wald"
)

interaction_joint_test <- tibble(
  test = "Overall alcohol-category-by-sex interaction",
  n = nrow(des_int$variables),
  numerator_df = 2,
  denominator_df = degf(des_int),
  p = as.numeric(interaction_joint$p),
  significance = sig_stars(as.numeric(interaction_joint$p)),
  p_method = "Design-based joint Wald F test"
)

write_csv(
  interaction_joint_test,
  file.path(dir_tables, "interaction_joint_pvalue.csv")
)

############################################################
# 21 — SPEARMAN CORRELATION MATRIX (EXPLORATORY)
############################################################

corr_df <- final_analysis_sheet %>%
  select(
    drinks_per_day,
    drinking_frequency_rank,
    sleep_hours,
    age_years,
    bmi,
    waist_circumference,
    poverty_ratio,
    depression_score_no_sleep
  )

corr_labels <- c(
  drinks_per_day = "Drinks per drinking day",
  drinking_frequency_rank = "Drinking frequency (rank)",
  sleep_hours = "Sleep duration",
  age_years = "Age",
  bmi = "BMI",
  waist_circumference = "Waist circumference",
  poverty_ratio = "Poverty ratio",
  depression_score_no_sleep = "PHQ score (sleep item excluded)"
)

spearman_one <- function(x, y, x_name, y_name) {
  ok <- complete.cases(x, y)
  n_pair <- sum(ok)
  if (n_pair < 3) {
    return(tibble(variable_1 = x_name, variable_2 = y_name,
                  n = n_pair, rho = NA_real_, p = NA_real_))
  }
  tst <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  tibble(
    variable_1 = x_name,
    variable_2 = y_name,
    n = n_pair,
    rho = unname(tst$estimate),
    p = tst$p.value
  )
}

format_p <- function(p) {
  ifelse(is.na(p), "p = NA",
         ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p)))
}

corr_long <- bind_rows(lapply(names(corr_df), function(v1) {
  bind_rows(lapply(names(corr_df), function(v2) {
    spearman_one(corr_df[[v1]], corr_df[[v2]], v1, v2)
  }))
})) %>%
  mutate(
    significance = case_when(
      variable_1 == variable_2 ~ "Not applicable",
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ "NS"
    ),
    label = if_else(variable_1 == variable_2, "1.00",
                    paste0(sprintf("%.2f", rho), significance))
  )

write_csv(corr_long, file.path(dir_tables, "spearman_correlation_long.csv"))

rho_matrix <- corr_long %>%
  select(variable_1, variable_2, rho) %>%
  pivot_wider(names_from = variable_2, values_from = rho)
p_matrix <- corr_long %>%
  select(variable_1, variable_2, p) %>%
  pivot_wider(names_from = variable_2, values_from = p)

write_csv(rho_matrix, file.path(dir_tables, "spearman_rho_matrix.csv"))
write_csv(p_matrix, file.path(dir_tables, "spearman_p_matrix.csv"))

corr_plot_df <- corr_long %>%
  mutate(
    label_1 = factor(corr_labels[variable_1], levels = rev(unname(corr_labels))),
    label_2 = factor(corr_labels[variable_2], levels = unname(corr_labels))
  )

p_corr <- ggplot(corr_plot_df, aes(x = label_2, y = label_1, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = label), size = 3.3, fontface = "bold") +
  scale_fill_gradient2(
    low = "#3B4CC0", mid = "white", high = "#B40426",
    midpoint = 0, limits = c(-1, 1), name = "Spearman rho"
  ) +
  coord_equal() +
  labs(
    title = "Spearman correlations among alcohol, sleep, and covariates",
    subtitle = "Pairwise-complete, unweighted exploratory analysis",
    x = NULL, y = NULL,
    caption = "* p < 0.05; ** p < 0.01; *** p < 0.001. PHQ sleep item excluded."
  ) +
  theme_minimal(base_size = 11, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    plot.margin = margin(12, 12, 12, 12)
  )

ggsave(file.path(dir_figures, "figure_spearman_correlation.png"), p_corr,
       width = 10, height = 8.5, dpi = 600, bg = "white")
ggsave(file.path(dir_figures, "figure_spearman_correlation.pdf"), p_corr,
       width = 10, height = 8.5,
       device = if (capabilities("cairo")) cairo_pdf else "pdf", bg = "white")

############################################################
# 22 — WEIGHTED PREVALENCE + FIGURES
############################################################

des_prev_plot <- update(des_prev, sleep_yes = ifelse(sleep_trouble == "Yes", 1, 0))

prev_sleep_by_alcohol <- svyby(
  ~sleep_yes,
  ~alcohol_category,
  des_prev_plot,
  svymean,
  vartype = "se",
  na.rm = TRUE
)

p_prev_overall <- tryCatch(
  svychisq(~sleep_trouble + alcohol_category, des_prev, statistic = "F")$p.value,
  error = function(e) NA_real_
)

prev_df <- as.data.frame(prev_sleep_by_alcohol) %>%
  rename(prevalence = sleep_yes, se = se) %>%
  mutate(
    prevalence = prevalence * 100,
    se = se * 100,
    critical = qt(0.975, df = degf(des_prev_plot)),
    CI_low = pmax(0, prevalence - critical * se),
    CI_high = pmin(100, prevalence + critical * se),
    p_overall = p_prev_overall,
    p_method = "Design-adjusted Rao-Scott F test",
    alcohol_category = factor(
      alcohol_category,
      levels = c("Non-drinker", "Drinkers (non-binge)", "Binge drinker"),
      labels = c(
        "Lifetime abstainers",
        "Ever drinkers without past-year binge drinking",
        "Ever drinkers with past-year binge drinking"
      )
    )
  )

write_csv(prev_df, file.path(dir_tables, "sleep_trouble_prevalence_by_alcohol.csv"))

# Pairwise comparisons of trouble-sleeping odds among all three alcohol groups.
pairwise_from_m0 <- function(comparison, contrast_values) {
  b <- coef(m0)
  V <- vcov(m0)
  L <- setNames(rep(0, length(b)), names(b))
  L[names(contrast_values)] <- contrast_values
  log_or <- sum(L * b)
  se_log_or <- sqrt(as.numeric(t(L) %*% V %*% L))
  df_test <- degf(des_m0)
  statistic <- log_or / se_log_or
  p_value <- 2 * pt(abs(statistic), df = df_test, lower.tail = FALSE)
  critical <- qt(0.975, df = df_test)
  tibble(
    comparison = comparison,
    n = nrow(des_m0$variables),
    log_OR = log_or,
    se_log_OR = se_log_or,
    OR = exp(log_or),
    CI_low = exp(log_or - critical * se_log_or),
    CI_high = exp(log_or + critical * se_log_or),
    test_statistic = statistic,
    df_test = df_test,
    p = p_value,
    significance = sig_stars(p_value),
    p_method = "Design-based Wald t test (df = PSU - strata)"
  )
}

prevalence_pairwise <- bind_rows(
  pairwise_from_m0(
    "Ever drinkers without past-year binge vs lifetime abstainers",
    c("alcohol_categoryDrinkers (non-binge)" = 1)
  ),
  pairwise_from_m0(
    "Ever drinkers with past-year binge vs lifetime abstainers",
    c("alcohol_categoryBinge drinker" = 1)
  ),
  pairwise_from_m0(
    "Ever drinkers with vs without past-year binge",
    c("alcohol_categoryDrinkers (non-binge)" = -1,
      "alcohol_categoryBinge drinker" = 1)
  )
)

write_csv(
  prevalence_pairwise,
  file.path(dir_tables, "sleep_trouble_prevalence_pairwise_pvalues.csv")
)

# Figure 1: weighted prevalence with design-based overall p-value
p1 <- ggplot(prev_df, aes(x = prevalence, y = alcohol_category, color = alcohol_category)) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.12, linewidth = 0.9) +
  geom_point(size = 3.6) +
  geom_text(aes(label = sprintf("%.1f%%", prevalence)), hjust = -0.35,
            color = "black", fontface = "bold", size = 4) +
  scale_color_manual(values = c("#7A7A7A", "#3F6FA8", "#B45F4A"), guide = "none") +
  scale_x_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0.03, 0.18))) +
  labs(
    title = "Survey-weighted prevalence of reported trouble sleeping",
    subtitle = paste("Overall association:", format_p(p_prev_overall)),
    x = "Survey-weighted prevalence (95% CI)",
    y = NULL,
    caption = "NHANES 2017–2018; estimates account for strata, clusters, and MEC examination weights."
  ) +
  theme_classic(base_size = 12, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    plot.margin = margin(12, 12, 12, 12)
  )

ggsave(
  filename = file.path(dir_figures, "figure1_sleep_trouble_by_alcohol.png"),
  plot = p1,
  width = 8,
  height = 5,
  dpi = 600,
  bg = "white"
)
ggsave(file.path(dir_figures, "figure1_sleep_trouble_by_alcohol.pdf"), p1,
       width = 8.5, height = 5.5,
       device = if (capabilities("cairo")) cairo_pdf else "pdf", bg = "white")

# Figure 2: forest plot of primary alcohol ORs
forest_df <- primary_results %>%
  filter(term != "(Intercept)", grepl("alcohol_category", term)) %>%
  mutate(
    exposure = case_when(
      str_detect(term, "Drinkers \\(non-binge\\)") ~ "Without past-year binge drinking",
      str_detect(term, "Binge drinker") ~ "With past-year binge drinking",
      TRUE ~ clean_level_name(term)
    ),
    model = factor(model, levels = paste("Model", 0:4)),
    significance_group = factor(if_else(p < 0.05, "p < 0.05", "p >= 0.05"),
                                levels = c("p < 0.05", "p >= 0.05")),
    result_label = paste0(
      sprintf("OR %.2f (%.2f–%.2f); ", OR, CI_low, CI_high), format_p(p)
    )
  )

p2 <- ggplot(forest_df, aes(x = OR, y = model, color = exposure, shape = significance_group)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey35", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.15,
                 position = position_dodge(width = 0.48), linewidth = 0.85) +
  geom_point(size = 3.2, stroke = 1, position = position_dodge(width = 0.48)) +
  geom_text(aes(x = CI_high, label = result_label), hjust = -0.04,
            color = "black", fontface = "bold", size = 3.2,
            position = position_dodge(width = 0.48), show.legend = FALSE) +
  scale_x_log10(breaks = c(0.5, 1, 2, 4),
                expand = expansion(mult = c(0.05, 0.62))) +
  scale_color_manual(values = c(
    "Without past-year binge drinking" = "#3F6FA8",
    "With past-year binge drinking" = "#B45F4A"
  )) +
  scale_shape_manual(values = c("p < 0.05" = 16, "p >= 0.05" = 1)) +
  labs(
    title = "Alcohol-use pattern and reported trouble sleeping",
    subtitle = "Survey-weighted logistic regression; lifetime abstainers are the reference group",
    x = "Odds ratio (log scale)",
    y = NULL,
    color = NULL,
    shape = NULL,
    caption = "Points are odds ratios; horizontal lines are design-based 95% confidence intervals."
  ) +
  theme_classic(base_size = 12, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", color = "black"),
    legend.position = "bottom",
    legend.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    plot.margin = margin(12, 12, 12, 12)
  )

ggsave(
  filename = file.path(dir_figures, "figure2_forest_primary_models.png"),
  plot = p2,
  width = 10,
  height = 6,
  dpi = 600,
  bg = "white"
)
ggsave(file.path(dir_figures, "figure2_forest_primary_models.pdf"), p2,
       width = 12, height = 7,
       device = if (capabilities("cairo")) cairo_pdf else "pdf", bg = "white")

# Figure 3: survey-weighted alcohol-use distribution (descriptive; no p-value needed)
critical_alcohol <- qt(0.975, df = degf(des_alcohol))
alcohol_plot_df <- prop_alcohol %>%
  mutate(
    estimate = 100 * estimate,
    se = 100 * se,
    CI_low = pmax(0, estimate - critical_alcohol * se),
    CI_high = pmin(100, estimate + critical_alcohol * se),
    category = factor(
      category,
      levels = c("Non-drinker", "Drinkers (non-binge)", "Binge drinker"),
      labels = c(
        "Lifetime abstainers",
        "Ever drinkers without past-year binge drinking",
        "Ever drinkers with past-year binge drinking"
      )
    )
  )

p3 <- ggplot(alcohol_plot_df, aes(x = category, y = estimate, fill = category)) +
  geom_col(width = 0.66) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.14, linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", estimate)), vjust = -0.6,
            fontface = "bold", size = 4) +
  scale_fill_manual(values = c("#7A7A7A", "#3F6FA8", "#B45F4A"), guide = "none") +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Survey-weighted distribution of alcohol-use patterns",
    x = NULL, y = "Survey-weighted percentage (95% CI)",
    caption = "Descriptive distribution; inferential p-values are not applicable."
  ) +
  theme_classic(base_size = 12, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", color = "black"),
    axis.text.x = element_text(angle = 15, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    plot.margin = margin(12, 12, 12, 12)
  )

ggsave(file.path(dir_figures, "figure3_weighted_alcohol_distribution.png"), p3,
       width = 9, height = 6, dpi = 600, bg = "white")
ggsave(file.path(dir_figures, "figure3_weighted_alcohol_distribution.pdf"), p3,
       width = 9, height = 6,
       device = if (capabilities("cairo")) cairo_pdf else "pdf", bg = "white")

# Figure 4: guaranteed participant-flow/complete-case figure
flow_plot_df <- flow_counts %>%
  mutate(
    step_label = recode(
      step,
      N_adults = "Adults aged 18 years or older",
      N_positive_MEC_weight = "Positive MEC examination weight",
      N_non_missing_outcome = "Nonmissing trouble-sleeping outcome",
      N_non_missing_exposure = "Nonmissing outcome and alcohol pattern",
      N_complete_model0 = "Complete cases: Model 0",
      N_complete_model1 = "Complete cases: Model 1",
      N_complete_model2 = "Complete cases: Model 2",
      N_complete_model3 = "Complete cases: Model 3",
      N_complete_model4 = "Complete cases: Model 4"
    ),
    step_label = factor(step_label, levels = rev(step_label))
  )

p4 <- ggplot(flow_plot_df, aes(x = n, y = step_label)) +
  geom_col(fill = "#3F6FA8", width = 0.68) +
  geom_text(aes(label = scales::comma(n)), hjust = -0.15,
            fontface = "bold", size = 4) +
  scale_x_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.13))) +
  labs(
    title = "Analytic-sample flow across sequential models",
    x = "Participants, n", y = NULL,
    caption = "Survey models include participants with positive WTMEC2YR and complete model-specific data."
  ) +
  theme_classic(base_size = 12, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    plot.margin = margin(12, 12, 12, 12)
  )

ggsave(file.path(dir_figures, "figure4_participant_flow.png"), p4,
       width = 10, height = 7, dpi = 600, bg = "white")
ggsave(file.path(dir_figures, "figure4_participant_flow.pdf"), p4,
       width = 10, height = 7,
       device = if (capabilities("cairo")) cairo_pdf else "pdf", bg = "white")

# Optional diagram version of the participant flow
if (requireNamespace("DiagrammeR", quietly = TRUE) &&
    requireNamespace("htmlwidgets", quietly = TRUE)) {
  
  library(DiagrammeR)
  
  flow_list <- setNames(flow_counts$n, flow_counts$step)
  
  dot <- sprintf(
    "
    digraph flow {
      graph [rankdir=TB, layout=dot, splines=ortho, nodesep=0.35, ranksep=0.45]
      node [shape=box, style='rounded,filled', fillcolor=white, color=gray30, fontname=Helvetica, fontsize=11]
      edge [color=gray40]

      a [label='NHANES %s\\nAdults (age >= %s)\\nN = %s']
      b [label='Non-missing outcome\\n(sleep_trouble)\\nN = %s']
      c [label='Non-missing exposure\\n(alcohol_category)\\nN = %s']
      d [label='Complete Model 0\\nN = %s']
      e [label='Complete Model 1\\nN = %s']
      f [label='Complete Model 2\\nN = %s']
      g [label='Complete Model 3\\nN = %s']
      h [label='Complete Model 4\\nN = %s']

      a -> b -> c -> d -> e -> f -> g -> h
    }
    ",
    CYCLE, AGE_MIN,
    flow_list[["N_adults"]],
    flow_list[["N_non_missing_outcome"]],
    flow_list[["N_non_missing_exposure"]],
    flow_list[["N_complete_model0"]],
    flow_list[["N_complete_model1"]],
    flow_list[["N_complete_model2"]],
    flow_list[["N_complete_model3"]],
    flow_list[["N_complete_model4"]]
  )
  
  gr <- grViz(dot)
  
  html_path <- file.path(dir_figures, "figure_optional_flowchart.html")
  htmlwidgets::saveWidget(gr, file = html_path, selfcontained = TRUE)
  
  if (requireNamespace("DiagrammeRsvg", quietly = TRUE) &&
      requireNamespace("rsvg", quietly = TRUE)) {
    library(DiagrammeRsvg)
    library(rsvg)
    svg <- export_svg(gr)
    rsvg_png(charToRaw(svg), file.path(dir_figures, "figure_optional_flowchart.png"))
  }
}

############################################################
# 23 — FINAL STATUS MESSAGE
############################################################

cat("\n================ FINAL STATUS ================\n")
cat("Analysis completed successfully.\n\n")

cat("Master sheet:\n")
cat(" - ", file.path(dir_step, "final_analysis_sheet.csv"), "\n", sep = "")
cat(" - ", file.path(dir_step, "final_analysis_sheet.rds"), "\n", sep = "")

cat("\nMain tables:\n")
cat(" - ", file.path(dir_logs, "qc_summary.txt"), "\n", sep = "")
cat(" - ", file.path(dir_tables, "cleaning_qc_valid_values.csv"), "\n", sep = "")
cat(" - ", file.path(dir_tables, "missing_table.csv"), "\n", sep = "")
cat(" - ", file.path(dir_tables, "flow_counts.csv"), "\n", sep = "")
cat(" - ", file.path(dir_tables, "weighted_alcohol_props.csv"), "\n", sep = "")
cat(" - ", file.path(dir_tables, "weighted_sleep_props.csv"), "\n", sep = "")
cat(" - ", file.path(dir_tables, "table1_continuous.csv"), "\n", sep = "")
cat(" - ", file.path(dir_tables, "table1_categorical.csv"), "\n", sep = "")
cat(" - ", file.path(dir_tables, "primary_sleep_trouble_models.csv"), "\n", sep = "")
cat(" - ", file.path(dir_tables, "primary_alcohol_overall_pvalues.csv"), "\n", sep = "")
cat(" - ", file.path(dir_tables, "sleep_hours_model.csv"), "\n", sep = "")
cat(" - ", file.path(dir_tables, "secondary_binary_models.csv"), "\n", sep = "")
cat(" - ", file.path(dir_tables, "dose_response_models.csv"), "\n", sep = "")
cat(" - ", file.path(dir_tables, "sensitivity_models.csv"), "\n", sep = "")
cat(" - ", file.path(dir_tables, "interaction_models.csv"), "\n", sep = "")
cat(" - ", file.path(dir_tables, "interaction_joint_pvalue.csv"), "\n", sep = "")
cat(" - ", file.path(dir_tables, "spearman_correlation_long.csv"), "\n", sep = "")
cat(" - ", file.path(dir_tables, "spearman_rho_matrix.csv"), "\n", sep = "")
cat(" - ", file.path(dir_tables, "spearman_p_matrix.csv"), "\n", sep = "")
cat(" - ", file.path(dir_tables, "sleep_trouble_prevalence_by_alcohol.csv"), "\n", sep = "")
cat(" - ", file.path(dir_tables, "sleep_trouble_prevalence_pairwise_pvalues.csv"), "\n", sep = "")

cat("\nFigures:\n")
cat(" - ", file.path(dir_figures, "figure1_sleep_trouble_by_alcohol.png"), "\n", sep = "")
cat(" - ", file.path(dir_figures, "figure2_forest_primary_models.png"), "\n", sep = "")
cat(" - ", file.path(dir_figures, "figure3_weighted_alcohol_distribution.png"), "\n", sep = "")
cat(" - ", file.path(dir_figures, "figure4_participant_flow.png"), "\n", sep = "")
cat(" - ", file.path(dir_figures, "figure_spearman_correlation.png"), "\n", sep = "")

cat("==============================================\n")


####Additoional check 
sum(!is.na(final_analysis_sheet$alcohol_category))
sum(!is.na(final_analysis_sheet$sleep_trouble))
sum(!is.na(final_analysis_sheet$sleep_trouble) & !is.na(final_analysis_sheet$alcohol_category))


final_analysis_sheet %>%
  filter(!is.na(alcohol_category), is.na(sleep_trouble)) %>%
  select(SEQN, alcohol_category, sleep_trouble)

############################################################
# FINAL READABLE PUBLICATION FIGURES
# NHANES 2017–2018: alcohol-use pattern and sleep
#
# Run this AFTER the corrected analysis has created the CSV
# result tables. It does not refit or change any model.
#
# Figure order follows first citation in the Results:
#   Figure 1: participant flow
#   Figure 2: weighted prevalence
#   Figure 3: primary forest plot
#   Figure 4: exploratory Spearman correlations
############################################################

rm(list = ls())

required_packages <- c(
  "dplyr", "ggplot2", "patchwork", "readr",
  "scales", "stringr", "tibble", "tidyr"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(scales)
  library(stringr)
  library(tibble)
  library(tidyr)
})

############################################################
# 1. FOLDERS
############################################################

default_base_dir <- "E:/Neural Vista/projct alcohol and sleep"
base_dir <- Sys.getenv(
  "SLEEP_ALCOHOL_PROJECT_DIR",
  unset = default_base_dir
)

tables_dir <- Sys.getenv(
  "SLEEP_ALCOHOL_TABLES_DIR",
  unset = file.path(base_dir, "outputs", "tables")
)

figures_dir <- Sys.getenv(
  "SLEEP_ALCOHOL_FIGURES_DIR",
  unset = file.path(
    base_dir, "outputs", "figures", "publication_final"
  )
)

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  "flow_counts.csv",
  "sleep_trouble_prevalence_by_alcohol.csv",
  "primary_sleep_trouble_models.csv",
  "primary_alcohol_overall_pvalues.csv",
  "spearman_correlation_long.csv"
)

missing_files <- required_files[
  !file.exists(file.path(tables_dir, required_files))
]

if (length(missing_files) > 0) {
  stop(
    "Missing result file(s) in ", tables_dir, ":\n",
    paste0(" - ", missing_files, collapse = "\n")
  )
}

############################################################
# 2. PUBLICATION STYLE
############################################################

base_font <- "Times New Roman"

# Colour-blind-friendly and print-readable palette.
colour_abstainer <- "#666666"
colour_no_binge <- "#0072B2"
colour_binge <- "#D55E00"
colour_text <- "#111111"
colour_grid <- "#D8D8D8"
colour_band <- "#F2F5F7"

p_text <- function(p, include_p = FALSE) {
  answer <- dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
  if (include_p) paste0("p = ", answer) else answer
}

theme_publication <- function(base_size = 11.5) {
  theme_classic(base_size = base_size, base_family = base_font) +
    theme(
      text = element_text(colour = colour_text),
      plot.title = element_text(
        face = "bold", size = base_size + 2.5,
        margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        face = "bold", size = base_size,
        margin = margin(b = 9)
      ),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = colour_text),
      panel.border = element_rect(
        colour = "black", fill = NA, linewidth = 0.55
      ),
      plot.background = element_rect(
        colour = "black", fill = "white", linewidth = 0.55
      ),
      plot.margin = margin(10, 12, 10, 12),
      legend.position = "none"
    )
}

save_figure <- function(plot_object, stem, width, height) {
  png_path <- file.path(figures_dir, paste0(stem, ".png"))
  pdf_path <- file.path(figures_dir, paste0(stem, ".pdf"))
  
  ggsave(
    filename = png_path,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )
  
  pdf_device <- if (capabilities("cairo")) cairo_pdf else "pdf"
  
  ggsave(
    filename = pdf_path,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    device = pdf_device,
    bg = "white",
    limitsize = FALSE
  )
  
  message("Saved: ", png_path)
  message("Saved: ", pdf_path)
}

############################################################
# 3. FIGURE 1 — PARTICIPANT FLOW
############################################################

flow_raw <- read_csv(
  file.path(tables_dir, "flow_counts.csv"),
  show_col_types = FALSE
)

flow_lookup <- setNames(flow_raw$n, flow_raw$step)

needed_flow <- c(
  "N_adults", "N_positive_MEC_weight",
  "N_non_missing_outcome", "N_non_missing_exposure",
  "N_complete_model2", "N_complete_model3",
  "N_complete_model4"
)

if (!all(needed_flow %in% names(flow_lookup))) {
  stop("flow_counts.csv does not contain all required rows.")
}

flow_data <- tribble(
  ~stage, ~n,
  "Adults aged 18 years or older", flow_lookup[["N_adults"]],
  "Positive MEC examination weight", flow_lookup[["N_positive_MEC_weight"]],
  "Nonmissing trouble-sleeping outcome\nafter weight restriction", flow_lookup[["N_non_missing_outcome"]],
  "Primary outcome and alcohol pattern\n(complete cases for Models 0–1)", flow_lookup[["N_non_missing_exposure"]],
  "Complete cases: Model 2", flow_lookup[["N_complete_model2"]],
  "Complete cases: Model 3", flow_lookup[["N_complete_model3"]],
  "Complete cases: Model 4", flow_lookup[["N_complete_model4"]]
) %>%
  mutate(
    stage = factor(stage, levels = rev(stage)),
    n_label = comma(n)
  )

figure1 <- ggplot(flow_data, aes(x = n, y = stage)) +
  geom_col(width = 0.64, fill = colour_no_binge) +
  geom_text(
    aes(label = n_label),
    hjust = -0.10,
    family = base_font,
    fontface = "bold",
    size = 4.0
  ) +
  scale_x_continuous(
    limits = c(0, 6500),
    breaks = c(0, 2000, 4000, 6000),
    labels = comma,
    expand = c(0, 0)
  ) +
  labs(
    title = "Analytic-sample flow across sequential models",
    x = "Participants, n",
    y = NULL
  ) +
  theme_publication(base_size = 11.5) +
  theme(
    axis.text.y = element_text(
      face = "bold", size = 10.0,
      lineheight = 0.95
    ),
    panel.grid.major.x = element_line(
      colour = colour_grid, linewidth = 0.38
    ),
    panel.grid.minor = element_blank()
  )

save_figure(
  figure1,
  "Figure_1_participant_flow",
  width = 7.4,
  height = 5.5
)

############################################################
# 4. FIGURE 2 — WEIGHTED PREVALENCE
############################################################

prevalence_data <- read_csv(
  file.path(tables_dir, "sleep_trouble_prevalence_by_alcohol.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    short_group = case_when(
      alcohol_category == "Lifetime abstainers" ~
        "Lifetime abstainers",
      str_detect(alcohol_category, "without") ~
        "Ever drinkers without\npast-year binge drinking",
      TRUE ~
        "Ever drinkers with\npast-year binge drinking"
    ),
    short_group = factor(
      short_group,
      levels = c(
        "Ever drinkers with\npast-year binge drinking",
        "Ever drinkers without\npast-year binge drinking",
        "Lifetime abstainers"
      )
    ),
    group_colour = case_when(
      alcohol_category == "Lifetime abstainers" ~ colour_abstainer,
      str_detect(alcohol_category, "without") ~ colour_no_binge,
      TRUE ~ colour_binge
    ),
    result_label = sprintf(
      "%.1f%% (%.1f–%.1f)", prevalence, CI_low, CI_high
    )
  )

overall_p <- unique(prevalence_data$p_overall)
if (length(overall_p) != 1) {
  stop("Expected exactly one overall Rao–Scott p-value.")
}

figure2 <- ggplot(prevalence_data, aes(y = short_group)) +
  geom_errorbarh(
    aes(xmin = CI_low, xmax = CI_high),
    height = 0.13,
    linewidth = 0.85,
    colour = prevalence_data$group_colour
  ) +
  geom_point(
    aes(x = prevalence),
    size = 4.0,
    colour = prevalence_data$group_colour
  ) +
  geom_text(
    aes(x = CI_high + 0.6, label = result_label),
    hjust = 0,
    family = base_font,
    fontface = "bold",
    size = 3.8
  ) +
  scale_x_continuous(
    limits = c(10, 43),
    breaks = c(10, 20, 30, 40),
    labels = label_percent(scale = 1, accuracy = 1),
    expand = c(0, 0)
  ) +
  labs(
    title = "Survey-weighted prevalence of reported trouble sleeping",
    subtitle = paste0(
      "Overall Rao–Scott ", p_text(overall_p, include_p = TRUE)
    ),
    x = "Survey-weighted prevalence (95% CI)",
    y = NULL
  ) +
  theme_publication(base_size = 11.5) +
  theme(
    axis.text.y = element_text(
      face = "bold", size = 10.0,
      lineheight = 0.95
    ),
    panel.grid.major.x = element_line(
      colour = colour_grid, linewidth = 0.38
    ),
    panel.grid.minor = element_blank()
  )

save_figure(
  figure2,
  "Figure_2_weighted_prevalence",
  width = 7.4,
  height = 4.4
)

############################################################
# 5. FIGURE 3 — PRIMARY FOREST PLOT
#
# The figure uses separate label, CI, and numeric panels.
# It is exported at manuscript width, so Word does not have
# to shrink a 13-inch graphic into a 6.5-inch text column.
############################################################

primary <- read_csv(
  file.path(tables_dir, "primary_sleep_trouble_models.csv"),
  show_col_types = FALSE
)

primary_overall <- read_csv(
  file.path(tables_dir, "primary_alcohol_overall_pvalues.csv"),
  show_col_types = FALSE
)

model_levels <- paste("Model", 0:4)

forest_data <- primary %>%
  filter(str_detect(term, "^alcohol_category")) %>%
  mutate(
    model_index = match(model, model_levels) - 1,
    exposure = case_when(
      str_detect(term, "Drinkers \\(non-binge\\)") ~
        "No past-year binge",
      str_detect(term, "Binge drinker") ~
        "Past-year binge",
      TRUE ~ term
    ),
    exposure_order = if_else(exposure == "No past-year binge", 0, 1),
    y = 14 - 3 * model_index - exposure_order,
    colour = if_else(
      exposure == "No past-year binge",
      colour_no_binge,
      colour_binge
    ),
    result_text = sprintf("%.2f (%.2f–%.2f)", OR, CI_low, CI_high),
    estimate_p = p_text(p)
  ) %>%
  left_join(
    primary_overall %>%
      select(model, p_overall_alcohol),
    by = "model"
  ) %>%
  mutate(
    overall_p_text = if_else(
      exposure_order == 0,
      p_text(p_overall_alcohol),
      ""
    )
  ) %>%
  arrange(model_index, exposure_order)

if (nrow(forest_data) != 10) {
  stop("Expected 10 alcohol estimates in the primary result table.")
}

model_rows <- forest_data %>%
  group_by(model, model_index, n) %>%
  summarise(y = mean(y), .groups = "drop") %>%
  mutate(model_label = paste0(str_replace(model, "Model ", "M"), "\nn = ", comma(n)))

band_data <- tibble(
  model_index = 0:4,
  ymin = 12.55 - 3 * model_index,
  ymax = 14.45 - 3 * model_index,
  fill = if_else(model_index %% 2 == 0, colour_band, "white")
)

y_limits <- c(0.35, 15.45)

label_panel <- ggplot() +
  geom_rect(
    data = band_data,
    aes(xmin = 0, xmax = 2.65, ymin = ymin, ymax = ymax),
    fill = band_data$fill,
    colour = NA
  ) +
  geom_text(
    data = model_rows,
    aes(x = 0.02, y = y, label = model_label),
    hjust = 0,
    family = base_font,
    fontface = "bold",
    size = 3.25,
    lineheight = 0.92
  ) +
  geom_point(
    data = forest_data,
    aes(x = 0.82, y = y),
    shape = 15,
    size = 2.8,
    colour = forest_data$colour
  ) +
  geom_text(
    data = forest_data,
    aes(x = 0.98, y = y, label = exposure),
    hjust = 0,
    family = base_font,
    size = 3.15
  ) +
  annotate(
    "text", x = 0.02, y = 15.12, label = "Model",
    hjust = 0, family = base_font, fontface = "bold", size = 3.3
  ) +
  annotate(
    "text", x = 0.82, y = 15.12, label = "Alcohol-use group",
    hjust = 0, family = base_font, fontface = "bold", size = 3.3
  ) +
  scale_x_continuous(limits = c(0, 2.65), expand = c(0, 0)) +
  scale_y_continuous(limits = y_limits, expand = c(0, 0)) +
  theme_void(base_family = base_font) +
  theme(plot.margin = margin(5, 2, 5, 5))

ci_panel <- ggplot(forest_data, aes(y = y)) +
  geom_rect(
    data = band_data,
    aes(xmin = 0.48, xmax = 4.5, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = band_data$fill,
    colour = NA
  ) +
  geom_vline(
    xintercept = 1,
    linetype = 2,
    colour = "#555555",
    linewidth = 0.55
  ) +
  geom_errorbarh(
    aes(xmin = CI_low, xmax = CI_high),
    height = 0.14,
    linewidth = 0.75,
    colour = forest_data$colour
  ) +
  geom_point(
    aes(x = OR),
    size = 3.0,
    colour = forest_data$colour
  ) +
  scale_x_log10(
    limits = c(0.48, 4.5),
    breaks = c(0.5, 1, 2, 4),
    labels = c("0.5", "1", "2", "4"),
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = y_limits, expand = c(0, 0)) +
  labs(x = "Odds ratio (log scale)", y = NULL) +
  theme_classic(base_size = 10.5, base_family = base_font) +
  theme(
    axis.title.x = element_text(face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.grid.major.x = element_line(
      colour = colour_grid, linewidth = 0.35
    ),
    plot.margin = margin(5, 3, 5, 3)
  )

numeric_panel <- ggplot() +
  geom_rect(
    data = band_data,
    aes(xmin = 0, xmax = 4.35, ymin = ymin, ymax = ymax),
    fill = band_data$fill,
    colour = NA
  ) +
  geom_text(
    data = forest_data,
    aes(x = 0.02, y = y, label = result_text),
    hjust = 0,
    family = base_font,
    size = 3.10
  ) +
  geom_text(
    data = forest_data,
    aes(x = 2.72, y = y, label = estimate_p),
    hjust = 0,
    family = base_font,
    fontface = "bold",
    size = 3.10
  ) +
  geom_text(
    data = forest_data,
    aes(x = 3.60, y = y, label = overall_p_text),
    hjust = 0,
    family = base_font,
    fontface = "bold",
    size = 3.10
  ) +
  annotate(
    "text", x = 0.02, y = 15.12, label = "OR (95% CI)",
    hjust = 0, family = base_font, fontface = "bold", size = 3.25
  ) +
  annotate(
    "text", x = 2.72, y = 15.12, label = "p",
    hjust = 0, family = base_font, fontface = "bold", size = 3.25
  ) +
  annotate(
    "text", x = 3.60, y = 15.12, label = "Overall p",
    hjust = 0, family = base_font, fontface = "bold", size = 3.25
  ) +
  scale_x_continuous(limits = c(0, 4.35), expand = c(0, 0)) +
  scale_y_continuous(limits = y_limits, expand = c(0, 0)) +
  theme_void(base_family = base_font) +
  theme(plot.margin = margin(5, 5, 5, 2))

figure3 <- (
  label_panel | ci_panel | numeric_panel
) +
  plot_layout(widths = c(1.25, 1.35, 1.70)) +
  plot_annotation(
    title = "Alcohol-use pattern and reported trouble sleeping",
    subtitle = paste0(
      "Survey-weighted logistic regression; ",
      "lifetime abstainers are the reference group"
    ),
    theme = theme(
      text = element_text(family = base_font, colour = colour_text),
      plot.title = element_text(
        face = "bold", size = 14.0,
        margin = margin(b = 3)
      ),
      plot.subtitle = element_text(
        face = "bold", size = 10.5,
        margin = margin(b = 8)
      ),
      plot.background = element_rect(
        colour = "black", fill = "white", linewidth = 0.55
      ),
      plot.margin = margin(10, 10, 8, 10)
    )
  )

save_figure(
  figure3,
  "Figure_3_primary_forest",
  width = 7.5,
  height = 6.9
)

############################################################
# 6. FIGURE 4 — SPEARMAN CORRELATIONS
#
# Only the lower triangle is shown. Short axis labels and a
# manuscript-width export keep every coefficient readable.
############################################################

corr_long <- read_csv(
  file.path(tables_dir, "spearman_correlation_long.csv"),
  show_col_types = FALSE
)

variable_order <- c(
  "drinks_per_day",
  "drinking_frequency_rank",
  "sleep_hours",
  "age_years",
  "bmi",
  "waist_circumference",
  "poverty_ratio",
  "depression_score_no_sleep"
)

variable_labels <- c(
  drinks_per_day = "Drinks/day",
  drinking_frequency_rank = "Drinking frequency",
  sleep_hours = "Sleep duration",
  age_years = "Age",
  bmi = "BMI",
  waist_circumference = "Waist circumference",
  poverty_ratio = "Poverty ratio",
  depression_score_no_sleep = "PHQ score\n(sleep item excluded)"
)

corr_data <- corr_long %>%
  mutate(
    row_index = match(variable_1, variable_order),
    col_index = match(variable_2, variable_order)
  ) %>%
  filter(row_index > col_index) %>%
  mutate(
    stars = case_when(
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ ""
    ),
    label = paste0(sprintf("%.2f", rho), stars),
    label_colour = if_else(abs(rho) >= 0.60, "white", "black"),
    row_label = factor(
      variable_labels[variable_1],
      levels = variable_labels[variable_order]
    ),
    col_label = factor(
      variable_labels[variable_2],
      levels = variable_labels[variable_order]
    )
  )

figure4 <- ggplot(
  corr_data,
  aes(x = col_label, y = row_label, fill = rho)
) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(
    aes(label = label, colour = label_colour),
    family = base_font,
    fontface = "bold",
    size = 3.55,
    show.legend = FALSE
  ) +
  scale_colour_identity() +
  scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0,
    limits = c(-1, 1),
    breaks = c(-1, -0.5, 0, 0.5, 1),
    name = "Spearman ρ"
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  coord_fixed() +
  labs(
    title = "Spearman correlations among alcohol, sleep, and covariates",
    subtitle = "Pairwise-complete, unweighted exploratory analysis",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 10.5, base_family = base_font) +
  theme(
    plot.title = element_text(
      face = "bold", size = 13.5,
      margin = margin(b = 3)
    ),
    plot.subtitle = element_text(
      face = "bold", size = 10.5,
      margin = margin(b = 8)
    ),
    axis.text.x = element_text(
      angle = 40, hjust = 1, vjust = 1,
      face = "bold", colour = colour_text,
      size = 9.0, lineheight = 0.95
    ),
    axis.text.y = element_text(
      face = "bold", colour = colour_text,
      size = 9.0, lineheight = 0.95
    ),
    panel.grid = element_blank(),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9.0),
    plot.background = element_rect(
      colour = "black", fill = "white", linewidth = 0.55
    ),
    plot.margin = margin(10, 12, 12, 12)
  )

save_figure(
  figure4,
  "Figure_4_spearman_correlation",
  width = 7.4,
  height = 6.5
)

############################################################
# 7. FIGURE MANIFEST
############################################################

figure_manifest <- tribble(
  ~figure, ~file_stem, ~place_after,
  1, "Figure_1_participant_flow",
  "Study population paragraph reporting model-specific sample sizes",
  2, "Figure_2_weighted_prevalence",
  "Weighted prevalence paragraph",
  3, "Figure_3_primary_forest",
  "Primary analysis paragraph",
  4, "Figure_4_spearman_correlation",
  "Exploratory Spearman-correlation paragraph"
)

write_csv(
  figure_manifest,
  file.path(figures_dir, "figure_manifest.csv")
)

cat("\nFinal publication figures created in:\n", figures_dir, "\n\n")
print(figure_manifest)

cat(
  "\nImportant: insert each image at 6.5–7.0 inches wide in Word.\n",
  "Do not add a second caption inside the image. Use the manuscript caption below it.\n",
  sep = ""
)
###done###