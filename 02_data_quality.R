dev <- read.csv(
  "data/development.csv",
  na.strings = c("", "NA"),
  stringsAsFactors = FALSE
)

str(dev) #displays the structure of the dataset: its dimensions, column names, data types and some example values.

# c() combines items into a vector that contain column names, not the measurements themselves.

clinical_vars <- c(
  "age_m12_years",
  "sex",
  "disease_duration_m12_years",
  "edss_m12",
  "relapses_m0_m12",
  "dmt_class_m12",
  "dmt_duration_m12_months",
  "bmi_m12_kg_m2",
  "egfr_m12_ml_min_1_73m2"
)

mri_vars <- c(
  "brain_parenchymal_fraction_m12",
  "brain_volume_change_pct_m0_m12",
  "t2_lesion_volume_m12_ml",
  "cord_area_m12_mm2",
  "cord_area_change_pct_m0_m12"
)

blood_vars <- c(
  "snfl_m12_pg_ml",
  "sgfap_m12_pg_ml",
  "snfl_log2_change_m0_m12",
  "sgfap_log2_change_m0_m12"
)

predictor_vars <- c(clinical_vars, mri_vars, blood_vars)

length(predictor_vars)


# Now we can use these vectors to select the corresponding measurements:
predictor_data <- dev[predictor_vars]


# sum the number of columns that has NA in the predictor data, then create a summary
missing_counts <- colSums(is.na(predictor_data))

missing_table <- data.frame(
  variable = names(missing_counts),
  n_missing = as.integer(missing_counts),
  percent_missing = 100 * as.integer(missing_counts) / nrow(dev)
)

missing_table


# count participants who have at least one missing predictor
# complete.cases() identifies rows with no missing entries

sum(!complete.cases(predictor_data))


# summary() reports the minimum, quartiles, median, mean, maximum and missing count for numeric columns.
summary(dev[c(
  "age_m12_years",
  "edss_m12",
  "snfl_m12_pg_ml",
  "sgfap_m12_pg_ml",
  "t2_lesion_volume_m12_ml",
  "brain_volume_change_pct_m0_m12"
)])


# breaks = 25 requests approximately 25 intervals.
# These distributions should have a longer tail towards larger values, called right skew because mean>median.
hist(
  dev$snfl_m12_pg_ml,
  breaks = 25,
  main = "Distribution of month-12 sNfL",
  xlab = "sNfL (pg/mL)",
  col = "lightblue",
  border = "white"
)


hist(
  dev$t2_lesion_volume_m12_ml,
  breaks = 25,
  main = "Distribution of month-12 T2 lesion volume",
  xlab = "T2 lesion volume (mL)",
  col = "lightgreen",
  border = "white"
)



events <- subset(dev, pira_event == 1) #selects only the data that has pira=1

timing_checks <- c(
  onset_matches_analysis_time = all(
    events$pira_onset_study_month ==
      events$followup_months + 12
  ),
  
  confirmation_is_12_months_later = all(
    events$pira_confirmation_study_month ==
      events$pira_onset_study_month + 12
  ),
  
  confirmation_within_followup = all(
    events$pira_confirmation_study_month <=
      events$last_followup_study_month
  )
)

# The reasoning is:
# 1. Analytical time starts at study month 12. An event 6 months after that landmark belongs at study month 18.
# 2. The synthetic dataset assigns confirmation 12 months after onset.
# 3. Confirmation must occur while follow-up is available.

timing_checks



non_events <- subset(dev, pira_event == 0)

all(
  non_events$followup_months ==
    non_events$adjudicable_through_study_month - 12
)


# Creates a folder inside this path
dir.create("results", showWarnings = FALSE)

#creates a file, name it, insert path
write.csv(
  missing_table,
  "results/lesson02_missing_predictors.csv",
  row.names = FALSE
)

write.csv(
  data.frame(check = names(timing_checks), passed = unname(timing_checks)),
  "results/lesson02_timing_checks.csv",
  row.names = FALSE
)
