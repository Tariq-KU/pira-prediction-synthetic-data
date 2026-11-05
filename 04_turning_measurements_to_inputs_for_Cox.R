# add package "reciept" which keeps preprocessing instructions together. It prepares predictor measurements for modelling.

install.packages("recipes")
library(recipes)
library(survival)

dev <- read.csv("data/development.csv")

# recipe(): Write down the preprocessing instructions.
# prep(): Learn the required quantities from training data.
# bake(): Apply the learned instructions to a dataset.

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
  "log2_t2_lesion_m12",
  "cord_area_m12_mm2",
  "cord_area_change_pct_m0_m12"
)

blood_vars <- c(
  "log2_snfl_m12",
  "log2_sgfap_m12",
  "snfl_log2_change_m0_m12",
  "sgfap_log2_change_m0_m12"
)

predictor_groups <- list(    #stores the four collections of predictor names
  clinical = clinical_vars,
  clinical_mri = c(clinical_vars, mri_vars),
  clinical_blood = c(clinical_vars, blood_vars),
  combined = c(clinical_vars, mri_vars, blood_vars)
)

sapply(predictor_groups, length)  #counts the variables in each collection.



# creates a function that converts the sNfL, sGFAP and brain volume into the log2 scale and assign sex and DMT as categories
make_features <- function(d) {
  
  stopifnot(
    all(d$snfl_m12_pg_ml > 0, na.rm = TRUE),
    all(d$sgfap_m12_pg_ml > 0, na.rm = TRUE),
    all(d$t2_lesion_volume_m12_ml >= 0, na.rm = TRUE)
  )
  
  d$log2_snfl_m12 <- log2(d$snfl_m12_pg_ml)
  
  d$log2_sgfap_m12 <- log2(d$sgfap_m12_pg_ml)
  
  d$log2_t2_lesion_m12 <- log2(
    1 + d$t2_lesion_volume_m12_ml
  )
  
  d$sex <- factor(
    d$sex,
    levels = c("female", "male")
  )
  
  d$dmt_class_m12 <- factor(
    d$dmt_class_m12,
    levels = c(
      "anti_cd20",
      "other_high_efficacy",
      "moderate_oral",
      "injectable"
    )
  )
  
  stopifnot(
    !anyNA(d$sex),
    !anyNA(d$dmt_class_m12)
  )
  
  return(d)
}

dev_features <- make_features(dev)


# Create a temporary training and held-out split

set.seed(2026)  # makes the selection reproducible.

# seq_len(nrow(dev_features)): creates the row numbers.
#sample(): randomly selects 80% of them.

training_rows <- sample(
  seq_len(nrow(dev_features)),
  size = floor(0.80 * nrow(dev_features))
)

demo_train <- dev_features[training_rows, ]

demo_check <- dev_features[-training_rows, ]

c(
  training_participants = nrow(demo_train),
  heldout_participants = nrow(demo_check)
)


# Select the combined model’s predictors
# We will demonstrate preprocessing with the combined model because it contains all predictor types.

selected_vars <- predictor_groups$combined

train_x_raw <- demo_train[selected_vars] # here we select the training data that we just seperated and only apply the columns that are predictors for the combined model

check_x_raw <- demo_check[selected_vars]

dim(train_x_raw)
dim(check_x_raw)

# Here, “raw” means before imputation, scaling and numerical category coding. The fixed transformations from Chunk 3 have already happened.


# Write the preprocessing recipe (we are making a function here)

# The pipe symbol, |>, passes the result of one operation into the next.
# 1- recipe(~ ., ...) means “use all columns in this predictor-only dataset"
# 2- step_impute_median() will learn a median for each numerical predictor and use it for missing values.
# 3- step_normalize() uses: (x-training mean)/training standard deviation
# 4- step_dummy() converts factors into numerical columns. With one_hot = FALSE, a variable with k categories produces k-1 columns.

make_preprocessor <- function(training_predictors) {
  
  recipe(~ ., data = training_predictors) |>
    step_impute_median(all_numeric_predictors()) |>
    step_normalize(all_numeric_predictors()) |>
    step_dummy(
      all_nominal_predictors(),
      one_hot = FALSE
    )
}

recipe_plan <- make_preprocessor(train_x_raw)

recipe_plan


# Learn the preprocessing quantities from training participants
trained_recipe <- prep(
  recipe_plan,
  training = train_x_raw
)

# prep() calculates the quantities required by the recipe using the data supplied through training. It stores those quantities inside trained_recipe
# Compares between the training raw data and the processed data to learn by:
# 1. Calculates training medians.
# 2. Fills missing training values.
# 3. Calculates means and standard deviations from the resulting training values.
# 4. Records the category coding.


# You can inspect what it learned through the following where 1 refers to the first recipe step (median imputation), 2 refer to scaling and 3 for category coding.
tidy(trained_recipe, number = 1)
tidy(trained_recipe, number = 2)
tidy(trained_recipe, number = 3)



# Apply the same learned rules to both datasets
train_x <- bake(
trained_recipe,
new_data = train_x_raw
)

# here we apply the learned rules from the training data only to the validation data. That's why we didn't have prep for check
check_x <- bake(
  trained_recipe,
  new_data = check_x_raw
)


# just double checking on the data
names(train_x)
dim(train_x)
dim(check_x)


check_x  # this now has 3 columns for DMT because we had total 4. So 0,0,0 means the 4th DMT, 1,0,0 means he is on the 1st DMT and so on

# To see a table visualiztion of what I meant above
dummy_columns <- c(
  "sex_male",
  "dmt_class_m12_other_high_efficacy",
  "dmt_class_m12_moderate_oral",
  "dmt_class_m12_injectable"
)

head(train_x[dummy_columns])


# Check missingness, scaling and column consistency. We should have zero since we filled them
c(
  missing_training_values = sum(is.na(train_x)),
  missing_heldout_values = sum(is.na(check_x))
)

identical(names(train_x), names(check_x))

# Sanity check: inspect one scaled variable. We should find training mean approximately 0 and training standard deviation approximately 1.
# While the Held-out mean and standard deviation generally different from 0 and 1.
c(
  training_age_mean = mean(train_x$age_m12_years),
  training_age_sd = sd(train_x$age_m12_years),
  heldout_age_mean = mean(check_x$age_m12_years),
  heldout_age_sd = sd(check_x$age_m12_years)
)


# Create the numerical matrices and corresponding survival outcomes

# as.matrix() creates the numerical format that glmnet will use.
# x_train: Prepared predictor measurements for training participants.
# y_train: Their follow-up times and event indicators. (6+, 18, 24+ etc)
x_train <- as.matrix(train_x)

x_check <- as.matrix(check_x)

y_train <- with(
  demo_train,
  Surv(
    time = followup_months,
    event = pira_event,
    type = "right"
  )
)

y_check <- with(
  demo_check,
  Surv(
    time = followup_months,
    event = pira_event,
    type = "right"
  )
)


# Finally
dir.create("results", showWarnings = FALSE)

saveRDS(
  list(
    training_rows = training_rows,
    predictor_groups = predictor_groups,
    trained_recipe = trained_recipe,
    x_train = x_train,
    y_train = y_train,
    x_check = x_check,
    y_check = y_check
  ),
  "results/lesson04_preprocessing_demo.rds"
)
