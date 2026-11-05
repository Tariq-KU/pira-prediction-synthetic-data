# This part will generate predictions for the reserved validation participants and compare their censoring-adjusted Brier scores.
# The following parts will address uncertainty in that comparison, discrimination and calibration.


# We now apply the original models saved in Lesson 5 to participants who were excluded from their development.

# Each model brings its existing:
#   - Predictor definitions and preprocessing rules.
#   - Selected ridge penalty.
#   - Fitted coefficients.
#   - Baseline survival estimated from the development data.

# This package estimates the censoring distribution.
install.packages("prodlim")

library(survival)
library(recipes)
library(glmnet)
library(prodlim)

lesson5 <- readRDS(
  "results/lesson05_ridge_models.rds"
)

models <- lesson5$models

dev <- read.csv("data/development.csv")
val <- read.csv("data/validation.csv")

# the following function (from lesson 4) applies our prespecified transformations and categorical definitions.
# It does not estimate new imputation or scaling values.
# Now, we have the log2 for brain volume and blood markers
val_features <- make_features(val)

horizon <- 24
# This means 24 months after the month-12 landmark, corresponding to study month 36.


# Check that the cohorts contain different participants:
stopifnot(
  anyDuplicated(val$participant_id) == 0L,
  
  length(
    intersect(dev$participant_id, val$participant_id)
  ) == 0L,
  
  all(val$pira_event %in% c(0, 1)),
  
  all(val$followup_months > 0 &
        val$followup_months <= horizon)
)

# Now inspect the horizon outcomes:
table(
  val$pira_by_24m,
  useNA = "ifany"
)


# 2. Define a function that produces absolute risk predictions

# In Lesson 6, we used the linear predictor, eta, to calculate the C-index.
# For the Brier score, we need probabilities between zero and one.
# We therefore obtain predicted survival at month 24 and subtract it from one.

predict_risk_at <- function(model_bundle,
                            new_data,
                            horizon = 24) {
  
  # Apply preprocessing learned during development
  new_x <- bake(
    model_bundle$recipe,
    new_data = new_data[model_bundle$variables],
    composition = "matrix"
  )
  
  # Calculate survival using the development model
  survival_curves <- survival::survfit(
    model_bundle$fit,
    s = model_bundle$selected_lambda,
    x = model_bundle$x_train,
    y = model_bundle$y_train,
    newx = new_x,
    stype = 2,
    ctype = 1
  )
  
  # The training x and y allow the survival prediction method to estimate the baseline survival from the development data
  # newx supplies the validation participants’ predictors
  
  # Convert survival probability into event probability
  risk <- 1 - as.numeric(
    summary(
      survival_curves,
      times = horizon
    )$surv
  )
  
  stopifnot(
    length(risk) == nrow(new_data),
    all(is.finite(risk)),
    all(risk >= 0 & risk <= 1)
  )
  
  risk
}

#------------------------------------------------------------------------------
# 3. Generate predictions from all four models

model_names <- names(models)

validation_predictions <- data.frame(
  participant_id = val_features$participant_id
)

for (model_name in model_names) {
  
  validation_predictions[[model_name]] <-
    predict_risk_at(
      model_bundle = models[[model_name]],
      new_data = val_features,
      horizon = horizon
    )
}

head(validation_predictions)

# a value of 0.18 means that the model assigns an estimated 18% probability of PIRA onset during the 24-month prediction window.
# Each model is evaluated on the same participants.

#------------------------------------------------------------------------------

# 4. Add the simple benchmark
# Our benchmark assigns every validation participant the same risk, estimated from the development cohort’s Kaplan–Meier curve.

# Why? This asks whether individual predictor information improves on a simple constant prediction.

km_development <- survfit(
  Surv(followup_months, pira_event) ~ 1,
  data = dev
)

benchmark_risk_24 <- as.numeric(
  1 - summary(
    km_development,
    times = horizon
  )$surv
)

validation_predictions$benchmark <-
  benchmark_risk_24

benchmark_risk_24

# For the original development dataset, this is approximately 0.1528, or 15.28%. This repeats the benchmark calculation from Lesson 3


#------------------------------------------------------------------------------
# 5. Identify known outcomes and estimate censoring weights

# Now, we will create a function that does 3 things:
    # 1- Identify observed events, known non-events and unknown outcomes.
    # 2- Estimate the probability of remaining uncensored long enough.
    # 3- Assign a weight to each usable outcome.

# Note: a person who is event-free through month 24 has a known outcome. We therefore include non-events censored exactly at month 24.
make_ipcw_outcome <- function(data, horizon = 24) {
  
  time <- data$followup_months
  event <- data$pira_event
  
  # Identify outcomes known at the prediction horizon
  event_by_horizon <- event == 1 & time <= horizon
  
  known_non_event <- !event_by_horizon &
    time >= horizon
  
  known <- event_by_horizon | known_non_event
  
  outcome <- rep(NA_real_, nrow(data))
  outcome[event_by_horizon] <- 1
  outcome[known_non_event] <- 0
  
  # Estimate the censoring distribution
  censor_fit <- prodlim::prodlim(
    Hist(followup_months, pira_event) ~ 1,
    data = data,
    reverse = TRUE,
    type = "surv"
  )
  
  # Look up probability of remaining observable
  # through each relevant assessment time
  evaluation_times <- pmin(time, horizon)
  
  positions <- prodlim::sindex(
    jump.times = censor_fit$time,
    eval.times = evaluation_times,
    strict = TRUE
  )
  
  G_before <- c(
    1,
    as.numeric(censor_fit$surv)
  )[positions + 1L]
  
  stopifnot(
    all(is.finite(G_before[known])),
    all(G_before[known] > 0)
  )
  
  # Unknown horizon outcomes have no direct error contribution
  weight <- numeric(nrow(data))
  weight[known] <- 1 / G_before[known]
  
  data.frame(
    outcome = outcome,
    known = known,
    observation_probability = G_before,
    weight = weight
  )
}

#------------------------------------------------------------------------------
# Now, we apply the function:
ipcw_outcome <- make_ipcw_outcome(
  data = val_features,
  horizon = horizon
)

table(
  ipcw_outcome$outcome,
  useNA = "ifany"
)

summary(
  ipcw_outcome$weight[ipcw_outcome$known]
)

#------------------------------------------------------------------------------

# 6. Calculate the censoring-adjusted Brier scores

# we create the adjusted brier score calculation function for application later:
brier_ipcw <- function(predicted_risk, outcome_data) {
  
  stopifnot(
    length(predicted_risk) == nrow(outcome_data),
    all(is.finite(predicted_risk)),
    all(predicted_risk >= 0 & predicted_risk <= 1)
  )
  
  known <- outcome_data$known
  
  squared_error <- (
    predicted_risk[known] -
      outcome_data$outcome[known]
  )^2
  
  weighted_error <- outcome_data$weight[known] *
    squared_error
  
  sum(weighted_error) / nrow(outcome_data)
}

#------------------------------------------------------------------------------

# Apply the calculation to all four models and the benchmark:

prediction_names <- c(
  model_names,
  "benchmark"
)

brier_values <- setNames(
  numeric(length(prediction_names)),
  prediction_names
)

for (prediction_name in prediction_names) {
  
  brier_values[prediction_name] <- brier_ipcw(
    predicted_risk =
      validation_predictions[[prediction_name]],
    
    outcome_data = ipcw_outcome
  )
}

brier_results <- data.frame(
  model = names(brier_values),
  Brier_24m = as.numeric(brier_values),
  row.names = NULL
)

brier_results


#------------------------------------------------------------------------------
# 7. Calculate the primary comparison
# We will define improvement as: Brier improvement=Brier(clinical+MRI)-Brier(combined)
# Because lower Brier scores are better, if the result is + then blood addition is good

brier_improvement <- unname(
  brier_values["clinical_mri"] -
    brier_values["combined"]
)

primary_comparison <- data.frame(
  comparison = "Clinical-MRI minus combined",
  Brier_improvement = brier_improvement
)

primary_comparison
# This gives improvement of 0.002042948. How do we evaluate how helpful is this really????

# That is a reduction in mean squared probability error
# A positive point estimate alone is not enough to establish a reliable improvement
# We still need a confidence interval for this paired comparison.

#------------------------------------------------------------------------------

# 8. Save this stage of the evaluation

dir.create(
"results",
showWarnings = FALSE
)

saveRDS(
  list(
    horizon = horizon,
    validation_features = val_features,
    predictions = validation_predictions,
    ipcw_outcome = ipcw_outcome,
    benchmark_risk = benchmark_risk_24,
    brier_results = brier_results,
    primary_comparison = primary_comparison,
    session_info = sessionInfo()
  ),
  "results/lesson07_part1_validation.rds"
)

write.csv(
  validation_predictions,
  "results/lesson07_validation_predictions.csv",
  row.names = FALSE
)

write.csv(
  brier_results,
  "results/lesson07_brier_scores.csv",
  row.names = FALSE
)


# #------------------------------------------------------------------------------

####

# Plotting the ROC curve:
lesson7a <- readRDS("results/lesson07_part1_validation.rds")

roc_predictions <- lesson7a$predictions
roc_outcomes <- lesson7a$ipcw_outcome

# Check that the rows describe the same participants in the same order.
stopifnot(
  identical(
    as.character(roc_predictions$participant_id),
    as.character(lesson7a$validation_features$participant_id)
  ),
  nrow(roc_predictions) == nrow(roc_outcomes),
  lesson7a$horizon == 24
)


#------------------------------------------------------------------------------
# 2. Calculate the points that form each ROC curve.

# At each prediction threshold, we calculate:
  # - True-positive rate: the weighted proportion of PIRA cases classified as positive.
  # - False-positive rate: the weighted proportion of non-cases classified as positive.


make_weighted_roc <- function(predicted_risk, outcome_data) {
  
  # Keep participants whose outcome is known at the prediction horizon.
  known <- outcome_data$known
  
  p <- predicted_risk[known]
  y <- outcome_data$outcome[known]
  w <- outcome_data$weight[known]
  
  # Check the inputs before calculating the curve.
  stopifnot(
    all(is.finite(p)),
    all(p >= 0 & p <= 1),
    all(is.finite(w)),
    all(w > 0),
    all(y %in% c(0, 1)),
    any(y == 1),
    any(y == 0)
  )
  
  # Start above every prediction: initially nobody is classified positive.
  # Then progressively lower the threshold.
  thresholds <- c(
    Inf,
    sort(unique(p), decreasing = TRUE)
  )
  
  # Denominators: total case weight and total non-case weight.
  case_weight <- sum(w[y == 1])
  control_weight <- sum(w[y == 0])
  
  tpr <- vapply(thresholds, function(cutoff) {
    sum(w[y == 1 & p >= cutoff]) / case_weight
  }, numeric(1))
  
  fpr <- vapply(thresholds, function(cutoff) {
    sum(w[y == 0 & p >= cutoff]) / control_weight
  }, numeric(1))
  
  # Area under the curve using the trapezoidal rule.
  auc <- sum(
    diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2
  )
  
  list(
    fpr = fpr,
    tpr = tpr,
    auc = auc
  )
}

#------------------------------------------------------------------------------
# Now apply that function to our two models:

roc_clinical_mri <- make_weighted_roc(
  predicted_risk = roc_predictions$clinical_mri,
  outcome_data = roc_outcomes
)

roc_combined <- make_weighted_roc(
  predicted_risk = roc_predictions$combined,
  outcome_data = roc_outcomes
)

data.frame(
  model = c("Clinical-MRI", "Combined"),
  AUC_24m = c(
    roc_clinical_mri$auc,
    roc_combined$auc
  )
)

#------------------------------------------------------------------------------
# 3. Draw the curves together.


draw_roc <- function() {
  
  plot(
    roc_clinical_mri$fpr,
    roc_clinical_mri$tpr,
    type = "l",
    col = "#0072B2",
    lwd = 2,
    xlim = c(0, 1),
    ylim = c(0, 1),
    xaxs = "i",
    yaxs = "i",
    asp = 1,
    xlab = "False-positive rate (1 - specificity)",
    ylab = "True-positive rate (sensitivity)",
    main = "24-month PIRA ROC curves"
  )
  
  # Reference diagonal: no discrimination.
  abline(
    a = 0,
    b = 1,
    col = "grey65",
    lty = 3
  )
  
  # Add the combined model to the same plot.
  lines(
    roc_combined$fpr,
    roc_combined$tpr,
    col = "#D55E00",
    lwd = 2,
    lty = 2
  )
  
  legend(
    "bottomright",
    legend = c(
      sprintf("Clinical-MRI: AUC = %.3f", roc_clinical_mri$auc),
      sprintf("Combined: AUC = %.3f", roc_combined$auc)
    ),
    col = c("#0072B2", "#D55E00"),
    lwd = 2,
    lty = c(1, 2),
    bg = "white",
    cex = 0.85
  )
}

# Display in the RStudio Plots pane or your R Markdown output.
draw_roc()

#------------------------------------------------------------------------------
# TO save a high resulotion for the report:
dir.create("figures", showWarnings = FALSE)

png(
  filename = "figures/PIRA_ROC_24m.png",
  width = 2000,
  height = 2000,
  res = 300
)

draw_roc()

dev.off()