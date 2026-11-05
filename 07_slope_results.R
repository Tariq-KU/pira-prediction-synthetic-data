
# finish calibration using only the slope and its 95% confidence interval

# The purpose is to assess whether the predictions spread risks too widely or too narrowly. We will use IPCW-weighted logistic regression as a diagnostic; the original Cox models remain unchanged.

#------------------------------------------------------------------------------
# 1. Load the saved validation predictions

library(prodlim)

lesson7a <- readRDS(
  "results/lesson07_part1_validation.rds"
)

val_features <- lesson7a$validation_features
validation_predictions <- lesson7a$predictions
horizon <- lesson7a$horizon

model_names <- c(
  "clinical",
  "clinical_mri",
  "clinical_blood",
  "combined"
)

stopifnot(
  nrow(val_features) ==
    nrow(validation_predictions),
  
  identical(
    as.character(val_features$participant_id),
    as.character(validation_predictions$participant_id)
  ),
  
  all(model_names %in%
        names(validation_predictions))
)

# We reuse two functions from earlier parts:
#    - make_ipcw_outcome() — creates the outcome information and censoring weights.
#    - draw_rows_within_centres() — samples validation participants within centres.


# Run if restarted R
for (function_name in c(
  "make_ipcw_outcome",
  "draw_rows_within_centres"
)) {
  
  if (!exists(function_name, mode = "function")) {
    stop(
      "Rerun the function-definition chunk for: ",
      function_name
    )
  }
}

#------------------------------------------------------------------------------

# 2. Define the calibration-slope calculation

fit_calibration_slope <- function(
    predicted_risk,
    outcome_data
) {
  
  stopifnot(
    length(predicted_risk) ==
      nrow(outcome_data),
    
    all(is.finite(predicted_risk)),
    
    all(predicted_risk > 0 &
          predicted_risk < 1)
  )
  
  known <- outcome_data$known
  
  calibration_data <- data.frame(
    outcome =
      outcome_data$outcome[known],
    
    logit_prediction =
      qlogis(predicted_risk[known]),
    
    ipcw_weight =
      outcome_data$weight[known]
  )
  
  fit <- glm(
    outcome ~ logit_prediction,
    
    data = calibration_data,
    
    weights = ipcw_weight,
    
    family = quasibinomial(link = "logit")
  )
  
  if (
    !fit$converged ||
    any(!is.finite(coef(fit)))
  ) {
    stop(
      "The calibration regression did not fit reliably."
    )
  }
  
  unname(
    coef(fit)["logit_prediction"]
  )
}

#  The key operations are:
#   - qlogis() converts the predicted probabilities to log odds.
#   - glm() fits the diagnostic relationship, including an intercept.
#   - weights supplies the censoring adjustment.
#   - coef(fit)["logit_prediction"] extracts the slope.

# quasibinomial(link = "logit") provides the weighted logistic fit. We will obtain uncertainty through bootstrapping rather than relying on its default standard errors.


#------------------------------------------------------------------------------

# 3. Calculate the original slope estimates
ipcw_outcome <- make_ipcw_outcome(
  data = val_features,
  horizon = horizon
)

slope_estimates <- setNames(
  numeric(length(model_names)),
  model_names
)

for (model_name in model_names) {
  
  slope_estimates[model_name] <-
    fit_calibration_slope(
      predicted_risk =
        validation_predictions[[model_name]],
      
      outcome_data =
        ipcw_outcome
    )
}

slope_estimates

# These are the estimates from your original validation sample.
# The benchmark is excluded because everyone receives the same prediction, so its calibration slope cannot be estimated.


# 4. Estimate uncertainty through bootstrapping

# We resample validation participants, recalculate their censoring weights, and estimate the slope again. We keep their original predictions fixed.

B_calibration <- 2000L

slope_bootstrap <- matrix(
  NA_real_,
  nrow = B_calibration,
  ncol = length(model_names),
  
  dimnames = list(
    NULL,
    model_names
  )
)

set.seed(7301)

for (b in seq_len(B_calibration)) {
  
  rows <- draw_rows_within_centres(
    val_features
  )
  
  bootstrap_outcomes <- make_ipcw_outcome(
    data = val_features[
      rows, , drop = FALSE
    ],
    
    horizon = horizon
  )
  
  for (model_name in model_names) {
    
    slope_bootstrap[b, model_name] <-
      fit_calibration_slope(
        predicted_risk =
          validation_predictions[[model_name]][rows],
        
        outcome_data =
          bootstrap_outcomes
      )
  }
  
  if (b %% 200 == 0) {
    message(
      "Completed ",
      b,
      " of ",
      B_calibration,
      " repetitions"
    )
  }
}

# Every model uses the same sampled participants within a repetition.
# Only the diagnostic regressions are fitted again—the Cox models and ridge penalties stay fixed.


#------------------------------------------------------------------------------
# 5. Produce the final table

slope_limits <- apply(
  slope_bootstrap,
  MARGIN = 2,
  FUN = quantile,
  
  probs = c(0.025, 0.975),
  names = FALSE
)

slope_results <- data.frame(
  model = model_names,
  
  estimate =
    as.numeric(slope_estimates[model_names]),
  
  lower_95 =
    as.numeric(slope_limits[1, model_names]),
  
  upper_95 =
    as.numeric(slope_limits[2, model_names]),
  
  row.names = NULL
)

slope_results


# Slope near 1: Approximately appropriate spread of predictions
# Positive slope below 1: Predictions spread risks too widely
# Slope above 1: Predictions are too compressed

#------------------------------------------------------------------------------
# 6. Save the results

dir.create(
  "results",
  showWarnings = FALSE
)

saveRDS(
  list(
    horizon = horizon,
    
    method = "IPCW-weighted logistic calibration slope",
    
    slope_results = slope_results,
    
    slope_bootstrap = slope_bootstrap,
    
    repetitions = B_calibration,
    
    seed = 7301,
    
    session_info = sessionInfo()
  ),
  
  file = "results/lesson07_part3b_calibration_slopes.rds"
)

write.csv(
  slope_results,
  
  file = "results/lesson07_calibration_slopes.csv",
  
  row.names = FALSE
)