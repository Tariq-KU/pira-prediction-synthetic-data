
# We have calculated the models’ Brier scores.
# Now we want to understand how much those scores—and the difference between the models might vary if we evaluated them on another comparable sample of participants.
# Where do we get the comparable sample of participants???

# The combined model’s Brier score is lower by 0.004 in that validation sample. But a different sample might produce a smaller improvement, a larger improvement, or even a reversal.
# That is why we need a confidence interval around the difference, alongside the difference itself. Reporting uncertainty is an essential part of evaluating prediction models in an independent validation dataset.


# We will use bootstrap
# We will repeatedly draw samples with replacement from the validation participants.
# Each sample will contain the same number of participants as the original validation dataset, but some participants will appear more than once and others will be absent.
# For each, we calculate the diff in the Brier score between clinical+MRI and combined


# “Paired” means both models are evaluated on exactly the same sampled participants in each repetition

#------------------------------------------------------------------------------


# 1. Load the saved validation predictions

library(prodlim)

lesson7a <- readRDS(
  "results/lesson07_part1_validation.rds"
)

val_features <- lesson7a$validation_features
validation_predictions <- lesson7a$predictions
horizon <- lesson7a$horizon

prediction_names <- c(
  "clinical",
  "clinical_mri",
  "clinical_blood",
  "combined",
  "benchmark"
)


# We will rerun those functions:
    # make_ipcw_outcome() calculates the observed outcomes and censoring weights.
    # brier_ipcw() calculates a Brier score using those weights.

#------------------------------------------------------------------------------
# Next, check that the predictions and outcomes correspond to the same participants in the same order:
stopifnot(
  nrow(val_features) == nrow(validation_predictions),
  
  identical(
    as.character(val_features$participant_id),
    as.character(validation_predictions$participant_id)
  ),
  
  anyDuplicated(val_features$participant_id) == 0L,
  
  !anyNA(val_features$centre_id),
  
  all(
    prediction_names %in%
      names(validation_predictions)
  )
)


#------------------------------------------------------------------------------

# 2. Create one function that calculates all the scores

evaluate_brier <- function(data, predictions, horizon = 24) {
  
  stopifnot(
    nrow(data) == nrow(predictions),
    
    identical(
      as.character(data$participant_id),
      as.character(predictions$participant_id)
    )
  )
  
  # Estimate censoring weights for the supplied sample.
  outcome_data <- make_ipcw_outcome(
    data,
    horizon = horizon
  )
  
  # Create a named container for the scores.
  scores <- setNames(
    numeric(length(prediction_names)),
    prediction_names
  )
  
  # Calculate a Brier score for each set of predictions.
  for (model_name in prediction_names) {
    
    scores[model_name] <- brier_ipcw(
      predicted_risk = predictions[[model_name]],
      outcome_data = outcome_data
    )
  }
  
  # Return the scores and our primary comparison.
  c(
    scores,
    
    improvement = unname(
      scores["clinical_mri"] - scores["combined"]
    )
  )
}

# The function receives a dataset and its corresponding predictions.
# It returns six numbers:
    #five Brier scores
    # the primary between-model difference.

# Notice that make_ipcw_outcome() is inside this function. Therefore, every time we evaluate a new bootstrap sample, its censoring weights are recalculated.

#------------------------------------------------------------------------------
# Now apply the function to the original validation dataset:
original_scores <- evaluate_brier(
  data = val_features,
  predictions = validation_predictions,
  horizon = horizon
)

original_scores

#------------------------------------------------------------------------------

# 3. Define how we sample participants within centres

# Suppose one validation centre contributed 100 participants. In each repetition, we draw 100 participants with replacement from that centre. We do the same for every other centre.
# keeps each centre’s contribution to the validation sample unchanged while allowing participants to vary.

draw_rows_within_centres <- function(data) {
  
  sampled_rows <- integer(0)
  
  for (centre in unique(data$centre_id)) {
    
    # Find the rows belonging to this centre.
    centre_rows <- which(
      data$centre_id == centre
    )
    
    # Sample positions within those rows, with replacement.
    sampled_positions <- sample.int(
      n = length(centre_rows),
      size = length(centre_rows),
      replace = TRUE
    )
    
    # Add the selected participant rows to our sample.
    sampled_rows <- c(
      sampled_rows,
      centre_rows[sampled_positions]
    )
  }
  
  sampled_rows
}

#------------------------------------------------------------------------------

# 4. Try one repetition before running many

set.seed(7001)

rows <- draw_rows_within_centres(
  val_features
)

one_bootstrap_result <- evaluate_brier(
  data = val_features[rows, , drop = FALSE],
  
  predictions = validation_predictions[
    rows, , drop = FALSE
  ],
  
  horizon = horizon
)

one_bootstrap_result


# You can compare this repetition with the original results:
rbind(
  original = original_scores,
  one_bootstrap_sample = one_bootstrap_result
)


#------------------------------------------------------------------------------

# 5. Repeat the evaluation 2,000 times

# We will use 2,000 repetitions as a practical starting point for estimating the confidence interval endpoints.
# This is less computationally demanding than Lesson 6 because we are recalculating scores, without retraining the models.

B_validation <- 2000L

bootstrap_scores <- matrix(
  NA_real_,
  nrow = B_validation,
  ncol = length(original_scores),
  
  dimnames = list(
    NULL,
    names(original_scores)
  )
)

set.seed(7002)

for (b in seq_len(B_validation)) {
  
  rows <- draw_rows_within_centres(
    val_features
  )
  
  scores_b <- evaluate_brier(
    data = val_features[rows, , drop = FALSE],
    
    predictions = validation_predictions[
      rows, , drop = FALSE
    ],
    
    horizon = horizon
  )
  
  if (!all(is.finite(scores_b))) {
    stop(
      "A non-finite score occurred in repetition ",
      b
    )
  }
  
  bootstrap_scores[b, ] <- scores_b[
    colnames(bootstrap_scores)
  ]
  
  if (b %% 200 == 0) {
    message(
      "Completed ",
      b,
      " of ",
      B_validation,
      " repetitions"
    )
  }
}


head(bootstrap_scores)

dim(bootstrap_scores)


#------------------------------------------------------------------------------

# 6. Calculate the confidence intervals

# use an approximate 95% percentile bootstrap confidence interval.
# Imagine sorting the 2,000 improvements from smallest to largest. We take:
    # The 2.5th percentile as the lower endpoint.
    # The 97.5th percentile as the upper endpoint.

ci_limits <- apply(
  bootstrap_scores,
  MARGIN = 2,
  FUN = quantile,
  probs = c(0.025, 0.975),
  names = FALSE
)

bootstrap_se <- apply(
  bootstrap_scores,
  MARGIN = 2,
  FUN = sd
)

validation_ci <- data.frame(
  measure = names(original_scores),
  
  estimate = as.numeric(original_scores),
  
  lower_95 = as.numeric(
    ci_limits[1, ]
  ),
  
  upper_95 = as.numeric(
    ci_limits[2, ]
  ),
  
  bootstrap_SE = as.numeric(
    bootstrap_se
  ),
  
  row.names = NULL
)


# Now separate the individual scores from the primary comparison:
brier_ci <- validation_ci[
  validation_ci$measure != "improvement",
]

primary_comparison_ci <- validation_ci[
  validation_ci$measure == "improvement",
]

brier_ci

primary_comparison_ci


#------------------------------------------------------------------------------
# 7. Visualise the uncertainty

bootstrap_delta <- bootstrap_scores[
  , "improvement"
]

hist(
  bootstrap_delta,
  breaks = 30,
  
  xlim = range(
    c(
      bootstrap_delta,
      original_scores["improvement"],
      0
    )
  ),
  
  col = "grey85",
  border = "white",
  
  main = "Uncertainty in the Brier-score improvement",
  
  xlab = "Clinical-MRI Brier score minus combined Brier score"
)

abline(
  v = 0,
  col = "red",
  lty = 2,
  lwd = 2
)

abline(
  v = original_scores["improvement"],
  col = "blue",
  lwd = 2
)

abline(
  v = ci_limits[, "improvement"],
  col = "grey30",
  lty = 3,
  lwd = 2
)

legend(
  "topright",
  
  legend = c(
    "No difference",
    "Original estimate",
    "95% interval endpoints"
  ),
  
  col = c(
    "red",
    "blue",
    "grey30"
  ),
  
  lty = c(2, 1, 3),
  lwd = 2,
  bty = "n"
)


# Read the plot as follows:
#   - Values to the right of zero favour the combined model.
#   - Values to the left of zero favour the clinical–MRI model.
#   - A wider distribution indicates greater uncertainty about the improvement.
#   - The blue line marks the estimate from the original validation sample.

#------------------------------------------------------------------------------

# 8. Save the results:
dir.create(
  "results",
  showWarnings = FALSE
)

saveRDS(
  list(
    horizon = horizon,
    original_scores = original_scores,
    bootstrap_scores = bootstrap_scores,
    summary = validation_ci,
    primary_comparison = primary_comparison_ci,
    
    repetitions = B_validation,
    seed = 7002,
    
    resampling = paste(
      "Participants sampled within validation centres;",
      "predictions fixed;",
      "censoring weights re-estimated in each sample"
    ),
    
    session_info = sessionInfo()
  ),
  
  file = "results/lesson07_part2_paired_bootstrap.rds"
)

write.csv(
  validation_ci,
  
  file = "results/lesson07_validation_confidence_intervals.csv",
  
  row.names = FALSE
)