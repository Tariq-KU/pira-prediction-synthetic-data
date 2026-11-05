
# This is the final analysis step: independent-validation AUC at 24 months.
# We will produce:
  # Each model’s AUC with a 95% confidence interval.
  # The difference between the combined and clinical–MRI models, with its paired confidence interval.

#------------------------------------------------------------------------------
# 1. Load the validation predictions

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


# We need As in the preceding part, keep these earlier functions available:
    # make_ipcw_outcome()
    # draw_rows_within_centres()


#------------------------------------------------------------------------------

# 2. Define the AUC calculation

auc_ipcw <- function(
    predicted_risk,
    outcome_data
) {
  
  stopifnot(
    length(predicted_risk) ==
      nrow(outcome_data),
    
    all(is.finite(predicted_risk)),
    
    all(predicted_risk >= 0 &
          predicted_risk <= 1)
  )
  
  cases <- which(
    outcome_data$known &
      outcome_data$outcome == 1
  )
  
  noncases <- which(
    outcome_data$known &
      outcome_data$outcome == 0
  )
  
  if (
    length(cases) == 0L ||
    length(noncases) == 0L
  ) {
    stop(
      "AUC requires both cases and non-cases."
    )
  }
  
  # Compare every case with every non-case.
  pair_credit <- outer(
    predicted_risk[cases],
    predicted_risk[noncases],
    FUN = ">"
  ) + 0.5 * outer(
    predicted_risk[cases],
    predicted_risk[noncases],
    FUN = "=="
  )
  
  # Construct a censoring weight for each pair.
  pair_weights <- outer(
    outcome_data$weight[cases],
    outcome_data$weight[noncases],
    FUN = "*"
  )
  
  sum(pair_credit * pair_weights) /
    sum(pair_weights)
}



# 3. Calculate the original AUCs and their difference
# We will wrap the evaluation into a function so the same procedure can be repeated during bootstrapping.

evaluate_auc <- function(
    data,
    predictions,
    horizon = 24
) {
  
  stopifnot(
    nrow(data) == nrow(predictions),
    
    identical(
      as.character(data$participant_id),
      as.character(predictions$participant_id)
    )
  )
  
  outcome_data <- make_ipcw_outcome(
    data = data,
    horizon = horizon
  )
  
  scores <- setNames(
    numeric(length(model_names)),
    model_names
  )
  
  for (model_name in model_names) {
    
    scores[model_name] <- auc_ipcw(
      predicted_risk =
        predictions[[model_name]],
      
      outcome_data =
        outcome_data
    )
  }
  
  c(
    scores,
    
    gain_combined_vs_mri = unname(
      scores["combined"] -
        scores["clinical_mri"]
    )
  )
}

original_auc <- evaluate_auc(
  data = val_features,
  predictions = validation_predictions,
  horizon = horizon
)

original_auc


#------------------------------------------------------------------------------

# 4. Bootstrap the AUCs and the paired difference

# We resample participants within centres and evaluate all models on the same sampled participants.
# Censoring weights are recalculated each time.
B_auc <- 2000L

auc_bootstrap <- matrix(
  NA_real_,
  nrow = B_auc,
  ncol = length(original_auc),
  
  dimnames = list(
    NULL,
    names(original_auc)
  )
)

set.seed(7401)

for (b in seq_len(B_auc)) {
  
  rows <- draw_rows_within_centres(
    val_features
  )
  
  result_b <- evaluate_auc(
    data = val_features[
      rows, , drop = FALSE
    ],
    
    predictions = validation_predictions[
      rows, , drop = FALSE
    ],
    
    horizon = horizon
  )
  
  auc_bootstrap[b, ] <- result_b[
    colnames(auc_bootstrap)
  ]
  
  if (b %% 200 == 0) {
    message(
      "Completed ",
      b,
      " of ",
      B_auc,
      " AUC repetitions"
    )
  }
}

# This estimates uncertainty in validation performance. It does not retrain or tune the models.


#------------------------------------------------------------------------------

# 5. Create the results tables

auc_limits <- apply(
  auc_bootstrap,
  MARGIN = 2,
  FUN = quantile,
  
  probs = c(0.025, 0.975),
  names = FALSE
)

auc_summary <- data.frame(
  measure = names(original_auc),
  
  estimate =
    as.numeric(original_auc),
  
  lower_95 =
    as.numeric(auc_limits[1, ]),
  
  upper_95 =
    as.numeric(auc_limits[2, ]),
  
  row.names = NULL
)

auc_results <- auc_summary[
  auc_summary$measure %in% model_names,
]

auc_comparison <- auc_summary[
  auc_summary$measure == "gain_combined_vs_mri",
]

auc_results

auc_comparison


# auc_results contains each model’s AUC and confidence interval.
#auc_comparison contains the paired improvement from adding blood biomarkers

#------------------------------------------------------------------------------

# 6. Save the outputs

dir.create(
  "results",
  showWarnings = FALSE
)

saveRDS(
  list(
    horizon = horizon,
    
    method = "IPCW cumulative/dynamic AUC",
    
    auc_results = auc_results,
    
    primary_comparison = auc_comparison,
    
    original_auc = original_auc,
    
    auc_bootstrap = auc_bootstrap,
    
    repetitions = B_auc,
    
    seed = 7401,
    
    session_info = sessionInfo()
  ),
  
  file = "results/lesson07_part4_auc.rds"
)

write.csv(
  auc_summary,
  
  file = "results/lesson07_auc_results.csv",
  
  row.names = FALSE
)

#------------------------------------------------------------------------------

# Checking final outputs:
slope_results
auc_results
auc_comparison
