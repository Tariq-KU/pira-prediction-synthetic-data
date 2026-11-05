# In Lesson 5, we developed four models and selected their ridge penalties. Now we ask:
# How much better might these models appear to perform because we evaluate them on the same participants used to develop them?


# Suppose we evaluate our combined model on the original development participants and obtain a C-index of 0.76.
# Because those participants helped determine the preprocessing, penalty and coefficients, that score may overstate performance in new participants. We call this apparent performance.

# C-index examines pairs of participants whose outcomes allow us to compare their event ordering.

# We will use Harrell’s C-index:
#    0.5: no useful ordering beyond chance.
#    1.0: perfect ordering of evaluable pairs.
#    Higher values indicate better discrimination.

library(survival)
library(recipes)
library(glmnet)

lesson5 <- readRDS(
  "results/lesson05_ridge_models.rds"
)

models <- lesson5$models
predictor_groups <- lesson5$predictor_groups
lambda_grid <- lesson5$lambda_grid
fold_id <- lesson5$fold_id

dev <- read.csv("data/development.csv")
dev_features <- make_features(dev)



# We will use the created functions "make_features" & "make_preprocessor" & "fit_ridge_path" & "tune_ridge"
stopifnot(
  nrow(dev_features) == length(fold_id),
  anyDuplicated(dev_features$participant_id) == 0L
)

# define a function that calculates a model’s C-index on a supplied dataset:
score_c_index <- function(model_bundle, evaluation_data) {
  
  # Apply this model's existing preprocessing rules
  x_evaluation <- bake(
    model_bundle$recipe,
    new_data = evaluation_data[model_bundle$variables],
    composition = "matrix"
  )
  
  # Obtain one risk score per participant
  eta <- as.numeric(
    predict(
      model_bundle$fit,
      newx = x_evaluation,
      s = model_bundle$selected_lambda,
      type = "link"
    )
  )
  
  # Keep predictions and outcomes aligned
  score_data <- data.frame(
    time = evaluation_data$followup_months,
    event = evaluation_data$pira_event,
    eta = eta
  )
  
  result <- survival::concordance(
    Surv(time, event) ~ eta,
    data = score_data,
    reverse = TRUE,
    timewt = "n",
    na.action = na.fail
  )
  
  as.numeric(result$concordance)
}

# The above function does three things: applies the model’s preprocessing, obtains its predictions, and compares their ordering with the observed outcomes.


# Now, we calculate the apparent performance (this is the discrimination of the model on the development data)

model_names <- names(predictor_groups)

apparent_c <- setNames(
  numeric(length(model_names)),
  model_names
)

for (model_name in model_names) {
  
  apparent_c[model_name] <- score_c_index(
    model_bundle = models[[model_name]],
    evaluation_data = dev_features
  )
}

apparent_c  # this shows as the c-index across the 4 models



# Create one bootstrap sample
# A bootstrap sample is created by sampling participants with replacement.

# Because our participants come from several centres, we will sample within each centre, preserving the original number of participants contributed by each centre.

draw_bootstrap_rows <- function(data) {
  
  sampled_rows <- integer(0)
  
  for (centre in unique(data$centre_id)) {
    
    centre_rows <- which(data$centre_id == centre)
    
    selected_rows <- sample(
      centre_rows,
      size = length(centre_rows),
      replace = TRUE
    )
    
    sampled_rows <- c(sampled_rows, selected_rows)
  }
  
  sampled_rows
}

# The This function above returns row numbers. We use those row numbers to construct the bootstrap dataset:

set.seed(6001)

boot_rows <- draw_bootstrap_rows(dev_features)

boot_data <- dev_features[boot_rows, ]

nrow(boot_data)

length(unique(boot_rows)) #this should return lower rows because some are repeats thus not unique (sanity check)


# Keep repeated participants together during cross-validation

# We must repeat penalty selection within the bootstrap sample.
#However, suppose participant A appears three times. If one copy enters a checking fold while another enters its training set, the model has already seen the same participant.
# To prevent that, each copy retains its original fold assignment:

boot_folds <- fold_id[boot_rows]

table(
  fold = boot_folds,
  event = boot_data$pira_event
)

# Now, all copies of a given participant therefore remain together.


# 7.Develop a new combined model in this bootstrap sample

# Now is the time to develop a new model on the bootstrapped data.
# 1st,select its penalty using only the bootstrap data:
boot_tuning <- tune_ridge(
  data = boot_data,
  variables = predictor_groups$combined,
  fold_id = boot_folds,
  lambdas = lambda_grid
)

# Above, this repeats the five-fold process from Lesson 5. Within every round, preprocessing is learned from that round’s training participants.


# Next, fit the model using the full bootstrap sample:
boot_model <- fit_ridge_path(
  training_data = boot_data,
  variables = predictor_groups$combined,
  lambdas = lambda_grid
)

boot_model$selected_lambda <- boot_tuning$best_lambda

# Now, The bootstrap model has its own imputation and scaling values, selected penalty and fitted coefficients.


## 8. Evaluate the same bootstrap model in two places
# First, score it on the bootstrap sample that trained it:
c_on_bootstrap <- score_c_index(
  model_bundle = boot_model,
  evaluation_data = boot_data
)

# Then score that same model, with the same preprocessing and coefficients, on the original development dataset:
c_on_original <- score_c_index(
  model_bundle = boot_model,
  evaluation_data = dev_features
)

# Finally, Calculate the difference (C-index):
optimism_example <- c_on_bootstrap - c_on_original

data.frame(
  C_on_bootstrap = c_on_bootstrap,
  C_on_original = c_on_original,
  optimism = optimism_example
)


## 9. Put one complete bootstrap repetition into a function
# We have now performed the procedure manually for the combined model.
# The following function repeats it for all four models, using the same bootstrap participants and fold assignments for each model.

one_bootstrap <- function(b) {
  
  # Give each repetition a reproducible random seed
  set.seed(6000 + b)
  
  rows <- draw_bootstrap_rows(dev_features)
  
  bootstrap_data <- dev_features[rows, ]
  bootstrap_folds <- fold_id[rows]
  
  # Check that each checking fold contains events
  events_by_fold <- tapply(
    bootstrap_data$pira_event,
    bootstrap_folds,
    sum
  )
  
  stopifnot(
    length(events_by_fold) == 5L,
    all(events_by_fold > 0)
  )
  
  result_rows <- list()
  
  for (model_name in model_names) {
    
    # Repeat penalty selection
    tuned <- tune_ridge(
      data = bootstrap_data,
      variables = predictor_groups[[model_name]],
      fold_id = bootstrap_folds,
      lambdas = lambda_grid
    )
    
    # Refit using the whole bootstrap sample
    fitted <- fit_ridge_path(
      training_data = bootstrap_data,
      variables = predictor_groups[[model_name]],
      lambdas = lambda_grid
    )
    
    fitted$selected_lambda <- tuned$best_lambda
    
    # Evaluate the same fitted model in both datasets
    c_training <- score_c_index(
      fitted,
      bootstrap_data
    )
    
    c_original <- score_c_index(
      fitted,
      dev_features
    )
    
    result_rows[[model_name]] <- data.frame(
      bootstrap = b,
      model = model_name,
      selected_lambda = tuned$best_lambda,
      C_on_bootstrap = c_training,
      C_on_original = c_original,
      optimism = c_training - c_original,
      lambda_at_boundary =
        tuned$best_lambda %in% range(lambda_grid)
    )
  }
  
  do.call(rbind, result_rows)
}


# Test one repetition:
first_bootstrap <- one_bootstrap(1)

first_bootstrap


## 10. Repeat the procedure
# Start with 20 repetitions to check that everything runs and to understand the output:
# This process is repeated with 20 repetitions across all 4 models
B <- 500L

bootstrap_runs <- vector("list", B)

dir.create(
  "results",
  showWarnings = FALSE
)

for (b in seq_len(B)) {
  
  bootstrap_runs[[b]] <- one_bootstrap(b)
  
  if (b %% 10 == 0 || b == B) {
    
    message(
      "Completed bootstrap ",
      b,
      " of ",
      B
    )
    
    saveRDS(
      list(
        completed = b,
        target = B,
        runs = bootstrap_runs
      ),
      "results/lesson06_bootstrap_progress.rds"
    )
  }
}

bootstrap_results <- do.call(
  rbind,
  bootstrap_runs
)

rownames(bootstrap_results) <- NULL

head(bootstrap_results)

# Twenty repetitions are for learning and checking the workflow.
# For the analysis run, change B to 500L and rerun this chunk. Guidance recommends at least 500 repetitions as a practical target for estimating optimism; computational time will depend on your machine and the repeated tuning.



# Calculate the average optimism and corrected C-index
# This also calculates the monte carlo error
summary_rows <- list()

for (model_name in model_names) {
  
  model_results <- bootstrap_results[
    bootstrap_results$model == model_name,
  ]
  
  mean_optimism <- mean(
    model_results$optimism
  )
  
  corrected_c <- apparent_c[[model_name]] -
    mean_optimism
  
  monte_carlo_se <- sd(
    model_results$optimism
  ) / sqrt(nrow(model_results))
  
  summary_rows[[model_name]] <- data.frame(
    model = model_name,
    apparent_C = apparent_c[[model_name]],
    mean_optimism = mean_optimism,
    corrected_C = corrected_c,
    optimism_MCSE = monte_carlo_se,
    bootstrap_repetitions = nrow(model_results),
    boundary_selections = sum(
      model_results$lambda_at_boundary
    )
  )
}

internal_validation <- do.call(
  rbind,
  summary_rows
)

rownames(internal_validation) <- NULL

internal_validation


saveRDS(
  list(
    apparent_C = apparent_c,
    bootstrap_results = bootstrap_results,
    summary = internal_validation,
    resampling = "Participants sampled within centres",
    fold_handling = "Original fold retained for every copy",
    session_info = sessionInfo()
  ),
  "results/lesson06_internal_validation.rds"
)

write.csv(
  internal_validation,
  "results/lesson06_internal_validation_summary.csv",
  row.names = FALSE
)