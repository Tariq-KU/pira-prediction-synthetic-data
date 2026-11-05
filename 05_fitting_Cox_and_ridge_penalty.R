# we prepared the predictors so that the model could use them. We now need to learn the relationship between those predictors and the timing of PIRA.

# By the end of this lesson, you will have four fitted models, each with a penalty selected through cross-validation.
#You will also see how a Cox model produces an individual’s predicted 24-month risk.

install.packages("glmnet")

library(survival) # Represent censored outcomes and calculate survival predictions
library(recipes)  # Learn and apply preprocessing rules
library(glmnet)   # Fit ridge-penalised Cox models


dev <- read.csv("data/development.csv")
# for this one, we need the make_features function to be loaded from lesson 4
# remember, this is the function for the log2 conversion
dev_features <- make_features(dev)

nrow(dev_features)
table(dev_features$pira_event)


# Now for the ridge penalization, where we will try a range of strengths:

# This creates 60 values from 100 down to 0.0001, spaced evenly on a logarithmic scale.
lambda_grid <- 10^seq(
  from = 2,
  to = -4,
  length.out = 60
)

head(lambda_grid)
tail(lambda_grid)


# Divide the development participants into five folds (fold=group of participants)
#In each cross-validation round:
# 1. Four folds provide the training participants.
# 2. The remaining fold provides participants for checking the candidate models.
# 3. We repeat until every fold has served as the checking group.

# We will distribute event and non-event participants approximately evenly across the folds.

set.seed(2026)

K <- 5L

# fold_id stores one fold number for each participant.
# which(...) finds the rows belonging to the current event group.
# rep(...) creates a roughly balanced allocation to folds 1–5
# sample(...) shuffles that allocation.

# I don't exactly understand how this works

fold_id <- integer(nrow(dev_features))

for (event_group in c(0, 1)) {
  
  rows <- which(   # stores the row where we start from PIRA=0
    dev_features$pira_event == event_group
  )
  
  fold_id[rows] <- sample(
    rep(seq_len(K), length.out = length(rows))
  )
}

table(
  fold = fold_id,
  event = dev_features$pira_event
)

# We will use the same folds for all four models, so their tuning uses the same participant divisions.



# Build one reusable function for preprocessing and fitting
# The next function does three connected jobs:
# 1. Learns preprocessing from the supplied training participants.
# 2. Creates their predictor matrix and survival outcome.
# 3. Fits a ridge Cox model for every lambda in our grid.

fit_ridge_path <- function(training_data, variables, lambdas) {
  
  # Select only the prespecified predictors
  raw_x <- training_data[variables]
  
  # Learn preprocessing from these training participants
  trained_recipe <- prep(
    make_preprocessor(raw_x),
    training = raw_x
  )
  
  # Apply those rules and return a numeric matrix
  x <- bake(
    trained_recipe,
    new_data = raw_x,
    composition = "matrix"
  )
  
  # Create the censored survival outcome
  y <- with(
    training_data,
    Surv(followup_months, pira_event, type = "right")
  )
  
  # Fit the ridge Cox models across the lambda grid
  fit <- glmnet(
    x = x,
    y = y,
    family = "cox",
    alpha = 0,
    lambda = lambdas,
    standardize = FALSE,
    cox.ties = "breslow"
  )
  
  # Keep everything required to use the fitted model later
  list(
    recipe = trained_recipe,
    variables = variables,
    x_train = x,
    y_train = y,
    fit = fit
  )
}


# family = "cox"	Fit a Cox survival model
# alpha = 0	Use ridge penalisation
# lambda = lambdas	Fit the supplied penalty strengths
# standardize = FALSE	Avoid scaling again; our recipe already scales continuous predictors
# cox.ties = "breslow"	Use a specified method for events recorded at the same time


# Inspect one cross-validation round before automating it (temporarily hold out fold 1)
training_data <- dev_features[fold_id != 1, ]
checking_data <- dev_features[fold_id == 1, ]

# here we are applying the function we just created
example_fit <- fit_ridge_path(
  training_data = training_data,
  variables = predictor_groups$combined,
  lambdas = lambda_grid
)

# The model and preprocessing rules have now been learned using folds 2–5 with fold 1 held out

# Now, apply those same preprocessing rules to fold 1:

x_check <- bake(
  example_fit$recipe,
  new_data = checking_data[predictor_groups$combined],
  composition = "matrix"
)

# Why do we use bake() here instead of another prep()?
# Because fold 1 must receive the medians, means and standard deviations learned from folds 2–5.


y_check <- with(
  checking_data,
  Surv(followup_months, pira_event, type = "right")
)

dim(example_fit$x_train)
dim(x_check)

# Now calculate the checking participants’ linear predictors (what also is this??): 
eta_check <- predict(
  example_fit$fit,
  newx = x_check,
  s = lambda_grid,
  type = "link"
)

dim(eta_check)


# This returns a matrix:
# One row per checking participant.
# One column per lambda.
# type = "link" returns \(\eta\). For a Cox model, even type = "response" returns relative hazard, not an absolute event probability.



# Score the candidate penalties
# We need a score that can use event times and censoring. For tuning, we will use Cox partial-likelihood deviance.

fold_deviance <- glmnet::coxnet.deviance(
  pred = eta_check,
  y = y_check,
  std.weights = FALSE,
  cox.ties = "breslow"
)

head(fold_deviance)

# This produces one score for each lambda. Lower is better.


# Repeat the process across all five folds

tune_ridge <- function(data, variables, fold_id, lambdas) {
  
  fold_numbers <- sort(unique(fold_id))
  
  fold_loss <- matrix(
    NA_real_,
    nrow = length(fold_numbers),
    ncol = length(lambdas)
  )
  
  for (k in fold_numbers) {
    
    training_data <- data[fold_id != k, ]
    checking_data <- data[fold_id == k, ]
    
    bundle <- fit_ridge_path(
      training_data,
      variables,
      lambdas
    )
    
    x_check <- bake(
      bundle$recipe,
      new_data = checking_data[variables],
      composition = "matrix"
    )
    
    y_check <- with(
      checking_data,
      Surv(followup_months, pira_event, type = "right")
    )
    
    eta <- predict(
      bundle$fit,
      newx = x_check,
      s = lambdas,
      type = "link"
    )
    
    fold_loss[k, ] <- glmnet::coxnet.deviance(
      pred = eta,
      y = y_check,
      std.weights = FALSE,
      cox.ties = "breslow"
    )
  }
  
  cv_score <- colSums(fold_loss) / sum(data$pira_event)
  
  best_index <- which.min(cv_score)
  
  if (best_index %in% c(1L, length(lambdas))) {
    warning(
      "Selected lambda is at a grid boundary. Review and expand the grid."
    )
  }
  
  list(
    results = data.frame(
      lambda = lambdas,
      cv_deviance_per_event = cv_score
    ),
    best_lambda = lambdas[best_index],
    fold_loss = fold_loss
  )
}


combined_tuning <- tune_ridge(
  data = dev_features,
  variables = predictor_groups$combined,
  fold_id = fold_id,
  lambdas = lambda_grid
)

combined_tuning$best_lambda


#------------------------------------------------------------------------------
# Inspect its tuning curve, The red line marks the selected penalty.

draw_ridgee<- function() {
plot(
  log10(combined_tuning$results$lambda),
  combined_tuning$results$cv_deviance_per_event,
  type = "b",
  pch = 16,
  xlab = "log10(lambda): larger values mean stronger shrinkage",
  ylab = "Cross-validated deviance per event",
  main = "Choosing the ridge penalty"
)

abline(
  v = log10(combined_tuning$best_lambda),
  col = "red",
  lty = 2
)

}

draw_ridgee();

# TO save a high resulotion for the report:
dir.create("figures", showWarnings = FALSE)

png(
  filename = "figures/Ridge Penalty.png",
  width = 2000,
  height = 2000,
  res = 300
)

draw_ridgee()

dev.off()

#------------------------------------------------------------------------------


# We have tuned the combined model. Now tune the other three using the same folds

tuning <- list(
  combined = combined_tuning
)

other_models <- setdiff(
  names(predictor_groups),
  "combined"
)

for (model_name in other_models) {
  
  tuning[[model_name]] <- tune_ridge(
    data = dev_features,
    variables = predictor_groups[[model_name]],
    fold_id = fold_id,
    lambdas = lambda_grid
  )
}


# Each model gets its own selected lambda. Different predictor sets may require different amounts of shrinkage.


# Now refit each model using all development participants:

models <- list()

for (model_name in names(predictor_groups)) {
  
  bundle <- fit_ridge_path(
    training_data = dev_features,
    variables = predictor_groups[[model_name]],
    lambdas = lambda_grid
  )
  
  bundle$selected_lambda <- tuning[[model_name]]$best_lambda
  
  models[[model_name]] <- bundle
}

# Why fit again?
# The cross-validation models each learned from approximately 80% of the development participants. Their role was to select lambda.
# Having selected lambda, we can use all available development participants to estimate the model that will later be evaluated independently. Its preprocessing is also relearned from the full development cohort.
# We retain the fitted lambda path, but use the selected lambda whenever extracting coefficients or predictions.



# Inspect the resulting model setting (return the predictor columns and selected lambda for each model)
# So here, we are testing for all 4 models right? (yes, run "model_overview" below)
model_overview <- data.frame(
  model = names(models),
  
  predictor_columns = vapply(
    models,
    function(m) ncol(m$x_train),
    integer(1)
  ),
  
  selected_lambda = vapply(
    models,
    function(m) m$selected_lambda,
    numeric(1)
  ),
  
  row.names = NULL
)

model_overview


# Inspect the learned coefficients (find coefficinets and hazard ratio for the combined model)
# This produces the coeff for the combined model only.
# You can change the first line to inspect another model (models$clinical, models$clinical_mri, models$clinical_blood).
combined_model <- models$combined

beta <- as.matrix(
  coef(
    combined_model$fit,
    s = combined_model$selected_lambda
  )
)

coefficient_table <- data.frame(
  predictor = rownames(beta),
  coefficient = as.numeric(beta[, 1]),
  hazard_ratio = exp(as.numeric(beta[, 1])),
  row.names = NULL
)

coefficient_table



# Convert the fitted model into a 24-month risk prediction. The model needs two components to calculate absolute risk:
# 1. The participant’s predictor score.
# 2. The baseline survival estimated from the development outcomes.

# We will demonstrate the calculation using the first three development participants:
example_people <- dev_features[1:3, ]

example_x <- bake(
  combined_model$recipe,
  new_data = example_people[combined_model$variables],
  composition = "matrix"
)

# Calculate their survival curves:
example_survival <- survival::survfit(
  combined_model$fit,
  s = combined_model$selected_lambda,
  x = combined_model$x_train,
  y = combined_model$y_train,
  newx = example_x,
  stype = 2,
  ctype = 1
)

# x and y supply the original training data needed to estimate baseline survival
# newx supplies the predictors for the people receiving predictions.
# stype = 2 and ctype = 1 specify survival calculated from the cumulative hazard using the Breslow baseline calculation, consistent with our fitting choice

# Now, Extract survival at 24 months and subtract it from one:
survival_24 <- as.numeric(
  summary(example_survival, times = 24)$surv
)

example_predictions <- data.frame(
  participant_id = example_people$participant_id,
  predicted_pira_risk_24m = 1 - survival_24
)

example_predictions

# Note: These three predictions demonstrate the calculation. They are not validation results because these participants helped train the model.


# Save the fitted models and tuning decisions
dir.create(
  "results",
  showWarnings = FALSE
)

saveRDS(
  list(
    models = models,
    tuning = tuning,
    fold_id = fold_id,
    lambda_grid = lambda_grid,
    predictor_groups = predictor_groups,
    
    R_version = R.version.string,
    
    package_versions = vapply(
      c("survival", "recipes", "glmnet"),
      function(pkg) as.character(packageVersion(pkg)),
      character(1)
    )
  ),
  
  file = "results/lesson05_ridge_models.rds"
)

write.csv(
  coefficient_table,
  "results/lesson05_combined_coefficients.csv",
  row.names = FALSE
)

model_overview

head(coefficient_table)

example_predictions