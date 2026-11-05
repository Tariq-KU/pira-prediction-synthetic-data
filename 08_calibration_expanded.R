
# 1. Load the validation data and fixed predictions

library(survival)

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

prediction_names <- c(
  model_names,
  "benchmark"
)

stopifnot(
  nrow(val_features) ==
    nrow(validation_predictions),
  
  identical(
    as.character(val_features$participant_id),
    as.character(validation_predictions$participant_id)
  ),
  
  all(
    prediction_names %in%
      names(validation_predictions)
  )
)

#------------------------------------------------------------------------------

# 2. Create a function that estimates observed PIRA risk
# We will need the same calculation for the entire validation dataset and later for each risk group.

km_risk_at <- function(data, horizon = 24) {
  
  if (max(data$followup_months) < horizon) {
    stop(
      "This sample does not have follow-up reaching ",
      horizon,
      " months."
    )
  }
  
  km_fit <- survival::survfit(
    survival::Surv(
      followup_months,
      pira_event
    ) ~ 1,
    
    data = data,
    conf.int = 0.95,
    conf.type = "log-log"
  )
  
  km_at_horizon <- summary(
    km_fit,
    times = horizon,
    extend = FALSE
  )
  
  c(
    observed_risk =
      1 - as.numeric(km_at_horizon$surv),
    
    lower_95 =
      1 - as.numeric(km_at_horizon$upper),
    
    upper_95 =
      1 - as.numeric(km_at_horizon$lower)
  )
}

# The function first estimates the probability of remaining PIRA-free, then subtracts that probability from one.


# Now calculate the observed risk for the whole validation sample:
observed_overall <- km_risk_at(
  data = val_features,
  horizon = horizon
)

observed_overall
#------------------------------------------------------------------------------

# 3. Compare average predicted risk with observed risk

# For each model, we average its predictions across all validation participants, including those censored early.
# Why?? Because their predictions are available even when their eventual outcomes are not.

mean_predicted <- vapply(
  validation_predictions[prediction_names],
  FUN = mean,
  FUN.VALUE = numeric(1)
)

mean_calibration <- data.frame(
  model = prediction_names,
  
  predicted_risk_pct =
    100 * as.numeric(mean_predicted),
  
  observed_risk_pct =
    100 * unname(
      observed_overall["observed_risk"]
    ),
  
  difference_pp =
    100 * (
      as.numeric(mean_predicted) -
        unname(observed_overall["observed_risk"])
    ),
  
  row.names = NULL
)

print(
  mean_calibration,
  digits = 4,
  row.names = FALSE
)

#------------------------------------------------------------------------------

# How the calibration plot will work?
# We will now sort participants by a model’s predicted risk and divide them into approximately five equally sized groups.

# Within each group, we calculate:
  #  1. The average predicted risk.
  #  2. The Kaplan–Meier estimate of observed risk.
  #  3. A confidence interval for that observed risk.

#------------------------------------------------------------------------------

# 4. Calculate calibration within five risk groups

grouped_calibration <- function(
    data,
    predicted_risk,
    horizon = 24,
    groups = 5
) {
  
  stopifnot(
    length(predicted_risk) == nrow(data),
    all(is.finite(predicted_risk)),
    all(predicted_risk >= 0),
    all(predicted_risk <= 1)
  )
  
  # Define boundaries using the predicted-risk distribution.
  boundaries <- unique(
    as.numeric(
      quantile(
        predicted_risk,
        probs = seq(
          0,
          1,
          length.out = groups + 1
        )
      )
    )
  )
  
  if (length(boundaries) < 3) {
    stop(
      "There are too few distinct predictions ",
      "to create risk groups."
    )
  }
  
  # Assign participants to groups.
  risk_group <- cut(
    predicted_risk,
    breaks = boundaries,
    include.lowest = TRUE,
    labels = FALSE
  )
  
  results <- list()
  
  for (g in sort(unique(risk_group))) {
    
    rows <- which(risk_group == g)
    
    group_data <- data[
      rows, , drop = FALSE
    ]
    
    observed <- km_risk_at(
      data = group_data,
      horizon = horizon
    )
    
    results[[length(results) + 1L]] <- data.frame(
      group = g,
      
      participants = length(rows),
      
      observed_events = sum(
        group_data$pira_event == 1 &
          group_data$followup_months <= horizon
      ),
      
      early_censored = sum(
        group_data$pira_event == 0 &
          group_data$followup_months < horizon
      ),
      
      mean_predicted =
        mean(predicted_risk[rows]),
      
      observed_risk =
        unname(observed["observed_risk"]),
      
      lower_95 =
        unname(observed["lower_95"]),
      
      upper_95 =
        unname(observed["upper_95"])
    )
  }
  
  do.call(rbind, results)
}

#------------------------------------------------------------------------------
# 5. Apply this to the four models

calibration_groups <- setNames(
  vector(
    "list",
    length(model_names)
  ),
  model_names
)

for (model_name in model_names) {
  
  calibration_groups[[model_name]] <-
    grouped_calibration(
      data = val_features,
      
      predicted_risk =
        validation_predictions[[model_name]],
      
      horizon = horizon,
      
      groups = 5
    )
}

# Inspect the two models in our primary comparison:
calibration_groups$clinical_mri

calibration_groups$combined

# Group 1 contains the lowest predictions and group 5 the highest, assuming five groups were formed.


#------------------------------------------------------------------------------

# Draw the calibration plots

plot_grouped_calibration <- function(
    calibration_table,
    title
) {
  
  plot(
    x = 100 * calibration_table$mean_predicted,
    y = 100 * calibration_table$observed_risk,
    
    xlim = c(0, 100),
    ylim = c(0, 100),
    
    xlab = "Mean predicted PIRA risk (%)",
    ylab = "Observed PIRA risk (%)",
    
    main = title,
    
    pch = 19,
    col = "steelblue4"
  )
  
  # Perfect agreement.
  abline(
    a = 0,
    b = 1,
    lty = 2,
    col = "grey40"
  )
  
  # Draw confidence intervals where estimable.
  usable_ci <- with(
    calibration_table,
    is.finite(lower_95) &
      is.finite(upper_95) &
      upper_95 > lower_95
  )
  
  arrows(
    x0 = 100 * calibration_table$mean_predicted[usable_ci],
    y0 = 100 * calibration_table$lower_95[usable_ci],
    
    x1 = 100 * calibration_table$mean_predicted[usable_ci],
    y1 = 100 * calibration_table$upper_95[usable_ci],
    
    angle = 90,
    code = 3,
    length = 0.04,
    col = "steelblue4"
  )
}


# Display the 4 models:

model_labels <- c(
clinical = "Clinical",
clinical_mri = "Clinical + MRI",
clinical_blood = "Clinical + blood",
combined = "Clinical + MRI + blood"
)

old_par <- par(
  mfrow = c(2, 2),
  mar = c(4.5, 4.5, 3, 1)
)

for (model_name in model_names) {
  
  plot_grouped_calibration(
    calibration_table =
      calibration_groups[[model_name]],
    
    title =
      model_labels[[model_name]]
  )
}

par(old_par)


# The common axis limits make the panels directly comparable. You may need to enlarge the RStudio plot window.

#------------------------------------------------------------------------------

# 7. Save the calibration results:

dir.create(
"results",
showWarnings = FALSE
)

saveRDS(
  list(
    horizon = horizon,
    
    observed_overall = observed_overall,
    
    mean_calibration = mean_calibration,
    
    grouped_calibration = calibration_groups,
    
    requested_groups = 5,
    
    session_info = sessionInfo()
  ),
  
  file = "results/lesson07_part3a_calibration.rds"
)

write.csv(
  mean_calibration,
  
  file = "results/lesson07_mean_calibration.csv",
  
  row.names = FALSE
)


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

#### Calibration Slope & Smooth Calibration Curve

# 1. Load the data and check the required functions
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
  identical(
    as.character(val_features$participant_id),
    as.character(validation_predictions$participant_id)
  ),
  
  nrow(val_features) ==
    nrow(validation_predictions),
  
  all(model_names %in%
        names(validation_predictions))
)


# We will reuse:
  # make_ipcw_outcome() from Part 1.
  # draw_rows_within_centres() from Part 2.


# 2. Prepare the censoring weights
ipcw_outcome <- make_ipcw_outcome(
  data = val_features,
  horizon = horizon
)

table(ipcw_outcome$known)

# This reuses the outcome and censoring conventions established in Part 1.


# We now return to IPCW because the regression needs individual outcome records. In Part 3A, Kaplan–Meier directly estimated observed risk within groups.


#------------------------------------------------------------------------------
# 3. Define the diagnostic regression function
fit_calibration_diagnostics <- function(
    predicted_risk,
    outcome_data,
    risk_grid
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
  
  # Estimate the calibration slope.
  slope_fit <- glm(
    outcome ~ logit_prediction,
    
    data = calibration_data,
    
    weights = ipcw_weight,
    
    family = quasibinomial(link = "logit")
  )
  
  # Estimate a flexible calibration curve.
  curve_fit <- glm(
    outcome ~ splines::ns(
      logit_prediction,
      df = 3
    ),
    
    data = calibration_data,
    
    weights = ipcw_weight,
    
    family = quasibinomial(link = "logit")
  )
  
  stopifnot(
    slope_fit$converged,
    curve_fit$converged,
    
    all(is.finite(coef(slope_fit))),
    all(is.finite(coef(curve_fit)))
  )
  
  curve_estimates <- predict(
    curve_fit,
    
    newdata = data.frame(
      logit_prediction =
        qlogis(risk_grid)
    ),
    
    type = "response"
  )
  
  stopifnot(
    all(is.finite(curve_estimates))
  )
  
  list(
    slope = unname(
      coef(slope_fit)["logit_prediction"]
    ),
    
    joint_intercept = unname(
      coef(slope_fit)["(Intercept)"]
    ),
    
    curve = data.frame(
      predicted_risk = risk_grid,
      
      estimated_observed_risk =
        as.numeric(curve_estimates)
    )
  )
}

# This function will estimate both the slope and the smooth curve for one model.


# 4. Estimate the original slopes and curves

# We need a sequence of predicted risks at which to draw each curve. We will display the central 90% of each model’s predictions, from its 5th to 95th percentile, to avoid emphasising sparsely represented extremes.

risk_grids <- list()
original_calibration <- list()

for (model_name in model_names) {
  
  predictions <-
    validation_predictions[[model_name]]
  
  display_limits <- quantile(
    predictions,
    probs = c(0.05, 0.95),
    names = FALSE
  )
  
  risk_grids[[model_name]] <- seq(
    from = display_limits[1],
    to = display_limits[2],
    length.out = 100
  )
  
  original_calibration[[model_name]] <-
    fit_calibration_diagnostics(
      predicted_risk = predictions,
      
      outcome_data = ipcw_outcome,
      
      risk_grid = risk_grids[[model_name]]
    )
}

# Extract the slope estimates:
slope_estimates <- vapply(
  original_calibration,
  
  FUN = function(result) {
    result$slope
  },
  
  FUN.VALUE = numeric(1)
)

slope_estimates

#------------------------------------------------------------------------------

# 5. Bootstrap the slopes and curves

  # This follows the same validation-bootstrap logic you already learned.
   #For each repetition, we:
    #   1. Resample participants within centres.
    #   2. Keep their original predictions.
    #   3. Re-estimate censoring weights.
    #   4. Refit the diagnostic regressions.
    #   5. Store the slopes and curve estimates.

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

curve_bootstrap <- lapply(
  risk_grids,
  
  FUN = function(grid) {
    matrix(
      NA_real_,
      nrow = B_calibration,
      ncol = length(grid)
    )
  }
)

set.seed(7301)

for (b in seq_len(B_calibration)) {
  
  rows <- draw_rows_within_centres(
    val_features
  )
  
  bootstrap_data <- val_features[
    rows, , drop = FALSE
  ]
  
  bootstrap_outcomes <- make_ipcw_outcome(
    data = bootstrap_data,
    horizon = horizon
  )
  
  for (model_name in model_names) {
    
    result_b <- fit_calibration_diagnostics(
      predicted_risk =
        validation_predictions[[model_name]][rows],
      
      outcome_data =
        bootstrap_outcomes,
      
      risk_grid =
        risk_grids[[model_name]]
    )
    
    slope_bootstrap[b, model_name] <-
      result_b$slope
    
    curve_bootstrap[[model_name]][b, ] <-
      result_b$curve$estimated_observed_risk
  }
  
  if (b %% 100 == 0) {
    message(
      "Completed ",
      b,
      " of ",
      B_calibration,
      " calibration repetitions"
    )
  }
}

# This will take longer than bootstrapping Brier scores because we fit diagnostic regressions in every repetition.


#------------------------------------------------------------------------------

# 6. Create the calibration-slope results table
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
  
  bootstrap_SE =
    as.numeric(
      apply(
        slope_bootstrap,
        MARGIN = 2,
        FUN = sd
      )[model_names]
    ),
  
  row.names = NULL
)

slope_results

# The column meanings are the same as before, but the measured quantity is now the calibration slope.

#------------------------------------------------------------------------------
################
################
################

#This part can be excluded
# 7. Add uncertainty intervals to the smooth curves

# For every displayed risk value, we have 2,000 bootstrap estimates of the corresponding event probability.
# We take their percentiles separately at each risk value:
calibration_curves <- list()

for (model_name in model_names) {
  
  curve_limits <- apply(
    curve_bootstrap[[model_name]],
    
    MARGIN = 2,
    FUN = quantile,
    
    probs = c(0.025, 0.975),
    names = FALSE
  )
  
  curve_table <-
    original_calibration[[model_name]]$curve
  
  curve_table$lower_95 <-
    as.numeric(curve_limits[1, ])
  
  curve_table$upper_95 <-
    as.numeric(curve_limits[2, ])
  
  calibration_curves[[model_name]] <-
    curve_table
}

#------------------------------------------------------------------------------

# Plot the curves

plot_smooth_calibration <- function(
    curve_table,
    predictions,
    title
) {
  
  plot(
    0, 0,
    type = "n",
    
    xlim = c(0, 100),
    ylim = c(0, 100),
    
    xlab = "Predicted 24-month PIRA risk (%)",
    ylab = "Estimated observed PIRA risk (%)",
    
    main = title
  )
  
  # Pointwise confidence band.
  polygon(
    x = 100 * c(
      curve_table$predicted_risk,
      rev(curve_table$predicted_risk)
    ),
    
    y = 100 * c(
      curve_table$lower_95,
      rev(curve_table$upper_95)
    ),
    
    col = adjustcolor(
      "steelblue",
      alpha.f = 0.20
    ),
    
    border = NA
  )
  
  # Ideal agreement.
  abline(
    a = 0,
    b = 1,
    lty = 2,
    col = "grey40"
  )
  
  # Estimated calibration curve.
  lines(
    x = 100 * curve_table$predicted_risk,
    
    y = 100 *
      curve_table$estimated_observed_risk,
    
    col = "steelblue4",
    lwd = 2
  )
  
  # Show where participants' predictions lie.
  rug(
    100 * predictions,
    
    col = adjustcolor(
      "grey30",
      alpha.f = 0.35
    )
  )
}

###
# Display the 4 models:
model_labels <- c(
  clinical = "Clinical",
  clinical_mri = "Clinical + MRI",
  clinical_blood = "Clinical + blood",
  combined = "Clinical + MRI + blood"
)

old_par <- par(
  mfrow = c(2, 2),
  mar = c(4.5, 4.5, 3, 1)
)

for (model_name in model_names) {
  
  plot_smooth_calibration(
    curve_table =
      calibration_curves[[model_name]],
    
    predictions =
      validation_predictions[[model_name]],
    
    title =
      model_labels[[model_name]]
  )
}

par(old_par)

#------------------------------------------------------------------------------

# 9. Save the results:

dir.create(
  "results",
  showWarnings = FALSE
)

saveRDS(
  list(
    horizon = horizon,
    
    method = paste(
      "IPCW-weighted logistic calibration;",
      "natural cubic spline with df = 3"
    ),
    
    slope_results = slope_results,
    
    calibration_curves = calibration_curves,
    
    original_calibration = original_calibration,
    
    slope_bootstrap = slope_bootstrap,
    
    curve_bootstrap = curve_bootstrap,
    
    risk_grids = risk_grids,
    
    display_quantiles = c(0.05, 0.95),
    
    bootstrap_repetitions = B_calibration,
    
    seed = 7301,
    
    session_info = sessionInfo()
  ),
  
  file = "results/lesson07_part3b_calibration.rds"
)

write.csv(
  slope_results,
  
  file = "results/lesson07_calibration_slopes.csv",
  
  row.names = FALSE
)