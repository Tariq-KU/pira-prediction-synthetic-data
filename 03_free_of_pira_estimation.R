
 
library(survival) #to use Surv() to combines follow-up time and event status & survfit() to calculates the Kaplan–Meier estimate.

dev <- read.csv("data/development.csv")

stopifnot(  # stops the analysis if a condition fails. If all checks pass, it prints nothing.
  !anyNA(dev$followup_months),
  !anyNA(dev$pira_event),
  all(dev$followup_months > 0 & dev$followup_months <= 24), #this makes sure any PIRA developed is at max at study month 36 so we have enough time to confirm it
  all(dev$pira_event %in% c(0, 1)) #PIRA events should be either 0 or 1
)

table(dev$pira_event)

table(dev$outcome_status) # to distinguish early censoring from completing the prediction window


# Now we create the survival outcome to estimate the survival curve for the outcome of this
pira_survival <- Surv(
  time = dev$followup_months,
  event = dev$pira_event,
  type = "right" # for right censoring
)

head(pira_survival, 12) #shows us the first 12 only, there are a lot


# To compare these symbols with the original records, run this
head(
  dev[c(
    "participant_id",
    "followup_months",
    "pira_event",
    "outcome_status"
  )],
  12
)

# this just bascically means the ratio of observed PIRA to total participants
observed_event_fraction <- mean(dev$pira_event)

100 * observed_event_fraction


# Now we can calculate the Kaplan-Meier Curve
km_fit <- survfit(
  Surv(
    time = followup_months,
    event = pira_event,
    type = "right"
  ) ~ 1,
  data = dev,
  conf.int = 0.95,       # Requests 95% confidence limits.
  conf.type = "log-log"  # Selects a method for calculating those limits within the probability range.
)

km_fit


# Display the risk set behind the curve (the table that we calculate from in my book)
km_summary <- summary(
  km_fit,
  times = c(6, 12, 18, 24)
)

risk_set_table <- data.frame(
  months_after_landmark = km_summary$time,
  study_month = km_summary$time + 12,
  at_risk = km_summary$n.risk,
  events_since_previous_time = km_summary$n.event,  # retrieves the risk-set counts.
  censored_since_previous_time = km_summary$n.censor
)

risk_set_table

# Calculation of the survuval rate
# The first calculation produces the fraction avoiding an event at each event time:
conditional_survival <- 1 -
  risk_set_table$events_since_previous_time /
  risk_set_table$at_risk

# cumulative product which is required when we have censoring
survival_by_hand <- cumprod(conditional_survival)

manual_check <- data.frame(
  month = risk_set_table$months_after_landmark,
  conditional_survival = conditional_survival,
  survival_by_hand = survival_by_hand,
  survival_from_R = km_summary$surv
)

round(manual_check, 4)


# Convert event-free survival into PIRA risk. This gives us the S(t) for survival along with the 95% CI
km_estimates <- data.frame(
months_after_landmark = km_summary$time,
pira_free_percent = 100 * km_summary$surv,
pira_risk_percent = 100 * (1 - km_summary$surv),
risk_lower_95_percent = 100 * (1 - km_summary$upper),
risk_upper_95_percent = 100 * (1 - km_summary$lower)
)

round(km_estimates, 2)



# ------------------------------------------------------------------
# Plot the development-cohort Kaplan–Meier curve
# ------------------------------------------------------------------

# km_fit must be the original DEVELOPMENT-cohort Kaplan–Meier fit.

dev_km_at_horizon <- summary(
  km_fit,
  times = 24,
  extend = FALSE
)

draw_KM <- function() {
  
  previous_settings <- par(no.readonly = TRUE)
  on.exit(par(previous_settings))
  
  par(
    mar = c(5.0, 5.0, 3.0, 1.2),
    mgp = c(3.0, 0.8, 0),
    las = 1
  )
  
  plot(
    km_fit,
    conf.int = TRUE,
    conf.type = "log-log",
    mark.time = TRUE,
    pch = 3,
    mark.col = "#D55E00",    # Contrasting orange
    cex = 2,              # Larger censoring symbols
    col = "#236B85",
    lwd = 2,
    xlim = c(0, 24),
    ylim = c(0.75, 1),
    yaxt = "n",          # Add our own y-axis labels below
    xaxs = "i",
    xaxt = "n",
    xlab = "Months after the month-12 prediction landmark",
    ylab = "Estimated probability of remaining PIRA-free",
    main = "Kaplan–Meier curve: synthetic development cohort"
  )
  
  # Match the six-month assessment intervals.
  axis(
    side = 1,
    at = seq(0, 24, by = 6)
  )
  
  # Show PIRA-free probabilities from 0.75 to 1.00.
  axis(
    side = 2,
    at = seq(0.75, 1.00, by = 0.05),
    labels = sprintf("%.2f", seq(0.75, 1.00, by = 0.05)),
    las = 1
  )
  
  # Calculate annotations from the fitted curve.
  legend(
    "bottomleft",
    legend = c(
      paste0("Participants: ", sum(km_fit$n)),
      sprintf(
        "24-month PIRA risk: %.2f%%",
        100 * (1 - dev_km_at_horizon$surv)
      ),
      "Confidence limits: 95%; + marks: censoring"
    ),
    bty = "n",
    cex = 0.85
  )
}

# Display in the RStudio Plots pane.
draw_KM()


# ------------------------------------------------------------------
# Save versions matching the test-cohort figure
# ------------------------------------------------------------------

dir.create("figures", showWarnings = FALSE)

# Vector PDF: remains sharp when resized.
pdf(
  file = "figures/PIRA_KM_development_cohort.pdf",
  width = 6.8,
  height = 5.2,
  pointsize = 11,
  useDingbats = FALSE
)

draw_KM()
dev.off()

# High-resolution PNG for Word.
png(
  filename = "figures/PIRA_KM_development_cohort.png",
  width = 6.8,
  height = 4.3,
  units = "in",
  res = 600,
  pointsize = 11,
  bg = "white"
)

draw_KM()
dev.off()

# Show the folder containing the exported figures.
normalizePath("figures")
#------------------------------------------------------------------------------

# Create the simple benchmark for later lessons
# This benchmark must come from the development data and remain fixed during validation. Our main scientific comparison will still be clinical–MRI versus clinical–MRI plus blood biomarkers.
km_24 <- summary(km_fit, times = 24)

benchmark_risk_24 <- as.numeric(
  1 - km_24$surv
)

comparison <- data.frame(
  measure = c(
    "Observed event fraction",
    "Kaplan-Meier estimated 24-month risk"
  ),
  percent = 100 * c(
    observed_event_fraction,
    benchmark_risk_24
  )
)

comparison


# Save the results for later lessons
# saveRDS(): preserves an R object, including the fitted curve and its components. You can reload it in a later session

dir.create("results", showWarnings = FALSE)

write.csv(
  risk_set_table,
  "results/lesson03_risk_sets.csv",
  row.names = FALSE
)

write.csv(
  km_estimates,
  "results/lesson03_km_estimates.csv",
  row.names = FALSE
)

saveRDS(
  km_fit,
  "results/lesson03_km_fit.rds"
)

saveRDS(
  benchmark_risk_24,
  "results/lesson03_benchmark_risk_24.rds"
)

