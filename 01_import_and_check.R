# Import the development dataset
dev <- read.csv(
  "data/development.csv",
  na.strings = c("", "NA"),  # Ensures treating blank cells and "NA" text as missing data
  stringsAsFactors = FALSE   # Imports text as standard characters, not categorical factors (should be done by default in 4.0.0)
)


dim(dev) # returns the number of rows and columns
names(dev) # shows you the names of all columns
View(dev)  # opens a spreadsheet view of the imported table + shows visualization of data distribution across columns

#--------------------------

# Check that participants are not duplicated
sum(duplicated(dev$participant_id))
# The $ means “take this column from this table.”, duplicated() identifies IDs that have already appeared.

# Examines the outcome labels (table counts how often each value occurs)
table(dev$pira_event)

# This gives more detail about the people who had no outcome (where they full 36month with no PIRA? or where they censored before the full 36month with no PIRA)
table(dev$outcome_status)

sum(is.na(dev$pira_by_24m)) # these relate to the censored before month 36

