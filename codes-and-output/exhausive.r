



N <- 2^15



# Load necessary libraries
library(survival)
library(MASS)


# Load dataset
survival_data_adjusted <- read.csv("survival_data_adjusted.csv")

# Display first few rows to verify\head(survival_data_adjusted)

# Create survival object with entry time (left truncation)
surv_obj <- with(survival_data_adjusted, Surv(time = entry_time, 
                                              time2 = survival_time, 
                                              event = status))


head(survival_data_adjusted)



# Fit initial Cox model (full model with all predictors)
# Add the columns 'adjusted_entry_time', 'survival_time', and 'status' to the dataset
survival_data_adjusted <- survival_data_adjusted[, c(#"adjusted_entry_time", "survival_time", "status", 
                                                     names(survival_data_adjusted)[grepl("^X", names(survival_data_adjusted))])]

# Fit the Cox model
cox_min <- coxph(surv_obj ~ X3+X4, data = survival_data_adjusted)


library(ggplot2)
library(broom)

coef_df <- tidy(cox_min, conf.int = TRUE)

# # Plot
# ggplot(coef_df, aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high)) +
#   geom_pointrange() +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   coord_flip() +
#   labs(title = "Coefficient Estimates with 95% Confidence Intervals",
#        x = "Covariates", y = "Log Hazard Ratio (coef)") +
#   theme_minimal()



# Load necessary library
library(survival)
library(compiler)
compiler::enableJIT(3)
# Assuming 'survival_data_adjusted' is your dataset and 'surv_obj' is your Surv object
# Extract predictor variable names (e.g., those starting with "X")
predictor_names <- grep("^X", names(survival_data_adjusted), value = TRUE)


set.seed(123)  # For reproducibility







library(survival)
library(progress)

# Predictor variable names
predictor_names <- grep("^X", names(survival_data_adjusted), value = TRUE)

num_predictors <- length(predictor_names)
total_combinations <- 2^num_predictors

# Progress bar initialization
pb <- progress_bar$new(total = total_combinations, format = "[:bar] :percent eta: :eta")

# Preallocate storage
results <- vector("list", total_combinations)

# Main loop to iterate through all possible combinations
for (i in 0:(total_combinations - 1)) {
  pb$tick()

  # Convert integer to binary selection vector
  binary_vector <- as.integer(intToBits(i))[1:num_predictors]

  selected_vars <- predictor_names[which(binary_vector == 1)]

  if (length(selected_vars) == 0) {
    formula_str <- "surv_obj ~ 1"
  } else {
    formula_str <- paste("surv_obj ~", paste(selected_vars, collapse = " + "))
  }

  model <- tryCatch(
    coxph(as.formula(formula_str), data = survival_data_adjusted),
    error = function(e) NULL
  )

  aic_value <- if (!is.null(model)) extractAIC(model)[2] else NA

  # Store results
  results[[i + 1]] <- list(
    combination = paste(binary_vector, collapse = ""),
    formula = formula_str,
    AIC = aic_value
  )

  # Periodic garbage collection
  if ((i + 1) %% 100 == 0) {
    gc(verbose = FALSE)
  }
}

# Convert list to data frame
results_df <- do.call(rbind, lapply(results, data.frame, stringsAsFactors = FALSE))

# Clean and sort results
results_clean <- na.omit(results_df)
results_sorted <- results_clean[order(results_clean$AIC), ]

# Top 5 combinations
top_5 <- head(results_sorted, 5)
print(top_5)

# Save top 5 to CSV
write.csv(top_5, file = "top_5_cox_models_exhausive.csv", row.names = FALSE)
cat("Top 5 Cox model combinations saved to 'top_5_cox_models.csv'\n")