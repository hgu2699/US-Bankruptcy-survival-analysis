


N <- 2^15



# Load necessary libraries
library(survival)
library(MASS)


# Load dataset
survival_data_adjusted <- read.csv("survival_data_adjusted.csv")

# Display first few rows to verify\head(survival_data_adjusted)

# Create survival object with entry time (left truncation)
surv_obj <- with(survival_data_adjusted, Surv(time = adjusted_entry_time, 
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

# Initialize progress bar
pb <- progress_bar$new(total = N, format = "  [:bar] :percent eta: :eta")

# Initialize storage for results
results <- data.frame(
  combination = character(N),
  formula = character(N),
  AIC = numeric(N),
  stringsAsFactors = FALSE
)

set.seed(123)
gc_interval <- 100
for (i in 1:N) {
  # Update progress bar
  pb$tick()

  binary_vector <- rbinom(length(predictor_names), 1, 0.5)
  selected_vars <- predictor_names[which(binary_vector == 1)]

  if (length(selected_vars) == 0) {
    results$combination[i] <- paste(rep(0, length(predictor_names)), collapse = "")
    results$formula[i] <- "surv_obj ~ 1"
    results$AIC[i] <- NA
    next
  }

  formula_str <- paste("surv_obj ~", paste(selected_vars, collapse = " + "))
  model <- tryCatch(
    coxph(as.formula(formula_str), data = survival_data_adjusted),
    error = function(e) NULL
  )

  aic_value <- if (!is.null(model)) extractAIC(model)[2] else NA

  results$combination[i] <- paste(binary_vector, collapse = "")
  results$formula[i] <- formula_str
  results$AIC[i] <- aic_value

 if (i %% gc_interval == 0) {
    gc(verbose = FALSE)
  }
}

# After completion:
cat("\nCompleted feature selection.\n")









# Remove rows with NA AIC values
results_clean <- na.omit(results)

# Sort the results by AIC in ascending order
results_sorted <- results_clean[order(results_clean$AIC), ]

# Display the top 5 combinations
top_5 <- head(results_sorted, 5)
print(top_5)


# Save the top 5 combinations to a CSV file
write.csv(top_5, file = "top_5_cox_models.csv", row.names = FALSE)

# Confirmation message
cat("Top 5 Cox model combinations saved to 'top_5_cox_models.csv'\n")
