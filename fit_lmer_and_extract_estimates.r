#fit_lmer_and_extract_estimates function
#Version: 2.1.2
#Last_update:2024/10/10
#Programmer: David_P_Quevedo

fit_lmer_and_extract_estimates <- function(data, mode = "mode") {
  # Identify columns ending with _bcd and _bcf
  bcd_columns <- grep("_bcd$", names(data), value = TRUE)
  bcf_columns <- grep("_bcf$", names(data), value = TRUE)
  
  # Combine the two sets of columns
  all_columns <- c(bcd_columns, bcf_columns)
  
  # Create a dataframe to store results
  results <- data.frame(variable = character(),
                        predictor = character(),
                        estimate = numeric(),
                        stringsAsFactors = FALSE)
  
  # Determine the formula based on the mode
  if (mode == "individual") {
    formula_template <- "~ sexselF * popsubF * sex + (1 | generation) + (1 | line.unique)"
  } else if (mode == "populational") {
    formula_template <- "~ sexselF * popsubF + (1 | generation) + (1 | line.unique)"
  } else {
    stop("Invalid mode. Choose either 'individual' or 'populational'.")
  }
  
  # Loop through each variable and fit the model
  for (var in all_columns) {
    # Calculate counts of zeros and non-zeros
    zero_counts <- data %>%
      count(is_zero = data[[var]] == 0)
    
    # Print zero counts for debugging
    # print(zero_counts)
    
    # Calculate the proportion of non-zeros
    true_count <- zero_counts %>% filter(is_zero == TRUE) %>% pull(n)
    false_count <- zero_counts %>% filter(is_zero == FALSE) %>% pull(n)
    total_count <- sum(true_count, false_count)
    
    # Calculate proportion
    proportion <- false_count / total_count
    
    # Print the proportion for debugging
    print(paste("Variable:", var, " - Proportion of non-zeros:", proportion))
    
    # Check if the proportion of non-zeros is greater than 0.2
    if (!is.na(proportion) && proportion > 0.4) {
      # Fit the glmmTMB model with warnings suppressed and control options
      formula <- as.formula(paste(var, formula_template))
      model <- suppressWarnings(glmmTMB(formula, data = data, family = gaussian(), zi = ~1, 
                                         REML = FALSE, 
                                         control = glmmTMBControl(optimizer = "nlminb", 
                                                                   optCtrl = list(maxit = 1000))))
      
      # Extract the estimates from the conditional model
      estimates <- fixef(model)$cond
    } else {
      # Fit the linear mixed model
      formula <- as.formula(paste(var, formula_template))
      model <- lmer(formula, data = data, REML = FALSE)
      # Extract the estimates
      estimates <- fixef(model)
    }
    
    # Store relevant estimates in the results dataframe
    # Skip the intercept
    predictors <- names(estimates)[-1]  # Get predictor names, skipping intercept
    estimates <- estimates[-1]  # Estimate values, skipping intercept
    
    results <- rbind(results, data.frame(variable = var,
                                         predictor = predictors,
                                         estimate = estimates,
                                         stringsAsFactors = FALSE))
  }
  
  # Pivot to wide format with predictors as columns and estimates under them
  results_wide <- results %>%
    pivot_wider(names_from = predictor, values_from = estimate, values_fill = NA)
  results_df_str<-as.data.frame(results_wide)
  return(results_df_str)
}
