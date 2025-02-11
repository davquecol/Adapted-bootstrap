#Date:2024_10_04
#Programmer: David_Quevedo 
#Operative system: Windows 11 Pro.
#Reviewer1(rvw1): 
#Reviewer2(rvw2): 
#Minor revisions comments rvw1:
#Minor revisions comments rvw2:

#################################################################################
############### Sexual Selection and Population Structure Effects ###############
#################################################################################

# Calling all the own functions we are going to need 

# Setting working directory for functions and sourcing them

set_working_directory <<- function(Folder, Subfolder, Group, Type) {
    home_dir <- path.expand("~")
    home_dir <- gsub("\\\\", "/", home_dir)
    parts <- strsplit(home_dir, "/")[[1]]
    Local_disk <- paste0(parts[1], "/")
    Users <- parts[2]
    Users_value <- parts[3]
    sep <- "/"
    Environment <- paste0(Local_disk, sep, Users, sep, Users_value, sep, Folder, sep, 
                          Subfolder, sep, Group, sep, Type)
    setwd(Environment)
  }
  
set_working_directory("Dropbox", "david_quevedo_SexNet", "thesis", "functions")

source_all_functions <- function(directory) {
  # Get the list of all .R files in the specified directory
  r_files <- list.files(directory, pattern = "\\.r$", full.names = TRUE)
  
  # Source each file
  for (file in r_files) {
    source(file)
  }
}

directory<-getwd()
source_all_functions(directory)

# Calling all the required libraries

libraries_for_parallel<-c("doParallel","foreach","pbapply","parallel")
libraries_for_modelling<-c("lmerTest", "lme4", "glmmTMB")
libraries_for_data_structure<-c("tidyr","dplyr","purrr","tnet")

sapply(c(libraries_for_parallel,libraries_for_modelling,libraries_for_data_structure),load_or_install)

#STEP 1: Call all the data from repository, structure it and check for directed links in 
#order to correct them and make them equal (for undirected networks). 

processed_data_from_ETXT <- start_progress()

#@ Pre-matrix comparison
# Number of pairs with the same data for bcd:  2408 
# Number of pairs with different data for bcd:  0 
# Number of pairs with the same data for bcf:  2408 
# Number of pairs with different data for bcf:  0 
#@ Post-matrix comparison
# Number of pairs with the same data for bcd:  2408 
# Number of pairs with different data for bcd:  0 
# Number of pairs with the same data for bcf:  2408 
# Number of pairs with different data for bcf:  0 

# This comparison shows that Ethovision gave us an undirected matrix (same link weight for both nodes)
# and that there is not any error after purging based on threshold (same number of pairs in tnat 
# format and matrix format after using the 11 seconds of bcd threshold to purgue). 

min_bcd_value <- processed_data_from_ETXT$min_bcd_value
bcd_data <- processed_data_from_ETXT$bcd_data
bcf_data <- processed_data_from_ETXT$bcf_data
both_gaths_filtered <- processed_data_from_ETXT$both_gaths_filtered
alldata <- processed_data_from_ETXT$alldata
data_to_calculate_net_data <- processed_data_from_ETXT$data_to_calculate_net
data_to_calculate_net_data <- data_to_calculate_net_data[order(data_to_calculate_net_data$video),]
data_to_calculate_net_data$sexsel<-as.factor(data_to_calculate_net_data$sexsel)
data_to_calculate_net_data$popsub<-as.factor(data_to_calculate_net_data$popsub)
data_to_calculate_net_data$generation<-as.factor(data_to_calculate_net_data$generation)
data_to_calculate_net_data$line.unique<-as.factor(data_to_calculate_net_data$line.unique)
rownames(data_to_calculate_net_data)<-1:nrow(data_to_calculate_net_data)


#STEP 2: Calculate individual and populational network data and structure it 
#in a dataframe with the original videorecording data, processing the variables of interest 
#with lmer models and keeping the betas 


net_list <- split(data_to_calculate_net_data, data_to_calculate_net_data$video)
observed_data_list <- lapply(net_list, net_data_calculation)
observed_data<-do.call(rbind,observed_data_list)

ind_data_observed <- observed_data[, -which(names(observed_data) %in% 
                                    c("gm_bcd", "gm_bcf", "bi_bcd", "bi_bcf", 
                                      "density_bcd", "density_bcf", 
                                      "gm_divby_bi_bcd", "gm_divby_bi_bcf"))]
pop_data_observed <- observed_data[, -which(names(observed_data) %in% 
                                    c("degree_bcd", "degree_bcf", "betweenness_bcd", "betweenness_bcf", 
                                      "entropy_bcd", "entropy_bcf", 
                                      "strength_bcd", "strength_bcf","colour","id","sex",
									                   "velocity","totaldistmov.mm","movingtime","interpolation","sbnf"))]
pop_data_observed <- pop_data_observed[!duplicated(pop_data_observed$video), ]


# Fitting lmer or glmmTMB models and extracting estimates0. 
# I set a mode (individual or populational) depending on the data income.

results_df_individual <- fit_lmer_and_extract_estimates_2.0 (ind_data_observed, mode = "individual")
results_df_populational <- fit_lmer_and_extract_estimates_2.0(pop_data_observed, mode = "populational")

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

iterations = 100

# Obtaining betas for all the variables from randomisations

betas_from_randomisations <- random_netcal_model_betas_pro(iterations = iterations)

# Function to extract betas for a specified variable from each dataset

extract_variable_betas <- function(betas_data, variable_name) {
  
  # Initialize an empty list to store the extracted results
  betas_list <- list()
  
  # Loop through individual and populational results
  for (i in 1:length(betas_data$betas_for_indv)) {
    # Extract betas for the specified variable from the individual results
    indv_result <- betas_data$betas_for_indv[[i]]
    pop_result <- betas_data$betas_for_pop[[i]]
    
    # Filter for the specified variable from individual results
    indv_betas <- indv_result[indv_result$variable == variable_name, , drop = FALSE]
    
    # Filter for the specified variable from populational results
    pop_betas <- pop_result[pop_result$variable == variable_name, , drop = FALSE]

    # Check if both results have data
    if (nrow(indv_betas) > 0 && nrow(pop_betas) > 0) {
      # Create a combined data frame with separate columns for individual and populational
      combined_betas <- data.frame(iteration = i, indv_betas, pop_betas[-1])  # Remove duplicate 'variable' column
      colnames(combined_betas) <- c("iteration", paste(colnames(indv_betas)[-1], "_indv", sep = ""), 
                                      paste(colnames(pop_betas)[-1], "_pop", sep = ""))
      betas_list[[i]] <- combined_betas
    } else if (nrow(indv_betas) > 0) {
      # Only individual result is available
      indv_betas_with_iter <- cbind(iteration = i, indv_betas)
      colnames(indv_betas_with_iter)[-1] <- paste(colnames(indv_betas_with_iter)[-1], "_indv", sep = "")
      betas_list[[i]] <- indv_betas_with_iter
    } else if (nrow(pop_betas) > 0) {
      # Only populational result is available
      pop_betas_with_iter <- cbind(iteration = i, pop_betas)
      colnames(pop_betas_with_iter)[-1] <- paste(colnames(pop_betas_with_iter)[-1], "_pop", sep = "")
      betas_list[[i]] <- pop_betas_with_iter
    }
  }
  
  # Combine the results into a single data frame, filtering out any empty entries
  final_betas_df <- do.call(rbind, betas_list)
  
  # Remove any rows that are completely NA
  final_betas_df <- final_betas_df[rowSums(is.na(final_betas_df)) != ncol(final_betas_df), , drop = FALSE]

  return(final_betas_df)
}

# Function to extract betas for all variables iterating the previous function across variable names

betas_for_all <- function(betas_data) {
  vars_indv <- betas_data$betas_for_indv[[1]]$variable
  vars_pop <- betas_data$betas_for_pop[[1]]$variable
  all_vars <- unique(c(vars_indv, vars_pop))  # Combine and ensure uniqueness

  # Initialize a list to hold results for all variables
  results_list <- list()

  # Loop over all variables and extract betas
  for (var in all_vars) {
    # Extract betas for the current variable
    extracted_betas <- extract_variable_betas(betas_data, var)

    # Check if the extracted betas data frame is not empty
    if (nrow(extracted_betas) > 0) {
      # Create a dynamic variable name for the results list
      var_name <- paste0(var, "_betas")
      
      # Store the results in the list with the dynamic name
      results_list[[var_name]] <- extracted_betas
    }
  }

  # Return the list of data frames
  return(results_list)
}

betas_all <- betas_for_all(betas_from_randomisations)


betas_combined_indv = do.call(rbind, betas_from_randomisations$betas_for_indv)

#Generating a dataframe for each variable 

strength_bcd_betas<-betas_all$strength_bcd_betas
strength_bcf_betas<-betas_all$strength_bcf_betas
degree_bcd_betas<-betas_all$degree_bcd_betas
degree_bcf_betas<-betas_all$degree_bcf_betas
density_bcd_betas<-betas_all$density_bcd_betas
gm_bcd_betas<-betas_all$gm_bcd_betas
gm_bcf_betas<-betas_all$gm_bcf_betas
betweenness_bcd_betas<-betas_all$betweenness_bcd_betas
betweenness_bcf_betas<-betas_all$betweenness_bcf_betas

minp <- function(sim,mod){
  pval1 <- 2*((1 + sum(sim >= mod))/(length(sim)+1))
  pval2 <- 2*((1 + sum(sim <= mod))/(length(sim)+1))
  min(pval1, pval2)
}

# Function to calculate the p-value from shuffled estimates and an empirical estimate
calculate_p_value <- function(shuffled_estimates, empirical_estimate) {
  # Check if the inputs are valid
  if (!is.numeric(empirical_estimate) || length(empirical_estimate) != 1) {
    stop("empirical_estimate should be a single numeric value.")
  }
  if (!is.numeric(shuffled_estimates) || length(shuffled_estimates) < 1) {
    stop("shuffled_estimates should be a numeric vector of length greater than 0.")
  }
  
  # Calculate the two-sided p-value
  p_value <- mean(abs(shuffled_estimates) >= abs(empirical_estimate))
  
  return(p_value)
}
