#random_netcal_model_betas_pro function
#Version: 7.1.1
#Last_update:2024/10/15
#Programmer: David_P_Quevedo

random_netcal_model_betas_pro <- function(iterations) { #pro version (using clusters, parallel and monitoring with progress bar)
  
   start_time <- Sys.time()
  
  # Set up parallel processing
  num_cores <- detectCores() - 2  # My laptop has 22 cores, is better to work in pairs, so I leave 2 free
  cl <- makeCluster(num_cores)
  
  # Export necessary data and functions to each node
  clusterExport(cl, varlist = c("data_to_calculate_net_data", 
								                "shuffle_sexsel_popsub_sex_strat",
                                "net_data_calculation",
                                "fit_lmer_and_extract_estimates",
                                "load_or_install"))
  
  # Load necessary libraries on each worker node
  clusterEvalQ(cl, {
    options(repos = c(CRAN = "https://cloud.r-project.org"))  # Ensure CRAN mirror is set
    libraries_to_run<-c("dplyr","tidyr","tidyverse","igraph","tnet","purrr","lmerTest","car","glmmTMB")
    sapply(libraries_to_run,load_or_install)
  })

  # Initialize lists to store results
  betas_for_indv <- list()  # Use lists for dynamic sizing
  betas_for_pop <- list()    # Use lists for dynamic sizing
  
  # To track successful iterations
  successful_iterations <- 0
  max_attempts <- iterations * 1  # Maximum attempts to avoid infinite loop
  
  # Use pblapply for progress bar
  pboptions(type = "txt")
  results <- pblapply(1:max_attempts, function(i) {
    suppressWarnings({
      tryCatch({
        # Perform the swap operation
        df_new <- shuffle_sexsel_popsub_sex_strat(data_to_calculate_net_data)
        
        # Calculate the net data 
        net_list <- split(df_new, df_new$video)
        df_net <- lapply(net_list, net_data_calculation)
        df_net <- do.call(rbind, df_net)
        
        # Structure data for individual and populational variables 
        df_net_ind <- df_net %>% select(-c(gm_bcd, gm_bcf, bi_bcd, bi_bcf, 
                                             density_bcd, density_bcf, 
                                             gm_divby_bi_bcd, gm_divby_bi_bcf))
        df_net_pop <- df_net %>% select(-c(degree_bcd, degree_bcf, 
                                             betweenness_bcd, betweenness_bcf, 
                                             entropy_bcd, entropy_bcf, 
                                             strength_bcd, strength_bcf,
                                             colour, id, sex, velocity, 
                                             totaldistmov.mm, movingtime, 
                                             interpolation, sbnf)) %>%
                       distinct(video, .keep_all = TRUE)
        
        # Check for empty data frames and skip if empty
        if (nrow(df_net_ind) == 0 || nrow(df_net_pop) == 0) {
          return(NULL)  # Return NULL to skip this iteration
        }
        
        # Run the models and extract the betas
        results_for_individual <- suppressWarnings(fit_lmer_and_extract_estimates(df_net_ind, mode = "individual"))
        results_for_individual$iteration <- successful_iterations + 1
        
        results_for_populational <- suppressWarnings(fit_lmer_and_extract_estimates(df_net_pop, mode = "populational"))
        results_for_populational$iteration <- successful_iterations + 1
        
        successful_iterations <<- successful_iterations + 1  # Increment successful iteration count
        
        list(individual = results_for_individual, populational = results_for_populational)
        
      }, error = function(e) {
        message(paste("Error in iteration", i, ":", e$message))
        return(NULL)  # Return NULL in case of error
      })
    })
  }, cl = cl)

  # Stop the cluster
  stopCluster(cl)
  
  # Filter out NULL results
  results <- Filter(Negate(is.null), results)
  
  # Extract results into dynamic lists
  for (result in results) {
    if (!is.null(result$individual)) {
      betas_for_indv[[length(betas_for_indv) + 1]] <- result$individual
    }
    if (!is.null(result$populational)) {
      betas_for_pop[[length(betas_for_pop) + 1]] <- result$populational
    }
  }
  
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  cat("Elapsed time:", elapsed_time, "\n")
  
  return(list(betas_for_indv = betas_for_indv, betas_for_pop = betas_for_pop))
}


#End here

