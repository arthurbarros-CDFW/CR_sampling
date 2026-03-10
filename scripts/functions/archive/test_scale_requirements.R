#' @title test_scale_requirements
#'
#' @description fill
#' 
#' @param N
#' 
#' @param theta 
#' 
#' @param k_iterations 
#' 
#' @param tag_rate
#' 
#' @param survey "hatchery" or "trib"
#' 
#' @param scales_n
#' 
#' @return fill
#' 

test_scale_requirements<-function(N,
                                  theta,
                                  k_iterations,
                                  scale_iterations,
                                  tag_rate,
                                  survey,
                                  scales_n_seq,
                                  n_replicates,
                                  CI,
                                  #distribution flags
                                  N_dist = "fixed",
                                  N_min = NULL, N_max = NULL,
                                  theta_dist = "fixed",
                                  theta_shape1 = 1, theta_shape2 = 1,
                                  theta_min = 0.1, theta_max = 0.5,
                                  tag_rate_dist = "fixed",
                                  tag_rate_shape1 = 1, tag_rate_shape2 = 1,
                                  age_probs_dist = "fixed",
                                  age_dirichlet_alpha = c(5, 5, 5, 5),
                                  hatchery_props_dist = "fixed",
                                  hatchery_dirichlet_alpha = c(56, 42, 2)
                                  ){
  #create results output to fill
  results_grid <- expand.grid(
    scales_n = scales_n_seq,
    replicate = 1:n_replicates,
    age = 2:5,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      estimated_prop_natural = NA,
      true_prop_natural = NA,
      bias = NA,
      CI=CI,
      upper_CI_natural = NA,
      lower_CI_natural = NA,
      ci_width = NA,
      coverage = NA,
      sim_time = NA,
      age_present_in_sample = FALSE,
      age_present_in_pop = FALSE,
      warning_message = NA_character_
    )
  
  #run for loops over scales_n_seq and n_replicates
  for(r in 1:n_replicates){
    #simulate population
    pop<-sim_pop(N=N,tag_rate=tag_rate)
    for(s in seq_along(scales_n_seq)){
      #run simulation with current scales_n
      time_start<-Sys.time()
      sim_result<-tryCatch({
        sim_wrapper(
          pop=pop,
          theta=theta,
          k_iterations=k_iterations,
          scale_iterations=scale_iterations,
          tag_rate=tag_rate,
          survey="trib",
          scales_n=scales_n_seq[s],
          CI=CI
        )
      },error=function(e){
        message(paste("Error in simulation:", e$message))
        return(NULL)
      })
      
      time_end<-Sys.time()
      sim_time<-time_end-time_start
      
      #skip if simulation failed
      if (is.null(sim_result)) {
        warning_msg <- "sim_wrapper failed"
        for (a in 2:5) {
          idx <- which(results_grid$scales_n == scales_n_seq[s] &
                         results_grid$replicate == r & 
                         results_grid$age == a)
          if (length(idx) > 0) {
            results_grid$warning_message[idx] <- warning_msg
            results_grid$sim_time[idx] <- sim_time
          }
        }
        print(paste("FAILED: scales_n =", scales_n_seq[s],
                    "boot iteration =", r))
        next
      }
      
      est_summary_stats=sim_result$est_summary_stats
      true_summary_stats=sim_result$true_summary_stats
      
      #store results
      for(a in 2:5){
        
        #get row index for values to fill
        idx<-which(results_grid$scales_n==scales_n_seq[s] &
                     results_grid$replicate==r & 
                     results_grid$age==a)
        
        if (length(idx) == 0) next
        
        est_row <- which(est_summary_stats$age == a)
        
        if (length(est_row) > 0) {
          results_grid$estimated_prop_natural[idx] <- est_summary_stats$mean_proportion_natural[est_summary_stats$age == a]
          results_grid$true_prop_natural[idx] <- true_summary_stats$true_proportion_natural[true_summary_stats$age == a]
          results_grid$bias[idx] <- results_grid$estimated_prop_natural[idx] - results_grid$true_prop_natural[idx]
          
          results_grid$lower_CI_natural[idx] <- est_summary_stats$lower_CI_natural[est_summary_stats$age == a]
          results_grid$upper_CI_natural[idx] <- est_summary_stats$upper_CI_natural[est_summary_stats$age == a]
          #CI width
          results_grid$ci_width[idx] <- est_summary_stats$upper_CI_natural[est_summary_stats$age == a] - 
            est_summary_stats$lower_CI_natural[est_summary_stats$age == a]
          
          #does CI contain true value?
          results_grid$coverage[idx] <- 
            est_summary_stats$lower_CI_natural[est_summary_stats$age == a] <= true_summary_stats$true_proportion_natural[true_summary_stats$age == a] &
            est_summary_stats$upper_CI_natural[est_summary_stats$age == a] >= true_summary_stats$true_proportion_natural[true_summary_stats$age == a]
          
          results_grid$sim_time[idx] <- sim_time
          
          results_grid$age_present_in_sample[idx] <- TRUE
        } else {
          results_grid$estimated_prop_natural[idx] <- NA
          results_grid$age_present_in_sample[idx] <- FALSE
          results_grid$warning_message[idx] <- "No scale samples for this age class"
          
          results_grid$bias[idx] <- -results_grid$true_prop_natural[idx]  # Bias = 0 - true
          results_grid$lower_CI_natural[idx] <- 0
          results_grid$upper_CI_natural[idx] <- 0
          results_grid$ci_width[idx] <- 0
          results_grid$coverage[idx] <- FALSE
          }
        
        
      }
      print(paste("finished scales_n = ",scales_n_seq[s],
                  " boot iteration = ",r,
                  " time = ",sim_time,"seconds"))
    }
  }
  
  #calculate summary stats across replicates
  summary_stats <- results_grid %>%
    group_by(scales_n, age) %>%
    summarize(
      mean_bias = mean(bias, na.rm = TRUE),
      sd_bias = sd(bias, na.rm = TRUE),
      rmse = sqrt(mean(bias^2, na.rm = TRUE)),
      mean_ci_width = mean(ci_width, na.rm = TRUE),
      coverage_rate = mean(coverage, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(list(
    "raw_results" = results_grid,
    "summary" = summary_stats
  ))
}
