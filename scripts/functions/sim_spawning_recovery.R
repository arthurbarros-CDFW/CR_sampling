#' @title sim_spawning_recovery
#'
#' @description simulate recoveries of a simulated adult SRFC spawning-grounds 
#' population produced in the sim_pop(). 
#' Recoveries at a hatchery are 'complete' so we know how many tags are there, 
#' but for fish spawning in the river there is some amount of tagged fish
#' that go unrecovered. 
#' This script also returns estimates of 'k',the estimate of unrecovered tags 
#' available for each tag recovered by running est_unrecovered_tags() 
#' 
#' @param sim_population data frame of a simulated population of returning SRFC
#' 
#' @param theta probability of recovering a given fish, or the sampling fraction.
#' 
#' @param iterations number of draws from negative-binomial distribution to 
#' estimate k, where k~NB(1,theta).
#' 
#' @param tag_rate 
#' 
#' @return a list with 1) 'recovered_fish' a dataframe of recovered/sampled fish 
#' from the simulated population and 2) 'hatchery_ages' a dataframe of estimates 
#' of number of unrecovered tags for each recovered tags.

sim_spawning_recovery <- function(sim_population,
                                  theta=0.2,
                                  tag_rate=0.25,
                                  iterations=1000) {
  
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Please install data.table for this version: install.packages('data.table')")
  }
  
  #turn sim_population into a data.table for speed
  dt_pop <- data.table::as.data.table(sim_population)
  
  #recovered fish from population
  n_total <- nrow(dt_pop)
  n_recover <- round(n_total * theta)
  recover_idx <- sample.int(n_total, n_recover, replace = FALSE)
  n_recoveries <- dt_pop[recover_idx]
  
  #get tagged fish from recoveries
  recovered_tagged <- n_recoveries[has_cwt == TRUE]
  
  #if statement to throw empty value back if no recovered tags
  if (nrow(recovered_tagged) > 0) {
    #count by age
    counts <- recovered_tagged[, .(N = .N), by = age]
    
    #vectorized simulation
    N_vals <- counts$N #counts for each age group
    n_ages <- length(N_vals) #number of unique ages
    
    #generate all random draws at once
    #generate matrix of k_iteration columns, n_ages rows
    #fill with k draws
    unrecovered_tags <- matrix(
      stats::rnbinom(n = n_ages * iterations, size = N_vals, prob = theta),
      nrow = n_ages,
      ncol = iterations
    )
    
    #add original recovered tags to k draws to get total number of tags in pop
    total_tags <- unrecovered_tags + N_vals 
    
    #hatchery tag rate expansion for estimate of number of hatchery origin fish in pop
    total_ho <- total_tags / tag_rate
    
    #create result using data.table melting
    result_dt <- data.table::data.table(
      age = rep(counts$age, iterations),
      N = rep(N_vals, iterations),
      total_tags = as.vector(total_tags),
      total_hatchery = as.vector(total_ho),
      k_iteration = rep(1:iterations, each = n_ages)
    )
    
    hatchery_ages <- as.data.frame(result_dt)
    
  } else {
    hatchery_ages <- data.frame(
      age = numeric(0),
      N = numeric(0),
      total_tags = numeric(0),
      total_hatchery = numeric(0),
      k_iteration = numeric(0)
    )
  }
  
  survey_results<-list(
    recovered_fish = as.data.frame(n_recoveries),
    hatchery_ages = hatchery_ages
  )
  
  return(survey_results)
  
}
