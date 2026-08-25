#' @title sim_spawning_recovery
#'
#' @description simulate recoveries of a simulated adult SRFC spawning-grounds 
#' population produced in the sim_pop(). 
#' Recoveries at a hatchery are 'complete' so we know how many tags are there, 
#' but for fish spawning in the river there is some amount of tagged fish
#' that go unrecovered. 
#' 
#' @param sim_population data frame of a simulated population of returning SRFC
#' 
#' @param theta probability of recovering a given fish, or the sampling fraction.
#' 
#' @param iterations number of draws from negative-binomial distribution to 
#' estimate k, where k~NB(1,theta).
#' 
#' @return a 'recovered_fish' a dataframe of recovered/sampled fish 
#' from the simulated population

sim_spawning_recovery <- function(sim_population,
                                  theta=0.2) {
  
  #require data table package
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
  
  return(n_recoveries)
  
}
