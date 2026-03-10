#' @title sim_pop
#'
#' @description Simulate returning adult SRFC population for a given site
#' 
#' @param N true population size. In a hatchery this should be known, at a spawning
#' this is estimated.
#' 
#' @param h_prop true proportion of hatchery-origin fish
#' 
#' @param tag_rate cwt tagging rate of hatchery-origin fish
#' 
#' @param ages range of potential returning adult ages
#'
#' @param age_probs probability seed for ages. Note: should I add variation 
#' between natural and hatchery-origin fish?
#'
#' @param hatcheries list of potential hatchery sources
#' 
#' @param hatchery_props probability/proportion of hatchery sources
#'
#' @return a simulated population of adult returns with fish ID, origin,
#' age, hatchery source if applicable, and if tag present

sim_pop<-function(N,
                  h_prop=0.75, 
                  tag_rate, 
                  ages=c(2,3,4,5), 
                  age_probs=c(0.15, 0.60, 0.20, 0.05), 
                  hatcheries=c("Coleman", "Feather River", "Nimbus"), 
                  hatchery_props=c(0.56, 0.42, 0.02),
                  #distribution flags
                  N_dist = c("fixed", "poisson", "uniform"),
                  tag_rate_dist = c("fixed", "beta"),
                  age_probs_dist = c("fixed", "dirichlet"),
                  hatchery_props_dist = c("fixed", "dirichlet"),
                  #distribution parameters
                  N_min = NULL, N_max = NULL,  # for uniform
                  tag_rate_shape1 = 1, tag_rate_shape2 = 1,  # for beta
                  age_dirichlet_alpha = rep(5, 4),  # concentration for ages
                  hatchery_dirichlet_alpha = c(56, 42, 2)){ 
  
  if (N_dist[1] == "poisson") {
    N <- rpois(1, N)
  } else if (N_dist[1] == "uniform") {
    N <- round(runif(1, N_min, N_max))
  }
  #else N is fixed and uses provided value
  
  #handle age_probs distribution
  if (age_probs_dist[1] == "dirichlet") {
    age_probs <- gtools::rdirichlet(1, age_dirichlet_alpha)[1, ]
  }
  
  #handle hatchery_prop distribution
  if (hatchery_props_dist[1] == "dirichlet") {
    hatchery_props <- gtools::rdirichlet(1, hatchery_dirichlet_alpha)[1, ]
  }
  
  #create simulated population
  sim_population <- tibble(
    fish_id = 1:N, #assign each fish an "ID"
    
    #fish origin
    origin = sample(c("hatchery", "natural"), 
                    size = N, 
                    replace = TRUE, 
                    prob = c(h_prop, 1 - h_prop)),
    
    #fish age
    age = sample(ages, 
                 size = N, 
                 replace = TRUE, 
                 prob = age_probs),
    
    #fish hatchery source (where applicable)
    hatchery_source = ifelse(origin == "hatchery",
                             sample(hatcheries, 
                                    size = sum(origin == "hatchery"), 
                                    replace = TRUE, 
                                    prob = hatchery_props),
                             NA)
  )
  sim_population <- sim_population %>%
    mutate(
      #assign cwt to hatchery origin fish
      has_cwt = ifelse(origin == "hatchery",
                       rbinom(n = n(), size = 1, prob = tag_rate),
                       FALSE))
  
  return(sim_population)
}
