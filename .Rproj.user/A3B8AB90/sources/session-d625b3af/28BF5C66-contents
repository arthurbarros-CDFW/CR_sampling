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
#' @return a simulated population of adult returns with fish ID, origin,
#' age, hatchery source if applicable, and if tag present

sim_pop<-function(N,
                  h_prop=0.75, 
                  tag_rate=0.25, 
                  ages=c(2,3,4), 
                  age_probs=c(0.15, 0.60, 0.25), 
                  random_age_probs = TRUE,
                  age_var_origin = TRUE,
                  age_dirichlet_concentration = 100){ 
  
  #set up a Dirichlet distribution function to add some random noise
  generate_age_probs <- function(base_probs, concentration) {
    #use Dirichlet distribution: parameters = base_probs * concentration
    dirichlet_params <- base_probs * concentration
    #ensure all parameters > 0
    dirichlet_params <- pmax(dirichlet_params, 0.01)
    #generate random probabilities that sum to 1
    rdirichlet(1, dirichlet_params)[1, ]
  }
  
  #set up age probabilities based on options
  if (random_age_probs && age_var_origin) {
    #generate different random age distributions for hatchery and natural fish
    age_probs_hatchery <- generate_age_probs(age_probs, age_dirichlet_concentration)
    age_probs_natural <- generate_age_probs(age_probs, age_dirichlet_concentration)
    
  } else if (random_age_probs) {
    #use a single age distribution for all fish
    age_probs_random <- generate_age_probs(age_probs, age_dirichlet_concentration)
    age_probs_hatchery <- age_probs_random
    age_probs_natural <- age_probs_random
    
  } else {
    #fixed age probabilities
    age_probs_hatchery <- age_probs
    age_probs_natural <- age_probs
  }
  
  #create fish data frame
  pop_origin <- data.frame(
    fish_id=1:N,
    origin = ifelse(rbinom(N, 1, h_prop) == 1, "hatchery", "natural")
  )

  sim_population = pop_origin%>%
    mutate(age=ifelse(
      origin == "hatchery",
      sample(ages, size = n(), replace = TRUE, prob = age_probs_hatchery),
      sample(ages, size = n(), replace = TRUE, prob = age_probs_natural)
    ),
    has_cwt = ifelse(
      origin == "hatchery",
      rbinom(n = n(), size = 1, prob = tag_rate),
      0
    )
    )

  return(sim_population)
}
