#' @title sim_k_draws
#'
#' @description 
#' 
#' @param sim_population 
#' 
#' @param theta 
#' 
#' @param iterations 
#' 
#' @param tag_rate 
#' 
#' @return 

tag_k_draws <- function(tag_data,
                        theta=0.2,
                        iterations=1000,
                        freq,
                        f_prod){
  
  #generate all random draws at once
  #generate matrix of k_iteration columns, n_ages rows
  #fill with k draws
  unrecovered_tags <- matrix(
    stats::rnbinom(n = iterations, size = freq, prob = theta),
    nrow = 1,
    ncol = iterations
  )
  
  total_tags<-unrecovered_tags+freq
  
  total_hatchery <- total_tags / f_prod
  
  result_dt <- data.table::data.table(
    release_group = rep(tag_data$release_group, iterations),
    brood_year=rep(tag_data$brood_year),
    tag_age=rep(tag_data$tag_age),
    return_year=rep(tag_data$return_year),
    location=rep(tag_data$location),
    theta=rep(theta),
    tag_rate=rep(tag_data$tag_rate),
    tag_freq = rep(tag_data$frequency, iterations),
    total_tags = as.vector(total_tags),
    total_hatchery = as.vector(total_hatchery),
    k_iteration = rep(1:iterations)
  )
  return(result_dt)
}
