#' @title est_unrecovered_tags
#'
#' @description for each tag recovered from a spawning ground estimates 
#' the number of tags present but unrecovered (k) using a negative binomial 
#' draw k~NB(1,theta) from E. Chen 2023. This is only run for escapement to 
#' tributary spawning grounds.
#' Runs 1000 iterations to characterize uncertainty from sampling.
#' 
#' @param n_recovered number of CWT tags recovered in spawning ground survey.
#' 
#' @param theta probability of recovering a given fish, or the sampling fraction.
#' 
#' @param iterations number of draws from negative-binomial distribution to 
#' estimate k, where k~NB(1,theta).
#' 
#' @return a list with mean_k, simulated k draws dataframe, median_k, and 95% CI

est_unrecovered_tags<-function(n_recovered, 
                               theta,
                               iterations=1000){
  k_sim <- rnbinom(n = iterations, size = 1, prob = theta)
  mean_k <- mean(k_sim)
  return(list(
    mean_k = mean_k,
    k_sim = k_sim,
    median_k = median(k_sim),
    ci_95 = quantile(k_sim, probs = c(0.025, 0.975))
  ))
}
