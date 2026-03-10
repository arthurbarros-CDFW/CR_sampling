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
#' from the simulated population and 2) 'k_expansions' a dataframe of estimates 
#' of number of unrecovered tags for each recovered tags.

sim_spawning_recovery<-function(sim_population,
                                theta=0.2,
                                tag_rate=0.25,
                                iterations=1000){
  
  n_recoveries=sample_n(tbl=sim_population,
                        size=nrow(sim_population)*theta,
                        replace=F)
  
  recovered_tagged_fish <- n_recoveries %>% 
    filter(has_cwt == TRUE)
  recovered_untagged_fish<-n_recoveries%>%
    filter(has_cwt == FALSE)
  
  tag_estimates<-list()
  k_sims<-list()
  #for loop to estimate k for each recovered tag
  for (i in 1:nrow(recovered_tagged_fish)) {
    estimates <- est_unrecovered_tags(n_recovered = 1, theta = theta,iterations)
    tag_estimates[[i]] <- list(
      fish_id = recovered_tagged_fish$fish_id[i],
      origin = recovered_tagged_fish$origin[i],
      age = recovered_tagged_fish$age[i],
      hatchery_source = recovered_tagged_fish$hatchery_source[i],
      mean_k = estimates$mean_k,
      median_k = estimates$median_k,
      ci_lower = estimates$ci_95[1],
      ci_upper = estimates$ci_95[2]
    )
    k_sims[[i]]<-list(
      fish_id=recovered_tagged_fish$fish_id[i],
      k_sim=estimates$k_sim
    )
  }
  
  tag_estimates_df <- bind_rows(tag_estimates)
  tag_estimates_df <- tag_estimates_df %>%
    mutate(
      sum_k = mean_k + 1  # recovered tag + unrecovered tags
    )
  
  #get estimate of all tagged fish and their ages
  N_tagged_est <- sum(tag_estimates_df$sum_k) #est of all tagged fish in population
  
  #group and sum to get age specific estimates of number of tagged fish in population
  N_tagged_age_est<-tag_estimates_df%>%
    group_by(age)%>%
    summarise(tagged_age_total=sum(sum_k))%>%
    mutate(tagged_age_prop=tagged_age_total/sum(tagged_age_total))
  
  sum(N_tagged_age_est$tagged_age_total)
  
  #est of all hatchery origin fish (tagged and untagged) by age
  estimated_hatchery_ages <- N_tagged_age_est%>%
    mutate(hatchery_age_total=tagged_age_total/tag_rate)
  
  #estimate untagged hatchery origin fish
  estimated_hatchery_ages<-estimated_hatchery_ages%>%
    mutate(untagged_hatchery_ages=hatchery_age_total-tagged_age_total)
  
  results<-list("recovered_fish"=data.frame(n_recoveries),
                "k_expansions"=data.frame(tag_estimates_df),
                "hatchery_ages"=data.frame(estimated_hatchery_ages),
                "k_sims"=(k_sims))
  
  return(results)
}
