#' @title sim_wrapper
#'
#' @description fill
#' 
#' @param N
#' 
#' @param theta 
#' 
#' @param iterations 
#' 
#' @param tag_rate
#' 
#' @param survey "hatchery" or "trib"
#' 
#' @param scales_n
#' 
#' @return fill

sim_wrapper<-function(N,
                      theta,
                      iterations,
                      tag_rate,
                      survey,
                      scales_n){
  
  #simulate population
  pop<-sim_pop(N=N,tag_rate=tag_rate)
  
  #get true population totals
  total_tagged_fish<-pop %>% 
    filter(has_cwt == TRUE)
  total_untagged_fish<-pop%>%
    filter(has_cwt==FALSE)
  
  # if the survey is at a spawning site we need to:
  # 1) simulate the number of spawning fish recovered
  # 2) estimate the number of unrecovered tags
  # however if the survey is at a hatchery there is no need,
  # as all fish are "recovered" and we know the population of hatchery-origin fish
  
  if(survey=="trib"){
    #run simulated spawning tributary recovery
    spawning_results<-sim_spawning_recovery(pop,theta,tag_rate,iterations)
    untagged_total<-N-sum(spawning_results$hatchery_ages$tagged_age_total)
    
    #get dataframe of all untagged fish recovered
    untagged_recovered<-spawning_results$recovered_fish%>%
      filter(has_cwt==0)
    
    #get dataframe of all tagged fish ages with k-expansions
    tagged_ages<-spawning_results$hatchery_ages
    
  } else {
    print("ensure survey is set to either 'hatchery' or 'trib'")
  }
  
  #get scales subsample from untagged fish
  scale_samples<-sample_n(untagged_recovered,scales_n,replace = F)
  
  untagged_scale_ages<-scale_samples%>%
    group_by(age)%>%
    summarise(sample_untagged_count=n())%>%
    mutate(untagged_age_prop=sample_untagged_count/sum(sample_untagged_count))%>%
    mutate(untagged_age_total=untagged_age_prop*untagged_total)
  
  age_counts<-untagged_scale_ages%>%
    left_join(tagged_ages,by="age")%>%
    mutate(natural_origin_ages=untagged_age_total-untagged_hatchery_ages,
           natural_proportions=natural_origin_ages/sum(natural_origin_ages))
  
  #run sim_error
  accuracy_metrics<-sim_error(age_counts,pop)
  
  results<-list("accuracy_metrics"=accuracy_metrics,
                "age_counts"=age_counts)
}

