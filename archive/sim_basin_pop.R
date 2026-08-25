#' @title sim_location_pop
#'
#' @description Simulate returning adult SRFC population for a given site
#' 
#' @param data a cleaned and formatted data set of survey based estimates 
#' containing distributions of age-class abundances for each spawning 
#' ground/hatchery and each year, as well as tag rates. Data formatted in the 
#' clean_recovery_estimates.R script.
#' 
#' @param iter index of which estimate iteration to pull from, this is used 
#' when bootstrapping across available data for uncertainty.
#' 
#' @param year the target year for which population data will be pulled.
#' 
#' @return 
sim_basin_pop<-function(data,iter,year){ 
  
  d<-data%>%filter(iteration==iter,return_year==year)
  
  #get location props for year
  loc_props<-d%>%
    group_by(location)%>%
    summarise(total_count=sum(count))%>%
    mutate(loc_proportion=total_count/sum(total_count))
  
  expand_fish_efficient <- function(dt) {
    dt_nonzero <- dt[count > 0]
    
    fish_list <- lapply(1:nrow(dt_nonzero), function(i) {
      row <- dt_nonzero[i]
      data.table(
        return_year = rep(row$return_year, row$count),
        location = rep(row$location, row$count),
        age = rep(row$age, row$count),
        tag_rate = rep(row$tag_rate, row$count),
        origin = rep(row$origin, row$count)
      )
    })
    rbindlist(fish_list)
  }
  
  iter_pop<-expand_fish_efficient(d)
  iter_pop<-iter_pop%>%
    mutate(tag_rate=round(tag_rate,2))
  
  loc_list<-
  
  #we need to standardize the tag_rate for year x location combos,
  iter_tag_rates<-unique(select(iter_pop,age,tag_rate,origin,location))
  iter_tag_rates<-iter_tag_rates%>%
    group_by(location,origin)%>%
    summarize(tag_rate=mean(tag_rate))
  
  iter_tag_rates$year=year
  
  iter_pop<-select(iter_pop,-tag_rate)%>%
    left_join(iter_tag_rates)
  
  iter_pop<-iter_pop%>%
    mutate(has_cwt=ifelse(
      origin == "hatchery",
      rbinom(n = n(), size = 1, prob = tag_rate),
      0
    ))
  
  iter_pop<-select(iter_pop,-tag_rate)
  
  iter_pop<-iter_pop%>%
    mutate(fish_id=1:nrow(iter_pop))%>%
    filter(return_year==year)
  
  results<-list("iter_pop"=iter_pop,
                "tag_rates"=iter_tag_rates)
  return(results)
}
