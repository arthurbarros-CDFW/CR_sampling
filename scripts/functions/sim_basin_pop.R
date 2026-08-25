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
  d<-d%>%mutate(brood_year=return_year-age)
  
  d<-d%>%rename(tag_age=age)
  
  tag_list<-unique(select(d,brood_year,tag_rate,location,tag_age,return_year))
  
  #get location props for year
  loc_props<-d%>%
    group_by(location,brood_year)%>%
    summarise(count=sum(count))
  loc_props <- loc_props %>%
    group_by(brood_year) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup()
  loc_props<-loc_props%>%
    select(-count)
  
  expand_fish_efficient <- function(dt) {
    dt_nonzero <- dt[count > 0]
    
    fish_list <- lapply(1:nrow(dt_nonzero), function(i) {
      row <- dt_nonzero[i]
      data.table(
        return_year = rep(row$return_year, row$count),
        location = rep(row$location, row$count),
        tag_age = rep(row$tag_age, row$count),
        tag_rate = rep(row$tag_rate, row$count),
        origin = rep(row$origin, row$count)
      )
    })
    rbindlist(fish_list)
  }
  
  iter_pop<-expand_fish_efficient(d)

  #loc tag rates?
  loc_tag_rates<-unique(select(iter_pop,tag_age,tag_rate,origin,location,return_year))
  loc_tag_rates<-loc_tag_rates%>%
    mutate(brood_year=return_year-tag_age,
           brood_location=location)%>%
    filter(origin=="hatchery")%>%
    group_by(brood_year,brood_location)%>%
    summarise(tag_rate=mean(tag_rate))%>%
    select(brood_year,brood_location,tag_rate)
  
  #add brood location and year to iter_pop
  iter_pop<-iter_pop%>%
    mutate(brood_year=return_year-tag_age)

  assign_brood_location <- function(row, props_df) {
    if (runif(1) < 0.9) {
      return(row$location)
    } else {
      #get location proportions
      year_props <- props_df[props_df$brood_year == row$brood_year, ]
      return(sample(year_props$location, size = 1, prob = year_props$proportion))
    }
  }
  
  iter_pop$brood_location <- apply(iter_pop, 1, function(row) {
    assign_brood_location(as.list(row), loc_props)
  })
  
  iter_pop<-iter_pop%>%
    select(-tag_rate)%>%
    left_join(loc_tag_rates)
  iter_pop$tag_rate[is.na(iter_pop$tag_rate)] <- .25
  
  iter_pop<-iter_pop%>%
    mutate(has_cwt=ifelse(
      origin == "hatchery",
      rbinom(n = n(), size = 1, prob = tag_rate),
      0
    ))
  
  iter_pop<-iter_pop%>%
    mutate(fish_id=1:nrow(iter_pop))%>%
    filter(return_year==year)
  
  #add CWT UID
  iter_pop<-iter_pop%>%
    mutate(release_group=ifelse(has_cwt==1,paste(brood_location,brood_year),NA),
           tag_rate=ifelse(origin=="hatchery",tag_rate,NA))%>%
    mutate(tag_age=ifelse(has_cwt==0,NA,tag_age))
  
  results<-list("iter_pop"=iter_pop)
  return(results)
}
