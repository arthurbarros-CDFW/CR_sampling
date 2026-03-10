#' @title sim_location_pop
#'
#' @description Simulate returning adult SRFC population for a given site
#' 
sim_location_pop<-function(data,iter,year){ 
  
  d<-data%>%filter(iteration==iter)
  
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
    mutate(has_cwt=ifelse(
      origin == "hatchery",
      rbinom(n = n(), size = 1, prob = tag_rate),
      0
    ))
  
  iter_pop<-select(iter_pop,-tag_rate)
  
  iter_pop<-iter_pop%>%
    mutate(fish_id=1:nrow(iter_pop))%>%
    filter(return_year==year)
  
  return(iter_pop)
}
