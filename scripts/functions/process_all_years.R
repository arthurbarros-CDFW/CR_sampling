#function for stripping yearly sim data and processing all year data
process_all_years <- function(sim_results,true_values="TRUE") {
  map_dfr(names(sim_results), function(year) {
    #get regions
    year<-as.character(year)
    regions <- names(sim_results[[year]])
    
    if(true_values==TRUE){
      #get age_summary_stats
      map_dfr(regions, function(region) {
        if (!is.null(sim_results[[year]][[region]]$true_summary_stats)) {
          sim_results[[year]][[region]]$true_summary_stats %>%
            mutate(
              region = region,
              year = year
            )
        }
      })
    }else{
      #get age_summary_stats
      map_dfr(regions, function(region) {
        if (!is.null(sim_results[[year]][[region]]$age_summary_stats)) {
          sim_results[[year]][[region]]$age_summary_stats %>%
            mutate(
              region = region,
              year = year
            )
        }
      })
    }
  })
}