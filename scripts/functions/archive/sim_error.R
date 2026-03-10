#' @title sim_error
#'
#' @description filler
#' 
#' @param age_counts
#' 
#' @param pop
#' 
#' @return fill

sim_error<-function(summary_stats,pop){
  
  true_natural_counts <- pop %>%
    filter(origin == "Natural") %>%
    group_by(age) %>%
    summarise(true_natural_count = n(), .groups = 'drop')%>%
    mutate(true_proportions=true_natural_count/sum(true_natural_count))
  
  accuracy_df <- true_natural_counts %>%
    left_join(age_counts, by = "age")
  
  accuracy_metrics <- accuracy_df %>%
    mutate(
      #absolute error
      abs_error = natural_proportions - true_proportions,
      
      #relative error (percentage)
      rel_error = ifelse(true_count > 0, 
                         abs_error / true_proportions * 100, 
                                 NA),
      
      #absolute percentage error
      abs_perc_error = ifelse(true_count > 0, 
                                      abs(abs_error) / true_proportions * 100, 
                                      NA))
  
  accuracy_metrics<-select(accuracy_metrics,age,true_count,true_proportions,
                           natural_origin_ages,natural_proportions,abs_error)
  
  return(accuracy_metrics)
}
