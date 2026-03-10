#' @title summary_stats
#'
#' @description fill
#' 
#' @param N
#' 
#' @param theta 
#' 
#' @param k_iterations 
#' 
#' @param tag_rate
#' 
#' @param survey "hatchery" or "trib"
#' 
#' @param scales_n
#' 
#' @return fill

sim_summary_stats<-function(data,
                            target_moe){
  
  d_age<-data$age_results%>%
    mutate(meets_target_moe_prop=ifelse(moe_prop<target_moe,T,F))
  d_age<-d_age%>%
    mutate(meets_target_moe_count=ifelse((moe_count/est_count_natural)<target_moe,T,F))
  
  d_totals<-data$total_results%>%
    mutate(meets_target_moe_prop=ifelse(moe_prop<target_moe,T,F))
  
  #calculate summary stats across replicates
  age_summary_stats <- d_age %>%
    group_by(scales_n, age) %>%
    summarize(
      mean_bias = mean(bias, na.rm = TRUE),
      sd_bias = sd(bias, na.rm = TRUE),
      rmse = sqrt(mean(bias^2, na.rm = TRUE)),
      mean_ci_width_prop = mean(ci_width_prop, na.rm = TRUE),
      coverage_rate = mean(coverage, na.rm = TRUE),
      mean_moe_prop=mean(moe_prop, na.rm = TRUE),
      mean_moe_count=mean(moe_count,na.rm=TRUE),
      sd_moe_prop=sd(moe_prop, na.rm = TRUE),
      sd_moe_count=sd(moe_count,na.rm=TRUE),
      pct_meeting_target_moe_prop = mean(meets_target_moe_prop) * 100,
      pct_meeting_target_moe_count = mean(meets_target_moe_count) * 100,
      .groups = "drop"
    )
  
  totals_summary_stats <- d_totals %>%
    group_by(scales_n) %>%
    summarize(
      mean_bias = mean(bias, na.rm = TRUE),
      sd_bias = sd(bias, na.rm = TRUE),
      rmse = sqrt(mean(bias^2, na.rm = TRUE)),
      mean_ci_width_prop = mean(ci_width_prop, na.rm = TRUE),
      coverage_rate = mean(coverage, na.rm = TRUE),
      mean_moe_prop=mean(moe_prop, na.rm = TRUE),
      mean_moe_count=mean(moe_count,na.rm=TRUE),
      sd_moe_prop=sd(moe_prop, na.rm = TRUE),
      sd_moe_count=sd(moe_count,na.rm=TRUE),
      pct_meeting_target_moe_prop = mean(meets_target_moe_prop) * 100,
      .groups = "drop"
    )
  
  return(list(
    "age_summary_stats"=age_summary_stats,
    "totals_summary_stats"=totals_summary_stats
  ))
}
