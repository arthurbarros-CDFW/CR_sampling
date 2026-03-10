#' @title sim_wrapper
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

sim_wrapper<-function(pop,
                      theta,
                      k_iterations,
                      scale_iterations,
                      tag_rate,
                      survey,
                      scales_n,
                      CI){
  
  #currently we are only focusing on tribs as they have the most complexity/uncertainty
  #run simulated spawning tributary recovery
  spawning_results<-sim_spawning_recovery(pop,theta,tag_rate,k_iterations)
    
  #get dataframe of all untagged fish recovered
  #these are the fish we can get scales from
  untagged_recovered <- spawning_results$recovered_fish %>%
      filter(has_cwt == 0)
  
  #total number of untagged recoveries
  untagged_total=nrow(untagged_recovered)
  
  #get scales subsample from untagged fish based on scales_n input
  scale_samples<-sample_n(untagged_recovered,scales_n,replace = F)
  
  #convert scale samples to data.table for speed
  scale_samples_dt <- data.table::as.data.table(scale_samples)

  #hatchery origin ages for each k_iteration
  hatchery_ages_list <- spawning_results$hatchery_ages
  
  #split up hatchery ages by k_iteration before loop for speed
  hatchery_ages_split <- split(hatchery_ages_list, hatchery_ages_list$k_iteration)

  #make blank boot results vector
  boot_results <- vector("list", scale_iterations)
  
  time_start<-Sys.time()
  for(j in 1:scale_iterations){
    #bootstrap using data.table sampling index methods for speed
    sample_idx <- sample.int(scales_n, scales_n, replace = TRUE)
    scales_boot <- scale_samples_dt[sample_idx]
    
    #get tagged ages for this iteration using same k_iteration
    tagged_ages_boot <- as.data.table(hatchery_ages_split[[j]])
    setDT(tagged_ages_boot)  #convert to data.table
    
    #get total est of tagged fish in pop
    total_tagged_boot <- sum(tagged_ages_boot$total_tags, na.rm = TRUE)
    #get total est of untagged fish in pop
    untagged_total_boot <- nrow(pop) - total_tagged_boot
    
    #for untagged scale boot sample
    #group by age and get count
    #then calculate age proportion
    #then calculate total estimate of untagged fish in that age for population 
    untagged_scale_ages <- scales_boot[
      , .(sample_untagged_count = .N), by = age
    ][
      , untagged_age_prop := sample_untagged_count / sum(sample_untagged_count)
    ][
      , untagged_age_total := untagged_age_prop * untagged_total_boot
    ]
    
    #merge untagged_scale_ages and tagged_ages_boot
    #fast merge and calculation
    age_counts <- merge(
      untagged_scale_ages, 
      tagged_ages_boot, 
      by = "age", 
      all.x = TRUE
    )
    
    #calculate estimate of untagged hatchery origin fish in pop at each age
    age_counts[, untagged_hatchery := total_hatchery - total_tags]
    
    #calculate estimate of natural origin fish in pop at each age
    age_counts[, total_natural := untagged_age_total - untagged_hatchery]
    
    age_counts<-age_counts%>%
      mutate(total_natural=ifelse(total_natural<0,0,total_natural))
    
    #calculate estimate of proportion of origin fish in pop at each age
    age_counts[, natural_proportions := total_natural/sum(total_natural)]
    
    #save age results in boot_results list
    boot_results[[j]] <- as.data.frame(age_counts)
  }
  
  time_end<-Sys.time()
  
  boot_time<-time_end-time_start #~10-15 seconds, variable for 1000 iterations
  
  #estimate uncertainty in natural age counts and proportions
  #combine iterations into one data table
  all_boot_dt <- rbindlist(boot_results, idcol = "iteration")
  
  #calculate summary statistics by age class
  est_summary_stats <- all_boot_dt[
    !is.na(natural_proportions),
    .(
      est_count_hatchery=mean(total_hatchery),
      est_count_natural=mean(total_natural),
      mean_proportion_natural = mean(natural_proportions, na.rm = TRUE),
      median_proportion_natural = median(natural_proportions, na.rm = TRUE),
      sd_proportion_natural = sd(natural_proportions, na.rm = TRUE),
      se_proportion_natural = sd(natural_proportions, na.rm = TRUE) / sqrt(.N),
      lower_CI_natural = quantile(natural_proportions, probs = (1-CI)/2, na.rm = TRUE),
      upper_CI_natural = quantile(natural_proportions, probs = 1-(1-CI)/2, na.rm = TRUE),
      n_iterations = .N
    ),
    by = age
  ]
  
  setorder(est_summary_stats, age)
  
  #set true population statistics
  true_summary_stats<-pop%>%
    group_by(age,origin)%>%
    summarize(true_count=n())
  
  true_summary_stats <- true_summary_stats %>%
    pivot_wider(
      id_cols = age,
      names_from = origin,
      values_from = c(true_count),
      names_prefix = "true_count_"
    )
  
  true_summary_stats<-true_summary_stats%>%
    mutate(
      total_hatchery = sum(true_count_hatchery),
      total_natural = sum(true_count_natural)
    )%>%
    mutate(
      prop_hatchery = true_count_hatchery / sum(total_hatchery),
      prop_natural = true_count_natural / sum(total_natural)
    )
  
  
  true_summary_stats <- pop %>%
    group_by(age, origin) %>%
    summarize(true_count = n()) %>%
    pivot_wider(
      id_cols = age,
      names_from = origin,
      values_from = c(true_count),
      names_prefix = "true_count_"
    )
  
  
  total_count_hatchery<-sum(true_summary_stats$true_count_hatchery)
  total_count_natural<-sum(true_summary_stats$true_count_natural)
  
  true_summary_stats <- true_summary_stats %>%
    mutate(
      true_proportion_hatchery = true_count_hatchery/total_count_hatchery,
      true_proportion_natural = true_count_natural/total_count_natural
    )
  
  
  results<-list("est_summary_stats"=est_summary_stats,
                "true_summary_stats"=true_summary_stats)
  return(results)
}

