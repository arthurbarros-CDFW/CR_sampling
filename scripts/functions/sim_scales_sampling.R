#' @title sim_scales_sampling
#'
#' @description simulate sampling of scales from recovered fish output of 
#' sim_spawning_recovery.R and produce estimates of natural-origin abundance
#' at each age-class, as well as uncertainty metrics.
#' 
#' @param pop the simulated population data frame produced by sim_pop().
#' 
#' @param spawning_results the output of the sim_spawning_recovery() function. 
#' Must be the list with two data frames spawning_recoveries$recovered_fish 
#' and spawning_recoveries$hatchery_ages.
#' 
#' @param iterations number of draws from from a negative-binomial distribution 
#' to estimate k. Default 1000. In this function also the number of iterations 
#' to bootstrap resampling of the scale subsample to characterize scale 
#' sampling uncertainty.
#' 
#' @param scales_n number of scale samples to collect and "age" from recovered 
#' fish.
#' 
#' @param CI confidence intervals used to calculate upper and lower confidence 
#' intervals of natural-origin age class estimates produced from 
#' iterations bootstrapping. Default 0.95.
#' 
#' @return 

sim_scales_sampling<-function(pop,
                              spawning_results,
                      iterations=1000,
                      scales_n,
                      CI=0.95){
    
  #get dataframe of all untagged fish recovered
  #these are the fish we can get scales from
  untagged_recovered <- spawning_results$recovered_fish %>%
      filter(has_cwt == 0)
  
  #total number of untagged recoveries
  untagged_total=nrow(untagged_recovered)
  
  time_start<-Sys.time()
  
  #make blank boot results vector
  boot_results <- vector("list", iterations)
  
  if(untagged_total>scales_n){
    #get scales subsample from untagged fish based on scales_n input
    scale_samples<-sample_n(untagged_recovered,scales_n,replace = F)
  } else {
    scale_samples<-untagged_recovered
    scales_n=nrow(scale_samples)
  }
    
  #convert scale samples to data.table for speed
  scale_samples_dt <- data.table::as.data.table(scale_samples)
  
  #hatchery origin ages for each k_iteration
  hatchery_ages_list <- spawning_results$hatchery_ages
    
  #split up hatchery ages by k_iteration before loop for speed
  hatchery_ages_split <- split(hatchery_ages_list, hatchery_ages_list$k_iteration)
    
  #make blank boot results vector
  boot_results <- vector("list", iterations)
    
  for(j in 1:iterations){
    #bootstrap using data.table sampling index methods for speed
    sample_idx <- sample.int(scales_n, scales_n, replace = TRUE)
    scales_boot <- scale_samples_dt[sample_idx]
    
    if(length(hatchery_ages_split)==0){
      age<-c(2,3,4)
      N<-total_tags<-total_hatchery<-c(0,0,0)
      iterations=c(j,j,j)
      
      tagged_ages_boot<-data.table(age,N,total_tags,total_hatchery,iterations)
    }else{
      #get tagged ages for this iteration using same k_iteration
      tagged_ages_boot <- as.data.table(hatchery_ages_split[[j]])
    }
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
    
    #age list to compare
    all_ages <- data.table(age = 2:4)
    
    #merge with untagged_scale_ages to ensure all ages are present
    untagged_scale_ages_complete <- merge(
      all_ages, 
      untagged_scale_ages, 
      by = "age", 
      all.x = TRUE
    )
    
    #replace NAs with zeros for counts and appropriate values for derived columns
    untagged_scale_ages_complete[
      is.na(sample_untagged_count), 
      `:=`(
        sample_untagged_count = 0,
        untagged_age_prop = 0,
        untagged_age_total = 0
      )
    ]
      
    #merge untagged_scale_ages and tagged_ages_boot
    #fast merge and calculation
    age_counts <- merge(
      untagged_scale_ages_complete, 
      tagged_ages_boot, 
      by = "age", 
      all.x = TRUE
    )
    
    #replace any NAs in tagged_ages_boot columns with 0
    age_counts[
      is.na(total_tags), 
      `:=`(
        total_tags = 0,
        total_hatchery = 0
      )
    ]
    
    
    #calculate estimate of untagged hatchery origin fish in pop at each age
    age_counts[, untagged_hatchery := total_hatchery - total_tags]
    
    #calculate estimate of natural origin fish in pop at each age
    age_counts[, total_natural := untagged_age_total - untagged_hatchery]
    
    age_counts<-age_counts%>%
      mutate(total_natural=ifelse(total_natural<0,0,total_natural))
    
    #calculate estimate of proportion of origin fish in pop at each age
    total_natural_sum <- sum(age_counts$total_natural)
    suppressWarnings(if(total_natural_sum > 0) {
      age_counts[, natural_proportions := total_natural / total_natural_sum]
    } else {
      age_counts[, natural_proportions := 0]
    })
    
    #save age results in boot_results list
    boot_results[[j]] <- as.data.frame(age_counts)
  }
  
  time_end<-Sys.time()
  
  boot_time<-time_end-time_start #~10-15 seconds, variable for 1000 iterations
  
  #estimate uncertainty in natural age counts and proportions
  #combine iterations into one data table
  all_boot_dt <- rbindlist(boot_results, idcol = "iteration")
  
  #calculate summary statistics by age class
  age_summary_stats <- all_boot_dt[
    !is.na(natural_proportions),
    .(
      mean_count_hatchery=mean(total_hatchery),
      mean_count_natural=mean(total_natural),
      mean_proportion_natural = mean(natural_proportions, na.rm = TRUE),
      sd_proportion_natural = sd(natural_proportions, na.rm = TRUE),
      se_proportion_natural = sd(natural_proportions, na.rm = TRUE) / sqrt(.N),
      lower_CI_prop = quantile(natural_proportions, probs = (1-CI)/2, na.rm = TRUE),
      upper_CI_prop = quantile(natural_proportions, probs = 1-(1-CI)/2, na.rm = TRUE),
      lower_CI_count = quantile(total_natural, probs = (1-CI)/2, na.rm = TRUE),
      upper_CI_count = quantile(total_natural, probs = 1-(1-CI)/2, na.rm = TRUE),
      n_iterations = .N
    ),
    by = age
  ]
  
  setorder(age_summary_stats, age)
  
  #calculate summary statistics by age class
  iteration_overall <- all_boot_dt[
    !is.na(natural_proportions),
    .(
      total_hatchery_iter = sum(total_hatchery),
      total_natural_iter = sum(total_natural)
    ),
    by = iteration
  ][
    , overall_prop_natural := total_natural_iter / (total_hatchery_iter + total_natural_iter)
  ]
  
  #set true population statistics
  true_summary_stats <- pop %>%
    group_by(age, origin) %>%
    summarize(true_count = n(), .groups = "drop") %>%
    complete(age, origin = c("hatchery", "natural"), fill = list(true_count = 0)) %>%
    pivot_wider(
      id_cols = age,
      names_from = origin,
      values_from = true_count,
      names_prefix = "true_count_"
    )
  
  total_count_hatchery=sum(true_summary_stats$true_count_hatchery)
  total_count_natural=sum(true_summary_stats$true_count_natural)
  
  true_summary_stats <- true_summary_stats %>%
    mutate(
      true_proportion_hatchery = if(total_count_hatchery > 0){
        true_count_hatchery / sum(true_count_hatchery)
      }else{0},
      true_proportion_natural = if(total_count_natural > 0){
        true_count_natural / sum(true_count_natural)
      }else{0}
    )
  
  
  
  results<-list("age_summary_stats"=age_summary_stats,
                "true_summary_stats"=true_summary_stats,
                "scales_collected"=scales_n)
  return(results)
}

