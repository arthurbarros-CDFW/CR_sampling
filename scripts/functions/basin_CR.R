#' @title basin_CR
#'
#' @description 
#' 
#' @param 
#' 
#' @param 
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


basin_CR<-function(basin_N,
                   spawning_results,
                   iterations=1000,
                   CI=0.95){
  
  #convert scale samples to data.table for speed
  scale_samples_dt <- data.table::as.data.table(spawning_results$basin_scales)
  
  scales_n<-nrow(scale_samples_dt)
  
  #hatchery origin ages for each k_iteration
  hatchery_ages<- spawning_results$hatchery_ages
  
  #split up hatchery ages by k_iteration before loop for speed
  hatchery_ages_split <- split(hatchery_ages, hatchery_ages$k_iteration)
  
  #make blank boot results vector
  time_start<-Sys.time()
  boot_results <- vector("list", iterations)
  
  #bootstrapping scale samples to estimate uncertainty
  for(i in 1:iterations){
    #bootstrap using data.table sampling index methods for speed
    sample_idx <- sample.int(scales_n, scales_n, replace = TRUE)
    scales_boot <- scale_samples_dt[sample_idx]
    
    if(length(hatchery_ages_split)==0){
      age<-c(2,3,4)
      N<-total_tags<-total_hatchery<-c(0,0,0)
      k_iteration=c(i,i,i)
      tagged_ages_boot<-data.table(tag_age,N,total_tags,total_hatchery,k_iteration)
    }else{
      #get tagged ages for this iteration using same k_iteration
      tagged_ages_boot <- as.data.table(hatchery_ages_split[[i]])
    }
    setDT(tagged_ages_boot)  #convert to data.table
    
    #estimate total of tagged fish in basin pop
    total_tagged_boot<-sum(tagged_ages_boot$total_tags,na.rm=T)
    
    #estimate of total untagged fish in basin pop
    total_untagged_boot<-basin_N-total_tagged_boot
    
    #for untagged scale boot sample
    #group by age and get count
    #then calculate age proportion
    #then calculate total estimate of untagged fish in that age for population 
    untagged_scale_ages <- scales_boot[
      , .(sample_untagged_count = .N), by = scale_age
    ][
      , untagged_age_prop := sample_untagged_count / sum(sample_untagged_count)
    ][
      , untagged_age_total := untagged_age_prop * total_untagged_boot
    ]
    
    #age list to ensure no missing ages
    all_ages <- data.table(scale_age = 2:4)
    
    #merge with untagged_scale_ages to ensure all ages are present
    #this was more important when just estimating for smaller populations, but
    #still useful
    untagged_scale_ages_complete <- merge(
      all_ages, 
      untagged_scale_ages, 
      by = "scale_age", 
      all.x = TRUE
    )
    
    #replace NAs with zeros
    untagged_scale_ages_complete[
      is.na(sample_untagged_count), 
      `:=`(
        sample_untagged_count = 0,
        untagged_age_prop = 0,
        untagged_age_total = 0
      )
    ]
    
    untagged_scale_ages_complete<-untagged_scale_ages_complete%>%
      mutate(age=scale_age)
    
    tagged_ages_boot<-tagged_ages_boot%>%
      mutate(age=tag_age)
    
    #merge untagged_scale_ages and tagged_ages_boot
    #fast merge and calculation
    age_counts <- merge(
      untagged_scale_ages_complete, 
      tagged_ages_boot, 
      by = "age", 
      all.x = TRUE
    )
    
    #replace any NAs in age_counts columns with 0
    age_counts[
      is.na(total_tags), 
      `:=`(
        total_tags = 0,
        total_hatchery = 0
      )
    ]
    
    #estimate untagged hatchery origin fish in pop at each age
    age_counts[,untagged_hatchery:=total_hatchery-total_tags]
    
    #calculate estimate of natural origin fish in pop at each age
    age_counts[, total_natural := untagged_age_total - untagged_hatchery]
    
    boot_results[[i]] <- as.data.frame(age_counts)
  }
  time_end<-Sys.time()
  boot_time<-time_end-time_start
  
  #estimate uncertainty in estimates
  
  #combine iterations into one data table
  all_boot_dt <- rbindlist(boot_results, idcol = "iteration")
  
  #calculate summary statistics by age class
  age_summary_stats <- all_boot_dt%>%
    group_by(age)%>%
    summarize(
      mean_hatchery=mean(total_hatchery),
      lower_CI_hatchery = quantile(total_hatchery, probs = (1-CI)/2, na.rm = TRUE),
      upper_CI_hatchery = quantile(total_hatchery, probs = 1-(1-CI)/2, na.rm = TRUE),
      mean_natural=mean(total_natural),
      lower_CI_natural = quantile(total_natural, probs = (1-CI)/2, na.rm = TRUE),
      upper_CI_natural = quantile(total_natural, probs = 1-(1-CI)/2, na.rm = TRUE),
      n_iterations = iterations
    )
  
  age_summary_stats_long <- age_summary_stats %>%
    pivot_longer(
      cols = c(mean_hatchery, lower_CI_hatchery, upper_CI_hatchery,
               mean_natural, lower_CI_natural, upper_CI_natural),
      names_to = c(".value", "origin"),
      names_pattern = "(mean|lower_CI|upper_CI)_(hatchery|natural)"
    )
  
  age_summary_stats_long<-age_summary_stats_long%>%
    rename(estimate=mean)
  
  setorder(age_summary_stats_long, age)
  
  return(age_summary_stats_long)
}
