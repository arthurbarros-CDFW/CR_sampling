rm( list = ls()) #clear env
#data simulations
library(tidyverse)
library(data.table)
library(ggplot2)
library(gtools)

#Load functions
sapply(list.files("scripts/functions", pattern = "\\.R$", full.names = TRUE), source)

#set.seed(123)


##########################
#location runs
##########################

scales_n_seq =seq(200,2000,100)#number of scales to age from untagged recoveries
#h_prop_dist = seq(0.5,0.9,0.1)

k_iterations=100 #
scale_iterations=100
n_replicates=100 #times to repeat total simulation
CI=0.95
target_moe=0.10
tag_rate=0.25
trib_data<-readRDS("data/all_returns_clean.Rds")
theta_data<-readRDS("data/location_year_theta.Rds")
iter=1

year_list<-2010:2020

nopooling_iteration_results<-list()

results_list<-list()

for( i in 1:n_replicates){
  iter<-i
  iter_start<-Sys.time()
  
  target_year=sample(year_list,1)
  pop<-sim_location_pop(trib_data,iter,target_year)
  loc_list<-unique(pop$location)
  
  est_summary_stats<-data.frame()
  true_summary_stats<-data.frame()
  for(l in 1:length(loc_list)){
    loc_start<-Sys.time()
    target_loc<-loc_list[l]
    loc_pop<-pop%>%filter(location==target_loc)
    target_theta<-theta_data%>%
      filter(location==target_loc,
             return_year==target_year)
    if(nrow(target_theta)==0){target_theta=0.2} else{
      target_theta<-round(target_theta$theta,2)
    }
    
    spawning_results<-sim_spawning_recovery(loc_pop,
                                            target_theta,tag_rate,
                                            k_iterations)
    
    scale_iters<-list()
    actual_scales_n_run <- c()
    
    for(s in 1: length(scales_n_seq)){
      
      untagged_recovered <- spawning_results$recovered_fish %>%
        filter(has_cwt == 0)
      
      #total number of untagged recoveries
      untagged_total=nrow(untagged_recovered)
      
      # Determine actual sample size to use
      actual_scales_n <- scales_n_seq[s]
      if(untagged_total < scales_n_seq[s]) {
        actual_scales_n <- untagged_total
      }
      
      # Check if we've already run this actual sample size
      if(actual_scales_n %in% actual_scales_n_run) {
        message(paste("Skipping location", target_loc, "- actual_scales_n =", 
                      actual_scales_n, "already run"))
        next
      }
      
      # Record that we're running this sample size
      actual_scales_n_run <- c(actual_scales_n_run, actual_scales_n)
      
      # Print message if we're using a different sample size than intended
      if(untagged_total < scales_n_seq[s]) {
        message(paste("Location", target_loc, "- untagged_total =", untagged_total,
                      "is less than desired scales_n =", scales_n_seq[s],
                      ". Using actual_scales_n =", actual_scales_n))
      }
      
      sim_result<-tryCatch({
        sim_scales_sampling(
          pop=loc_pop,
          spawning_results,
          theta=target_theta,
          k_iterations=k_iterations,
          scale_iterations=scale_iterations,
          tag_rate=tag_rate,
          survey="trib",
          scales_n=scales_n_seq[s],
          CI=CI
        )
      },error=function(e){
        message(paste("Error in simulation:", e$message))
        return(NULL)
      })
      
      #skip if simulation failed
      if (is.null(sim_result)) {
        warning_msg <- "sim_wrapper failed"
        for (a in 2:4) {
          idx <- which(age_grid$scales_n == scales_n_seq[s] &
                         age_grid$replicate == r & 
                         age_grid$age == a)
          if (length(idx) > 0) {
            age_grid$warning_message[idx] <- warning_msg
            age_grid$sim_time[idx] <- sim_time
          }
        }
        print(paste("FAILED: scales_n =", scales_n_seq[s],
                    "boot iteration =", r))
        next
      }
      
      sim_result$age_summary_stats$location=target_loc
      sim_result$age_summary_stats$scales_collected=sim_result$scales_collected
      sim_result$true_summary_stats$location=target_loc
      
      est_summary_stats<-est_summary_stats%>%
        rbind(sim_result$age_summary_stats)
      true_summary_stats<-true_summary_stats%>%
        rbind(sim_result$true_summary_stats)
      
    }
    
    true_summary_stats<-unique(true_summary_stats)
    loc_end<-Sys.time()
    loc_time<-loc_end-loc_start
    print(paste("location:",target_loc," time:",loc_time," iter:",i))
    
  }
  iter_end<-Sys.time()
  iter_time=iter_end-iter_start
  iter_results<-list("est_summary_stats"=est_summary_stats,
                     "true_summary_Stats"=true_summary_stats)
  results_list[[i]]=iter_results
  print(paste("ITERATION",iter,"DONE, time:",round(iter_time,2)))
}
saveRDS(results_list,"outputs/location_based_sims.Rds")

results_list<-readRDS("outputs/location_based_sims.Rds")

#analyse
# Combine all iterations into one dataframe
  all_results <- bind_rows(lapply(1:length(results_list), function(i) {
    results_list[[i]]$est_summary_stats %>%
      mutate(iteration = i)
  }))
  
  all_results <- all_results %>%
    mutate(scales_collected = floor(scales_collected / 100) * 100)
  
  all_results[, moe_count := (upper_CI_count - lower_CI_count)/2]
  all_results[, moe_prop := (upper_CI_prop - lower_CI_prop)/2]
  
  all_results[, rel_moe_count := moe_count / est_count_natural * 100]
  all_results[, rel_moe_prop := moe_prop / mean_proportion_natural * 100]
  
  all_results[, cv_natural := sd_proportion_natural / mean_proportion_natural]
  
  
locations<-c(unique(select(all_results,location)))
for(l in 1:length(locations$location)){
  target_loc<-locations$location[l]
  d<-all_results%>%filter(location==target_loc,
                             !scales_collected<200)
  
  ggplot(d, aes(x = factor(scales_collected), y = (rel_moe_count), color = factor(age))) +
    geom_boxplot()+
    labs(title = target_loc,
         x = "Scales Collected",
         y = "Count % Margin of Error",
         color = "Age") +
    #scale_y_continuous(limits=c(0,200),
     #                  breaks=seq(from=0,to=200,by=20))+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  
  ggsave(paste("outputs/locations/count_moe_",target_loc,".png",sep=""))
  
  ggplot(d, aes(x = factor(scales_collected), y = (moe_prop*100), color = factor(age))) +
    geom_boxplot()+
    labs(title = target_loc,
         x = "Scales Collected",
         y = "Proportion % Margin of Error",
         color = "Age") +
    #scale_y_continuous(limits=c(0,200),
    #                   breaks=seq(from=0,to=200,by=20))+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  
  ggsave(paste("outputs/locations/prop_moe_",target_loc,".png",sep=""))
}

