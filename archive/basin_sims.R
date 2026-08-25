rm( list = ls()) #clear env
#data simulations
library(tidyverse)
library(data.table)
library(ggplot2)
library(gtools)

#Load functions
sapply(list.files("scripts/functions/simulation archive", pattern = "\\.R$", full.names = TRUE), source)

set.seed(67) #for reproducibility

##########################
#Without pooling
##########################
#The below should give us a list of all the estimates for each loc x year x
#spawning size
#we can use those age count estimates and join them to get basin wide
#this is our "do n_scales at each location approach"

iterations=1000 #
scales_n_seq =seq(500,500,100)
boot_replicates=100 #times to repeat total simulation
CI=0.90
h_prop=.75

#read in data
trib_data<-readRDS("data/all_returns_clean.Rds")
theta_data<-readRDS("data/location_year_theta.Rds")

year_list<-2010:2020

iter_results<-list()

for(i in 1:boot_replicates){
  iter_start<-Sys.time()
  
  #iteration population
  iter_year<-sample(year_list,1)
  pop<-sim_basin_pop(trib_data,i,iter_year)
  tag_rates<-pop$tag_rates
  
  pop<-pop$iter_pop
  loc_list<-unique(pop$location)#pull list of sampled locations
  
  true_pop_ages<-pop%>%
    filter(origin=="natural")%>%
    group_by(age)%>%
    summarise(true_abundance=n())
  
  #spawning recoveries for each location
  recoveries_list<-list()
  for(l in 1:length(loc_list)){ #loop recoveries sim for each location
    iter_loc<-loc_list[l]
    loc_pop<-pop%>%filter(location==iter_loc) #get pop data for locaiton
    target_theta<-theta_data%>% #get theta for location x year
      filter(location==iter_loc,
             return_year==iter_year)
    if(nrow(target_theta)==0){target_theta=0.2} else{ #if no theta, set default
      target_theta<-round(target_theta$theta,2)
    }
    
    #get spawning recoveries
    spawning_recoveries<-sim_spawning_recovery(loc_pop,
                                               theta=target_theta,
                                               iterations=iterations)
    recoveries_list[[l]]<-spawning_recoveries
  }
  names(recoveries_list)<-loc_list
  
  scales_n_results <- list()
  
  for(scales_idx in 1:length(scales_n_seq)) {
    current_scales_n <- scales_n_seq[scales_idx]
    
    location_results<-list()
    scales_collected<-list()
    for(l in 1:length(loc_list)){
      target_loc<-loc_list[[l]]
      loc_pop<-pop%>%filter(location==target_loc)
      
      d<-recoveries_list[[l]]
      
      sim_result<-tryCatch({
        sim_scales_sampling(
          pop=loc_pop,
          spawning_results = d,
          iterations=iterations,
          scales_n=current_scales_n,
          CI=.95
        )
      },error=function(e){
        message(paste("Error in simulation:", e$message))
        return(NULL)
      })
      
      scales_collected[[l]]<-sim_result$scales_collected
      l_results <- list_rbind(sim_result$boot_results) %>%
        select(k_iteration, age, total_natural) %>%
        pivot_wider(
          id_cols = k_iteration,
          names_from = age,
          values_from = total_natural
        )
      location_results[[l]]=l_results
    }
    
    names(location_results)<-names(scales_collected)<-loc_list
    
    basin_results <- bind_rows(location_results, .id = "location") %>%
      group_by(k_iteration) %>%
      summarise(
        age_2 = sum(`2`, na.rm = TRUE),
        age_3 = sum(`3`, na.rm = TRUE),
        age_4 = sum(`4`, na.rm = TRUE)
      ) %>%
      ungroup()
    
    basin_summary <- basin_results %>%
      pivot_longer(
        cols = starts_with("age_"),
        names_to = "age",
        values_to = "value"
      ) %>%
      group_by(age) %>%
      summarise(
        mean_natural_count = mean(value, na.rm = TRUE),
        lower_CI_count = quantile(value, (1-CI)/2, na.rm = TRUE),
        upper_CI_count = quantile(value, 1-(1-CI)/2, na.rm = TRUE)
      ) %>%
      ungroup()
    
    basin_summary <- basin_summary %>%
      mutate(age = as.numeric(str_extract(age, "\\d+")))
    
    basin_summary$scales_collected=sum(unlist(scales_collected))
    basin_summary$target_location_scales=current_scales_n
    
    basin_summary<-basin_summary%>%left_join(true_pop_ages)
    
    scales_n_results[[scales_idx]]<-basin_summary
    
  }
  names(scales_n_results)=scales_n_seq
  
  iter_end<-Sys.time()
  iter_time<-iter_end-iter_start
  iter_results[[i]]<-list(
    "scales_results" = scales_n_results,
    "iter_time"=iter_time)
  
  print(paste("ITERATION",i,"DONE, time:",round(iter_time,2)))
}
saveRDS(iter_results,"outputs/basin_results_nopooling.Rds")



##########################
#with pooling
##########################
trib_data<-readRDS("data/all_returns_clean.Rds")
theta_data<-readRDS("data/location_year_theta.Rds")

year_list<-2010:2020

boot_replicates=500 #times to repeat total simulation
scales_n_seq =seq(2000,2000,1)
iter_results<-list()
pooling_iteration_results<-list()
iterations=100

for(i in 1:boot_replicates){
  iter_start<-Sys.time()
  
  #sample population
  iter_year<-sample(year_list,1)
  pop<-sim_basin_pop(trib_data,i,iter_year)
  sim_tag_rates<-pop$tag_rates
  pop<-pop$iter_pop
  loc_list<-unique(pop$location)#pull list of sampled locations
  
  #spawning recoveries for each location
  recoveries_list<-list()
  for(l in 1:length(loc_list)){ #loop recoveries sim for each location
    iter_loc<-loc_list[l]
    loc_pop<-pop%>%filter(location==iter_loc) #get pop data for locaiton
    
    target_theta<-theta_data%>% #get theta for location x year
      filter(location==iter_loc,
             return_year==iter_year)
    if(nrow(target_theta)==0){target_theta=0.2} else{ #if no theta, set default
      target_theta<-round(target_theta$theta,2)
    }
    
    target_tag_rate<-sim_tag_rates%>%
      filter(location==iter_loc,
             year==iter_year,
             origin=="hatchery")
    target_tag_rate<-target_tag_rate$tag_rate
    if(is.null(target_tag_rate)){target_tag_rate=0.25}
    
    #get spawning recoveries
    spawning_recoveries<-sim_spawning_recovery(loc_pop,
                                               theta=target_theta,
                                               tag_rate=target_tag_rate,
                                               iterations=iterations)
    spawning_recoveries$recovered_fish$theta=target_theta
    recoveries_list[[l]]<-spawning_recoveries
  }
  names(recoveries_list)<-loc_list
  
  all_recovered_fish<-all_hatchery_estimates<-data.frame()
  
  for(l in 1:length(recoveries_list)){
    all_recovered_fish<-all_recovered_fish%>%
      rbind(recoveries_list[[l]]$recovered_fish)
    hatcheries<-recoveries_list[[l]]$hatchery_ages
    if(nrow(hatcheries)==0){
      hatcheries<-data.frame(
        "age"=c(2,3,4),
        "N"=c(0,0,0),
        "total_tags"=c(0,0,0),
        "total_hatchery"=c(0,0,0),
        "k_iteration"=c(1,1,1)
      )
      hatcheries$location=loc_list[l]
    }else{
      hatcheries$location=loc_list[l]
    }
    all_hatchery_estimates<-all_hatchery_estimates%>%rbind(hatcheries)
  }
  
  setDT(all_hatchery_estimates)
  basin_hatchery_estimates <- all_hatchery_estimates[, .(
    N = sum(N),
    total_tags = sum(total_tags),
    total_hatchery = sum(total_hatchery)
  ), by = .(age, k_iteration)]
  
  #get locations proportional escapement
  #for this we want to estimate the contribution of untagged fish from each location
  #we can average the estimates from k-iterations
  basin_tagged<-all_hatchery_estimates%>%
    group_by(k_iteration)%>%
    summarise(total_tagged=sum(total_tags))
  basin_tagged=mean(basin_tagged$total_tagged)
  
  basin_untagged=nrow(pop)-basin_tagged
  
  loc_tagged<-all_hatchery_estimates%>%
    group_by(location,k_iteration)%>%
    summarise(total_tagged=sum(total_tags))
  loc_tagged<-loc_tagged%>%
    group_by(location)%>%
    summarise(total_tagged=mean(total_tagged))
  
  location_esc_prop <- pop %>%
    group_by(location) %>%
    summarize(esc = n())
  
  loc_untagged_prop<-loc_tagged%>%
    left_join(location_esc_prop)%>%
    mutate(total_untagged=esc-total_tagged)%>%
    mutate(untagged_prop=total_untagged/basin_untagged)
  
  #pre-alocate scale_iters list
  scale_iters <- vector("list", length(scales_n_seq))
  names(scale_iters) <- scales_n_seq
  
  
  for(scales_idx in 1:length(scales_n_seq)) {
    
    #get locational scales
    current_scales_n <- scales_n_seq[scales_idx]
    
    #use round_to_sum to get scales per location
    location_esc <- loc_untagged_prop %>%
      mutate(target_scales = round_to_sum(untagged_prop, current_scales_n))
    
    #find amount of scales actually available per location
    location_availability <- all_recovered_fish %>%
      filter(has_cwt == 0) %>%
      group_by(location) %>%
      summarise(available = n())
    
    location_esc <- location_esc %>%
      left_join(location_availability, by = "location") %>%
      mutate(
        available = ifelse(is.na(available), 0, available),
        #initially take min(allocation, available)
        actual_scales = pmin(target_scales, available),
        shortfall = target_scales - actual_scales
      )
    
    total_shortfall <- sum(location_esc$shortfall)
    
    untagged_recoveries<-all_recovered_fish%>%
      filter(has_cwt==0)
    loc_recovered_fish=data.frame()
    for(l in 1:length(loc_list)){
      target_loc<-loc_list[l]
      ld<-location_esc%>%filter(location==target_loc)
      
      available_scales<-untagged_recoveries%>%
        filter(location==target_loc)
      
      if(nrow(available_scales)>ld$target_scales){
        #get scales subsample from untagged fish based on scales_n input
        loc_scales<-sample_n(available_scales,ld$actual_scales,replace = F)
      } else {
        loc_scales <- available_scales  #take all available
      }
      loc_recovered_fish<-loc_recovered_fish%>%rbind(loc_scales)
    }
    
    actual_total_scales <- nrow(loc_recovered_fish)
    
    basin_spawning_results<-list("recovered_fish"=loc_recovered_fish,
                                 "hatchery_ages"=basin_hatchery_estimates)
    
    basin_sim_result<-sim_scales_sampling(
      pop=pop,
      spawning_results = basin_spawning_results,
      iterations=iterations,
      scales_n=current_scales_n,0.95)
    
    basin_sim_result$target_scales_n<-current_scales_n
    basin_sim_result$actual_scales_n<-actual_total_scales
    scale_iters[[scales_idx]]<-basin_sim_result
  }
  names(scale_iters)<-scales_n_seq
  iter_end<-Sys.time()
  print(paste("boot",i,"DONE, time:",round(iter_end-iter_start,2)))
  
  pooling_iteration_results[[i]]<-list(
    "scale_iters"=scale_iters)
}
saveRDS(pooling_iteration_results,"outputs/basin_results_wpooling_2000scales.Rds")

##########################
#with pooling no weighting
##########################
trib_data<-readRDS("data/all_returns_clean.Rds")
theta_data<-readRDS("data/location_year_theta.Rds")

year_list<-2010:2020

boot_replicates=100 #times to repeat total simulation
scales_n_seq =seq(500,2000,100)
iter_results<-list()
pooling_iteration_results<-list()
iterations=1000

for(i in 1:boot_replicates){
  iter_start<-Sys.time()
  
  #sample population
  iter_year<-sample(year_list,1)
  pop<-sim_basin_pop(trib_data,i,iter_year)
  sim_tag_rates<-pop$tag_rates
  pop<-pop$iter_pop
  loc_list<-unique(pop$location)#pull list of sampled locations
  
  #spawning recoveries for each location
  recoveries_list<-list()
  for(l in 1:length(loc_list)){ #loop recoveries sim for each location
    iter_loc<-loc_list[l]
    loc_pop<-pop%>%filter(location==iter_loc) #get pop data for locaiton
    
    target_theta<-theta_data%>% #get theta for location x year
      filter(location==iter_loc,
             return_year==iter_year)
    if(nrow(target_theta)==0){target_theta=0.2} else{ #if no theta, set default
      target_theta<-round(target_theta$theta,2)
    }
    
    target_tag_rate<-sim_tag_rates%>%
      filter(location==iter_loc,
             year==iter_year,
             origin=="hatchery")
    target_tag_rate<-target_tag_rate$tag_rate
    if(is.null(target_tag_rate)){target_tag_rate=0.25}
    
    #get spawning recoveries
    spawning_recoveries<-sim_spawning_recovery(loc_pop,
                                               theta=target_theta,
                                               tag_rate=target_tag_rate,
                                               iterations=iterations)
    spawning_recoveries$recovered_fish$theta=target_theta
    
    recoveries_list[[l]]<-spawning_recoveries
  }
  names(recoveries_list)<-loc_list
  
  all_recovered_fish<-all_hatchery_estimates<-data.frame()
  
  for(l in 1:length(recoveries_list)){
    all_recovered_fish<-all_recovered_fish%>%
      rbind(recoveries_list[[l]]$recovered_fish)
    all_hatchery_estimates<-all_hatchery_estimates%>%
      rbind(recoveries_list[[l]]$hatchery_ages)       
  }
  
  setDT(all_hatchery_estimates)
  basin_hatchery_estimates <- all_hatchery_estimates[, .(
    N = sum(N),
    total_tags = sum(total_tags),
    total_hatchery = sum(total_hatchery)
  ), by = .(age, k_iteration)]
  
  #pre-alocate scale_iters list
  scale_iters <- vector("list", length(scales_n_seq))
  names(scale_iters) <- scales_n_seq
  basin_spawning_results<-list("recovered_fish"=all_recovered_fish,
                               "hatchery_ages"=basin_hatchery_estimates)
  
  for(scales_idx in 1:length(scales_n_seq)) {
    current_scales_n=scales_n_seq[[scales_idx]]
    
    basin_sim_result<-sim_scales_sampling(
      pop=pop,
      spawning_results = basin_spawning_results,
      iterations=iterations,
      scales_n=current_scales_n)
    
    basin_sim_result$target_scales_n<-current_scales_n
    scale_iters[[scales_idx]]<-basin_sim_result
  }
  names(scale_iters)<-scales_n_seq
  print(paste(i, round(Sys.time()-iter_start,4)))
  
  pooling_iteration_results[[i]]<-list(
    "scale_iters"=scale_iters)
}
#saveRDS(pooling_iteration_results,"outputs/basin_results_wpooling_noweighting.Rds")

##########################
#With pooling tag_rate 100%
##########################

trib_data<-readRDS("data/all_returns_clean.Rds")
theta_data<-readRDS("data/location_year_theta.Rds")

year_list<-2010:2020

boot_replicates=100 #times to repeat total simulation
scales_n_seq =seq(100,2000,100)
iter_results<-list()
pooling_iteration_results<-list()
iterations=100

for(i in 1:boot_replicates){
  iter_start<-Sys.time()

  #sample population
  iter_year<-sample(year_list,1)
  pop<-sim_basin_pop(trib_data,i,iter_year)
  pop<-pop$iter_pop
  pop<-pop%>%
    mutate(has_cwt=ifelse(origin=="hatchery",1,0))
  loc_list<-unique(pop$location)#pull list of sampled locations
  
  #spawning recoveries for each location
  recoveries_list<-list()
  for(l in 1:length(loc_list)){ #loop recoveries sim for each location
    iter_loc<-loc_list[l]
    loc_pop<-pop%>%filter(location==iter_loc) #get pop data for locaiton
    
    target_theta<-theta_data%>% #get theta for location x year
      filter(location==iter_loc,
             return_year==iter_year)
    if(nrow(target_theta)==0){target_theta=0.2} else{ #if no theta, set default
      target_theta<-round(target_theta$theta,2)
    }
    
    target_tag_rate<-1
    
    #get spawning recoveries
    spawning_recoveries<-sim_spawning_recovery(loc_pop,
                                               theta=target_theta,
                                               tag_rate=target_tag_rate,
                                               iterations=iterations)
    spawning_recoveries$recovered_fish$theta=target_theta
    recoveries_list[[l]]<-spawning_recoveries
  }
  names(recoveries_list)<-loc_list
  
  all_recovered_fish<-all_hatchery_estimates<-data.frame()
  
  for(l in 1:length(recoveries_list)){
    all_recovered_fish<-all_recovered_fish%>%
      rbind(recoveries_list[[l]]$recovered_fish)
    hatcheries<-recoveries_list[[l]]$hatchery_ages
    if(nrow(hatcheries)==0){
      hatcheries<-data.frame(
        "age"=c(2,3,4),
        "N"=c(0,0,0),
        "total_tags"=c(0,0,0),
        "total_hatchery"=c(0,0,0),
        "k_iteration"=c(1,1,1)
      )
      hatcheries$location=loc_list[l]
    }else{
      hatcheries$location=loc_list[l]
    }
    all_hatchery_estimates<-all_hatchery_estimates%>%rbind(hatcheries)
  }
  
  setDT(all_hatchery_estimates)
  basin_hatchery_estimates <- all_hatchery_estimates[, .(
    N = sum(N),
    total_tags = sum(total_tags),
    total_hatchery = sum(total_hatchery)
  ), by = .(age, k_iteration)]
  
  #get locations proportional escapement
  #for this we want to estimate the contribution of untagged fish from each location
  #we can average the estimates from k-iterations
  basin_tagged<-all_hatchery_estimates%>%
    group_by(k_iteration)%>%
    summarise(total_tagged=sum(total_tags))
  basin_tagged=mean(basin_tagged$total_tagged)
  
  basin_untagged=nrow(pop)-basin_tagged
  
  loc_tagged<-all_hatchery_estimates%>%
    group_by(location,k_iteration)%>%
    summarise(total_tagged=sum(total_tags))
  loc_tagged<-loc_tagged%>%
    group_by(location)%>%
    summarise(total_tagged=mean(total_tagged))
  
  location_esc_prop <- pop %>%
    group_by(location) %>%
    summarize(esc = n())
  
  loc_untagged_prop<-loc_tagged%>%
    left_join(location_esc_prop)%>%
    mutate(total_untagged=esc-total_tagged)%>%
    mutate(untagged_prop=total_untagged/basin_untagged)
  
  #pre-alocate scale_iters list
  scale_iters <- vector("list", length(scales_n_seq))
  names(scale_iters) <- scales_n_seq
  
  for(scales_idx in 1:length(scales_n_seq)) {
    
    #get location scales
    current_scales_n <- scales_n_seq[scales_idx]
    
    #use round_to_sum to get scales per location
    location_contribution <- loc_untagged_prop %>%
      mutate(target_scales = round_to_sum(untagged_prop, current_scales_n))
    
    #find amount of scales actually available per location
    location_availability <- all_recovered_fish %>%
      filter(has_cwt == 0) %>%
      group_by(location) %>%
      summarise(available = n())
    
    location_contribution <- location_contribution %>%
      left_join(location_availability, by = "location") %>%
      mutate(
        available = ifelse(is.na(available), 0, available),
        # Initially take min(allocation, available)
        actual_scales = pmin(target_scales, available),
        shortfall = target_scales - actual_scales
      )
    
    total_shortfall <- sum(location_contribution$shortfall)
    
    untagged_recoveries<-all_recovered_fish%>%
      filter(has_cwt==0)
    loc_recovered_fish=data.frame()
    for(l in 1:length(loc_list)){
      target_loc<-loc_list[l]
      ld<-location_contribution%>%filter(location==target_loc)
      ld<-ld%>%
        mutate(actual_scales=ifelse(actual_scales<0,0,actual_scales))
      
      available_scales<-untagged_recoveries%>%
        filter(location==target_loc)
      
      if(nrow(available_scales)>ld$target_scales){
        #get scales subsample from untagged fish based on scales_n input
        loc_scales<-sample_n(available_scales,ld$actual_scales,replace = F)
      } else {
        loc_scales <- available_scales  #take all available
      }
      loc_recovered_fish<-loc_recovered_fish%>%rbind(loc_scales)
    }
    
    actual_total_scales <- nrow(loc_recovered_fish)
    
    basin_spawning_results<-list("recovered_fish"=loc_recovered_fish,
                                 "hatchery_ages"=basin_hatchery_estimates)
    
    basin_sim_result<-sim_scales_sampling(
      pop=pop,
      spawning_results = basin_spawning_results,
      iterations=iterations,
      scales_n=current_scales_n)
    
    basin_sim_result$target_scales_n<-current_scales_n
    basin_sim_result$actual_scales_n<-actual_total_scales
    scale_iters[[scales_idx]]<-basin_sim_result
  }
  names(scale_iters)<-scales_n_seq
  
  cat(red(paste("\niteration",i,"DONE, time:",round(Sys.time()-iter_start,4),"\n")))
  
  pooling_iteration_results[[i]]<-list(
    "scale_iters"=scale_iters)
}

saveRDS(pooling_iteration_results,"outputs/basin_results_wpooling_tagrate100_boot100.Rds")

#########################
#TESTING, trying to figure out why tag_rate =1 results in systematic bias?
###########################
results<-data.frame()
for(i in 1:100){
  #sample population
  iter_year<-sample(year_list,1)
  pop<-sim_basin_pop(trib_data,i,iter_year)
  pop<-pop$iter_pop
  pop<-pop%>%
    mutate(has_cwt=ifelse(origin=="hatchery",1,0))
  loc_list<-unique(pop$location)#pull list of sampled locations
  
  #spawning recoveries for each location
  recoveries_list<-list()
  for(l in 1:length(loc_list)){ #loop recoveries sim for each location
    iter_loc<-loc_list[l]
    loc_pop<-pop%>%filter(location==iter_loc) #get pop data for locaiton
    
    target_theta<-theta_data%>% #get theta for location x year
      filter(location==iter_loc,
             return_year==iter_year)
    if(nrow(target_theta)==0){target_theta=0.2} else{ #if no theta, set default
      target_theta<-round(target_theta$theta,2)
    }
    
    target_tag_rate<-1
    
    #get spawning recoveries
    spawning_recoveries<-sim_spawning_recovery(loc_pop,
                                               theta=target_theta,
                                               tag_rate=target_tag_rate,
                                               iterations=iterations)
    spawning_recoveries$recovered_fish$theta=target_theta
    recoveries_list[[l]]<-spawning_recoveries
  }
  names(recoveries_list)<-loc_list
  
  #test location recoveries and hatchery estimates
  true_loc_hatchery<-pop%>%
    filter(origin=="hatchery")%>%
    group_by(location,age)%>%
    summarise(true_count=n())
  
  location_hatchery_estimates<-data.frame()
  for(l in 1:length(recoveries_list)){
    d<-recoveries_list[[l]]$hatchery_ages
    if(nrow(d)>0){
      d$location=names(recoveries_list[l])
      d$theta<-unique(recoveries_list[[l]]$recovered_fish$theta)
      location_hatchery_estimates<-location_hatchery_estimates%>%
        rbind(d)
    }
  }
  location_hatchery_estimates<-location_hatchery_estimates%>%
    left_join(true_loc_hatchery)%>%
    mutate(bias=total_hatchery-true_count)%>%
    mutate(rel_bias=bias/true_count)
  
  location_hatchery_estimates<-location_hatchery_estimates%>%
    group_by(location,age)%>%
    summarise(mean_bias=mean(rel_bias))
  results<-results%>%rbind(location_hatchery_estimates)
  print(i)
}
ggplot(results)+
  geom_boxplot(aes(x=factor(age),y=mean_bias))+
  facet_grid(.~location)
test<-results%>%
  group_by(location,age)%>%
  summarize(mean_bias=mean(mean_bias))
            