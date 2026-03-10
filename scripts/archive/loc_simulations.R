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
#STEP 1: for each year and location, get spawning results
#produce list of all recovered fish available
#produce list of hatchery origin estimates
#produce list of sim_results
#perform yearly basin wide estimate by summing individual locations
##########################
#The below should give us a list of all the estimates for each loc x year
#we can use those age count estimates and join them to get basin wide
#this is our "do n_scales at each location approach"

tag_rate=0.25 #fixed CFM rate
k_iterations=10 #
scale_iterations=10 # needs to be same as k_iterations (consider just one iter value)
survey="trib" #hatchery or trib
scales_n_seq =seq(1000,1000,1)
iterations=100 #times to repeat total simulation
CI=0.95
h_pct=75
h_prop=.75
target_moe=0.20

trib_data<-readRDS("data/all_returns_clean.Rds")
theta_data<-readRDS("data/location_year_theta.Rds")
iter=1

year_list<-2010:2020

nopooling_iteration_results<-list()

for( i in 1:iterations){
  iter<-i
  iter_start<-Sys.time()
  
  recovered_fish<-sim_results<-hatchery_estimates<-basin_estimates<-vector("list", length(year_list))
  
  names(recovered_fish) <-names(hatchery_estimates)<-names(sim_results)<-names(basin_estimates)<- year_list
  
  for(y in 1:length(year_list)){
    target_year=year_list[y]
    pop<-sim_location_pop(trib_data,iter,target_year)
    loc_list<-unique(pop$location)
    
    recovered_fish[[y]] <- vector("list", length(loc_list))
    names(recovered_fish[[y]]) <- loc_list
    
    hatchery_estimates[[y]] <- vector("list", length(loc_list))
    names(hatchery_estimates[[y]]) <- loc_list
    
    sim_results[[y]]<- vector("list", length(loc_list))
    names(sim_results[[y]]) <- loc_list
    
    for(l in 1:length(loc_list)){
      time_start<-Sys.time()
      
      target_loc<-loc_list[l]
      loc_pop<-pop%>%filter(location==target_loc)
      target_theta<-theta_data%>%
        filter(location==target_loc,
               return_year==target_year)
      if(nrow(target_theta)==0){target_theta=0.2} else{
        target_theta<-round(target_theta$theta,2)
      }
      
      #gather some pop sim info
      n_ho<-nrow(loc_pop%>%filter(origin=="hatchery"))
      n_no<-nrow(loc_pop%>%filter(origin=="natural"))
      n_tags<-nrow(loc_pop%>%filter(has_cwt==1))
      
      spawning_results<-sim_spawning_recovery(loc_pop,
                                              target_theta,tag_rate,
                                              k_iterations)
      
      sim_result<-tryCatch({
        sim_scales_sampling(
          pop=loc_pop,
          spawning_results = spawning_results,
          theta=target_theta,
          k_iterations=k_iterations,
          scale_iterations=scale_iterations,
          tag_rate=tag_rate,
          survey="trib",
          scales_n=scales_n_seq,
          CI=CI
        )
      },error=function(e){
        message(paste("Error in simulation:", e$message))
        return(NULL)
      })
      time_end<-Sys.time()
      sim_time<-time_end-time_start
      
      pop_totals<-list("n_ho"=n_ho,
                       "n_no"=n_no,
                       "n_tags"=n_tags)
      
      sim_result$pop_totals<-pop_totals
      sim_result$age_summary_stats$scales_collected<-sim_result$scales_collected
      
      recovered_fish[[y]][[l]] <- spawning_results$recovered_fish
      hatchery_estimates[[y]][[l]] <- spawning_results$hatchery_ages
      sim_results[[y]][[l]]<-sim_result
      
      time_end<-Sys.time()
      sim_time<-time_end-time_start
      print(paste(iter,target_loc,target_year,round(sim_time,2)))
    }
  }
  
  saveRDS(sim_results,paste("outputs/pooling/no_pooling_iter",iter,".Rds",sep="_"))
  
  all_sim_data <- process_all_years(sim_results,true_values = FALSE)
  all_true_data<- process_all_years(sim_results,true_values = TRUE)
  
  basin_sim_ages <- all_sim_data %>%
    group_by(year, age) %>%
    summarise(
      total_est_count_hatchery = sum(est_count_hatchery, na.rm = TRUE),
      total_est_count_natural = sum(est_count_natural, na.rm = TRUE),
      .groups = "drop"
    )
  basin_scales_collected<-unique(select(all_sim_data,year,region,scales_collected))
  basin_scales_total<-basin_scales_collected%>%
    group_by(year)%>%
    summarise(scales_collected=sum(scales_collected))
  
  basin_true_ages <- all_true_data %>%
    group_by(year, age) %>%
    summarise(
      total_true_count_hatchery = sum(true_count_hatchery, na.rm = TRUE),
      total_true_count_natural = sum(true_count_natural, na.rm = TRUE),
      .groups = "drop"
    )
  
  iter_end<-Sys.time()
  iter_time<-iter_end-iter_start
  
  nopooling_iteration_results[[i]]<-list(
    "basin_sim_ages"=basin_sim_ages,
    "basin_scales_total"=basin_scales_total,
    "basin_true_ages"=basin_true_ages,
    "all_sim_data"=all_sim_data,
    "all_true_data"=all_true_data,
    "iter_time"=iter_time)
  
  print(paste("ITERATION",iter,"DONE, time:",round(iter_time,2)))
}

saveRDS(nopooling_iteration_results,"outputs/pooling/results_no_pooling.Rds")

##########################
#STEP 2: Produce basin pooling estimates
##########################
#we can use the recovered fish list from each loc x year
#and randomly select from among all the loc recovered fish in a year
#to produce our basin estimates
tag_rate=0.25 #fixed CFM rate
k_iterations=10 #
scale_iterations=10 # needs to be same as k_iterations (consider just one iter value)
survey="trib" #hatchery or trib
scales_n_seq =seq(1000,5000,1000)
iterations=100 #times to repeat total simulation
CI=0.95
h_pct=75
h_prop=.75
target_moe=0.20

pooling_iteration_results<-list()

for( i in 1:iterations){
  iter<-i
  iter_start<-Sys.time()
  
  n_years <- length(year_list)
  recovered_fish <- vector("list", n_years)
  hatchery_estimates <- vector("list", n_years)
  basin_sim_results <- vector("list", n_years)
  
  names(recovered_fish) <-names(hatchery_estimates)<-names(basin_sim_results)<-names(basin_estimates)<- year_list
  
  for(y in 1:length(year_list)){
    target_year=year_list[y]
    pop<-sim_location_pop(trib_data,iter,target_year)
    
    loc_list<-unique(pop$location)
    n_locs <- length(loc_list)
    sim_results<-vector("list", length(loc_list))
    names(sim_results) <- loc_list
    
    #pre-allocate location-level lists
    recovered_fish[[y]] <- vector("list", n_locs)
    hatchery_estimates[[y]] <- vector("list", n_locs)
    names(recovered_fish[[y]]) <- names(hatchery_estimates[[y]]) <- loc_list

    #gather some pop sim info
    pop_dt <- as.data.table(pop)
    n_ho <- pop_dt[origin == "hatchery", .N]
    n_no <- pop_dt[origin == "natural", .N]
    n_tags <- pop_dt[has_cwt == 1, .N]
    
    pop_totals<-list("n_ho"=n_ho,
                     "n_no"=n_no,
                     "n_tags"=n_tags)
    
    theta_year <- theta_data[theta_data$return_year == target_year, ]
    
    for(l in 1:length(loc_list)){
      time_start<-Sys.time()
      
      target_loc<-loc_list[l]
      loc_pop <- pop[pop$location == target_loc, ]
      
      target_theta <- theta_year[theta_year$location == target_loc, "theta"]
      if(nrow(target_theta) == 0) {
        target_theta <- 0.2
      } else {
        target_theta <- round(target_theta, 2)
      }
      target_theta<-as.numeric(target_theta[[1]])
      
      spawning_results<-sim_spawning_recovery(loc_pop,
                                              target_theta,tag_rate,
                                              k_iterations)
      
      recovered_fish[[y]][[l]] <- spawning_results$recovered_fish
      hatchery_estimates[[y]][[l]] <- spawning_results$hatchery_ages
      

      time_end<-Sys.time()
      sim_time<-time_end-time_start
      print(paste(iter,target_loc,target_year,round(sim_time,2)))
    }
    
    #take recovered fish and hatchery estimates from all regions for target year
    all_recovered_fish<-rbindlist(recovered_fish[[y]])
    all_hatchery_estimates<-rbindlist(hatchery_estimates[[y]])
    
    setDT(all_hatchery_estimates)
    basin_hatchery_estimates <- all_hatchery_estimates[, .(
      N = sum(N),
      total_tags = sum(total_tags),
      total_hatchery = sum(total_hatchery)
    ), by = .(age, k_iteration)]
    
    basin_spawning_results<-list("recovered_fish"=all_recovered_fish,
                              "hatchery_ages"=basin_hatchery_estimates)
    
    #pre-alocate scale_iters list
    scale_iters <- vector("list", length(scales_n_seq))
    names(scale_iters) <- scales_n_seq
    
    for(s in 1:length(scales_n_seq)){
      basin_sim_result<-tryCatch({
        sim_scales_sampling(
          pop=pop,
          spawning_results = basin_spawning_results,
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
      
      basin_sim_result$pop_totals<-pop_totals
      basin_sim_result$age_summary_stats$scales_collected<-basin_sim_result$scales_collected
      scale_iters[[s]]<-basin_sim_result
    }
    names(scale_iters)<-scales_n_seq
    
    time_end<-Sys.time()
    sim_time<-time_end-time_start
    
    basin_sim_results[[y]]<-scale_iters
  }
  
  saveRDS(basin_sim_results,paste("outputs/pooling/pooling_iter",iter,".Rds",sep="_"))
  
  basin_sim_ages_list <- vector("list", n_years)
  basin_sim_scales_list <- vector("list", n_years)
  basin_true_ages_list <- vector("list", n_years)
  
  for(y in 1:length(n_years)) {
    
    if(!is.null(basin_sim_results[[y]])) {
        
      sim_dt_list <- vector("list", length(scales_n_seq))
      scales_dt_list <- vector("list", length(scales_n_seq))
      
      for(s in seq_along(scales_n_seq)){
      #get age summary stats and add year column
        year_sim_data <- as.data.table(basin_sim_results[[y]][[s]]$age_summary_stats)
        year_sim_data[, year := year_list[y]]
        sim_dt_list[[s]] <- year_sim_data
      
        d_sim<-d_sim%>%rbind(year_sim_data)
        
        scales_dt_list[[s]] <- data.table(
          scales_collected = unique(basin_sim_results[[y]][[s]]$scales_collected),
          year = year_list[y]
        )
      }
      basin_sim_ages_list[[y]] <- rbindlist(sim_dt_list, use.names = TRUE, fill = TRUE)
      basin_sim_scales_list[[y]] <- rbindlist(scales_dt_list, use.names = TRUE)
      
      year_true_data <- as.data.table(basin_sim_results[[y]][[1]]$true_summary_stats)
      year_true_data[, year := year_list[y]]
      basin_true_ages_list[[y]] <- year_true_data
    }
    }
  
  #combine all years into one dataframe
  basin_sim_ages <- rbindlist(basin_sim_ages_list, use.names = TRUE, fill = TRUE)[
    , .(year, age, est_count_hatchery, est_count_natural, scales_collected)
  ]
  
  basin_scales_total <- rbindlist(basin_sim_scales_list, use.names = TRUE)
  
  basin_true_ages <- rbindlist(basin_true_ages_list, use.names = TRUE)[
    , .(year, age, true_count_hatchery, true_count_natural)
  ]
    
  iter_end<-Sys.time()
  iter_time<-iter_end-iter_start
  
  pooling_iteration_results[[i]]<-list(
    "basin_sim_ages"=basin_sim_ages,
    "basin_scales_total"=basin_scales_total,
    "basin_true_ages"=basin_true_ages,
    "basin_sim_results"=basin_sim_results,
    "iter_time"=iter_time)
  
  print(paste("ITERATION",iter,"DONE, time:",round(iter_time,2)))
}

saveRDS(pooling_iteration_results,"outputs/pooling/results_with_pooling.Rds")

#plotting
ggplot(pooling_iteration_results[[3]]$basin_sim_ages,
       aes(x=factor(year),y=est_count_natural,colour = factor(age)))+
  geom_point()

ggplot(pooling_iteration_results[[3]]$basin_true_ages,
       aes(x=factor(year),y=true_count_natural,colour = factor(age)))+
  geom_point()

#more viz

#1) non-pooling first
nopooling_iteration_results<-readRDS("outputs/pooling/results_no_pooling.Rds")

sim_list <- basin_scales<-list()

for(i in 1:length(nopooling_iteration_results)) {
  #for each year combine all iterations
  iteration_est <- nopooling_iteration_results[[i]]$basin_sim_ages
  iteration_true<-nopooling_iteration_results[[i]]$basin_true_ages
  iteration_est$iter<-iteration_true$iter<-i
  sim_list[[i]] <- iteration_est%>%left_join(iteration_true)
  basin_scales[[i]]<-nopooling_iteration_results[[i]]$basin_scales_total
}

#combine all iterations
all_est <- rbindlist(sim_list)
all_scales<-rbindlist(basin_scales)

#get bias estimates
all_est<-all_est%>%
  mutate(iter_relative_bias=(total_est_count_natural-total_true_count_natural)/total_true_count_natural*100)
  
no_pooling_summarize_all<-all_est%>%
  group_by(year,age)%>%
  summarize(mean_bias=mean(iter_relative_bias),
            lower_CI_bias = quantile(iter_relative_bias, probs = (1-CI)/2, na.rm = TRUE),
            upper_CI_bias = quantile(iter_relative_bias, probs = 1-(1-CI)/2, na.rm = TRUE))


#1) pooling next
pooling_iteration_results<-readRDS("outputs/pooling/results_with_pooling.Rds")

sim_list <- basin_scales<-list()

for(i in 1:length(pooling_iteration_results)) {
  #for each year combine all iterations
  iteration_est <- pooling_iteration_results[[i]]$basin_sim_ages
  iteration_true<-pooling_iteration_results[[i]]$basin_true_ages
  iteration_est$iter<-iteration_true$iter<-i
  sim_list[[i]] <- iteration_est%>%left_join(iteration_true)
  basin_scales[[i]]<-pooling_iteration_results[[i]]$basin_scales_total
}

#combine all iterations
all_est <- rbindlist(sim_list)
all_scales<-rbindlist(basin_scales)

#get bias estimates
all_est<-all_est%>%
  mutate(iter_relative_bias=(est_count_natural-true_count_natural)/true_count_natural*100)

pooling_summarize_all<-all_est%>%
  group_by(year,age)%>%
  summarize(mean_bias=mean(iter_relative_bias),
            lower_CI_bias = quantile(iter_relative_bias, probs = (1-CI)/2, na.rm = TRUE),
            upper_CI_bias = quantile(iter_relative_bias, probs = 1-(1-CI)/2, na.rm = TRUE))
