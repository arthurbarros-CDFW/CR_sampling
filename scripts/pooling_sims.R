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
k_iterations=100 #
scale_iterations=100 # needs to be same as k_iterations (consider just one iter value)
survey="trib" #hatchery or trib
scales_n_seq =seq(500,500,1)
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
  
  target_year=sample(year_list,1)
  pop<-sim_location_pop(trib_data,iter,target_year)
  loc_list<-unique(pop$location)
    
  recovered_fish <- vector("list", length(loc_list))
  names(recovered_fish) <- loc_list
    
  hatchery_estimates <- vector("list", length(loc_list))
  names(hatchery_estimates) <- loc_list
    
  sim_results<- vector("list", length(loc_list))
  names(sim_results) <- loc_list
    
  spawning_list<-list()
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
      
    spawning_results<-sim_spawning_recovery(loc_pop,
                                            target_theta,tag_rate,
                                            k_iterations)
    spawning_list[[l]]<-spawning_results
  }
  names(spawning_list)<-loc_list
  
  boot_results<-list()

  for(s in 1:scale_iterations){
    scale_iteration_results<-list()
    scales_collected<-list()
    for(l in 1:length(loc_list)){
      target_loc<-loc_list[[l]]
      loc_pop<-pop%>%filter(location==target_loc)
      sim_result<-tryCatch({
        sim_scales_sampling(
          pop=loc_pop,
          spawning_results = spawning_list[[target_loc]],
          theta=target_theta,
          k_iterations=k_iterations,
          scale_iterations=1,#here we set to 1 and do our scale resampling every iteration?
          tag_rate=tag_rate,
          survey="trib",
          scales_n=scales_n_seq,
          CI=CI
        )
      },error=function(e){
        message(paste("Error in simulation:", e$message))
        return(NULL)
      })
      
      scale_iteration_results[[l]]<-sim_result
      scales_collected[[l]]<-sim_result$scales_collected
    }
    
    names(scale_iteration_results)<-loc_list
    
    regions<-loc_list
    
    basin_sim_ages<-map_dfr(regions, function(region) {
      if (!is.null(scale_iteration_results[[region]]$age_summary_stats)) {
        scale_iteration_results[[region]]$age_summary_stats %>%
          mutate(
            region = region
          )
      }
    })
    
    basin_sim_ages <- basin_sim_ages %>%
      group_by(age) %>%
      summarise(
        est_count_hatchery = sum(est_count_hatchery, na.rm = TRUE),
        est_count_natural = sum(est_count_natural, na.rm = TRUE),
        .groups = "drop"
      )
    
    basin_true_ages<-map_dfr(regions, function(region) {
      if (!is.null(scale_iteration_results[[region]]$true_summary_stats)) {
        scale_iteration_results[[region]]$true_summary_stats %>%
          mutate(
            region = region
          )
      }
    })
    
    basin_true_ages <- basin_true_ages %>%
      group_by(age) %>%
      summarise(
        true_count_hatchery = sum(true_count_hatchery, na.rm = TRUE),
        true_count_natural = sum(true_count_natural, na.rm = TRUE),
        .groups = "drop"
      )
    
    basin_sim_ages$est_total_natural<-sum(basin_sim_ages$est_count_natural)
    basin_sim_ages$est_total_hatchery<-sum(basin_sim_ages$est_count_hatchery)
    
    basin_true_ages$true_total_natural<-sum(basin_true_ages$true_count_natural)
    basin_true_ages$true_total_hatchery<-sum(basin_true_ages$true_count_hatchery)
    
    basin_sim_ages<-basin_sim_ages%>%
      mutate(est_proportion_natural=est_count_natural/est_total_natural,
             est_proportion_hatchery=est_count_hatchery/est_total_hatchery)
    basin_true_ages<-basin_true_ages%>%
      mutate(true_proportion_natural=true_count_natural/true_total_natural,
             true_proportion_hatchery=true_count_hatchery/true_total_hatchery)
    
    scales_collected<-sum(unlist(scales_collected))
    
    
    
    scale_results<-list("basin_sim_ages"=basin_sim_ages,
                        "basin_true_ages"=basin_true_ages,
                        "scales_collected"=scales_collected,
                        "scale_iteration_results"=scale_iteration_results)
    boot_results[[s]]<-scale_results
  }
  saveRDS(boot_results,paste("outputs/pooling/no_pooling_iter",iter,".Rds",sep="_"))

  all_basin_sim_ages <- imap_dfr(boot_results, function(x, idx) {
    x$basin_sim_ages %>%
      mutate(
        scale_iter = idx
      )
  })
  
  all_basin_true_ages <- imap_dfr(boot_results, function(x, idx) {
    x$basin_true_ages %>%
      mutate(
        scale_iter = idx
      )
  })
  
  all_basin_true_ages<-unique(select(all_basin_true_ages,-scale_iter))
  
  sim_summary <- data.table(all_basin_sim_ages)[
    !is.na(est_proportion_natural),
    .(
      est_count_hatchery=mean(est_count_hatchery),
      est_count_natural=mean(est_count_natural),
      mean_proportion_natural = mean(est_proportion_natural, na.rm = TRUE),
      median_proportion_natural = median(est_proportion_natural, na.rm = TRUE),
      sd_proportion_natural = sd(est_proportion_natural, na.rm = TRUE),
      se_proportion_natural = sd(est_proportion_natural, na.rm = TRUE) / sqrt(.N),
      lower_CI_prop = quantile(est_proportion_natural, probs = (1-CI)/2, na.rm = TRUE),
      upper_CI_prop = quantile(est_proportion_natural, probs = 1-(1-CI)/2, na.rm = TRUE),
      lower_CI_count = quantile(est_count_natural, probs = (1-CI)/2, na.rm = TRUE),
      upper_CI_count = quantile(est_count_natural, probs = 1-(1-CI)/2, na.rm = TRUE),
      n_iterations = .N
    ),
    by = age
  ]
  
  sim_summary<-sim_summary%>%
    mutate(moe_prop=(upper_CI_prop-lower_CI_prop)/2)
  
  iter_end<-Sys.time()
  iter_time<-iter_end-iter_start
  nopooling_iteration_results[[i]]<-list(
    "basin_sim_summary"=sim_summary,
    "basin_scales_collected"=scales_collected,
    "basin_true_ages"=basin_true_ages,
    "iter_time"=iter_time)
  
  print(paste("ITERATION",iter,"DONE, time:",round(iter_time,2)))
}

saveRDS(nopooling_iteration_results,"outputs/pooling/results_no_pooling_test.Rds")

##########################
#STEP 2: Produce basin pooling estimates
##########################
#we can use the recovered fish list from each loc x year
#and randomly select from among all the loc recovered fish in a year
#to produce our basin estimates
tag_rate=0.25 #fixed CFM rate
k_iterations=100 #
scale_iterations=100 # needs to be same as k_iterations (consider just one iter value)
survey="trib" #hatchery or trib
scales_n_seq =seq(500,5000,500)
iterations=100 #times to repeat total simulation
CI=0.95
h_pct=75
h_prop=.75
target_moe=0.20
trib_data<-readRDS("data/all_returns_clean.Rds")
theta_data<-readRDS("data/location_year_theta.Rds")
iter=1

year_list<-2010:2020

pooling_iteration_results<-list()

for( i in 1:iterations){
  iter<-i
  iter_start<-Sys.time()
  
  target_year=sample(year_list,1)
  pop<-sim_location_pop(trib_data,iter,target_year)
  loc_list<-unique(pop$location)
  theta_year <- theta_data[theta_data$return_year == target_year, ]
  
  #gather some pop sim info
  pop_dt <- as.data.table(pop)
  n_ho <- pop_dt[origin == "hatchery", .N]
  n_no <- pop_dt[origin == "natural", .N]
  n_tags <- pop_dt[has_cwt == 1, .N]
  
  pop_totals<-list("n_ho"=n_ho,
                   "n_no"=n_no,
                   "n_tags"=n_tags)
  
  #get list of all recovered fish across locations
  sim_results<-vector("list", length(loc_list))
  names(sim_results) <- loc_list
  n_locs <- length(loc_list)
  recovered_fish<-hatchery_estimates<- basin_sim_results<- vector("list", n_locs)
  names(recovered_fish) <-names(hatchery_estimates)<-names(basin_sim_results)<- loc_list
  for(l in 1:length(loc_list)){
    time_start<-Sys.time()
    target_loc<-loc_list[l]
    loc_pop<-pop%>%filter(location==target_loc)
      
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
      
    recovered_fish[[l]] <- spawning_results$recovered_fish
    hatchery_estimates[[l]] <- spawning_results$hatchery_ages
    hatchery_estimates[[l]]<-hatchery_estimates[[l]]%>%
      mutate(location=target_loc,
             theta=target_theta)
    recovered_fish[[l]]<-recovered_fish[[l]]%>%
      mutate(theta=target_theta)
    
    time_end<-Sys.time()
    sim_time<-time_end-time_start
    #print(paste(iter,target_loc,target_year,round(sim_time,2)))
  }
  #take recovered fish and hatchery estimates from all regions for target year
  all_recovered_fish<-rbindlist(recovered_fish)
  all_hatchery_estimates<-rbindlist(hatchery_estimates)
    
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
    time_start<-Sys.time()  
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
      basin_sim_result$age_summary_stats<-basin_sim_result$age_summary_stats%>%
        mutate(moe_prop=(upper_CI_prop-lower_CI_prop)/2)
      
      scale_iters[[s]]<-basin_sim_result
      time_end<-Sys.time()
      scale_time<-time_end-time_start
  }
  print(paste("scale boots:",round(scale_time,2)))
  names(scale_iters)<-scales_n_seq
    
  saveRDS(scale_iters,paste("outputs/pooling/pooling_iter",iter,".Rds",sep="_"))
  
  iter_end<-Sys.time()
  iter_time<-iter_end-iter_start
  
  pooling_iteration_results[[i]]<-list(
    "scale_iters"=scale_iters,
    "iter_time"=iter_time)
  
  print(paste("ITERATION",iter,"DONE, time:",round(iter_time,2)))
}

saveRDS(pooling_iteration_results,"outputs/pooling/results_with_pooling.Rds")

#visualize while running pooling sim
est_test<-list()
true_test<-list()
for(r in 1:length(scale_iters)){
  est_test[[r]]<-scale_iters[[r]]$age_summary_stats
  true_test[[r]]<-scale_iters[[r]]$true_summary_stats
}
est_test<-rbindlist(est_test)
true_test<-rbindlist(true_test)
true_test<-unique(true_test)
ggplot(est_test,aes(x=factor(age),
                                  y=(est_count_natural),
                                  color=factor(age)))+
  #geom_bar(stat="identity")+
  geom_point(size=3)+
  geom_errorbar(aes(ymin = lower_CI_count, ymax = upper_CI_count), 
                width = 0.2) +
  geom_point(data=true_test,aes(x=factor(age),
                            y=true_count_natural),size=3,color="black",shape=9)+
  theme_bw()+
  labs(title="basin wide age count estimates",
       subtitle = paste("scales collected: ",unique(test$scales_collected)),
       x = "age", 
       y = "estimate of natural origin count",
       color="age")+
  facet_grid(.~scales_collected)
ggsave(file="outputs/basin_count_estimates_nscales.png",width=10,height=6)
