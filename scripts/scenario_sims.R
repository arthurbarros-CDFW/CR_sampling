rm( list = ls()) #clear env
#data simulations
library(tidyverse)
library(data.table)
library(ggplot2)
library(gtools)

#Load functions
sapply(list.files("scripts/functions", pattern = "\\.R$", full.names = TRUE), source)

set.seed(67) #for reproducibility

##########################
#scale sample size run
##########################
#simulate population
N=100000 #true population size
theta=0.2 #probability of recovering a given fish or sampling fraction
tag_rate=0.25 #fixed CFM rate
iterations=1000 #k and scale iterations
scales_n_seq =seq(400,2000,100)
boot_replicates=100 #times to repeat total simulation
CI=0.95
h_prop=.75

sim_results<-list()
for(i in 1:boot_replicates){
  iter_start<-Sys.time()
  #simulate population
  sim1<-sim_pop(N=N,
                tag_rate=tag_rate,
                h_prop = h_prop)
  
  spawning_recoveries<-sim_spawning_recovery(sim1,
                                             theta=theta,
                                             tag_rate=tag_rate,
                                             iterations = iterations)
  scale_iters<-list()
  for(s in 1:length(scales_n_seq)){
    scale_results<-sim_scales_sampling(pop=sim1,
                                       spawning_results=spawning_recoveries,
                                       scales_n=scales_n_seq[s],
                                       iterations = iterations)
    sim_stats<-sim_summary_stats(scale_results)
    scale_iters[[s]]=list(
      "scale_results"=scale_results,
      "sim_stats"=sim_stats
    )
  }
  names(scale_iters)=scales_n_seq
  
  sim_results[[i]]=scale_iters
  iter_end<-Sys.time()
  iter_time<-iter_end-iter_start
  print(paste("Scale scenarios, iteration ",i," finished, time: ",round(iter_time,4)))
}
saveRDS(sim_results,"outputs/scale_sim_results_N100K.Rds")

##########################
#theta size run
##########################
#simulate population
N=10000 #true population size
theta_list=seq(0.1,0.9,.1) #probability of recovering a given fish or sampling fraction
tag_rate=0.25 #fixed CFM rate
iterations=1000 #k and scale iterations
scales_n_seq =500
boot_replicates=100 #times to repeat total simulation
CI=0.95
h_prop=.75

theta_results<-list()
for(i in 1:boot_replicates){
  iter_start<-Sys.time()
  #simulate population
  sim1<-sim_pop(N=N,
                tag_rate=tag_rate,
                h_prop = h_prop)
  
  theta_iters<-list()
  for(t in 1:length(theta_list)){
    spawning_recoveries<-sim_spawning_recovery(sim1,
                                               theta=theta_list[t],
                                               tag_rate=tag_rate,
                                               iterations = iterations)
    
    scale_results<-sim_scales_sampling(pop=sim1,
                                       spawning_results=spawning_recoveries,
                                       scales_n=scales_n_seq[1],
                                       iterations = iterations)
    sim_stats<-sim_summary_stats(scale_results)
    theta_iters[[t]]=list(
      "scale_results"=scale_results,
      "sim_stats"=sim_stats
    )
  }
  names(theta_iters)=theta_list
  
  theta_results[[i]]=theta_iters
  iter_end<-Sys.time()
  iter_time<-iter_end-iter_start
  print(paste("Theta scenarios, iteration ",i," finished, time: ",round(iter_time,4)))
}
saveRDS(theta_results,"outputs/theta_sim_results.Rds")
test<-readRDS("outputs/theta_sim_results.Rds")
##########################
#false negatives run
##########################
#simulate population
N=10000 #true population size
theta=0.2 #probability of recovering a given fish or sampling fraction
tag_rate=0.25 #fixed CFM rate
iterations=1000 #k and scale iterations
scales_n_seq =seq(400,2000,100)
boot_replicates=100 #times to repeat total simulation
CI=0.95
h_prop=.95

fnegative_results<-list()
for(i in 1:boot_replicates){
  iter_start<-Sys.time()
  #simulate population
  sim1<-sim_pop(N=N,
                tag_rate=tag_rate,
                h_prop = h_prop)
  
  spawning_recoveries<-sim_spawning_recovery(sim1,
                                             theta=theta,
                                             tag_rate=tag_rate,
                                             iterations = iterations)
  scale_iters<-list()
  for(s in 1:length(scales_n_seq)){
    scale_results<-sim_scales_sampling(pop=sim1,
                                       spawning_results=spawning_recoveries,
                                       scales_n=scales_n_seq[s],
                                       iterations = iterations)
    sim_stats<-sim_summary_stats(scale_results)
    scale_iters[[s]]=list(
      "scale_results"=scale_results,
      "sim_stats"=sim_stats
    )
  }
  names(scale_iters)=scales_n_seq
  
  fnegative_results[[i]]=scale_iters
  iter_end<-Sys.time()
  iter_time<-iter_end-iter_start
  print(paste("False negatives, iteration ",i," finished, time: ",round(iter_time,4)))
}
saveRDS(fnegative_results,"outputs/fnegative_sim_results.Rds")

##########################
#false positive run
##########################
#simulate population
N=10000 #true population size
theta=0.2 #probability of recovering a given fish or sampling fraction
tag_rate=0.25 #fixed CFM rate
iterations=1000 #k and scale iterations
scales_n_seq =seq(400,2000,100)
boot_replicates=100 #times to repeat total simulation
CI=0.95
h_prop=1

fpositive_results<-list()
for(i in 1:boot_replicates){
  iter_start<-Sys.time()
  #simulate population
  sim1<-sim_pop(N=N,
                tag_rate=tag_rate,
                h_prop = h_prop)
  
  spawning_recoveries<-sim_spawning_recovery(sim1,
                                             theta=theta,
                                             tag_rate=tag_rate,
                                             iterations = iterations)
  scale_iters<-list()
  for(s in 1:length(scales_n_seq)){
    scale_results<-sim_scales_sampling(pop=sim1,
                                       spawning_results=spawning_recoveries,
                                       scales_n=scales_n_seq[s],
                                       iterations = iterations)
    sim_stats<-sim_summary_stats(scale_results)
    scale_iters[[s]]=list(
      "scale_results"=scale_results,
      "sim_stats"=sim_stats
    )
  }
  names(scale_iters)=scales_n_seq
  
  fpositive_results[[i]]=scale_iters
  iter_end<-Sys.time()
  iter_time<-iter_end-iter_start
  print(paste("False positives, iteration ",i," finished, time: ",round(iter_time,4)))
}
saveRDS(fpositive_results,"outputs/fpositive_sim_results.Rds")

