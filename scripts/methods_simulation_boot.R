rm( list = ls()) #clear env
#basin methods simulation
#this script aims to simulate basin populations and replicate methods
# to produce a basin-wide cohort reconstruction
library(tidyverse,quietly = "true")
library(data.table)
library(ggplot2)
library(gtools)
library(crayon)

options(dplyr.summarise.inform = FALSE)  #suppress annoying dplyr messages

#Load functions
sapply(list.files("scripts/functions", pattern = "\\.R$", full.names = TRUE), source)

set.seed(67) #for reproducibility

#read in data
trib_data<-readRDS("data/all_returns_clean.Rds")
theta_data<-readRDS("data/location_year_theta.Rds")

################################
#regular pooling
###############################

year_list<-2010:2020
boot_replicates<-100
iterations=1000
boot_scales<-seq(1000,10000,1000)
iter_results<-list()
pooling_iteration_results<-list()

for(i in 1:boot_replicates){
  iter_start<-Sys.time()
  ###################################
  #simulate population
  ###################################
  iter_year=sample(year_list,1)
  #generate basin-wide population
  #now updated with brood_year and brood_location
  basin_pop<-sim_basin_pop(trib_data,i,iter_year)$iter_pop
  
  #grab total basin pop
  basin_N=nrow(basin_pop)
  
  loc_list<-unique(basin_pop$location)
  
  ###################################
  #Methods Step 1: Escapement Surveys
  #spawning recoveries for each location
  #this should be the data collected and prepared by survey teams
  ###################################
  loc_recoveries<-list()
  for(l in 1:length(loc_list)){
    iter_loc<-loc_list[[l]]
    loc_pop<-basin_pop%>%filter(location==iter_loc) #get pop data for location
    
    target_theta<-theta_data%>% #get theta for location x year
      filter(location==iter_loc,
             return_year==iter_year)
    if(nrow(target_theta)==0){target_theta=0.2} else{ #if no theta, set default
      target_theta<-round(target_theta$theta,2)
    }
    
    recovered_fish<-sim_spawning_recovery(loc_pop,target_theta)
    
    loc_recoveries[[l]]<-recovered_fish
  }
  names(loc_recoveries)=loc_list
  
  #for each location produce the required data that a coordinator would need
  loc_data<-list()
  for(l in 1:length(loc_list)){
    iter_loc<-loc_list[[l]]
    d<-loc_recoveries[[l]]
    N_loc<-basin_pop%>% #estimate N_loc escapement estimate (this would be done by the survey group and have uncertainty in it)
      filter(location==iter_loc)%>%
      nrow()
    n_nocwt=d%>%filter(has_cwt==0) #number of recovered fish without cwt
    n_hascwt=d%>%filter(has_cwt==1) #number of recovered fish with cwt
    
    cwt_uid_freq<-n_hascwt%>% #for each cwt group, provide frequency
      group_by(location,release_group,tag_rate,return_year,
               brood_year,tag_age,brood_location)%>%
      summarise(frequency=n())
    theta=nrow(d)/N_loc #calculate theta
    
    ans<-list("N_loc"=N_loc,
              "n_recovered"=nrow(d),
              "theta"=theta,
              "F_prod"=cwt_uid_freq)
    loc_data[[l]]=ans
  }
  names(loc_data)=loc_list
  
  ###################################
  #Methods Step 2: Use k-draws to estimate N_nocwt
  #this process could be done by either a coordinator or surveys themselves
  ###################################
  #for each location run iterated k-draws to estimate N_nocwt and uncertainty
  loc_kdraws<-list()
  loc_results<-list()
  survey_year=iter_year #set survey year for surveys with no cwt recovered
  for(l in 1:length(loc_list)){
    d<-loc_data[[l]]
    d_fprod<-d$F_prod
    
    uid_list<-unique(d_fprod$release_group)
    uid_results<-data.frame()
    if(length(uid_list)>1){
      for(u in 1:length(uid_list)){
        tag_data<-d_fprod%>%filter(release_group==uid_list[u])
        tag_outputs<-tag_k_draws(tag_data=tag_data,
                                 theta=d$theta,
                                 iterations=iterations,
                                 freq=tag_data$frequency,
                                 f_prod=tag_data$tag_rate)
        
        uid_results<-uid_results%>%rbind(tag_outputs)
      }
    }else{
      uid_results <- data.table::data.table(
        cwt = rep(NA,iterations),
        brood_year=rep(NA,iterations),
        tag_age=rep(NA,iterations),
        return_year=rep(survey_year,iterations),
        location=rep(loc_list[l],iterations),
        theta=rep(d$theta,iterations),
        tag_rate=rep(NA,iterations),
        tag_freq =rep(0,iterations),
        total_tags = rep(0,iterations),
        total_hatchery = rep(0,iterations),
        k_iteration = rep(0,iterations)
      )
    }
    
    loc_kdraws[[l]]<-uid_results
    
    location_estimates<-uid_results%>%
      group_by(k_iteration,location,return_year,theta)%>%
      summarise(recovered_tags=sum(tag_freq),
                mean_tag_rate=round(mean(tag_rate),2),
                total_tags=sum(total_tags),
                total_hatchery=sum(total_hatchery))
    location_estimates<-location_estimates%>%
      group_by(location,return_year,mean_tag_rate,
               recovered_tags)%>%
      summarise(total_tags=round(mean(total_tags)),
                total_hatchery=round(mean(total_hatchery)))
    location_estimates<-location_estimates%>%
      mutate(N_loc=d$N_loc)%>% #add location total escapement
      mutate(N_nocwt=N_loc-total_tags) #estimate N untagged
    
    loc_results<-loc_results%>%rbind(location_estimates)
  }
  
  names(loc_kdraws)=loc_list #save k-draws for later use in uncertainty
  
  scale_results<-list()
  
  for(s in 1:length(boot_scales)){
    
    scale_boot_start<-Sys.time()
    target_scales=boot_scales[s]
    
    #estimate N_nocwt proportion from each survey
    basin_nocwt<-sum(loc_results$N_nocwt)
    loc_results<-loc_results%>%
      mutate(nocwt_prop=N_nocwt/basin_nocwt)
    
    #use round_to_sum function to properly round scales allotment
    loc_results$target_scales=round_to_sum(loc_results$nocwt_prop,target_scales)
    
    #get number of nocwt scales available at each location
    scales_available<-data.frame()
    for(l in loc_list){
      d<-loc_recoveries[[l]]
      d<-d%>%
        filter(has_cwt==0)
      available=data.frame(location=l,
                           nocwt_scales_available=nrow(d))
      scales_available<-scales_available%>%rbind(available)
    }
    
    loc_results<-loc_results%>%
      left_join(scales_available)
    
    #if we are short on available scales, deal with that
    loc_results<-loc_results%>%
      mutate(
        nocwt_scales_available=ifelse(is.na(nocwt_scales_available),0
                                      ,nocwt_scales_available),
        #set actual_scales to either target_scales or available, which ever is smallest
        actual_scales=pmin(target_scales,nocwt_scales_available),
        shortfall=target_scales-actual_scales
      )
    
    total_shortfall=sum(loc_results$shortfall)
    
    ###################################
    #Methods Step 3: Sample scales from locations
    #this is done randomly by the survey, but simulated here
    ###################################
    
    loc_scales<-vector(mode = "list", length = length(loc_list))
    names(loc_scales)<-loc_list
    for(l in loc_list){
      d<-loc_recoveries[[l]]
      scales_n=loc_results%>%
        filter(location==l)
      d<-d%>%
        filter(has_cwt==0)
      d_sample<-sample_n(d,scales_n$actual_scales,replace = F)
      loc_scales[[l]]=d_sample
    }
    
    basin_scales<-data.frame()
    for(l in loc_list){
      basin_scales<-basin_scales%>%rbind(loc_scales[[l]])
    }
    
    #now we assume scales are sent in for aging, and we will just get the results as basin_scales
    
    ###################################
    #Methods Step 4: get basin hachery estimates
    #done by coordinator
    ###################################
    
    #first take loc_kdraws and get basin wide estimate of hatchery origin fish for each iteration
    loc_hatcheries<-data.frame()
    for(l in loc_list){
      d<-loc_kdraws[[l]]
      d<-d%>%
        group_by(return_year,tag_age,location,k_iteration,theta)%>%
        summarise(tag_freq=sum(tag_freq),
                  total_tags=sum(total_tags),
                  total_hatchery=sum(total_hatchery))
      loc_hatcheries<-loc_hatcheries%>%rbind(d)
    }
    
    #sometimes locations have no hatchery fish so remove here
    loc_hatcheries<-loc_hatcheries%>%filter(!is.na(tag_age))
    
    #sum location based hatchery estimates from k-draws into basin-wide
    basin_hatchery_estimates<-loc_hatcheries%>%
      group_by(return_year,tag_age,k_iteration)%>%
      summarise(tag_freq = sum(tag_freq),
                total_tags = sum(total_tags),
                total_hatchery = sum(total_hatchery))
    
    spawning_results=list("basin_scales"=basin_scales, #irl this will be our aging results
                          "hatchery_ages"=basin_hatchery_estimates)
    
    cr_results<-basin_CR(basin_N=basin_N,
                         spawning_results = spawning_results,
                         iterations = iterations,
                         scales_n = target_scales,
                         CI=0.9)
    
    ###################################
    #compare results to true values
    ###################################
    
    true_values<-basin_pop%>%
      group_by(tag_age,origin)%>%
      summarize(abundance=n())%>%
      pivot_wider(names_from=origin,
                  values_from=abundance)
    
    compare<-cr_results%>%left_join(true_values)
    compare<-compare%>%
      mutate(count_accuracy=abs(mean_natural-natural)/natural)
    
    #save outputs
    scale_results[[s]]<-list(
      "comparison"=compare,
      "loc_results"=loc_results
    )
    
    scale_boot_end<-Sys.time()
    cat(red(paste("\nscale boot",boot_scales[s],"DONE, time:",round(scale_boot_end-scale_boot_start,2),"\n")))
  }
  
  names(scale_results)=boot_scales
  pooling_iteration_results[[i]]<-scale_results
  iter_end<-Sys.time()
  cat(red(paste("\niteration",i,"DONE, time:",round(iter_end-iter_start,2),"\n")))
}
saveRDS(pooling_iteration_results,"outputs/methods_comparisons.Rds")

pooling_iteration_results<-readRDS("outputs/methods_comparisons.Rds")
#plotting

target_accuracy=.150
pooling_estimates<-data.frame()
for(i in 1:length(pooling_iteration_results)){
  d<-pooling_iteration_results[[i]]
  for(s in 1:length(d)){
    ds<-d[[s]]$comparison
    ds$n_scales<-as.numeric(names(d[s]))
    
    pooling_estimates<-pooling_estimates%>%rbind(ds)
  }
}

stats<-pooling_estimates%>%
  group_by(age,n_scales)%>%
  summarise(
    pct_meeting_accuracy=mean(meets_target_accuracy)*100
  )

p1<-ggplot(data=stats,
           aes(x=(n_scales),y=pct_meeting_accuracy))+
  geom_point( size=3, aes(color=factor(age)))+
  geom_hline(yintercept=90,color="red")+
  scale_x_continuous(limits=c(1000,10000),
                     breaks=seq(from=1000,to=10000,by=1000))+
  theme_bw()+
  labs(x = "basin scale samples", y = "% of estimates that were within +/- 15% of true value",color="Age",linetype="Age")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none")+
  facet_grid(factor(age)~.)
p1
ggsave(p1,file="outputs/count_accuracy.png",width=10,height=6)



################################
#regular pooling low samples 100-1000
###############################

year_list<-2010:2020
boot_replicates<-100
iterations=100
boot_scales<-seq(100,1000,100)
iter_results<-list()
pooling_iteration_results<-list()

for(i in 1:boot_replicates){
  iter_start<-Sys.time()
  ###################################
  #simulate population
  ###################################
  iter_year=sample(year_list,1)
  #generate basin-wide population
  #now updated with brood_year and brood_location
  basin_pop<-sim_basin_pop(trib_data,i,iter_year)$iter_pop
  
  #grab total basin pop
  basin_N=nrow(basin_pop)
  
  loc_list<-unique(basin_pop$location)
  
  ###################################
  #Methods Step 1: Escapement Surveys
  #spawning recoveries for each location
  #this should be the data collected and prepared by survey teams
  ###################################
  loc_recoveries<-list()
  for(l in 1:length(loc_list)){
    iter_loc<-loc_list[[l]]
    loc_pop<-basin_pop%>%filter(location==iter_loc) #get pop data for location
    
    target_theta<-theta_data%>% #get theta for location x year
      filter(location==iter_loc,
             return_year==iter_year)
    if(nrow(target_theta)==0){target_theta=0.2} else{ #if no theta, set default
      target_theta<-round(target_theta$theta,2)
    }
    
    recovered_fish<-sim_spawning_recovery(loc_pop,target_theta)
    
    loc_recoveries[[l]]<-recovered_fish
  }
  names(loc_recoveries)=loc_list
  
  #for each location produce the required data that a coordinator would need
  loc_data<-list()
  for(l in 1:length(loc_list)){
    iter_loc<-loc_list[[l]]
    d<-loc_recoveries[[l]]
    N_loc<-basin_pop%>% #estimate N_loc escapement estimate (this would be done by the survey group and have uncertainty in it)
      filter(location==iter_loc)%>%
      nrow()
    n_nocwt=d%>%filter(has_cwt==0) #number of recovered fish without cwt
    n_hascwt=d%>%filter(has_cwt==1) #number of recovered fish with cwt
    
    cwt_uid_freq<-n_hascwt%>% #for each cwt group, provide frequency
      group_by(location,release_group,tag_rate,return_year,
               brood_year,tag_age,brood_location)%>%
      summarise(frequency=n())
    theta=nrow(d)/N_loc #calculate theta
    
    ans<-list("N_loc"=N_loc,
              "n_recovered"=nrow(d),
              "theta"=theta,
              "F_prod"=cwt_uid_freq)
    loc_data[[l]]=ans
  }
  names(loc_data)=loc_list
  
  ###################################
  #Methods Step 2: Use k-draws to estimate N_nocwt
  #this process could be done by either a coordinator or surveys themselves
  ###################################
  #for each location run iterated k-draws to estimate N_nocwt and uncertainty
  loc_kdraws<-list()
  loc_results<-list()
  survey_year=iter_year #set survey year for surveys with no cwt recovered
  for(l in 1:length(loc_list)){
    d<-loc_data[[l]]
    d_fprod<-d$F_prod
    
    uid_list<-unique(d_fprod$release_group)
    uid_results<-data.frame()
    if(length(uid_list)>1){
      for(u in 1:length(uid_list)){
        tag_data<-d_fprod%>%filter(release_group==uid_list[u])
        tag_outputs<-tag_k_draws(tag_data=tag_data,
                                 theta=d$theta,
                                 iterations=iterations,
                                 freq=tag_data$frequency,
                                 f_prod=tag_data$tag_rate)
        
        uid_results<-uid_results%>%rbind(tag_outputs)
      }
    }else{
      uid_results <- data.table::data.table(
        cwt = rep(NA,iterations),
        brood_year=rep(NA,iterations),
        tag_age=rep(NA,iterations),
        return_year=rep(survey_year,iterations),
        location=rep(loc_list[l],iterations),
        theta=rep(d$theta,iterations),
        tag_rate=rep(NA,iterations),
        tag_freq =rep(0,iterations),
        total_tags = rep(0,iterations),
        total_hatchery = rep(0,iterations),
        k_iteration = rep(0,iterations)
      )
    }
    
    loc_kdraws[[l]]<-uid_results
    
    location_estimates<-uid_results%>%
      group_by(k_iteration,location,return_year,theta)%>%
      summarise(recovered_tags=sum(tag_freq),
                mean_tag_rate=round(mean(tag_rate),2),
                total_tags=sum(total_tags),
                total_hatchery=sum(total_hatchery))
    location_estimates<-location_estimates%>%
      group_by(location,return_year,mean_tag_rate,
               recovered_tags)%>%
      summarise(total_tags=round(mean(total_tags)),
                total_hatchery=round(mean(total_hatchery)))
    location_estimates<-location_estimates%>%
      mutate(N_loc=d$N_loc)%>% #add location total escapement
      mutate(N_nocwt=N_loc-total_tags) #estimate N untagged
    
    loc_results<-loc_results%>%rbind(location_estimates)
  }
  
  names(loc_kdraws)=loc_list #save k-draws for later use in uncertainty
  
  scale_results<-list()
  
  for(s in 1:length(boot_scales)){
    
    scale_boot_start<-Sys.time()
    target_scales=boot_scales[s]
    
    #estimate N_nocwt proportion from each survey
    basin_nocwt<-sum(loc_results$N_nocwt)
    loc_results<-loc_results%>%
      mutate(nocwt_prop=N_nocwt/basin_nocwt)
    
    #use round_to_sum function to properly round scales allotment
    loc_results$target_scales=round_to_sum(loc_results$nocwt_prop,target_scales)
    
    #get number of nocwt scales available at each location
    scales_available<-data.frame()
    for(l in loc_list){
      d<-loc_recoveries[[l]]
      d<-d%>%
        filter(has_cwt==0)
      available=data.frame(location=l,
                           nocwt_scales_available=nrow(d))
      scales_available<-scales_available%>%rbind(available)
    }
    
    loc_results<-loc_results%>%
      left_join(scales_available)
    
    #if we are short on available scales, deal with that
    loc_results<-loc_results%>%
      mutate(
        nocwt_scales_available=ifelse(is.na(nocwt_scales_available),0
                                      ,nocwt_scales_available),
        #set actual_scales to either target_scales or available, which ever is smallest
        actual_scales=pmin(target_scales,nocwt_scales_available),
        shortfall=target_scales-actual_scales
      )
    
    total_shortfall=sum(loc_results$shortfall)
    
    ###################################
    #Methods Step 3: Sample scales from locations
    #this is done randomly by the survey, but simulated here
    ###################################
    
    loc_scales<-vector(mode = "list", length = length(loc_list))
    names(loc_scales)<-loc_list
    for(l in loc_list){
      d<-loc_recoveries[[l]]
      scales_n=loc_results%>%
        filter(location==l)
      d<-d%>%
        filter(has_cwt==0)
      d_sample<-sample_n(d,scales_n$actual_scales,replace = F)
      loc_scales[[l]]=d_sample
    }
    
    basin_scales<-data.frame()
    for(l in loc_list){
      basin_scales<-basin_scales%>%rbind(loc_scales[[l]])
    }
    
    #now we assume scales are sent in for aging, and we will just get the results as basin_scales
    
    ###################################
    #Methods Step 4: get basin hachery estimates
    #done by coordinator
    ###################################
    
    #first take loc_kdraws and get basin wide estimate of hatchery origin fish for each iteration
    loc_hatcheries<-data.frame()
    for(l in loc_list){
      d<-loc_kdraws[[l]]
      d<-d%>%
        group_by(return_year,tag_age,location,k_iteration,theta)%>%
        summarise(tag_freq=sum(tag_freq),
                  total_tags=sum(total_tags),
                  total_hatchery=sum(total_hatchery))
      loc_hatcheries<-loc_hatcheries%>%rbind(d)
    }
    
    #sometimes locations have no hatchery fish so remove here
    loc_hatcheries<-loc_hatcheries%>%filter(!is.na(tag_age))
    
    #sum location based hatchery estimates from k-draws into basin-wide
    basin_hatchery_estimates<-loc_hatcheries%>%
      group_by(return_year,tag_age,k_iteration)%>%
      summarise(tag_freq = sum(tag_freq),
                total_tags = sum(total_tags),
                total_hatchery = sum(total_hatchery))
    
    spawning_results=list("basin_scales"=basin_scales, #irl this will be our aging results
                          "hatchery_ages"=basin_hatchery_estimates)
    
    cr_results<-basin_CR(basin_N=basin_N,
                         spawning_results = spawning_results,
                         iterations = iterations,
                         scales_n = target_scales,
                         CI=0.9)
    
    ###################################
    #compare results to true values
    ###################################
    
    true_values<-basin_pop%>%
      group_by(tag_age,origin)%>%
      summarize(abundance=n())%>%
      pivot_wider(names_from=origin,
                  values_from=abundance)
    
    compare<-cr_results%>%left_join(true_values)
    compare<-compare%>%
      mutate(count_accuracy=abs(mean_natural-natural)/natural)
    
    #save outputs
    scale_results[[s]]<-list(
      "comparison"=compare,
      "loc_results"=loc_results
    )
    
    scale_boot_end<-Sys.time()
    cat(red(paste("\nscale boot",boot_scales[s],"DONE, time:",round(scale_boot_end-scale_boot_start,2),"\n")))
  }
  
  names(scale_results)=boot_scales
  pooling_iteration_results[[i]]<-scale_results
  iter_end<-Sys.time()
  cat(red(paste("\niteration",i,"DONE, time:",round(iter_end-iter_start,2),"\n")))
}
saveRDS(pooling_iteration_results,"outputs/methods_comparisons_lowsamples.Rds")

pooling_iteration_results<-readRDS("outputs/methods_comparisons_lowsamples.Rds")
#plotting

target_accuracy=.150
pooling_estimates<-data.frame()
for(i in 1:length(pooling_iteration_results)){
  d<-pooling_iteration_results[[i]]
  for(s in 1:length(d)){
    ds<-d[[s]]$comparison
    ds$n_scales<-as.numeric(names(d[s]))
    ds$age=ds$tag_age
    pooling_estimates<-pooling_estimates%>%rbind(ds)
  }
}

stats<-pooling_estimates%>%
  mutate(meets_target_accuracy=ifelse(count_accuracy<=target_accuracy,T,F))%>%
  group_by(age,n_scales)%>%
  summarise(
    pct_meeting_accuracy=mean(meets_target_accuracy)*100
  )

p1<-ggplot(data=stats,
           aes(x=(n_scales),y=pct_meeting_accuracy))+
  geom_point( size=3, aes(color=factor(age)))+
  geom_hline(yintercept=90,color="red")+
  scale_x_continuous(limits=c(100,1000),
                     breaks=seq(from=100,to=1000,by=100))+
  theme_bw()+
  labs(x = "basin scale samples", y = "% of estimates that were within +/- 15% of true value",color="Age",linetype="Age")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none")+
  facet_grid(factor(age)~.)
p1
ggsave(p1,file="outputs/count_accuracy_lowsamples.png",width=10,height=6)

################################
#regular pooling impact of theta
###############################

year_list<-2010:2020
boot_replicates<-100
iterations=100
boot_scales<-500
theta_seq<-seq(0.1,1,.1)
iter_results<-list()
pooling_iteration_results<-list()

for(i in 1:boot_replicates){
  theta_results<-list()
  for(t in 1:length(theta_seq)){
    t_inc<-theta_seq[[t]]
    iter_start<-Sys.time()
    ###################################
    #simulate population
    ###################################
    iter_year=sample(year_list,1)
    #generate basin-wide population
    #now updated with brood_year and brood_location
    basin_pop<-sim_basin_pop(trib_data,i,iter_year)$iter_pop
    
    #grab total basin pop
    basin_N=nrow(basin_pop)
    
    loc_list<-unique(basin_pop$location)
    
    ###################################
    #Methods Step 1: Escapement Surveys
    #spawning recoveries for each location
    #this should be the data collected and prepared by survey teams
    ###################################
    loc_recoveries<-list()
    for(l in 1:length(loc_list)){
      iter_loc<-loc_list[[l]]
      loc_pop<-basin_pop%>%filter(location==iter_loc) #get pop data for location
      #target_theta<-theta_data%>% #get theta for location x year
      #  filter(location==iter_loc,
      #         return_year==iter_year)
      #target_theta<-target_theta%>%
      #  mutate(theta=ifelse(theta>=1,theta,theta+t_inc))
      #target_theta<-target_theta%>%
      #  mutate(theta=ifelse(theta>1,1,theta))
      #if(nrow(target_theta)==0){target_theta=0.2+t_inc} else{ #if no theta, set default
      #  target_theta<-round(target_theta$theta,2)
      #}
      #hatcheries<-c("NFH","FRH","Battle")
      target_theta<-t_inc
      
      
      recovered_fish<-sim_spawning_recovery(loc_pop,target_theta)
      
      loc_recoveries[[l]]<-recovered_fish
    }
    names(loc_recoveries)=loc_list
    
    #for each location produce the required data that a coordinator would need
    loc_data<-list()
    for(l in 1:length(loc_list)){
      iter_loc<-loc_list[[l]]
      d<-loc_recoveries[[l]]
      N_loc<-basin_pop%>% #estimate N_loc escapement estimate (this would be done by the survey group and have uncertainty in it)
        filter(location==iter_loc)%>%
        nrow()
      n_nocwt=d%>%filter(has_cwt==0) #number of recovered fish without cwt
      n_hascwt=d%>%filter(has_cwt==1) #number of recovered fish with cwt
      
      cwt_uid_freq<-n_hascwt%>% #for each cwt group, provide frequency
        group_by(location,release_group,tag_rate,return_year,
                 brood_year,tag_age,brood_location)%>%
        summarise(frequency=n())
      theta=nrow(d)/N_loc #calculate theta
      
      ans<-list("N_loc"=N_loc,
                "n_recovered"=nrow(d),
                "theta"=theta,
                "F_prod"=cwt_uid_freq)
      loc_data[[l]]=ans
    }
    names(loc_data)=loc_list
    
    ###################################
    #Methods Step 2: Use k-draws to estimate N_nocwt
    #this process could be done by either a coordinator or surveys themselves
    ###################################
    #for each location run iterated k-draws to estimate N_nocwt and uncertainty
    loc_kdraws<-list()
    loc_results<-list()
    survey_year=iter_year #set survey year for surveys with no cwt recovered
    for(l in 1:length(loc_list)){
      d<-loc_data[[l]]
      d_fprod<-d$F_prod
      
      uid_list<-unique(d_fprod$release_group)
      uid_results<-data.frame()
      if(length(uid_list)>1){
        for(u in 1:length(uid_list)){
          tag_data<-d_fprod%>%filter(release_group==uid_list[u])
          tag_outputs<-tag_k_draws(tag_data=tag_data,
                                   theta=d$theta,
                                   iterations=iterations,
                                   freq=tag_data$frequency,
                                   f_prod=tag_data$tag_rate)
          
          uid_results<-uid_results%>%rbind(tag_outputs)
        }
      }else{
        uid_results <- data.table::data.table(
          cwt = rep(NA,iterations),
          brood_year=rep(NA,iterations),
          tag_age=rep(NA,iterations),
          return_year=rep(survey_year,iterations),
          location=rep(loc_list[l],iterations),
          theta=rep(d$theta,iterations),
          tag_rate=rep(NA,iterations),
          tag_freq =rep(0,iterations),
          total_tags = rep(0,iterations),
          total_hatchery = rep(0,iterations),
          k_iteration = rep(0,iterations)
        )
      }
      
      loc_kdraws[[l]]<-uid_results
      
      location_estimates<-uid_results%>%
        group_by(k_iteration,location,return_year,theta)%>%
        summarise(recovered_tags=sum(tag_freq),
                  mean_tag_rate=round(mean(tag_rate),2),
                  total_tags=sum(total_tags),
                  total_hatchery=sum(total_hatchery))
      location_estimates<-location_estimates%>%
        group_by(location,return_year,mean_tag_rate,
                 recovered_tags)%>%
        summarise(total_tags=round(mean(total_tags)),
                  total_hatchery=round(mean(total_hatchery)))
      location_estimates<-location_estimates%>%
        mutate(N_loc=d$N_loc)%>% #add location total escapement
        mutate(N_nocwt=N_loc-total_tags) #estimate N untagged
      
      loc_results<-loc_results%>%rbind(location_estimates)
    }
    
    names(loc_kdraws)=loc_list #save k-draws for later use in uncertainty
    
    scale_results<-list()
    
    for(s in 1:length(boot_scales)){
      
      scale_boot_start<-Sys.time()
      target_scales=boot_scales[s]
      
      #estimate N_nocwt proportion from each survey
      basin_nocwt<-sum(loc_results$N_nocwt)
      loc_results<-loc_results%>%
        mutate(nocwt_prop=N_nocwt/basin_nocwt)
      
      #use round_to_sum function to properly round scales allotment
      loc_results$target_scales=round_to_sum(loc_results$nocwt_prop,target_scales)
      
      #get number of nocwt scales available at each location
      scales_available<-data.frame()
      for(l in loc_list){
        d<-loc_recoveries[[l]]
        d<-d%>%
          filter(has_cwt==0)
        available=data.frame(location=l,
                             nocwt_scales_available=nrow(d))
        scales_available<-scales_available%>%rbind(available)
      }
      
      loc_results<-loc_results%>%
        left_join(scales_available)
      
      #if we are short on available scales, deal with that
      loc_results<-loc_results%>%
        mutate(
          nocwt_scales_available=ifelse(is.na(nocwt_scales_available),0
                                        ,nocwt_scales_available),
          #set actual_scales to either target_scales or available, which ever is smallest
          actual_scales=pmin(target_scales,nocwt_scales_available),
          shortfall=target_scales-actual_scales
        )
      
      total_shortfall=sum(loc_results$shortfall)
      
      ###################################
      #Methods Step 3: Sample scales from locations
      #this is done randomly by the survey, but simulated here
      ###################################
      
      loc_scales<-vector(mode = "list", length = length(loc_list))
      names(loc_scales)<-loc_list
      for(l in loc_list){
        d<-loc_recoveries[[l]]
        scales_n=loc_results%>%
          filter(location==l)
        d<-d%>%
          filter(has_cwt==0)
        d_sample<-sample_n(d,scales_n$actual_scales,replace = F)
        loc_scales[[l]]=d_sample
      }
      
      basin_scales<-data.frame()
      for(l in loc_list){
        basin_scales<-basin_scales%>%rbind(loc_scales[[l]])
      }
      
      #now we assume scales are sent in for aging, and we will just get the results as basin_scales
      
      ###################################
      #Methods Step 4: get basin hachery estimates
      #done by coordinator
      ###################################
      
      #first take loc_kdraws and get basin wide estimate of hatchery origin fish for each iteration
      loc_hatcheries<-data.frame()
      for(l in loc_list){
        d<-loc_kdraws[[l]]
        d<-d%>%
          group_by(return_year,tag_age,location,k_iteration,theta)%>%
          summarise(tag_freq=sum(tag_freq),
                    total_tags=sum(total_tags),
                    total_hatchery=sum(total_hatchery))
        loc_hatcheries<-loc_hatcheries%>%rbind(d)
      }
      
      #sometimes locations have no hatchery fish so remove here
      loc_hatcheries<-loc_hatcheries%>%filter(!is.na(tag_age))
      
      #sum location based hatchery estimates from k-draws into basin-wide
      basin_hatchery_estimates<-loc_hatcheries%>%
        group_by(return_year,tag_age,k_iteration)%>%
        summarise(tag_freq = sum(tag_freq),
                  total_tags = sum(total_tags),
                  total_hatchery = sum(total_hatchery))
      
      basin_scales<-basin_scales%>%
        mutate(scale_age=return_year-brood_year)
      
      spawning_results=list("basin_scales"=basin_scales, #irl this will be our aging results
                            "hatchery_ages"=basin_hatchery_estimates)
      
      cr_results<-basin_CR(basin_N=basin_N,
                           spawning_results = spawning_results,
                           iterations = iterations,
                           CI=0.9)
      
      ###################################
      #compare results to true values
      ###################################
      
      true_values<-basin_pop%>%
        mutate(age=return_year-brood_year)%>%
        group_by(age,origin)%>%
        summarize(true_val=n())
      
      compare<-cr_results%>%left_join(true_values)
      compare<-compare%>%
        mutate(count_accuracy=abs(estimate-true_val)/true_val)
      
      compare$theta=t_inc
      
      #save outputs
      scale_results[[s]]<-list(
        "comparison"=compare,
        "loc_results"=loc_results
      )
      names(scale_results)=boot_scales
      scale_boot_end<-Sys.time()
      cat(red(paste("\ntheta boot",t,"DONE, time:",round(scale_boot_end-scale_boot_start,2),"\n")))
    }
    
    theta_results[[t]]<-scale_results
    
  }
  names(theta_results)<-theta_seq
  pooling_iteration_results[[i]]<-theta_results
  iter_end<-Sys.time()
  cat(red(paste("\niteration",i,"DONE, time:",round(iter_end-iter_start,2),"\n")))
}
saveRDS(pooling_iteration_results,"outputs/methods_comparisons_inc_theta.Rds")

################################
#regular pooling impact of tag_rate
###############################

year_list<-2010:2020
boot_replicates<-100
iterations=100
boot_scales<-500
tag_rate<-seq(0.1,1,.1)
iter_results<-list()
pooling_iteration_results<-list()

for(i in 1:boot_replicates){
  tag_results<-list()
  for(t in 1:length(tag_rate)){
    t_inc<-tag_rate[[t]]
    iter_start<-Sys.time()
    ###################################
    #simulate population
    ###################################
    iter_year=sample(year_list,1)
    #generate basin-wide population
    #now updated with brood_year and brood_location
    trib_data<-trib_data%>%
      mutate(tag_rate=ifelse(origin=="hatchery",t_inc,NA))
    basin_pop<-sim_basin_pop(trib_data,i,iter_year)$iter_pop
    
    basin_pop<-basin_pop%>%
      mutate(tag_rate=ifelse(origin=="hatchery",t_inc,NA))
    
    #grab total basin pop
    basin_N=nrow(basin_pop)
    
    loc_list<-unique(basin_pop$location)
    
    ###################################
    #Methods Step 1: Escapement Surveys
    #spawning recoveries for each location
    #this should be the data collected and prepared by survey teams
    ###################################
    loc_recoveries<-list()
    for(l in 1:length(loc_list)){
      iter_loc<-loc_list[[l]]
      loc_pop<-basin_pop%>%filter(location==iter_loc) #get pop data for location
      target_theta<-theta_data%>% #get theta for location x year
        filter(location==iter_loc,
               return_year==iter_year)
      if(nrow(target_theta)==0){target_theta=0.2} else{ #if no theta, set default
        target_theta<-round(target_theta$theta,2)
      }
      
      recovered_fish<-sim_spawning_recovery(loc_pop,target_theta)
      
      loc_recoveries[[l]]<-recovered_fish
    }
    names(loc_recoveries)=loc_list
    
    #for each location produce the required data that a coordinator would need
    loc_data<-list()
    for(l in 1:length(loc_list)){
      iter_loc<-loc_list[[l]]
      d<-loc_recoveries[[l]]
      N_loc<-basin_pop%>% #estimate N_loc escapement estimate (this would be done by the survey group and have uncertainty in it)
        filter(location==iter_loc)%>%
        nrow()
      n_nocwt=d%>%filter(has_cwt==0) #number of recovered fish without cwt
      n_hascwt=d%>%filter(has_cwt==1) #number of recovered fish with cwt
      
      cwt_uid_freq<-n_hascwt%>% #for each cwt group, provide frequency
        group_by(location,release_group,tag_rate,return_year,
                 brood_year,tag_age,brood_location)%>%
        summarise(frequency=n())
      theta=nrow(d)/N_loc #calculate theta
      
      ans<-list("N_loc"=N_loc,
                "n_recovered"=nrow(d),
                "theta"=theta,
                "F_prod"=cwt_uid_freq)
      loc_data[[l]]=ans
    }
    names(loc_data)=loc_list
    
    ###################################
    #Methods Step 2: Use k-draws to estimate N_nocwt
    #this process could be done by either a coordinator or surveys themselves
    ###################################
    #for each location run iterated k-draws to estimate N_nocwt and uncertainty
    loc_kdraws<-list()
    loc_results<-list()
    survey_year=iter_year #set survey year for surveys with no cwt recovered
    for(l in 1:length(loc_list)){
      d<-loc_data[[l]]
      d_fprod<-d$F_prod
      
      uid_list<-unique(d_fprod$release_group)
      uid_results<-data.frame()
      if(length(uid_list)>1){
        for(u in 1:length(uid_list)){
          tag_data<-d_fprod%>%filter(release_group==uid_list[u])
          tag_outputs<-tag_k_draws(tag_data=tag_data,
                                   theta=d$theta,
                                   iterations=iterations,
                                   freq=tag_data$frequency,
                                   f_prod=tag_data$tag_rate)
          
          uid_results<-uid_results%>%rbind(tag_outputs)
        }
      }else{
        uid_results <- data.table::data.table(
          cwt = rep(NA,iterations),
          brood_year=rep(NA,iterations),
          tag_age=rep(NA,iterations),
          return_year=rep(survey_year,iterations),
          location=rep(loc_list[l],iterations),
          theta=rep(d$theta,iterations),
          tag_rate=rep(NA,iterations),
          tag_freq =rep(0,iterations),
          total_tags = rep(0,iterations),
          total_hatchery = rep(0,iterations),
          k_iteration = rep(0,iterations)
        )
      }
      
      loc_kdraws[[l]]<-uid_results
      
      location_estimates<-uid_results%>%
        group_by(k_iteration,location,return_year,theta)%>%
        summarise(recovered_tags=sum(tag_freq),
                  mean_tag_rate=round(mean(tag_rate),2),
                  total_tags=sum(total_tags),
                  total_hatchery=sum(total_hatchery))
      location_estimates<-location_estimates%>%
        group_by(location,return_year,mean_tag_rate,
                 recovered_tags)%>%
        summarise(total_tags=round(mean(total_tags)),
                  total_hatchery=round(mean(total_hatchery)))
      location_estimates<-location_estimates%>%
        mutate(N_loc=d$N_loc)%>% #add location total escapement
        mutate(N_nocwt=N_loc-total_tags) #estimate N untagged
      
      loc_results<-loc_results%>%rbind(location_estimates)
    }
    
    names(loc_kdraws)=loc_list #save k-draws for later use in uncertainty
    
    scale_results<-list()
    
    for(s in 1:length(boot_scales)){
      
      scale_boot_start<-Sys.time()
      target_scales=boot_scales[s]
      
      #estimate N_nocwt proportion from each survey
      basin_nocwt<-sum(loc_results$N_nocwt)
      loc_results<-loc_results%>%
        mutate(nocwt_prop=N_nocwt/basin_nocwt)
      
      #use round_to_sum function to properly round scales allotment
      loc_results$target_scales=round_to_sum(loc_results$nocwt_prop,target_scales)
      
      #get number of nocwt scales available at each location
      scales_available<-data.frame()
      for(l in loc_list){
        d<-loc_recoveries[[l]]
        d<-d%>%
          filter(has_cwt==0)
        available=data.frame(location=l,
                             nocwt_scales_available=nrow(d))
        scales_available<-scales_available%>%rbind(available)
      }
      
      loc_results<-loc_results%>%
        left_join(scales_available)
      
      #if we are short on available scales, deal with that
      loc_results<-loc_results%>%
        mutate(
          nocwt_scales_available=ifelse(is.na(nocwt_scales_available),0
                                        ,nocwt_scales_available),
          #set actual_scales to either target_scales or available, which ever is smallest
          actual_scales=pmin(target_scales,nocwt_scales_available),
          shortfall=target_scales-actual_scales
        )
      
      total_shortfall=sum(loc_results$shortfall)
      
      ###################################
      #Methods Step 3: Sample scales from locations
      #this is done randomly by the survey, but simulated here
      ###################################
      
      loc_scales<-vector(mode = "list", length = length(loc_list))
      names(loc_scales)<-loc_list
      for(l in loc_list){
        d<-loc_recoveries[[l]]
        scales_n=loc_results%>%
          filter(location==l)
        d<-d%>%
          filter(has_cwt==0)
        d_sample<-sample_n(d,scales_n$actual_scales,replace = F)
        loc_scales[[l]]=d_sample
      }
      
      basin_scales<-data.frame()
      for(l in loc_list){
        basin_scales<-basin_scales%>%rbind(loc_scales[[l]])
      }
      basin_scales<-basin_scales%>%mutate(scale_age=return_year-brood_year)
      
      #now we assume scales are sent in for aging, and we will just get the results as basin_scales
      
      ###################################
      #Methods Step 4: get basin hachery estimates
      #done by coordinator
      ###################################
      
      #first take loc_kdraws and get basin wide estimate of hatchery origin fish for each iteration
      loc_hatcheries<-data.frame()
      for(l in loc_list){
        d<-loc_kdraws[[l]]
        d<-d%>%
          group_by(return_year,tag_age,location,k_iteration,theta)%>%
          summarise(tag_freq=sum(tag_freq),
                    total_tags=sum(total_tags),
                    total_hatchery=sum(total_hatchery))
        loc_hatcheries<-loc_hatcheries%>%rbind(d)
      }
      
      #sometimes locations have no hatchery fish so remove here
      loc_hatcheries<-loc_hatcheries%>%filter(!is.na(tag_age))
      
      #sum location based hatchery estimates from k-draws into basin-wide
      basin_hatchery_estimates<-loc_hatcheries%>%
        group_by(return_year,tag_age,k_iteration)%>%
        summarise(tag_freq = sum(tag_freq),
                  total_tags = sum(total_tags),
                  total_hatchery = sum(total_hatchery))
      
      spawning_results=list("basin_scales"=basin_scales, #irl this will be our aging results
                            "hatchery_ages"=basin_hatchery_estimates)
      
      cr_results<-basin_CR(basin_N=basin_N,
                           spawning_results = spawning_results,
                           iterations = iterations,
                           CI=0.9)
      
      ###################################
      #compare results to true values
      ###################################
      
      true_values<-basin_pop%>%
        mutate(age=return_year-brood_year)%>%
        group_by(age,origin)%>%
        summarize(true_val=n())
      
      compare<-cr_results%>%left_join(true_values)
      compare<-compare%>%
        mutate(count_accuracy=abs(estimate-true_val)/true_val)
      
      compare$tag_rate=t_inc
      
      #save outputs
      scale_results[[s]]<-list(
        "comparison"=compare,
        "loc_results"=loc_results
      )
      names(scale_results)=boot_scales
      scale_boot_end<-Sys.time()
      cat(red(paste("\ntag_rate boot",t,"DONE, time:",round(scale_boot_end-scale_boot_start,2),"\n")))
    }
    
    tag_results[[t]]<-scale_results
    
  }
  names(tag_results)<-tag_rate
  pooling_iteration_results[[i]]<-tag_results
  iter_end<-Sys.time()
  cat(red(paste("\niteration",i,"DONE, time:",round(iter_end-iter_start,2),"\n")))
}
saveRDS(pooling_iteration_results,"outputs/methods_comparisons_inc_tagrate.Rds")

pooling_iteration_results<-readRDS("outputs/methods_comparisons_inc_tagrate.Rds")
#plotting

################################
#with 100% tag rate
###############################
year_list<-2010:2020
boot_replicates<-500
iterations=1000
boot_scales<-seq(100,2000,100)
iter_results<-list()
pooling_iteration_results<-list()

for(i in 1:boot_replicates){
  iter_start<-Sys.time()
  ###################################
  #simulate population
  ###################################
  iter_year=sample(year_list,1)
  #generate basin-wide population
  #now updated with brood_year and brood_location
  basin_pop<-sim_basin_pop(trib_data,i,iter_year)$iter_pop
  
  basin_pop<-basin_pop%>%
    mutate(has_cwt=ifelse(origin=="hatchery",1,0),
           tag_rate=ifelse(origin=="hatchery",1,0))%>%
    mutate(cwt_uid=ifelse(origin=="hatchery",paste(brood_location,brood_year,sep=" "),cwt_uid))
  
  #grab total basin pop
  basin_N=nrow(basin_pop)
  
  loc_list<-unique(basin_pop$location)
  
  ###################################
  #Methods Step 1: Escapement Surveys
  #spawning recoveries for each location
  #this should be the data collected and prepared by survey teams
  ###################################
  loc_recoveries<-list()
  for(l in 1:length(loc_list)){
    iter_loc<-loc_list[[l]]
    loc_pop<-basin_pop%>%filter(location==iter_loc) #get pop data for location
    
    target_theta<-theta_data%>% #get theta for location x year
      filter(location==iter_loc,
             return_year==iter_year)
    if(nrow(target_theta)==0){target_theta=0.2} else{ #if no theta, set default
      target_theta<-round(target_theta$theta,2)
    }
    
    recovered_fish<-sim_spawning_recovery(loc_pop,target_theta)
    
    loc_recoveries[[l]]<-recovered_fish
  }
  names(loc_recoveries)=loc_list
  
  #for each location produce the required data that a coordinator would need
  loc_data<-list()
  for(l in 1:length(loc_list)){
    iter_loc<-loc_list[[l]]
    d<-loc_recoveries[[l]]
    N_loc<-basin_pop%>% #estimate N_loc escapement estimate (this would be done by the survey group and have uncertainty in it)
      filter(location==iter_loc)%>%
      nrow()
    n_nocwt=d%>%filter(has_cwt==0) #number of recovered fish without cwt
    n_hascwt=d%>%filter(has_cwt==1) #number of recovered fish with cwt
    
    cwt_uid_freq<-n_hascwt%>% #for each cwt group, provide frequency
      group_by(location,cwt_uid,tag_rate,return_year,
               brood_year,age,brood_location)%>%
      summarise(frequency=n())
    theta=nrow(d)/N_loc #calculate theta
    
    ans<-list("N_loc"=N_loc,
              "n_recovered"=nrow(d),
              "theta"=theta,
              "F_prod"=cwt_uid_freq)
    loc_data[[l]]=ans
  }
  names(loc_data)=loc_list
  
  ###################################
  #Methods Step 2: Use k-draws to estimate N_nocwt
  #this process could be done by either a coordinator or surveys themselves
  ###################################
  #for each location run iterated k-draws to estimate N_nocwt and uncertainty
  loc_kdraws<-list()
  loc_results<-list()
  survey_year=iter_year #set survey year for surveys with no cwt recovered
  for(l in 1:length(loc_list)){
    d<-loc_data[[l]]
    d_fprod<-d$F_prod
    
    uid_list<-unique(d_fprod$cwt_uid)
    uid_results<-data.frame()
    if(length(uid_list)>1){
      for(u in 1:length(uid_list)){
        tag_data<-d_fprod%>%filter(cwt_uid==uid_list[u])
        tag_outputs<-sim_k_draws(tag_data=tag_data,
                                 theta=d$theta,
                                 iterations=iterations,
                                 freq=tag_data$frequency,
                                 f_prod=1)
        
        uid_results<-uid_results%>%rbind(tag_outputs)
      }
    }else{
      uid_results <- data.table::data.table(
        cwt = rep(NA,iterations),
        brood_year=rep(NA,iterations),
        age=rep(NA,iterations),
        return_year=rep(survey_year,iterations),
        location=rep(loc_list[l],iterations),
        theta=rep(d$theta,iterations),
        tag_rate=rep(NA,iterations),
        tag_freq =rep(0,iterations),
        total_tags = rep(0,iterations),
        total_hatchery = rep(0,iterations),
        k_iteration = rep(0,iterations)
      )
    }
    
    loc_kdraws[[l]]<-uid_results
    
    location_estimates<-uid_results%>%
      group_by(k_iteration,location,return_year,theta)%>%
      summarise(recovered_tags=sum(tag_freq),
                mean_tag_rate=round(mean(tag_rate),2),
                total_tags=sum(total_tags),
                total_hatchery=sum(total_hatchery))
    location_estimates<-location_estimates%>%
      group_by(location,return_year,mean_tag_rate,
               recovered_tags)%>%
      summarise(total_tags=round(mean(total_tags)),
                total_hatchery=round(mean(total_hatchery)))
    location_estimates<-location_estimates%>%
      mutate(N_loc=d$N_loc)%>% #add location total escapement
      mutate(N_nocwt=N_loc-total_tags) #estimate N untagged
    
    loc_results<-loc_results%>%rbind(location_estimates)
  }
  
  names(loc_kdraws)=loc_list #save k-draws for later use in uncertainty
  
  scale_results<-list()
  
  for(s in 1:length(boot_scales)){
    
    scale_boot_start<-Sys.time()
    target_scales=boot_scales[s]
    
    #estimate N_nocwt proportion from each survey
    basin_nocwt<-sum(loc_results$N_nocwt)
    loc_results<-loc_results%>%
      mutate(nocwt_prop=N_nocwt/basin_nocwt)
    
    #use round_to_sum function to properly round scales allotment
    loc_results$target_scales=round_to_sum(loc_results$nocwt_prop,target_scales)
    
    #get number of nocwt scales available at each location
    scales_available<-data.frame()
    for(l in loc_list){
      d<-loc_recoveries[[l]]
      d<-d%>%
        filter(has_cwt==0)
      available=data.frame(location=l,
                           nocwt_scales_available=nrow(d))
      scales_available<-scales_available%>%rbind(available)
    }
    
    loc_results<-loc_results%>%
      left_join(scales_available)
    
    #if we are short on available scales, deal with that
    loc_results<-loc_results%>%
      mutate(
        nocwt_scales_available=ifelse(is.na(nocwt_scales_available),0
                                      ,nocwt_scales_available),
        #set actual_scales to either target_scales or available, which ever is smallest
        actual_scales=pmin(target_scales,nocwt_scales_available),
        shortfall=target_scales-actual_scales
      )
    
    total_shortfall=sum(loc_results$shortfall)
    
    ###################################
    #Methods Step 3: Sample scales from locations
    #this is done randomly by the survey, but simulated here
    ###################################
    
    loc_scales<-vector(mode = "list", length = length(loc_list))
    names(loc_scales)<-loc_list
    for(l in loc_list){
      d<-loc_recoveries[[l]]
      scales_n=loc_results%>%
        filter(location==l)
      d<-d%>%
        filter(has_cwt==0)
      d_sample<-sample_n(d,scales_n$actual_scales,replace = F)
      loc_scales[[l]]=d_sample
    }
    
    basin_scales<-data.frame()
    for(l in loc_list){
      basin_scales<-basin_scales%>%rbind(loc_scales[[l]])
    }
    
    #now we assume scales are sent in for aging, and we will just get the results as basin_scales
    
    ###################################
    #Methods Step 4: get basin hachery estimates
    #done by coordinator
    ###################################
    
    #first take loc_kdraws and get basin wide estimate of hatchery origin fish for each iteration
    loc_hatcheries<-data.frame()
    for(l in loc_list){
      d<-loc_kdraws[[l]]
      d<-d%>%
        group_by(return_year,age,location,k_iteration,theta)%>%
        summarise(tag_freq=sum(tag_freq),
                  total_tags=sum(total_tags),
                  total_hatchery=sum(total_hatchery))
      loc_hatcheries<-loc_hatcheries%>%rbind(d)
    }
    
    #sometimes locations have no hatchery fish so remove here
    loc_hatcheries<-loc_hatcheries%>%filter(!is.na(age))
    
    #sum location based hatchery estimates from k-draws into basin-wide
    basin_hatchery_estimates<-loc_hatcheries%>%
      group_by(return_year,age,k_iteration)%>%
      summarise(tag_freq = sum(tag_freq),
                total_tags = sum(total_tags),
                total_hatchery = sum(total_hatchery))
    
    spawning_results=list("basin_scales"=basin_scales, #irl this will be our aging results
                          "hatchery_ages"=basin_hatchery_estimates)
    
    cr_results<-basin_CR(basin_N=basin_N,
                         spawning_results = spawning_results,
                         iterations = iterations,
                         scales_n = target_scales,
                         CI=0.9)
    
    ###################################
    #compare results to true values
    ###################################
    
    true_values<-basin_pop%>%
      group_by(age,origin)%>%
      summarize(abundance=n())%>%
      pivot_wider(names_from=origin,
                  values_from=abundance)
    
    compare<-cr_results%>%left_join(true_values)
    compare<-compare%>%
      mutate(count_accuracy=abs(mean_natural-natural)/natural)
    
    #save outputs
    scale_results[[s]]<-list(
      "comparison"=compare,
      "loc_results"=loc_results
    )
    
    scale_boot_end<-Sys.time()
    cat(red(paste("\nscale boot",boot_scales[s],"DONE, time:",round(scale_boot_end-scale_boot_start,2),"\n")))
  }
  
  names(scale_results)=boot_scales
  pooling_iteration_results[[i]]<-scale_results
  iter_end<-Sys.time()
  cat(red(paste("\niteration",i,"DONE, time:",round(iter_end-iter_start,2),"\n")))
}
saveRDS(pooling_iteration_results,"outputs/methods_comparisons_tag100.Rds")

pooling_iterations<-readRDS("outputs/methods_comparisons_tag100.Rds")
target_accuracy=.150
pooling_estimates<-data.frame()
for(i in 1:length(pooling_iterations)){
  d<-pooling_iterations[[i]]
  for(s in 1:length(d)){
    ds<-d[[s]]$comparison
    ds$n_scales<-as.numeric(names(d[s]))
    
    pooling_estimates<-pooling_estimates%>%rbind(ds)
  }
}

stats<-pooling_estimates%>%
  mutate(meets_target_accuracy=ifelse(count_accuracy<=target_accuracy,T,F))%>%
  group_by(age,n_scales)%>%
  summarise(
    pct_meeting_accuracy=mean(meets_target_accuracy)*100
  )

p2<-ggplot(data=stats,
           aes(x=(n_scales),y=pct_meeting_accuracy))+
  geom_point( size=3, aes(color=factor(age)))+
  geom_hline(yintercept=90,color="red")+
  scale_x_continuous(limits=c(100,2000),
                     breaks=seq(from=100,to=2000,by=100))+
  theme_bw()+
  labs(x = "basin scale samples", y = "% of estimates that were within +/- 15% of true value",color="Age",linetype="Age")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none")+
  facet_grid(factor(age)~.)
p2
ggsave(p2,file="outputs/count_accuracy_tagrate100.png",width=10,height=6)

###############################
#setting sample scale size with limiting location
###############################