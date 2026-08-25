rm( list = ls()) #clear env
#basin methods simulation
#this script aims to simulate basin populations and replicate methods
# to produce a basin-wide cohort reconstruction
library(tidyverse)
library(data.table)
library(ggplot2)
library(gtools)

#Load functions
sapply(list.files("scripts/functions", pattern = "\\.R$", full.names = TRUE), source)

set.seed(67) #for reproducibility

#read in data
trib_data<-readRDS("data/all_returns_clean.Rds")
theta_data<-readRDS("data/location_year_theta.Rds")

###################################
#simulate population
###################################
iter_year<-2020
i=1
target_scales=2000

#generate basin-wide population
#now updated with brood_year and brood_location
basin_pop<-sim_basin_pop(trib_data,i,iter_year)$iter_pop

#grab total basin pop
basin_N=nrow(basin_pop)

loc_list<-unique(basin_pop$location)

###################################
#Methods Step 1: Escapement Surveys
#simulate spawning recoveries for each location
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
  
  release_group_freq<-n_hascwt%>% #for each cwt group, provide frequency
    group_by(location,release_group,tag_rate,return_year,
             brood_year,tag_age,brood_location)%>%
    summarise(frequency=n())
  theta=nrow(d)/N_loc #estimate theta
  
  ans<-list("N_loc"=N_loc,
            "n_recovered"=nrow(d),
            "n_nocwt"=nrow(n_nocwt),
            "theta"=theta,
            "F_prod"=release_group_freq)
  loc_data[[l]]=ans
}
names(loc_data)=loc_list
saveRDS(loc_data,"outputs/sim_loc_data.rds")
saveRDS(loc_recoveries,"outputs/loc_recoveries.rds")

###################################
#Methods Step 2: Use k-draws to estimate N_nocwt
#this process could be done by either a coordinator or surveys themselves
###################################
loc_data<-readRDS("outputs/sim_loc_data.rds")
#for each location run iterated k-draws to estimate N_nocwt and uncertainty
loc_kdraws<-list()
loc_results<-list()
survey_year=iter_year #set survey year for surveys with no cwt recovered
iterations=1000
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
              total_tags=sum(total_tags),
              total_hatchery=sum(total_hatchery))
  location_estimates<-location_estimates%>%
    group_by(location,return_year,
             recovered_tags)%>%
    summarise(total_tags=round(mean(total_tags)),
              total_hatchery=round(mean(total_hatchery)))
  location_estimates<-location_estimates%>%
    mutate(N_loc=d$N_loc)%>% #add location total escapement
    mutate(N_untagged=N_loc-total_tags) #estimate N untagged
  
  loc_results<-loc_results%>%rbind(location_estimates)
}

names(loc_kdraws)=loc_list #save k-draws for later use in uncertainty
saveRDS(loc_kdraws,"outputs/loc_kdraws.rds")
saveRDS(loc_results,"outputs/loc_results.rds")

loc_kdraws<-readRDS("outputs/loc_kdraws.rds")

#estimate N_untagged proportion from each survey
basin_untagged<-sum(loc_results$N_untagged)
loc_results<-loc_results%>%
  mutate(untagged_prop=N_untagged/basin_untagged)

#use round_to_sum function to properly round scales allotment
loc_results$target_scales=round_to_sum(loc_results$untagged_prop,target_scales)

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
basin_scales<-basin_scales%>%
  mutate(scale_age=return_year-brood_year)

#now we assume scales are sent in for aging, 
#and we will just get the results as basin_scales

###################################
#Methods Step 4: get basin hatchery estimates
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
  summarize(abundance=n())%>%
  pivot_wider(names_from=origin,
              values_from=abundance)

compare<-cr_results%>%left_join(true_values)
compare<-compare%>%
  mutate(count_accuracy=abs(mean_natural-natural)/natural)
