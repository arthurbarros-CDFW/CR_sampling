#no pooling testing location bias
#iteration population

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
#Without pooling
##########################
#The below should give us a list of all the estimates for each loc x year x
#spawning size
#we can use those age count estimates and join them to get basin wide
#this is our "do n_scales at each location approach"

iterations=1000 #
scales_n_seq =500

CI=0.95
h_prop=.75

#read in data
trib_data<-readRDS("data/all_returns_clean.Rds")
theta_data<-readRDS("data/location_year_theta.Rds")

year_list<-2010:2020
iter_year<-2015
pop<-sim_basin_pop(trib_data,1,iter_year)
tag_rates<-pop$tag_rates

pop<-pop$iter_pop
loc_list<-unique(pop$location)#pull list of sampled locations

true_loc_ages<-pop%>%
  group_by(age,location,origin)%>%
  summarise(true_abundance=n())

true_basin_ages<-pop%>%
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

current_scales_n <- 500
  
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

results_summary <- location_results %>%
  #convert the list to a dataframe
  bind_rows(.id = "location") %>%
  #reshape from wide to long
  pivot_longer(
    cols = `2`:`4`,
    names_to = "age",
    values_to = "abundance"
  ) %>%
  group_by(location, age) %>%
  summarise(
    mean_abundance = mean(abundance, na.rm = TRUE),
    sd_abundance = sd(abundance, na.rm = TRUE),
    n = n(),
    ci_lower = quantile(abundance, probs = (1-CI)/2, 
                        na.rm = TRUE),
    ci_upper = quantile(abundance, probs = 1-(1-CI)/2, 
                        na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(age = factor(age, levels = c("2", "3", "4")))


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

year_theta<-theta_data%>%filter(return_year==iter_year)

for(l in 1: length(loc_list)){
  d<-results_summary%>%filter(location==loc_list[[l]])
  t<-true_loc_ages%>%filter(location==loc_list[[l]])
  t<-t%>%mutate(age=factor(age))
  theta<-year_theta%>%filter(location==loc_list[l])
  ggplot()+
    geom_point(data=d,aes(x=age,y=mean_abundance,color=age))+
    geom_errorbar(data=d,
                  aes(x=age,
                      ymin = ci_lower, 
                      ymax = ci_upper,
                      color=age,
                      linetype= age), 
                  width = 0.2)+
    geom_point(data=t,aes(x=age,y=true_abundance),color="black",shape=4)+
    theme_bw()+
    labs(x = "Age", y = "Estimate of Abundance",color="Age",linetype="Age",
         title=loc_list[[l]],
         subtitle = paste("theta:",round(theta$theta,4)))
  ggsave(paste("outputs/test_figures/loc_nopooling/",loc_list[[l]],iter_year,".png"))
}
ggplot()+
  geom_point(data=basin_summary,aes(x=factor(age),y=mean_natural_count,color=factor(age)))+
  geom_errorbar(data=basin_summary,
                aes(x=factor(age),
                    ymin = lower_CI_count, 
                    ymax = upper_CI_count,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2)+
  geom_point(data=true_basin_ages,aes(x=factor(age),y=true_abundance),color="black",shape=4)+
  theme_bw()+
  labs(x = "Age", y = "Estimate of Abundance",color="Age",linetype="Age",
       title = "Basin wide")
ggsave(paste("outputs/test_figures/loc_nopooling/basin",iter_year,".png"))
