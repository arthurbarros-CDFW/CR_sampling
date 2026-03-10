rm( list = ls()) #clear env
#pooling results analysis
library(tidyverse)
library(data.table)
library(ggplot2)
library(gtools)

#Load functions
sapply(list.files("scripts/functions", pattern = "\\.R$", full.names = TRUE), source)
CI=0.95
target_moe<-0.2

##################################
#no pooling
##################################

nopooling_iteration_results<-readRDS("outputs/pooling/results_no_pooling_test.Rds")

sim_list <- basin_scales<-list()

for(i in 1:length(nopooling_iteration_results)) {
  iteration_est <- as.data.table(nopooling_iteration_results[[i]]$basin_sim_summary)
  iteration_true <- as.data.table(nopooling_iteration_results[[i]]$basin_true_ages)
  
  iteration_est[, iter := i]
  iteration_true[, iter := i]
  
  #data.table join
  sim_list[[i]] <- iteration_est[iteration_true, on = .( age, iter), nomatch = 0]
  basin_scales[[i]] <- nopooling_iteration_results[[i]]$basin_scales_collected
}

#combine all iterations
all_est <- rbindlist(sim_list)

all_est<-all_est%>%
  mutate(meets_target_moe=ifelse(moe_prop<=target_moe,T,F),
         count_bias=(est_count_natural-true_count_natural)/true_count_natural)

p<-ggplot(all_est, aes(x = factor(age), y = 100*moe_prop, color = factor(age))) +
  geom_boxplot() +
  labs(subtitle = mean(as.numeric(basin_scales)),
       x = "year", 
       y = "MOE %",
       color = "Age") +
  theme_minimal()
p
ggsave(p,file="outputs/no_pooling.png",width=8,height=5)

p_count_bias<-ggplot(all_est, aes(x = factor(age), y = 100*count_bias, color = factor(age))) +
  geom_boxplot() +
  labs(x = "year", 
       y = "count bias %",
       color = "Age") +
  theme_bw()
p_count_bias
ggsave(p_count_bias,file="outputs/nopooling_count_bias.png",width=8,height=5)

#no_pooling regional examination
# Get all no_pooling RDS files
no_pooling_files <- list.files(path = "outputs/pooling", 
                               pattern = "no_pooling_iter_.*\\.Rds$", 
                               full.names = TRUE)

iter_numbers <- as.numeric(gsub(".*no_pooling_iter_(\\d+)_.*", "\\1", no_pooling_files))
no_pooling_files<-no_pooling_files[iter_numbers >= 1 & iter_numbers <= 10]

# Initialize lists to store data
all_location_age_stats <- list()
all_basin_true_ages <- list()

# Loop through each file and extract data
for(file in no_pooling_files) {
  
  # Extract iteration number from filename
  iter_num <- as.numeric(gsub(".*no_pooling_iter_(\\d+)_.*", "\\1", basename(file)))
  
  # Read the boot_results data
  boot_results <- readRDS(file)
  
  # For each scale iteration, extract data
  for(scale_iter in 1:length(boot_results)) {
    
    scale_results <- boot_results[[scale_iter]]
    
    # Get the scale_iteration_results which contains location-specific data
    loc_results <- scale_results$scale_iteration_results
    
    # For each location, extract age_summary_stats
    for(loc_name in names(loc_results)) {
      
      loc_data <- loc_results[[loc_name]]
      
      # Extract basin_true_ages
      if(!is.null(loc_data$true_summary_stats)) {
        temp_data <- loc_data$true_summary_stats %>%
          mutate(
            iteration = iter_num,
            scale_iteration = scale_iter,
            total_scales = scale_results$scales_collected
          )
        temp_data$location=loc_name
        all_basin_true_ages[[length(all_basin_true_ages) + 1]] <- temp_data
      }
      
      if(!is.null(loc_data$age_summary_stats)) {
        
        # Add metadata columns
        temp_data <- loc_data$age_summary_stats %>%
          mutate(
            iteration = iter_num,
            scale_iteration = scale_iter,
            location = loc_name,
            scales_collected_loc = loc_data$scales_collected,
            total_scales = scale_results$scales_collected
          )
        
        # Store in list
        all_location_age_stats[[length(all_location_age_stats) + 1]] <- temp_data
      }
    }
  }
}

# Combine all into data frames
all_location_stats_est <- bind_rows(all_location_age_stats)
all_location_stats_true <- bind_rows(all_basin_true_ages)

# Make basin_true_ages unique (as in your original code)
all_location_stats_true <- all_location_stats_true %>%
  distinct(across(-c(iteration, scale_iteration, total_scales)), .keep_all = TRUE)

summarise_est<-all_location_stats_est%>%
  group_by(age,location,iteration,scales_collected_loc)%>%
  summarise(mean_count_natural=mean(est_count_natural,),
            lower_CI_count = quantile(est_count_natural, probs = (1-CI)/2, na.rm = TRUE),
            upper_CI_count = quantile(est_count_natural, probs = 1-(1-CI)/2, na.rm = TRUE),
            new_proportion_natural = mean(mean_proportion_natural, na.rm = TRUE),
            lower_CI_prop = quantile(mean_proportion_natural, probs = (1-CI)/2, na.rm = TRUE),
            upper_CI_prop = quantile(mean_proportion_natural, probs = 1-(1-CI)/2, na.rm = TRUE))

for(i in 1:length(unique(summarise_est$iteration))){
  d_est<-summarise_est%>%filter(iteration==i)
  d_true<-all_location_stats_true%>%filter(iteration==i)
  
  d_true<-d_true%>%
    left_join(select(d_est,iteration,location,age,scales_collected_loc))
  
  d_est <- d_est %>%
    group_by(location) %>%
    mutate(location_with_scales = paste0(location, "\n(scales: ", scales_collected_loc, ")"))
  
  d_true <- d_true %>%
    group_by(location) %>%
    mutate(location_with_scales = paste0(location, "\n(scales: ", scales_collected_loc, ")"))
  
  p<-ggplot(d_est,aes(x=factor(age),
                      y=log10(mean_count_natural),
                      color=factor(age)))+
    geom_point(size=3)+
    geom_errorbar(aes(ymin = log10(lower_CI_count), ymax = log10(upper_CI_count)), 
                  width = 0.2) +
    geom_point(data=d_true,aes(x=factor(age),
                                  y=log10(true_count_natural)),size=3,color="black",shape=9)+
        theme_bw()+
    labs(title="simulated estimates",
         x = "age", 
         y = "log10(count)",
         color="age")+
    facet_grid(.~location_with_scales)
  ggsave(file=paste("outputs/pooling/no_pooling_figures/",i,"_log10.png",sep=""),width=10,height=6)

  p<-ggplot(d_est,aes(x=factor(age),
                      y=(mean_count_natural),
                      color=factor(age)))+
    geom_point(size=3)+
    geom_errorbar(aes(ymin = (lower_CI_count), ymax = (upper_CI_count)), 
                  width = 0.2) +
    geom_point(data=d_true,aes(x=factor(age),
                               y=(true_count_natural)),size=3,color="black",shape=9)+
    theme_bw()+
    labs(title="simulated estimates",
         x = "age", 
         y = "count",
         color="age")+
    facet_grid(.~location_with_scales)
  ggsave(file=paste("outputs/pooling/no_pooling_figures/",i,".png",sep=""),width=10,height=6)
  
}

##################################
#with pooling
##################################
pooling_iteration_results<-readRDS("outputs/pooling/results_with_pooling.Rds")

sim_list <- basin_scales<-list()

for(i in 1:length(pooling_iteration_results)) {
  iteration_results <- pooling_iteration_results[[i]]$scale_iters
  
  scales_list<-as.numeric(names(iteration_results))
  d_sim<-d_true<-list()
  for(n in 1:length(scales_list)){
    scales<-scales_list[[n]]
    d_sim[[n]]<-iteration_results[[n]]$age_summary_stats
    d_true[[n]]<-iteration_results[[n]]$true_summary_stats
  }
  iter_sim_results<-rbindlist(d_sim)
  iter_true_results<-rbindlist(d_true)
  iter_true_results<-unique(iter_true_results)
  
  #data.table join
  sim_list[[i]] <- iter_sim_results[iter_true_results, on = .( age), nomatch = 0]
}

#combine all iterations
all_est <- rbindlist(sim_list)

all_est<-all_est%>%
  mutate(meets_target_moe=ifelse(moe_prop<=target_moe,T,F))

all_est<-all_est%>%
  mutate(prop_bias=mean_proportion_natural-true_proportion_natural,
         count_bias=(est_count_natural-true_count_natural)/true_count_natural)

p_moe<-ggplot(all_est, aes(x = factor(age), y = moe_prop*100, color = factor(age))) +
  geom_boxplot() +
  labs(x = "year", 
       y = "MOE %",
       color = "Age") +
  theme_bw()+
  facet_grid(.~as.numeric(scales_collected))
p_moe
ggsave(p_moe,file="outputs/pooling_moe.png",width=8,height=5)

p_bias<-ggplot(all_est, aes(x = factor(age), y = 100*prop_bias, color = factor(age))) +
  geom_boxplot() +
  labs(x = "year", 
       y = "bias %",
       color = "Age") +
  theme_bw()+
  facet_grid(.~as.numeric(scales_collected))
p_bias

ggsave(p_bias,file="outputs/pooling_prop_bias.png",width=8,height=5)

p_count_bias<-ggplot(all_est, aes(x = factor(age), y = 100*count_bias, color = factor(age))) +
  geom_boxplot() +
  labs(x = "year", 
       y = "count bias %",
       color = "Age") +
  theme_bw()+
  facet_grid(.~as.numeric(scales_collected))
p_count_bias
ggsave(p_count_bias,file="outputs/pooling_count_bias.png",width=8,height=5)

all_est<-all_est%>%
  mutate(moe_count=(upper_CI_count-lower_CI_count)/2,
         moe_count_relative=(moe_count/est_count_natural))

p_countmoe<-ggplot(all_est, aes(x = factor(age),
                                y = moe_count_relative*100,
                                color = factor(age))) +
  geom_boxplot() +
  labs(x = "year", 
       y = "MOE %",
       color = "Age") +
  scale_y_continuous(limits=c(0,100),
                     breaks=seq(from=0,to=100,by=10))+
  theme_bw()+
  facet_grid(.~as.numeric(scales_collected))
p_countmoe
ggsave(p_countmoe,file="outputs/pooling_countmoe.png",width=8,height=5)
