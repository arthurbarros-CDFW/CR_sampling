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
#individual run
##########################
#simulate population
N=10000 #true population size
theta=0.2 #probability of recovering a given fish or sampling fraction
tag_rate=0.25 #fixed CFM rate
k_iterations=1000 #
scale_iterations=1000 # needs to be same as k_iterations (consider just one iter value)
survey="trib" #hatchery or trib
scales_n_seq =seq(400,2000,100)
n_replicates=100 #times to repeat total simulation
CI=0.95
h_pct=75
h_prop=.75
target_moe=0.20

sim_results<-test_scale_requirements(N,
                                 theta,
                                 k_iterations,
                                 scale_iterations,
                                 tag_rate,
                                 survey,
                                 scales_n_seq,
                                 n_replicates,
                                 CI,
                                 h_prop,
                                 target_moe)

saveRDS(sim_results,file=paste("outputs/hprop_",h_pct,"_results.rds",sep=""))
sim_results<-readRDS(paste("outputs/hprop_",h_pct,"_results.rds",sep=""))

summary_stats<-sim_summary_stats(sim_results,target_moe)

p<-plot_scale_sims(summary_stats,N,n_replicates)
ggsave(p$ages_plot,file=paste("outputs/hprop_",h_pct,"_ageplots.png",sep=""),width=10,height=6)


##########################
#scenario runs
##########################
#N_dist<-seq(from=5000,to=10000,by=2500)
N_dist<-10000
#tag_rate_dist<-seq(0.25,0.75,by=0.25)
tag_rate_dist<-.25
theta_list=seq(0.2,0.2,.1) #probability of recovering a given fish or sampling fraction
scales_n_seq =seq(500,1000,100)#number of scales to age from untagged recoveries
#h_prop_dist = seq(0.5,0.9,0.1)
h_prop_dist = 1
h_pct_dist = 1
k_iterations=1000 #
scale_iterations=1000
n_replicates=100 #times to repeat total simulation
CI=0.95
target_moe=0.10

results_list<-list()

seq_start<-Sys.time()
for(n in 1:length(N_dist)){
  N_current=N_dist[n]
  for(t in 1:length(tag_rate_dist)){
    tag_rate_current=tag_rate_dist[t]
    for(h in 1:length(h_prop_dist)){
      h_current=h_prop_dist
      for(p in 1:length(theta_list)){
        boot_start<-Sys.time()
        theta<-theta_list[p]
        h_current=h_prop_dist[h]
        
        results<-test_scale_requirements(N_current,
                                         theta,
                                         k_iterations,
                                         scale_iterations,
                                         tag_rate_current,
                                         survey,
                                         scales_n_seq,
                                         n_replicates,
                                         CI,
                                         h_prop=h_current)
        
        saveRDS(results,file=paste("outputs/scenario_data/",
                                   "N",N_current,"_",
                                   "tag_rate",tag_rate_current,"_",
                                   "theta",theta,
                                   ".rds",sep=""))
        
        results$age_results<-results$age_results%>%
          filter(scales_sampled>=scales_n)
        
        summary_stats<-sim_summary_stats(results,target_moe)
        
        #plot<-plot_scale_sims(summary_stats,N_current,n_replicates)
        
        #ggsave(plot$ages_plot,file=paste("outputs/scenario_figures/ages/",
        #                              "N",N_current,"_",
        #                              "tag_rate",tag_rate_current,"_",
        #                              "theta",theta,
        #                              ".png",sep=""),scale=2)
        
        #ggsave(plot$totals_plot,file=paste("outputs/scenario_figures/totals/",
        #                                "N",N_current,"_",
        #                                "tag_rate",tag_rate_current,"_",
        #                                "theta",theta,
        #                                ".png",sep=""),scale=2)
      #  
        results_list[[p]]<-results
        
        boot_end<-Sys.time()
        boot_time<-boot_end-boot_start
        
        print(paste("FINISHED:","N:",N_current,
                    "tag_rate:",tag_rate_current,
                    "theta:",theta,
                    "h_prop:",h_current,
                    " boot time:",round(boot_time,4),sep=" "))
      }
      
    }
  }
}
names(results_list)<-theta_list
seq_end<-Sys.time()
seq_time<-seq_end-seq_start
saveRDS(results_list,"outputs/false_positives_test.Rds")

########################
#Vizualize
########################
#1)n_scales
sim_results<-readRDS("outputs/hprop_75_results.rds")
summary_stats<-sim_summary_stats(sim_results,target_moe=.20)
ages_results<-summary_stats$age_summary_stats

#number of iterations within target moe
p<- ggplot(ages_results, aes(x = scales_n, y = pct_meeting_target_moe_prop, color = factor(age))) +
  geom_line() +
  geom_point() +
  scale_y_continuous(limits=c(0,100),
                     breaks=seq(from=0,to=100,by=10))+
  scale_x_continuous(limits=c(min(ages_results$scales_n),max(ages_results$scales_n)),
                     breaks=seq(from=0,to=max(ages_results$scales_n),by=100))+
  labs(title = paste("% iterations meeting target MOE (",target_moe*100,"%)",sep=""),
       x = "Number of scale samples", 
       y = "% iterations meeting target MOE",
       color = "Age") +
  theme_minimal()
ggsave(p,file="outputs/fig1.png",width=10,height=6)

#2) theta plots
theta_results<-readRDS("outputs/theta_results.Rds")

t_list<-list()
for(i in 1: length(theta_results)){
  d=theta_results[[i]]$age_results
  d$theta<-as.numeric(names(theta_results)[[i]])
  t_list[[i]]<-d
}
theta_data<-rbindlist(t_list)

#plot one set of scale n
target_moe=.20
d<-theta_data%>%filter(scales_n==500)
d<-d%>%
  mutate(meets_target_moe_prop=ifelse(moe_prop<target_moe,T,F))%>%
  mutate(meets_target_moe_count=ifelse((moe_count/est_count_natural)<target_moe,T,F))
d<-d%>%
  group_by(age,theta)%>%
  summarise(mean_bias = mean(bias, na.rm = TRUE),
            sd_bias = sd(bias, na.rm = TRUE),
            rmse = sqrt(mean(bias^2, na.rm = TRUE)),
            mean_ci_width_prop = mean(ci_width_prop, na.rm = TRUE),
            coverage_rate = mean(coverage, na.rm = TRUE),
            mean_moe_prop=mean(moe_prop, na.rm = TRUE),
            mean_moe_count=mean(moe_count,na.rm=TRUE),
            sd_moe_prop=sd(moe_prop, na.rm = TRUE),
            sd_moe_count=sd(moe_count,na.rm=TRUE),
            pct_meeting_target_moe_prop = mean(meets_target_moe_prop) * 100,
            pct_meeting_target_moe_count = mean(meets_target_moe_count) * 100,
            .groups = "drop")

p<- ggplot(d, aes(x = theta, y = pct_meeting_target_moe_prop, color = factor(age))) +
  geom_line() +
  geom_point() +
  scale_y_continuous(limits=c(0,100),
                     breaks=seq(from=0,to=100,by=10))+
  scale_x_continuous(limits=c(0,1),
                     breaks=seq(from=0,to=1,by=.1))+
  labs(title = paste("% iterations meeting target MOE (",target_moe*100,"%)",sep=""),
       subtitle = "n scales = 500",
       x = "Theta value", 
       y = "% iterations meeting target MOE",
       color = "Age") +
  theme_minimal()
p
ggsave(p,file="outputs/theta_plot.png",width=10,height=6)

p<-ggplot(location_year_theta)+
  geom_histogram(aes(x=theta))+
  labs(title = "Histogram of observed theta values",
       x = "Theta value", 
       y = "Frequency") +
  theme_minimal()
p
ggsave(p,file="outputs/theta_hist.png",width=10,height=6)


