#testing outside of quarto because its too slow
rm( list = ls()) #clear env
#data simulations
library(tidyverse)
library(data.table)
library(ggplot2)
library(gtools)
library(gridExtra)
#Load functions
sapply(list.files("scripts/functions", pattern = "\\.R$", full.names = TRUE), source)
target_moe=0.2
CI=0.95

#########################################
#approach 1: scales_n
########################################

scale_scenario<-readRDS("outputs/scale_sim_results.Rds")
scale_sims<-data.frame()
target_moe=0.2
CI=0.95
for(i in 1:length(scale_scenario)){
  d<-scale_scenario[[i]]
  scale_sizes<-data.frame()
  for(s in 1:length(d)){
    ds<-d[[s]]$sim_stats
    ds$scales_collected=d[[s]]$scale_results$scales_collected
    ds$scales_n=as.numeric(names(d)[[s]])
    ds$mean_proportion_natural=d[[s]]$scale_results$est_summary_stats$mean_proportion_natural
    ds$mean_count_natural=d[[s]]$scale_results$est_summary_stats$mean_count_natural
    ds$lower_CI_count=d[[s]]$scale_results$est_summary_stats$lower_CI_count
    ds$upper_CI_count=d[[s]]$scale_results$est_summary_stats$upper_CI_count
    ds$lower_CI_prop=d[[s]]$scale_results$est_summary_stats$lower_CI_prop
    ds$upper_CI_prop=d[[s]]$scale_results$est_summary_stats$upper_CI_prop
    scale_sizes<-scale_sizes%>%rbind(ds)
  }
  scale_sizes$iter=i
  scale_sims<-scale_sims%>%rbind(scale_sizes)
}

target_moe=.20
CI=.95

scale_sims<-scale_sims%>%
  mutate(meets_target_moe_prop=ifelse(prop_moe<=target_moe,T,F),
         meets_target_moe_count=ifelse(count_moe_relative<=target_moe,T,F))
scale_sims<-scale_sims%>%
  mutate(count_CI_width=upper_CI_count-lower_CI_count,
         prop_CI_width=upper_CI_prop-lower_CI_prop)

scale_stats<-scale_sims%>%
  group_by(age,scales_n)%>%
  summarize(
    mean_scales_collected=mean(scales_collected),
    mean_prop_bias = mean(prop_bias, na.rm = TRUE),
    mean_count_bias = mean(count_bias, na.rm = TRUE),
    lower_CI_prop_bias = quantile(prop_bias, probs = (1-CI)/2, 
                                  na.rm = TRUE),
    upper_CI_prop_bias = quantile(prop_bias, probs = 1-(1-CI)/2, 
                                  na.rm = TRUE),
    lower_CI_count_bias = quantile(count_bias, 
                                   probs = (1-CI)/2, na.rm = TRUE),
    upper_CI_count_bias = quantile(count_bias, 
                                   probs = 1-(1-CI)/2, na.rm = TRUE),
    mean_prop_moe=mean(prop_moe, na.rm = TRUE),
    lower_CI_prop_moe = quantile(prop_moe, 
                                   probs = (1-CI)/2, na.rm = TRUE),
    upper_CI_prop_moe = quantile(prop_moe, 
                                   probs = 1-(1-CI)/2, na.rm = TRUE),
    mean_count_moe=mean(count_moe_relative,na.rm=TRUE),
    sd_prop_moe=sd(prop_moe, na.rm = TRUE),
    sd_count_moe=sd(count_moe_relative,na.rm=TRUE),
    mean_count_CI_width=mean(count_CI_width),
    mean_prop_CI_width=mean(prop_CI_width),
    pct_meeting_target_prop_moe = mean(meets_target_moe_prop) * 100,
    pct_meeting_target_count_moe = mean(meets_target_moe_count) * 100,
    .groups = "drop"
  )

#at our sim theta, didn't get enough scales at scales_n >1600
scale_stats<-scale_stats%>%filter(scales_n<=1600)

#find smallest scale_n where 95% meets target moe
threshold_n <- scale_stats %>%
  group_by(age) %>%
  filter(pct_meeting_target_prop_moe >= 95) %>%
  slice_min(scales_n) %>%  # Get the smallest scales_n that meets the threshold
  ungroup()
threshold_value<-max(threshold_n$scales_n)

pbias<-ggplot(scale_stats, aes(x = scales_n,
                            y = mean_prop_bias,
                            color = factor(age),
                            linetype= factor(age))) +
  geom_point(aes(color=factor(age)))+
  geom_errorbar(data=scale_stats,
                aes(ymin = lower_CI_prop_bias, 
                    ymax = upper_CI_prop_bias,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2) +
  labs(title="Bias",
       x = "scale sample size", 
       y = "proportion estimate bias",
       color = "Age",
       linetype="Age") +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none")+
  facet_grid(factor(age)~.)
pbias
pmoe<-ggplot(scale_stats, aes(x = scales_n,
                                     y = mean_prop_moe,
                                     color = factor(age),
                                     linetype= factor(age))) +
  geom_point(aes(color=factor(age)))+
  geom_errorbar(data=scale_stats,
                aes(ymin = lower_CI_prop_moe, 
                    ymax = upper_CI_prop_moe,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2) +
  labs(title="MOE",
       x = "scale sample size", 
       y = "proportion estimate MOE",
       color = "Age",
       linetype="Age") +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none")+
  facet_grid(factor(age)~.)
pmoe
pmeetsMOE<-ggplot(scale_stats, aes(x = scales_n,
                              y = pct_meeting_target_prop_moe,
                              color = factor(age),
                              linetype= factor(age))) +
  geom_point(aes(color=factor(age)))+
  geom_line(aes(color=factor(age),linetype = factor(age)))+
  geom_vline(xintercept = threshold_value, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title="Meeting target MOE of 20%",
       x = "scale sample size", 
       y = "% meeting target MOE",
       color = "Age",
       linetype="Age") +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank())+
  facet_grid(factor(age)~.)
pmeetsMOE

coverage_sims<-data.frame()
for(i in 1:length(scale_scenario)){
  d<-scale_scenario[[i]]
  scale_sizes<-data.frame()
  for(s in 1:length(d)){
    ds<-d[[s]]$scale_results$est_summary_stats
    ds$scales_collected=d[[s]]$scale_results$scales_collected
    ds$true_proportion_natural=d[[s]]$scale_results$true_summary_stats$true_proportion_natural
    ds$true_count_natural=d[[s]]$scale_results$true_summary_stats$true_count_natural
    ds$scales_n=as.numeric(names(d)[[s]])
    scale_sizes<-scale_sizes%>%rbind(ds)
  }
  scale_sizes$iter=i
  coverage_sims<-coverage_sims%>%rbind(scale_sizes)
}
coverage_sims<-coverage_sims%>%filter(scales_collected<=1500)

coverage_sims<-coverage_sims%>%
  mutate(coverage_count=ifelse(true_count_natural>=lower_CI_count &
                                 true_count_natural<=upper_CI_count,1,0),
         coverage_proportion=ifelse(true_proportion_natural>=lower_CI_prop &
                                      true_proportion_natural<=upper_CI_prop,1,0))

coverage_stats<-coverage_sims%>%
  group_by(scales_collected,age)%>%
  summarise(pct_count_coverage=mean(coverage_count),
            pct_prop_coverage=mean(coverage_proportion))

pcover<-ggplot(coverage_stats, aes(x = scales_collected,
                               y = pct_prop_coverage,
                               color = factor(age),
                               linetype= factor(age))) +
  geom_line() +
  geom_point() +
  scale_y_continuous(limits=c(.5,1),
                     breaks=seq(from=.5,to=1,by=.10))+
  scale_x_continuous(limits=c(min(coverage_stats$scales_collected),
                              max(coverage_stats$scales_collected)),
                     breaks=seq(from=0,to=max(coverage_stats$scales_collected),by=500))+
  labs(subtitle ="Proportional estimates",
       x = "scales collected", 
       y = "proportion estimate coverage",
       color = "Age",
       linetype="Age") +
  theme_bw()+facet_grid(factor(age)~.)+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank())
pcover

p<-grid.arrange(pmoe,pbias,pcover, nrow = 1)
ggsave(p,file="outputs/approach1_scalen_metrics.png",width=10,height=6)

ggsave(pmeetsMOE,file="outputs/approach1_scalen_pmeetsMOE.png",width=6,height=6)

#########################################
#approach 1: theta
########################################

theta_scenario<-readRDS("outputs/theta_sim_results.Rds")
theta_sims<-data.frame()
for(i in 1:length(theta_scenario)){
  d<-theta_scenario[[i]]
  theta_sizes<-data.frame()
  for(s in 1:length(d)){
    ds<-d[[s]]$sim_stats
    ds$scales_collected=d[[s]]$scale_results$scales_collected
    ds$theta=as.numeric(names(d)[[s]])
    theta_sizes<-theta_sizes%>%rbind(ds)
  }
  theta_sizes$iter=i
  theta_sims<-theta_sims%>%rbind(theta_sizes)
}

target_moe=.20
CI=.95

theta_sims<-theta_sims%>%
  mutate(meets_target_moe_prop=ifelse(prop_moe<=target_moe,T,F),
         meets_target_moe_count=ifelse(count_moe_relative<=target_moe,T,F))

theta_stats<-theta_sims%>%
  group_by(age,theta)%>%
  summarize(
    mean_scales_collected=mean(scales_collected),
    mean_prop_bias = mean(prop_bias, na.rm = TRUE),
    mean_count_bias = mean(count_bias, na.rm = TRUE),
    lower_CI_prop_bias = quantile(prop_bias, probs = (1-CI)/2, 
                                  na.rm = TRUE),
    upper_CI_prop_bias = quantile(prop_bias, probs = 1-(1-CI)/2, 
                                  na.rm = TRUE),
    lower_CI_count_bias = quantile(count_bias, 
                                   probs = (1-CI)/2, na.rm = TRUE),
    upper_CI_count_bias = quantile(count_bias, 
                                   probs = 1-(1-CI)/2, na.rm = TRUE),
    mean_prop_moe=mean(prop_moe, na.rm = TRUE),
    lower_CI_prop_moe = quantile(prop_moe, probs = (1-CI)/2, 
                                  na.rm = TRUE),
    upper_CI_prop_moe = quantile(prop_moe, probs = 1-(1-CI)/2, 
                                  na.rm = TRUE),
    mean_count_moe=mean(count_moe_relative,na.rm=TRUE),
    sd_prop_moe=sd(prop_moe, na.rm = TRUE),
    sd_count_moe=sd(count_moe_relative,na.rm=TRUE),
    pct_meeting_target_prop_moe = mean(meets_target_moe_prop) * 100,
    pct_meeting_target_count_moe = mean(meets_target_moe_count) * 100,
    .groups = "drop"
  )

#find smallest theta where 95% meets target prop moe
threshold_n <- theta_stats %>%
  group_by(age) %>%
  filter(pct_meeting_target_prop_moe >= 95) %>%
  slice_min(theta) %>%  # Get the smallest scales_n that meets the threshold
  ungroup()
threshold_value<-max(threshold_n$theta)

pmoe<-ggplot(data=theta_stats, aes(x = factor(theta),
                              y = mean_prop_moe,
                              color = factor(age),
                              linetype= factor(age))) +
  geom_point(aes(color=factor(age)))+
  geom_errorbar(data=theta_stats,
                aes(ymin = lower_CI_prop_moe, 
                    ymax = upper_CI_prop_moe,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2) +
  labs(title="MOE",
       x = "theta", 
       y = "proportion estimate MOE",
       color = "Age",
       linetype="Age") +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none")+
  facet_grid(factor(age)~.)
pmoe

pbias<-ggplot(theta_stats, aes(x = factor(theta),
                               y = mean_prop_bias,
                               color = factor(age),
                               linetype= factor(age))) +
  geom_point(aes(color=factor(age)))+
  geom_errorbar(aes(ymin = lower_CI_prop_bias, 
                    ymax = upper_CI_prop_bias,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2) +
  labs(title="Bias",
       x = "theta", 
       y = "proportion estimate bias",
       color = "Age",
       linetype="Age") +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none")+
  facet_grid(factor(age)~.)
pbias

pmeetsMOE<-ggplot(theta_stats, aes(x = theta,
                        y = pct_meeting_target_prop_moe,
                        color = factor(age),
                        linetype= factor(age))) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = threshold_value, linetype = "dashed", color = "red", linewidth = 1) +
  scale_y_continuous(limits=c(0,100),
                     breaks=seq(from=0,to=100,by=10))+
  scale_x_continuous(limits=c(min(theta_stats$theta),max(theta_stats$theta)),
                     breaks=seq(from=0,to=max(theta_stats$theta),by=.1))+
  labs(title="Meeting target MOE of 20%",
       x = "theta", 
       y = "% iterations meeting target MOE",
       color = "Age",
       linetype="Age") +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank())+
  facet_grid(factor(age)~.)

pmeetsMOE

p<-grid.arrange(pmoe,pbias,pmeetsMOE, nrow = 1)
ggsave(p,file="outputs/approach1_theta_metrics.png",width=10,height=6)

#########################################
#approach 1: tag rates
########################################

tagrate_scenario<-readRDS("outputs/tagrate_sim_results.Rds")
tagrate_sims<-data.frame()
for(i in 1:length(tagrate_scenario)){
  d<-tagrate_scenario[[i]]
  tagrate_sizes<-data.frame()
  for(s in 1:length(d)){
    ds<-d[[s]]$sim_stats
    ds$scales_collected=d[[s]]$scale_results$scales_collected
    ds$tagrate=as.numeric(names(d)[[s]])
    tagrate_sizes<-tagrate_sizes%>%rbind(ds)
  }
  tagrate_sizes$iter=i
  tagrate_sims<-tagrate_sims%>%rbind(tagrate_sizes)
}

target_moe=.20
CI=.95

tagrate_sims<-tagrate_sims%>%
  mutate(meets_target_moe_prop=ifelse(prop_moe<=target_moe,T,F),
         meets_target_moe_count=ifelse(count_moe_relative<=target_moe,T,F))

tagrate_stats<-tagrate_sims%>%
  group_by(age,tagrate)%>%
  summarize(
    mean_scales_collected=mean(scales_collected),
    mean_prop_bias = mean(prop_bias, na.rm = TRUE),
    mean_count_bias = mean(count_bias, na.rm = TRUE),
    lower_CI_prop_bias = quantile(prop_bias, probs = (1-CI)/2, 
                                  na.rm = TRUE),
    upper_CI_prop_bias = quantile(prop_bias, probs = 1-(1-CI)/2, 
                                  na.rm = TRUE),
    lower_CI_count_bias = quantile(count_bias, 
                                   probs = (1-CI)/2, na.rm = TRUE),
    upper_CI_count_bias = quantile(count_bias, 
                                   probs = 1-(1-CI)/2, na.rm = TRUE),
    mean_prop_moe=mean(prop_moe, na.rm = TRUE),
    lower_CI_prop_moe = quantile(prop_moe, probs = (1-CI)/2, 
                                 na.rm = TRUE),
    upper_CI_prop_moe = quantile(prop_moe, probs = 1-(1-CI)/2, 
                                 na.rm = TRUE),
    mean_count_moe=mean(count_moe_relative,na.rm=TRUE),
    sd_prop_moe=sd(prop_moe, na.rm = TRUE),
    sd_count_moe=sd(count_moe_relative,na.rm=TRUE),
    pct_meeting_target_prop_moe = mean(meets_target_moe_prop) * 100,
    pct_meeting_target_count_moe = mean(meets_target_moe_count) * 100,
    .groups = "drop"
  )

#find smallest theta where 95% meets target prop moe
threshold_n <- tagrate_stats %>%
  group_by(age) %>%
  filter(pct_meeting_target_prop_moe >= 95) %>%
  slice_min(tagrate) %>%  # Get the smallest scales_n that meets the threshold
  ungroup()
threshold_value<-max(threshold_n$tagrate)

pmoe<-ggplot(data=tagrate_stats, aes(x = factor(tagrate),
                                   y = mean_prop_moe,
                                   color = factor(age),
                                   linetype= factor(age))) +
  geom_point(aes(color=factor(age)))+
  geom_errorbar(aes(ymin = lower_CI_prop_moe, 
                    ymax = upper_CI_prop_moe,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2) +
  labs(title="MOE",
       x = "tag rate", 
       y = "proportion estimate MOE",
       color = "Age",
       linetype="Age") +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none")+
  facet_grid(factor(age)~.)
pmoe

pbias<-ggplot(tagrate_stats, aes(x = factor(tagrate),
                               y = mean_prop_bias,
                               color = factor(age),
                               linetype= factor(age))) +
  geom_point(aes(color=factor(age)))+
  geom_errorbar(aes(ymin = lower_CI_prop_bias, 
                    ymax = upper_CI_prop_bias,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2) +
  labs(title="Bias",
       x = "tag rate", 
       y = "proportion estimate bias",
       color = "Age",
       linetype="Age") +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none")+
  facet_grid(factor(age)~.)
pbias

pmeetsMOE<-ggplot(data=tagrate_stats, aes(x = (tagrate),
                                   y = pct_meeting_target_prop_moe,
                                   color = factor(age),
                                   linetype= factor(age))) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = threshold_value, linetype = "dashed", color = "red", linewidth = 1) +
  scale_x_continuous(limits=c(.25,1),
                     breaks=seq(from=.25,to=1,by=.25))+
  labs(title="Meeting target MOE of 20%",
       x = "tag rate", 
       y = "% iterations meeting target MOE",
       color = "Age",
       linetype="Age") +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank())+
  facet_grid(factor(age)~.)
pmeetsMOE

p<-grid.arrange(pmoe,pbias,pmeetsMOE, nrow = 1)
ggsave(p,file="outputs/approach1_tagrate_metrics.png",width=10,height=6)


tagrate_coverage<-data.frame()
for(i in 1:length(tagrate_scenario)){
  d<-tagrate_scenario[[i]]
  tagrate_sizes<-data.frame()
  for(s in 1:length(d)){
    ds<-d[[s]]$scale_results$est_summary_stats
    ds$scales_collected=d[[s]]$scale_results$scales_collected
    ds$true_proportion_natural=d[[s]]$scale_results$true_summary_stats$true_proportion_natural
    ds$true_count_natural=d[[s]]$scale_results$true_summary_stats$true_count_natural
    ds$tagrate=as.numeric(names(d)[[s]])
    tagrate_sizes<-tagrate_sizes%>%rbind(ds)
  }
  tagrate_sizes$iter=i
  tagrate_coverage<-tagrate_coverage%>%rbind(tagrate_sizes)
}
#coverage_sims<-coverage_sims%>%filter(scales_collected<=1500)

tagrate_coverage<-tagrate_coverage%>%
  mutate(coverage_count=ifelse(true_count_natural>=lower_CI_count &
                                 true_count_natural<=upper_CI_count,1,0),
         coverage_proportion=ifelse(true_proportion_natural>=lower_CI_prop &
                                      true_proportion_natural<=upper_CI_prop,1,0))

tagrate_stats<-tagrate_coverage%>%
  group_by(tagrate,age)%>%
  summarise(pct_count_coverage=mean(coverage_count),
            pct_prop_coverage=mean(coverage_proportion))
p1<-ggplot(tagrate_stats, aes(x = tagrate,
                               y = pct_prop_coverage*100,
                               color = factor(age),
                               linetype= factor(age))) +
  geom_line() +
  geom_point() +
  scale_y_continuous(limits=c(50,100),
                     breaks=seq(from=50,to=100,by=10))+
  scale_x_continuous(limits=c(.25,1),
                     breaks=seq(from=.25,to=1,by=.25))+
  labs(subtitle ="Proportional estimates",
       x = "scales collected", 
       y = "% iterations") +
  theme_bw()
p1

#########################################
#approach 2: No pooling
########################################
nopooling<-readRDS("outputs/basin_results_nopooling.Rds")
target_moe=0.2
nopooling_estimates<-data.frame()
for(i in 1:length(nopooling)){
  d<-nopooling[[i]]
  for(s in 1:length(d$scales_results)){
    ds<-d$scales_results[[s]]
    ds$iteration=i
    nopooling_estimates<-nopooling_estimates%>%rbind(ds)
  }
}

nopooling_totals<-nopooling_estimates%>%
  group_by(iteration,target_location_scales)%>%
  summarise(total_age_est=sum(mean_natural_count),
            total_age_true=sum(true_abundance))

nopooling_estimates<-nopooling_estimates%>%
  left_join(nopooling_totals)

nopooling_estimates<-nopooling_estimates%>%
  mutate(mean_natural_prop=mean_natural_count/total_age_est,
         true_natural_prop=true_abundance/total_age_true,
         lower_CI_prop=lower_ci/total_age_est,
         upper_CI_prop=upper_ci/total_age_est,
         prop_moe=(upper_CI_prop-lower_CI_prop)/2,
         prob_bias=mean_natural_prop-true_natural_prop)

nopooling_estimates<-nopooling_estimates%>%
  mutate(coverage_prop=ifelse(true_natural_prop>=lower_CI_prop &
                                true_natural_prop<=upper_CI_prop,1,0),
         meets_target_moe_prop=ifelse(prop_moe<=target_moe,1,0))

nopooling_stats<-nopooling_estimates%>%
  group_by(age,target_location_scales)%>%
  summarize(
    mean_scales_collected=mean(scales_collected),
    lower_CI_prop_moe=quantile(prop_moe, probs = (1-CI)/2, 
                                na.rm = TRUE),
    upper_CI_prop_moe=quantile(prop_moe, probs = 1-(1-CI)/2, 
                                na.rm = TRUE),
    mean_natural_prop=mean(prop_moe,na.rm=TRUE),
    pct_meeting_target_moe_prop = mean(meets_target_moe_prop) * 100,
    pct_prop_coverage=mean(coverage_prop),
    mean_prop_bias=mean(prob_bias,na.rm=TRUE),
    lower_CI_prop_bias=quantile(prob_bias, probs = (1-CI)/2, 
                                 na.rm = TRUE),
    upper_CI_prop_bias=quantile(prob_bias, probs = 1-(1-CI)/2, 
                                 na.rm = TRUE),
    .groups = "drop"
  )


pmoe<-ggplot(data=nopooling_stats, aes(x = (target_location_scales),
                                   y = mean_natural_prop,
                                   color = factor(age),
                                   linetype= factor(age))) +
  geom_point(aes(color=factor(age)))+
  geom_errorbar(data=nopooling_stats,
                aes(ymin = lower_CI_prop_moe, 
                    ymax = upper_CI_prop_moe,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2) +
  scale_x_continuous(limits=c(min(nopooling_stats$target_location_scales),
                              max(nopooling_stats$target_location_scales)),
                     breaks=seq(from=0,
                                to=max(nopooling_stats$target_location_scales),
                                by=200))+
  labs(title="MOE",
       x = "location scale sample size", 
       y = "proportion estimate MOE",
       color = "Age",
       linetype="Age") +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none")+
  facet_grid(factor(age)~.)
pmoe

pbias<-ggplot(nopooling_stats, aes(x = (target_location_scales),
                               y = mean_prop_bias,
                               color = factor(age),
                               linetype= factor(age))) +
  geom_point(aes(color=factor(age)))+
  geom_errorbar(aes(ymin = lower_CI_prop_bias, 
                    ymax = upper_CI_prop_bias,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2) +
  labs(title="Bias",
       x = "location scale sample size", 
       y = "proportion estimate bias",
       color = "Age",
       linetype="Age") +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none")+
  facet_grid(factor(age)~.)
pbias

pcover<-ggplot(nopooling_stats, aes(x = target_location_scales,
                                   y = pct_prop_coverage*100,
                                   color = factor(age),
                                   linetype= factor(age))) +
  geom_line() +
  geom_point() +
    scale_x_continuous(limits=c(min(nopooling_stats$target_location_scales),
                              max(nopooling_stats$target_location_scales)),
                     breaks=seq(from=0,
                                to=max(nopooling_stats$target_location_scales),
                                by=200))+
  scale_y_continuous(limits=c(0,100),
                     breaks=seq(from=0,to=100,by=20))+
  labs(title="proportion estimate coverage",
       x = "location scale sample size", 
       y = "% iterations",
       color = "Age",
       linetype="Age") +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank())+
  facet_grid(factor(age)~.)

pcover

p<-grid.arrange(pmoe,pbias,pcover, nrow = 1)
ggsave(p,file="outputs/approach2_nopooling_metrics.png",width=10,height=6)

#########################################
#approach 2: No pooling location breakdown
########################################
trib_data <- readRDS("data/all_returns_clean.Rds")
theta_data <- readRDS("data/location_year_theta.Rds")

target_year <- 2018

#draw population data
pop <- sim_basin_pop(trib_data,
                     iter = 10,
                     year = target_year)$iter_pop
#create list of survey locations
loc_list <- unique(pop$location)
loc_list

true_natural<-pop%>%
  filter(origin=="natural")
true_pop_ages <- pop%>%
  filter(origin == "natural")%>%
  group_by(age)%>%
  summarise(true_proportion = n()/nrow(true_natural))

recoveries_list <- list()
for(l in 1:length(loc_list)){ #loop recoveries sim for each location
  iter_loc <- loc_list[l]
  loc_pop <- pop%>%filter(location == iter_loc) #get pop data for location
  
  target_theta <- theta_data%>% #get theta for location x year
    filter(location == iter_loc,
           return_year == target_year)
  
  if(nrow(target_theta) == 0)
  {target_theta = 0.2} #if no theta, set default
  else 
  {target_theta <- round(target_theta$theta,2)}
  
  #get spawning recoveries
  spawning_recoveries <- sim_spawning_recovery(loc_pop,
                                               theta = target_theta,
                                               iterations = 1000)
  recoveries_list[[l]] <- spawning_recoveries
}
names(recoveries_list) <- loc_list

location_results <- list()
scales_collected <- list()
for(l in 1:length(loc_list)){
  target_loc <- loc_list[[l]]
  loc_pop <- pop %>% filter(location == target_loc)
  
  d <- recoveries_list[[l]]
  
  sim_result <- sim_scales_sampling(
    pop = loc_pop,
    spawning_results = d,
    iterations = 1000,
    scales_n = 400,
    CI = 0.95
  )
  
  scales_collected[[l]] <- sim_result$scales_collected
  l_results <- list_rbind(sim_result$boot_results) %>%
    select(k_iteration, age, total_natural) %>%
    pivot_wider(
      id_cols = k_iteration,
      names_from = age,
      values_from = total_natural
    )
  location_results[[l]] <- l_results
}

names(location_results) <- names(scales_collected) <- loc_list

basin_results <- bind_rows(location_results, .id = "location") %>%
  group_by(k_iteration) %>%
  summarise(
    age_2 = sum(`2`, na.rm = TRUE),
    age_3 = sum(`3`, na.rm = TRUE),
    age_4 = sum(`4`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ungroup()

CI = 0.95
basin_summary <- basin_results %>%
  pivot_longer(
    cols = starts_with("age_"),
    names_to = "age",
    values_to = "value",
    names_prefix = "age_"
  ) 
basin_totals<-basin_summary%>%
  group_by(k_iteration)%>%
  summarise(total=sum(value))

basin_summary<-basin_summary%>%
  left_join(basin_totals)%>%
  mutate(prop=value/total)

basin_stats<-basin_summary%>%
  group_by(age) %>%
  summarise(
    mean_natural_prop = mean(prop, na.rm = TRUE),
    lower_CI_prop = quantile(prop, (1 - CI) / 2, na.rm = TRUE),
    upper_CI_prop = quantile(prop, 1 - (1 - CI) / 2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ungroup()

basin_stats<-basin_stats%>%
  mutate(age=as.numeric(age))%>%
  left_join(true_pop_ages)

pbasin<-ggplot(data=basin_stats,
               aes(x=factor(age),y=mean_natural_prop))+
  geom_point( size=3, aes(color=factor(age)))+
  geom_errorbar(aes(ymin = lower_CI_prop, ymax = upper_CI_prop,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2) +
  geom_point(aes(x=factor(age),y=true_proportion),shape=17)+
  theme_bw()+
  labs(title=target_year,
       subtitle=paste("scales collected:",sum(unlist(scales_collected))),
       x = "Age", y = "Est proportion")+
  guides(color = "none", fill = "none", linetype = "none")
pbasin

ggsave(pbasin,file="outputs/approach1_nopooling_2018iter10.png",width=10,height=6)


location_summary <- bind_rows(location_results, .id = "location") %>%
  group_by(k_iteration,location) %>%
  summarise(
    age_2 = sum(`2`, na.rm = TRUE),
    age_3 = sum(`3`, na.rm = TRUE),
    age_4 = sum(`4`, na.rm = TRUE)
  ) %>%
  ungroup()

location_summary_long <- location_summary %>%
  pivot_longer(
    cols = starts_with("age_"),
    names_to = "age",
    values_to = "count",
    names_prefix = "age_"
  ) %>%
  mutate(age = as.integer(age))

location_iter_totals<-location_summary_long%>%
  group_by(location,k_iteration)%>%
  summarise(total=sum(count))
location_summary_long<-location_summary_long%>%
  left_join(location_iter_totals)%>%
  mutate(prop=count/total)

summary_by_location_age <- location_summary_long %>%
  group_by(location, age) %>%
  summarise(
    mean_natural_prop = mean(prop, na.rm = TRUE),
    lower_CI_prop = quantile(prop, 0.025, na.rm = TRUE),
    upper_CI_prop = quantile(prop, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

true_loc_abundance<-pop%>%
  filter(origin=="natural")%>%
  group_by(age,location)%>%
  summarize(true_abundance=n())
true_loc_totals<-pop%>%
  filter(origin=="natural")%>%
  group_by(location)%>%
  summarize(total_abundance=n())
true_loc_props<-true_loc_abundance%>%left_join(true_loc_totals)%>%
  mutate(true_prop=true_abundance/total_abundance)

p1<-ggplot(data=summary_by_location_age%>%filter(location=="Other"),
           aes(x=factor(age),y=mean_natural_prop))+
  geom_point( size=3, aes(color=factor(age)))+
  geom_errorbar(aes(ymin = lower_CI_prop, ymax = upper_CI_prop,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2) +
  geom_point(data=true_loc_props%>%filter(location=="Other"),aes(x=factor(age),y=true_prop),shape=17)+
  theme_bw()+
  labs(title = paste(target_year,", Other, theta: ",
                     round(theta_data %>% 
                             filter(location == "Other",
                                    return_year==target_year) %>% 
                             pull(theta),2)),
       subtitle=paste("scales collected: ",
                      scales_collected[["Other"]]),
       x = "Age", y = "estimated proportion")+
  guides(color = "none", fill = "none", linetype = "none")
p1

p2<-ggplot(data=summary_by_location_age%>%filter(location=="NFH"),
           aes(x=factor(age),y=mean_natural_prop))+
  geom_point( size=3, aes(color=factor(age)))+
  geom_errorbar(aes(ymin = lower_CI_prop, ymax = upper_CI_prop,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2) +
  geom_point(data=true_loc_props%>%filter(location=="NFH"),aes(x=factor(age),y=true_prop),shape=17)+
  theme_bw()+
  labs(title = paste(target_year,", NFH, theta: ",
                     round(theta_data %>% 
                             filter(location == "NFH",
                                    return_year==target_year) %>% 
                             pull(theta),2)),
       subtitle=paste("scales collected: ",
                      scales_collected[["NFH"]]),
       x = "Age", y = "estimated proportion")+
  guides(color = "none", fill = "none", linetype = "none")
p2

p3<-ggplot(data=summary_by_location_age%>%filter(location=="Feather"),
           aes(x=factor(age),y=mean_natural_prop))+
  geom_point( size=3, aes(color=factor(age)))+
  geom_errorbar(aes(ymin = lower_CI_prop, ymax = upper_CI_prop,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2) +
  geom_point(data=true_loc_props%>%filter(location=="Feather"),aes(x=factor(age),y=true_prop),shape=17)+
  theme_bw()+
  labs(title = paste(target_year,", Feather, theta: ",
                     round(theta_data %>% 
                             filter(location == "Feather",
                                    return_year==target_year) %>% 
                             pull(theta),2)),
       subtitle=paste("scales collected: ",
                      scales_collected[["Feather"]]),
       x = "Age", y = "estimated proportion")+
  guides(color = "none", fill = "none", linetype = "none")
p3

p<-grid.arrange(p1,p2,p3, nrow = 1)
ggsave(p,file="outputs/approach1_nopooling_2018iter10_locations.png",width=10,height=6)

#examining feather river
feather_sim_result <- sim_scales_sampling(
  pop = pop%>%filter(location=="Feather"),
  spawning_results = recoveries_list[[4]],
  iterations = 1000,
  scales_n = 400,
  CI = 0.95
)

feather_ex<-feather_sim_result$boot_results[[1]]
feather_hatch<-recoveries_list$Feather$hatchery_ages%>%filter(k_iteration==1)


#########################################
#approach 2: recovery proportions
#########################################
location_recoveries<-data.frame()
for(l in 1:length(loc_list)){
  d<-nrow(recoveries_list[[loc_list[[l]]]]$recovered_fish)
  loc_pop=pop%>%filter(location==loc_list[[l]])
  loc_esc=nrow(loc_pop)
  esc_ratio=loc_esc/nrow(pop)
  result=data.frame("recoveries"=d,
                    "location"=loc_list[[l]],
                    "esc_ratio"=esc_ratio)
  location_recoveries<-location_recoveries%>%rbind(result)
}
total_recoveries<-sum(location_recoveries$recoveries)
location_recoveries<-location_recoveries%>%
  mutate(recovery_ratio=recoveries/total_recoveries)

p<-ggplot(data=location_recoveries)+
  geom_point(size=3, aes(x=location, y=esc_ratio, color="Escapement", shape="Escapement"))+
  geom_point(size=3, aes(x=location, y=recovery_ratio, color="Recoveries", shape="Recoveries"))+
  theme_bw()+
  labs(x = "locations", 
       y = "Proportion",
       color = "Type",
       shape = "Type")+
  scale_color_manual(values = c("Escapement" = "red", "Recoveries" = "black"))+
  scale_shape_manual(values = c("Escapement" = 16, "Recoveries" = 17))+
  theme(legend.position = "right",  # or "bottom", "top", "left"
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "white", color = "black"))

ggsave(p,file="outputs/approach2_recoveries_props.png",width=8,height=6)

#################################
#approach 2:with pooling
#####################################
pooling_iterations<-readRDS("outputs/basin_results_wpooling_lowsamples.Rds")
target_moe=.10
pooling_estimates<-data.frame()
for(i in 1:length(pooling_iterations)){
  d<-pooling_iterations[[i]]
  for(s in 1:length(d$scale_iters)){
    ds<-d$scale_iters[[s]]
    ds$est_summary_stats$scales_collected=ds$target_scales_n
    ds$est_summary_stats$true_count_natural=ds$true_summary_stats$true_count_natural
    ds$est_summary_stats$true_proportion_natural=ds$true_summary_stats$true_proportion_natural
    pooling_estimates<-pooling_estimates%>%rbind(ds$est_summary_stats)
  }
}

pooling_estimates<-pooling_estimates%>%
  mutate(moe_count=(upper_CI_count-lower_CI_count)/2,
         moe_count_relative=ifelse(true_count_natural == 0, NA, 
                                   moe_count / true_count_natural),
         moe_prop=(upper_CI_prop-lower_CI_prop)/2,
         meets_target_moe_prop=ifelse(moe_prop<=target_moe,T,F),
         meets_target_moe_count=ifelse(moe_count_relative<=target_moe,T,F),
         coverage_count=ifelse(true_count_natural>=lower_CI_count &
                                 true_count_natural<=upper_CI_count,1,0),
         coverage_proportion=ifelse(true_proportion_natural>=lower_CI_prop &
                                      true_proportion_natural<=upper_CI_prop,1,0),
         count_bias=mean_count_natural-true_count_natural,
         rel_count_bias=count_bias/true_count_natural,
         prop_bias=mean_proportion_natural-true_proportion_natural
         )


pooling_stats<-pooling_estimates%>%
  group_by(age,scales_collected)%>%
  summarize(
    mean_scales_collected=mean(scales_collected),
    mean_prop_moe=mean(moe_prop, na.rm = TRUE),
    lower_CI_prop_moe=quantile(moe_prop, probs = (1-CI)/2, 
                               na.rm = TRUE),
    upper_CI_prop_moe=quantile(moe_prop, probs = 1-(1-CI)/2, 
                               na.rm = TRUE),
    lower_CI_count_moe=quantile(moe_count_relative, probs = (1-CI)/2, 
                                na.rm = TRUE),
    upper_CI_count_moe=quantile(moe_count_relative, probs = 1-(1-CI)/2, 
                                na.rm = TRUE),
    mean_count_moe=mean(moe_count_relative,na.rm=TRUE),
    pct_meeting_target_prop_moe = mean(meets_target_moe_prop) * 100,
    pct_meeting_target_count_moe = mean(meets_target_moe_count) * 100,
    pct_count_coverage=mean(coverage_count),
    pct_prop_coverage=mean(coverage_proportion),
    mean_count_bias=mean(rel_count_bias,na.rm=TRUE),
    mean_prop_bias=mean(prop_bias,na.rm=T),
    lower_CI_prop_bias=quantile(prop_bias, probs = (1-CI)/2, 
                                na.rm = TRUE),
    upper_CI_prop_bias=quantile(prop_bias, probs = 1-(1-CI)/2, 
                                na.rm = TRUE),
    lower_CI_count_bias=quantile(rel_count_bias, probs = (1-CI)/2, 
                                 na.rm = TRUE),
    upper_CI_count_bias=quantile(rel_count_bias, probs = 1-(1-CI)/2, 
                                 na.rm = TRUE),
    .groups = "drop"
  )

p1<-ggplot(data=pooling_stats,
           aes(x=(scales_collected),y=mean_prop_moe))+
  geom_point( size=3, aes(color=factor(age)))+
  geom_errorbar(data=pooling_stats,
                aes(ymin = lower_CI_prop_moe, 
                    ymax = upper_CI_prop_moe,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2) +
  scale_x_continuous(limits=c(100,1000),
                     breaks=seq(from=100,to=1000,by=200))+
  theme_bw()+
  labs(x = "basin scale samples", y = "proportion estimate MOE",color="Age",linetype="Age")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none")+
  facet_grid(factor(age)~.)
p1
p2<-ggplot(data=pooling_stats,
           aes(x=(scales_collected),y=mean_prop_bias))+
  geom_point( size=3, aes(color=factor(age)))+
  geom_errorbar(data=pooling_stats,
                aes(ymin = lower_CI_prop_bias, 
                    ymax = upper_CI_prop_bias,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2) +
  theme_bw()+
  scale_x_continuous(limits=c(100,1000),
                     breaks=seq(from=100,to=1000,by=200))+
  labs(x = "basin scale samples", y = "proportion estimate bias",color="Age",linetype="Age")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none")+
  facet_grid(factor(age)~.)
p2
p3<-ggplot(data=pooling_stats,
           aes(x=(mean_scales_collected),y=pct_prop_coverage*100,
               color = factor(age),
               linetype= factor(age)))+
  geom_point()+
  geom_line()+
  theme_bw()+
  scale_y_continuous(limits=c(0,100),
                     breaks=seq(from=0,to=100,by=20))+
  labs(x = "basin scale samples", y = "proportion estimate coverage",color="Age",linetype="Age")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank())+
  scale_x_continuous(limits=c(100,1000),
                     breaks=seq(from=100,to=1000,by=200))+
  facet_grid(factor(age)~.)
p3
p<-grid.arrange(p1,p2,p3, nrow = 1)
ggsave(p,file="outputs/approach2_wpooling_metrics.png",width=10,height=6)

#################################
#approach 2:with pooling no weighting
#####################################
pooling_iterations<-readRDS("outputs/basin_results_wpooling_noweighting.Rds")

pooling_estimates<-data.frame()
for(i in 1:length(pooling_iterations)){
  d<-pooling_iterations[[i]]
  for(s in 1:length(d$scale_iters)){
    ds<-d$scale_iters[[s]]
    ds$est_summary_stats$scales_collected=ds$target_scales_n
    ds$est_summary_stats$true_count_natural=ds$true_summary_stats$true_count_natural
    ds$est_summary_stats$true_proportion_natural=ds$true_summary_stats$true_proportion_natural
    pooling_estimates<-pooling_estimates%>%rbind(ds$est_summary_stats)
  }
}

pooling_estimates<-pooling_estimates%>%
  mutate(moe_count=(upper_CI_count-lower_CI_count)/2,
         moe_count_relative=ifelse(true_count_natural == 0, NA, 
                                   moe_count / true_count_natural),
         moe_prop=(upper_CI_prop-lower_CI_prop)/2,
         meets_target_moe_prop=ifelse(moe_prop<=target_moe,T,F),
         meets_target_moe_count=ifelse(moe_count_relative<=target_moe,T,F),
         coverage_count=ifelse(true_count_natural>=lower_CI_count &
                                 true_count_natural<=upper_CI_count,1,0),
         coverage_proportion=ifelse(true_proportion_natural>=lower_CI_prop &
                                      true_proportion_natural<=upper_CI_prop,1,0),
         count_bias=mean_count_natural-true_count_natural,
         rel_count_bias=count_bias/true_count_natural,
         prop_bias=mean_proportion_natural-true_proportion_natural
  )


pooling_stats<-pooling_estimates%>%
  group_by(age,scales_collected)%>%
  summarize(
    mean_scales_collected=mean(scales_collected),
    mean_prop_moe=mean(moe_prop, na.rm = TRUE),
    lower_CI_prop_moe=quantile(moe_prop, probs = (1-CI)/2, 
                               na.rm = TRUE),
    upper_CI_prop_moe=quantile(moe_prop, probs = 1-(1-CI)/2, 
                               na.rm = TRUE),
    lower_CI_count_moe=quantile(moe_count_relative, probs = (1-CI)/2, 
                                na.rm = TRUE),
    upper_CI_count_moe=quantile(moe_count_relative, probs = 1-(1-CI)/2, 
                                na.rm = TRUE),
    mean_count_moe=mean(moe_count_relative,na.rm=TRUE),
    pct_meeting_target_prop_moe = mean(meets_target_moe_prop) * 100,
    pct_meeting_target_count_moe = mean(meets_target_moe_count) * 100,
    pct_count_coverage=mean(coverage_count),
    pct_prop_coverage=mean(coverage_proportion),
    mean_count_bias=mean(rel_count_bias,na.rm=TRUE),
    mean_prop_bias=mean(prop_bias,na.rm=T),
    lower_CI_prop_bias=quantile(prop_bias, probs = (1-CI)/2, 
                                na.rm = TRUE),
    upper_CI_prop_bias=quantile(prop_bias, probs = 1-(1-CI)/2, 
                                na.rm = TRUE),
    lower_CI_count_bias=quantile(rel_count_bias, probs = (1-CI)/2, 
                                 na.rm = TRUE),
    upper_CI_count_bias=quantile(rel_count_bias, probs = 1-(1-CI)/2, 
                                 na.rm = TRUE),
    .groups = "drop"
  )

p1<-ggplot(data=pooling_stats,
           aes(x=(scales_collected),y=mean_prop_moe))+
  geom_point( size=3, aes(color=factor(age)))+
  geom_errorbar(data=pooling_stats,
                aes(ymin = lower_CI_prop_moe, 
                    ymax = upper_CI_prop_moe,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2) +
  scale_x_continuous(limits=c(500,2000),
                     breaks=seq(from=500,to=2000,by=200))+
  theme_bw()+
  labs(x = "basin scale samples", y = "proportion estimate MOE",color="Age",linetype="Age")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none")+
  facet_grid(factor(age)~.)
p1
p2<-ggplot(data=pooling_stats,
           aes(x=(scales_collected),y=mean_prop_bias))+
  geom_point( size=3, aes(color=factor(age)))+
  geom_errorbar(data=pooling_stats,
                aes(ymin = lower_CI_prop_bias, 
                    ymax = upper_CI_prop_bias,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2) +
  scale_x_continuous(limits=c(500,2000),
                     breaks=seq(from=500,to=2000,by=200))+
  theme_bw()+
  labs(x = "basin scale samples", y = "proportion estimate bias",color="Age",linetype="Age")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none")+
  facet_grid(factor(age)~.)
p2
p3<-ggplot(data=pooling_stats,
           aes(x=(scales_collected),y=pct_prop_coverage*100,
               color = factor(age),
               linetype= factor(age)))+
  geom_point()+
  geom_line()+
  theme_bw()+
  scale_y_continuous(limits=c(0,100),
                     breaks=seq(from=0,to=100,by=20))+
  labs(x = "basin scale samples", y = "proportion estimate coverage",color="Age",linetype="Age")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank())+
  facet_grid(factor(age)~.)
p3
p<-grid.arrange(p1,p2,p3, nrow = 1)
ggsave(p,file="outputs/approach2_noweight_wpooling_metrics.png",width=10,height=6)


#########################################
#approach 2: with pooling breakdown
########################################
trib_data <- readRDS("data/all_returns_clean.Rds")
theta_data <- readRDS("data/location_year_theta.Rds")

target_year <- 2018

#draw population data
pop <- sim_basin_pop(trib_data,
                     iter = 10,
                     year = target_year)
sim_tag_rates<-pop$tag_rates
pop<-pop$iter_pop
#create list of survey locations
loc_list <- unique(pop$location)
loc_list
iterations=1000

true_natural<-pop%>%
  filter(origin=="natural")
true_pop_ages <- pop%>%
  filter(origin == "natural")%>%
  group_by(age)%>%
  summarise(true_proportion = n()/nrow(true_natural))

recoveries_list <- list()
#spawning recoveries for each location
recoveries_list<-list()
for(l in 1:length(loc_list)){ #loop recoveries sim for each location
  iter_loc<-loc_list[l]
  loc_pop<-pop%>%filter(location==iter_loc) #get pop data for locaiton
  
  target_theta<-theta_data%>% #get theta for location x year
    filter(location==iter_loc,
           return_year==2018)
  if(nrow(target_theta)==0){target_theta=0.2} else{ #if no theta, set default
    target_theta<-round(target_theta$theta,2)
  }
  
  target_tag_rate<-sim_tag_rates%>%
    filter(location==iter_loc,
           year==2018,
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

#get locations proportional escapement
basin_esc<-nrow(pop)
location_esc_prop<-pop%>%
  group_by(location)%>%
  summarize(esc=n())
location_esc_prop<-location_esc_prop%>%
  mutate(prop_esc=esc/basin_esc)

current_scales_n=400

#use round_to_sum to get scales per location
location_esc <- location_esc_prop %>%
  mutate(target_scales = round_to_sum(prop_esc, current_scales_n))

#find amount of scales actually available per location
location_availability <- all_recovered_fish %>%
  filter(has_cwt == 0) %>%
  group_by(location) %>%
  summarise(available = n())

location_esc <- location_esc %>%
  left_join(location_availability, by = "location") %>%
  mutate(
    available = ifelse(is.na(available), 0, available),
    # Initially take min(allocation, available)
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
  scales_n=current_scales_n)

basin_sim_result$target_scales_n<-current_scales_n



pbasin<-ggplot(data=basin_sim_result$est_summary_stats,
               aes(x=factor(age),y=mean_proportion_natural))+
  geom_point( size=3, aes(color=factor(age)))+
  geom_errorbar(aes(ymin = lower_CI_prop, ymax = upper_CI_prop,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2) +
  geom_point(data=basin_sim_result$true_summary_stats,
             aes(x=factor(age),y=true_proportion_natural),shape=17)+
  theme_bw()+
  labs(title=target_year,
       subtitle=paste("scales collected: ",basin_sim_result$scales_collected),
       x = "Age", y = "Est proportion")+
  guides(color = "none", fill = "none", linetype = "none")
pbasin

ggsave(pbasin,file="outputs/approach2_wpooling_2018iter10.png",width=10,height=6)

#########################################
#approach 2: with pooling breakdown no weighting
########################################
trib_data <- readRDS("data/all_returns_clean.Rds")
theta_data <- readRDS("data/location_year_theta.Rds")

target_year <- 2018

#draw population data
pop <- sim_basin_pop(trib_data,
                     iter = 1,
                     year = target_year)
sim_tag_rates<-pop$tag_rates
pop<-pop$iter_pop
#create list of survey locations
loc_list <- unique(pop$location)
loc_list
iterations=1000

true_natural<-pop%>%
  filter(origin=="natural")
true_pop_ages <- pop%>%
  filter(origin == "natural")%>%
  group_by(age)%>%
  summarise(true_proportion = n()/nrow(true_natural))

#spawning recoveries for each location
recoveries_list<-list()
for(l in 1:length(loc_list)){ #loop recoveries sim for each location
  iter_loc<-loc_list[l]
  loc_pop<-pop%>%filter(location==iter_loc) #get pop data for locaiton
  
  target_theta<-theta_data%>% #get theta for location x year
    filter(location==iter_loc,
           return_year==2018)
  if(nrow(target_theta)==0){target_theta=0.2} else{ #if no theta, set default
    target_theta<-round(target_theta$theta,2)
  }
  
  target_tag_rate<-sim_tag_rates%>%
    filter(location==iter_loc,
           year==2018,
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

basin_spawning_results<-list("recovered_fish"=all_recovered_fish,
                             "hatchery_ages"=basin_hatchery_estimates)

current_scales_n=400

basin_sim_result<-sim_scales_sampling(
  pop=pop,
  spawning_results = basin_spawning_results,
  iterations=iterations,
  scales_n=current_scales_n)

pbasin<-ggplot(data=basin_sim_result$est_summary_stats,
               aes(x=factor(age),y=mean_proportion_natural))+
  geom_point( size=3, aes(color=factor(age)))+
  geom_errorbar(aes(ymin = lower_CI_prop, ymax = upper_CI_prop,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2) +
  geom_point(data=basin_sim_result$true_summary_stats,
             aes(x=factor(age),y=true_proportion_natural),shape=17)+
  theme_bw()+
  labs(title=target_year,
       subtitle=paste("scales collected: ",basin_sim_result$scales_collected),
       x = "Age", y = "Est proportion")+
  guides(color = "none", fill = "none", linetype = "none")
pbasin

ggsave(pbasin,file="outputs/approach2_wpooling_2018iter10.png",width=10,height=6)
#################################
#with pooling lowsamples
#####################################
pooling_iterations<-readRDS("outputs/basin_results_wpooling_lowsamples.Rds")

pooling_estimates<-data.frame()
for(i in 1:length(pooling_iterations)){
  d<-pooling_iterations[[i]]
  for(s in 1:length(d$scale_iters)){
    ds<-d$scale_iters[[s]]
    ds$est_summary_stats$scales_collected=ds$target_scales_n
    ds$est_summary_stats$true_count_natural=ds$true_summary_stats$true_count_natural
    ds$est_summary_stats$true_proportion_natural=ds$true_summary_stats$true_proportion_natural
    pooling_estimates<-pooling_estimates%>%rbind(ds$est_summary_stats)
  }
}

pooling_estimates<-pooling_estimates%>%
  mutate(moe_count=(upper_CI_count-lower_CI_count)/2,
         moe_count_relative=ifelse(true_count_natural == 0, NA, 
                                   moe_count / true_count_natural),
         moe_prop=(upper_CI_prop-lower_CI_prop)/2,
         meets_target_moe_prop=ifelse(moe_prop<=target_moe,T,F),
         meets_target_moe_count=ifelse(moe_count_relative<=target_moe,T,F),
         coverage_count=ifelse(true_count_natural>=lower_CI_count &
                                 true_count_natural<=upper_CI_count,1,0),
         coverage_proportion=ifelse(true_proportion_natural>=lower_CI_prop &
                                      true_proportion_natural<=upper_CI_prop,1,0),
         count_bias=mean_count_natural-true_count_natural,
         rel_count_bias=count_bias/true_count_natural,
         prop_bias=mean_proportion_natural-true_proportion_natural
  )


pooling_stats<-pooling_estimates%>%
  group_by(age,scales_collected)%>%
  summarize(
    mean_scales_collected=mean(scales_collected),
    mean_prop_moe=mean(moe_prop, na.rm = TRUE),
    lower_CI_prop_moe=quantile(moe_prop, probs = (1-CI)/2, 
                               na.rm = TRUE),
    upper_CI_prop_moe=quantile(moe_prop, probs = 1-(1-CI)/2, 
                               na.rm = TRUE),
    lower_CI_count_moe=quantile(moe_count_relative, probs = (1-CI)/2, 
                                na.rm = TRUE),
    upper_CI_count_moe=quantile(moe_count_relative, probs = 1-(1-CI)/2, 
                                na.rm = TRUE),
    mean_count_moe=mean(moe_count_relative,na.rm=TRUE),
    pct_meeting_target_prop_moe = mean(meets_target_moe_prop) * 100,
    pct_meeting_target_count_moe = mean(meets_target_moe_count) * 100,
    pct_count_coverage=mean(coverage_count),
    pct_prop_coverage=mean(coverage_proportion),
    mean_count_bias=mean(rel_count_bias,na.rm=TRUE),
    mean_prop_bias=mean(prop_bias,na.rm=T),
    lower_CI_prop_bias=quantile(prop_bias, probs = (1-CI)/2, 
                                na.rm = TRUE),
    upper_CI_prop_bias=quantile(prop_bias, probs = 1-(1-CI)/2, 
                                na.rm = TRUE),
    lower_CI_count_bias=quantile(rel_count_bias, probs = (1-CI)/2, 
                                 na.rm = TRUE),
    upper_CI_count_bias=quantile(rel_count_bias, probs = 1-(1-CI)/2, 
                                 na.rm = TRUE),
    .groups = "drop"
  )

p1<-ggplot(data=pooling_stats,
           aes(x=(scales_collected),y=mean_prop_moe))+
  geom_point( size=3, aes(color=factor(age)))+
  geom_errorbar(data=pooling_stats,
                aes(ymin = lower_CI_prop_moe, 
                    ymax = upper_CI_prop_moe,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2) +
  scale_x_continuous(limits=c(500,2000),
                     breaks=seq(from=500,to=2000,by=200))+
  theme_bw()+
  labs(x = "basin scale samples", y = "proportion MOE",color="Age",linetype="Age")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none")+
  facet_grid(factor(age)~.)
p1
p2<-ggplot(data=pooling_stats,
           aes(x=(scales_collected),y=mean_prop_bias))+
  geom_point( size=3, aes(color=factor(age)))+
  geom_errorbar(data=pooling_stats,
                aes(ymin = lower_CI_prop_bias, 
                    ymax = upper_CI_prop_bias,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2) +
  scale_x_continuous(limits=c(500,2000),
                     breaks=seq(from=500,to=2000,by=200))+
  theme_bw()+
  labs(x = "basin scale samples", y = "proportional bias",color="Age",linetype="Age")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none")+
  facet_grid(factor(age)~.)
p2
p3<-ggplot(data=pooling_stats,
           aes(x=(scales_collected),y=pct_prop_coverage*100,
               color = factor(age),
               linetype= factor(age)))+
  geom_point()+
  geom_line()+
  theme_bw()+
  scale_x_continuous(limits=c(500,2000),
                     breaks=seq(from=500,to=2000,by=200))+
  scale_y_continuous(limits=c(0,100),
                     breaks=seq(from=0,to=100,by=20))+
  labs(x = "basin scale samples", y = "proportional coverage",color="Age",linetype="Age")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank())+
  facet_grid(factor(age)~.)
p3
p<-grid.arrange(p1,p2,p3, nrow = 1)
ggsave(p,file="outputs/approach2_wpooling_lowsamples_metrics.png",width=10,height=6)


##################
#tracking bias?
#####################
#Load functions
sapply(list.files("scripts/functions", pattern = "\\.R$", full.names = TRUE), source)

trib_data<-readRDS("data/all_returns_clean.Rds")
theta_data<-readRDS("data/location_year_theta.Rds")
year_list<-2010:2020
for(o in 1:10){
  
  target_iter=o
  target_year=year_list[o]
  
  sim_pop<-sim_location_pop(trib_data,
                         iter=target_iter,
                         year=target_year)
  
  sim2<-sim_pop$iter_pop
  sim_tag_rates<-sim_pop$tag_rates
  iterations=1000
  
  sim2_true_stats<-sim2%>%
    group_by(age,location,origin)%>%
    summarize(true_abundance=n())
  sim2_true_hatchery<-sim2_true_stats%>%filter(origin=="hatchery")
  sim2_true_natural<-sim2_true_stats%>%filter(origin=="natural")
  
  loc_list<-unique(sim2$location)#pull list of sampled locations
  
  sim2_true_hatchery$location <- factor(sim2_true_hatchery$location, levels = loc_list)
  sim2_true_natural$location <- factor(sim2_true_natural$location, levels = loc_list)
  
  #p1<-ggplot(sim2_true_hatchery,aes(x=factor(age),y=true_abundance))+
  #  geom_bar(stat="identity",aes(fill=factor(age)))+
  #  theme_bw()+
  #  labs(x = "Age", y = "hatchery abundance",fill="Age")+
  #  facet_grid(.~location,drop=F)
  
  #p2<-ggplot(sim2_true_natural,aes(x=factor(age),y=true_abundance))+
  #  geom_bar(stat="identity",aes(fill=factor(age)))+
  #  theme_bw()+
  #  labs(x = "Age", y = "natural abundance",fill="Age")+
  #  facet_grid(.~location,drop=F)
  #grid.arrange(p1,p2, ncol = 1)
  
  #spawning recoveries for each location
  recoveries_list<-list()
  for(l in 1:length(loc_list)){ #loop recoveries sim for each location
    iter_loc<-loc_list[l]
    loc_pop<-sim2%>%filter(location==iter_loc) #get pop data for locaiton
    target_theta<-theta_data%>% #get theta for location x year
      filter(location==iter_loc,
             return_year==target_year)
    if(nrow(target_theta)==0){target_theta=0.2} else{ #if no theta, set default
      target_theta<-round(target_theta$theta,2)
    }
    target_tag_rate<-sim_tag_rates%>%
      filter(location==iter_loc,
             year==target_year,
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
  
  basin_spawning_results<-list("recovered_fish"=all_recovered_fish,
                               "hatchery_ages"=basin_hatchery_estimates)
  
  ########################
  #look at total_hatchery estimates plots
  hatchery_estimates<-basin_spawning_results$hatchery_ages
  hatchery_totals<-sim2%>%
    filter(origin=="hatchery")%>%
    group_by(age)%>%
    summarise(true_count=n())
  p1<-ggplot(hatchery_estimates%>%filter(age==2),aes(x=total_hatchery))+
    geom_histogram()+
    facet_grid(.~age)+
    geom_vline(xintercept=hatchery_totals$true_count[1],color="red")+
    labs(title="distribution of hatchery estimates",
         x = "Estimated total hatchery", y = "Frequency")+
    theme_bw()
  p2<-ggplot(hatchery_estimates%>%filter(age==3),aes(x=total_hatchery))+
    geom_histogram()+
    facet_grid(.~age)+
    geom_vline(xintercept=hatchery_totals$true_count[2],color="red")+
    labs(title="distribution of hatchery estimates",
         x = "Estimated total hatchery", y = "Frequency")+
    theme_bw()
  p3<-ggplot(hatchery_estimates%>%filter(age==4),aes(x=total_hatchery))+
    geom_histogram()+
    facet_grid(.~age)+
    geom_vline(xintercept=hatchery_totals$true_count[3],color="red")+
    labs(title="distribution of hatchery estimates",
         x = "Estimated total hatchery", y = "Frequency")+
    theme_bw()
  
  p<-grid.arrange(p1,p2,p3, ncol = 1)
  ggsave(p,file=paste("outputs/test_figures/",target_year,target_iter,"total_hatchery.png",sep="_"))
  
  
  recovered_stats<-all_recovered_fish%>%
    group_by(location,age,origin)%>%
    summarize(recovered_count=n())
  
  recovered_stats$location<-factor(recovered_stats$location, levels = loc_list)
  
  recovered_hatchery<-recovered_stats%>%filter(origin=="hatchery")
  recovered_natural<-recovered_stats%>%filter(origin=="natural")
  
 # p1<-ggplot(recovered_hatchery,aes(x=factor(age),y=recovered_count))+
#    geom_bar(stat="identity",aes(fill=factor(age)))+
#   theme_bw()+
#    labs(x = "Age", y = "recovered hatchery",fill="Age")+
#    facet_grid(.~location,drop=F)
  
 # p2<-ggplot(recovered_natural,aes(x=factor(age),y=recovered_count))+
  #  geom_bar(stat="identity",aes(fill=factor(age)))+
   # theme_bw()+
  #  labs(x = "Age", y = "recovered natural",fill="Age")+
  #  facet_grid(.~location,drop=F)
  #grid.arrange(p1,p2, ncol = 1)
  
  #estimate proportion of untagged population from each location
  sim2_untagged<-sim2%>%
    filter(has_cwt==0)
  sim2_untagged<-sim2_untagged%>%
    group_by(location)%>%
    summarize(count=n())
  
  pop_untagged<-nrow(sim2%>%
                       filter(has_cwt==0))
  
  sim2_untagged$true_proportion=sim2_untagged$count/pop_untagged
  
  total_recovered_count<-nrow(all_recovered_fish)
  recovered_untagged<-all_recovered_fish%>%
    group_by(location)%>%
    summarise(recovered_count=n())
  recovered_untagged$recovery_proportion=recovered_untagged$recovered_count/total_recovered_count
  
  recovered_untagged<-recovered_untagged%>%
    left_join(sim2_untagged)
  
  p1<-ggplot()+
    geom_point(data=recovered_untagged,aes(x=location,y=true_proportion),color="blue")+
    geom_point(data=recovered_untagged,aes(x=location,y=recovery_proportion),color="red")+
    theme_bw()+
    labs(x = "Age", y = "proportion of total untagged population")
  
  basin_thetas<-unique(select(all_recovered_fish,location,theta))
  
  recovered_untagged<-recovered_untagged%>%
    left_join(basin_thetas)
  
  recovered_untagged<-recovered_untagged%>%
    mutate(weight=1/theta)
  
  #sample with weighting
  fish_level <- all_recovered_fish%>%
    filter(has_cwt==0)
  fish_level$weight <- 1 / fish_level$theta
  
  fish_level$prob <- fish_level$weight / sum(fish_level$weight)
  
  set.seed(123)  #for reproducibility
  scales_n=1000
  sample_indices <- sample(1:nrow(fish_level), 
                           size = scales_n, 
                           replace = FALSE, 
                           prob = fish_level$prob)
  scale_sample <- fish_level[sample_indices, ]
  
  #check to see if proportions mimic true proportions of untagged
  location_scales<-scale_sample%>%
    group_by(location)%>%
    summarize(proportion=n()/scales_n)
  
  p2<-ggplot()+
    geom_point(data=recovered_untagged,aes(x=location,y=true_proportion),color="blue")+
    geom_point(data=location_scales,aes(x=location,y=proportion),color="red")+
    theme_bw()+
    labs(x = "Age", y = "proportion of total untagged population")
  p<-grid.arrange(p1,p2, ncol = 1)
  ggsave(p,file=paste("outputs/test_figures/",target_year,target_iter,"untagged_proportions.png",sep="_"))
  
  
  ########################
  #weighted scale sampling
  
  scales_n_seq =seq(1000,5000,1000)
  scale_iters <- vector("list", length(scales_n_seq))
  names(scale_iters) <- scales_n_seq
  
  for(scales_idx in 1:length(scales_n_seq)) {
    current_scales_n <- scales_n_seq[scales_idx]
    basin_sim_result<-weighted_scales_sampling(
      pop=sim2,
      spawning_results = basin_spawning_results,
      iterations=iterations,
      scales_n=current_scales_n)
    scale_iters[[scales_idx]]<-basin_sim_result
  }
  names(scale_iters)<-scales_n_seq
  
  weighted_sampling<-data.frame()
  for(i in 1:length(scale_iters)){
    d<-scale_iters[[i]]$est_summary_stats
    d<-d%>%left_join(scale_iters[[i]]$true_summary_stats)
    d$scales_collected<-scale_iters[[i]]$scales_collected
    weighted_sampling<-weighted_sampling%>%rbind(d)
  }
  
  ggplot(weighted_sampling,aes(x=factor(age),y=mean_count_natural))+
    geom_point(aes(color=factor(age)))+
    geom_errorbar(data=weighted_sampling,
                  aes(ymin = lower_CI_count, 
                      ymax = upper_CI_count,
                      color=factor(age),
                      linetype= factor(age)), 
                  width = 0.2) +
    geom_point(aes(y=true_count_natural),shape=2)+
    facet_grid(.~scales_collected)+
    labs(title="weighted scale sampling",
         x = "Age", y = "Estimated count natural",color="Age",linetype="Age")+
    theme_bw()
  ggsave(paste("outputs/test_figures/",target_year,target_iter,"age_estimates.png",sep="_"))
  
}

########################
#unweighted scale sampling
########################
scales_n_seq =seq(1000,5000,1000)
scale_iters <- vector("list", length(scales_n_seq))
names(scale_iters) <- scales_n_seq

for(scales_idx in 1:length(scales_n_seq)) {
  current_scales_n <- scales_n_seq[scales_idx]
  basin_sim_result<-sim_scales_sampling(
    pop=sim2,
    spawning_results = basin_spawning_results,
    iterations=iterations,
    scales_n=current_scales_n)
  scale_iters[[scales_idx]]<-basin_sim_result
}
names(scale_iters)<-scales_n_seq

unweighted_sampling<-data.frame()
for(i in 1:length(scale_iters)){
  d<-scale_iters[[i]]$est_summary_stats
  d<-d%>%left_join(scale_iters[[i]]$true_summary_stats)
  d$scales_collected<-scale_iters[[i]]$scales_collected
  unweighted_sampling<-unweighted_sampling%>%rbind(d)
}

ggplot(unweighted_sampling,aes(x=factor(age),y=mean_count_natural))+
  geom_point(aes(color=factor(age)))+
  geom_errorbar(data=unweighted_sampling,
                aes(ymin = lower_CI_count, 
                    ymax = upper_CI_count,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2) +
  geom_point(aes(y=true_count_natural),shape=2)+
  facet_grid(.~scales_collected)+
  labs(title="unweighted scale sampling",
       subtitle = paste("year:", target_year),
       x = "Age", y = "Estimated count natural",color="Age",linetype="Age")+
  theme_bw()





#################################
#with pooling tag_rate 100
#####################################
pooling_iterations<-readRDS("outputs/basin_results_wpooling_tagrate100.Rds")
target_moe=0.20
CI=0.95
pooling_estimates<-data.frame()
for(i in 1:length(pooling_iterations)){
  d<-pooling_iterations[[i]]
  for(s in 1:length(d$scale_iters)){
    ds<-d$scale_iters[[s]]
    ds$est_summary_stats$scales_collected=ds$target_scales_n
    ds$est_summary_stats$true_count_natural=ds$true_summary_stats$true_count_natural
    ds$est_summary_stats$true_count_hatchery=ds$true_summary_stats$true_count_hatchery
    ds$est_summary_stats$true_proportion_natural=ds$true_summary_stats$true_proportion_natural
    pooling_estimates<-pooling_estimates%>%rbind(ds$est_summary_stats)
  }
}

pooling_estimates<-pooling_estimates%>%
  mutate(moe_count=(upper_CI_count-lower_CI_count)/2,
         moe_count_relative=ifelse(true_count_natural == 0, NA, 
                                   moe_count / true_count_natural),
         moe_prop=(upper_CI_prop-lower_CI_prop)/2,
         meets_target_moe_prop=ifelse(moe_prop<=target_moe,T,F),
         meets_target_moe_count=ifelse(moe_count_relative<=target_moe,T,F),
         coverage_count=ifelse(true_count_natural>=lower_CI_count &
                                 true_count_natural<=upper_CI_count,1,0),
         coverage_proportion=ifelse(true_proportion_natural>=lower_CI_prop &
                                      true_proportion_natural<=upper_CI_prop,1,0),
         count_bias=mean_count_natural-true_count_natural,
         hatchery_count_bias=mean_count_hatchery-true_count_hatchery,
         rel_count_bias=count_bias/true_count_natural,
         rel_hatchery_bias=hatchery_count_bias/true_count_hatchery,
         prop_bias=mean_proportion_natural-true_proportion_natural
  )


pooling_stats<-pooling_estimates%>%
  group_by(age,scales_collected)%>%
  summarize(
    mean_scales_collected=mean(scales_collected),
    mean_prop_moe=mean(moe_prop, na.rm = TRUE),
    lower_CI_prop_moe=quantile(moe_prop, probs = (1-CI)/2, 
                               na.rm = TRUE),
    upper_CI_prop_moe=quantile(moe_prop, probs = 1-(1-CI)/2, 
                               na.rm = TRUE),
    lower_CI_count_moe=quantile(moe_count_relative, probs = (1-CI)/2, 
                                na.rm = TRUE),
    upper_CI_count_moe=quantile(moe_count_relative, probs = 1-(1-CI)/2, 
                                na.rm = TRUE),
    mean_count_moe=mean(moe_count_relative,na.rm=TRUE),
    pct_meeting_target_prop_moe = mean(meets_target_moe_prop) * 100,
    pct_meeting_target_count_moe = mean(meets_target_moe_count) * 100,
    pct_count_coverage=mean(coverage_count),
    pct_prop_coverage=mean(coverage_proportion),
    mean_count_bias=mean(rel_count_bias,na.rm=TRUE),
    mean_hatchery_bias=mean(rel_hatchery_bias,na.rm=TRUE),
    mean_prop_bias=mean(prop_bias,na.rm=T),
    lower_CI_prop_bias=quantile(prop_bias, probs = (1-CI)/2, 
                                na.rm = TRUE),
    upper_CI_prop_bias=quantile(prop_bias, probs = 1-(1-CI)/2, 
                                na.rm = TRUE),
    lower_CI_count_bias=quantile(rel_count_bias, probs = (1-CI)/2, 
                                 na.rm = TRUE),
    upper_CI_count_bias=quantile(rel_count_bias, probs = 1-(1-CI)/2, 
                                 na.rm = TRUE),
    .groups = "drop"
  )

p1<-ggplot(data=pooling_stats,
           aes(x=(scales_collected),y=mean_prop_moe))+
  geom_point( size=3, aes(color=factor(age)))+
  geom_errorbar(data=pooling_stats,
                aes(ymin = lower_CI_prop_moe, 
                    ymax = upper_CI_prop_moe,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2) +
  scale_x_continuous(limits=c(100,1000),
                     breaks=seq(from=100,to=1000,by=200))+
  theme_bw()+
  labs(x = "basin scale samples", y = "proportion MOE",color="Age",linetype="Age")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none")+
  facet_grid(factor(age)~.)
p1
p2<-ggplot(data=pooling_stats,
           aes(x=(scales_collected),y=mean_prop_bias))+
  geom_point( size=3, aes(color=factor(age)))+
  geom_errorbar(data=pooling_stats,
                aes(ymin = lower_CI_prop_bias, 
                    ymax = upper_CI_prop_bias,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2) +
  scale_x_continuous(limits=c(100,1000),
                     breaks=seq(from=100,to=1000,by=200))+
  theme_bw()+
  labs(x = "basin scale samples", y = "proportional bias",color="Age",linetype="Age")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none")+
  facet_grid(factor(age)~.)
p2
p3<-ggplot(data=pooling_stats,
           aes(x=(scales_collected),y=pct_prop_coverage*100,
               color = factor(age),
               linetype= factor(age)))+
  geom_point()+
  geom_line()+
  theme_bw()+
  scale_x_continuous(limits=c(100,1000),
                     breaks=seq(from=100,to=1000,by=200))+
  scale_y_continuous(limits=c(0,100),
                     breaks=seq(from=0,to=100,by=20))+
  labs(x = "basin scale samples", y = "proportional coverage",color="Age",linetype="Age")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank())+
  facet_grid(factor(age)~.)
p3
p<-grid.arrange(p1,p2,p3, nrow = 1)
ggsave(p,file="outputs/approach2_wpooling_tagrate100.png",width=10,height=6)

p4<-ggplot(data=pooling_stats,
           aes(x=(scales_collected),y=mean_hatchery_bias))+
  geom_point( size=3, aes(color=factor(age)))+
  
  scale_x_continuous(limits=c(100,1000),
                     breaks=seq(from=100,to=1000,by=200))+
  theme_bw()+
  labs(x = "basin scale samples", y = "hatchery bias",color="Age",linetype="Age")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none")+
  facet_grid(factor(age)~.)
p4
