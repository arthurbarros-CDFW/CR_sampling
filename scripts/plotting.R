rm( list = ls()) #clear env
#basin methods simulation
#this script aims to simulate basin populations and replicate methods
# to produce a basin-wide cohort reconstruction
library(tidyverse,quietly = "true")
library(data.table)
library(ggplot2)
library(gtools)
library(crayon)

CI=.95

################################
#regular pooling
################################
pooling_iteration_results<-readRDS("outputs/methods_comparisons.Rds")
#plotting

pooling_estimates<-data.frame()
for(i in 1:length(pooling_iteration_results)){
  d<-pooling_iteration_results[[i]]
  for(s in 1:length(d)){
    ds<-d[[s]]$comparison
    ds$n_scales<-as.numeric(names(d[s]))
    
    pooling_estimates<-pooling_estimates%>%rbind(ds)
  }
}

pooling_estimates <- pooling_estimates %>%
  pivot_longer(
    cols = c(mean_hatchery, lower_CI_hatchery, upper_CI_hatchery,
             mean_natural, lower_CI_natural, upper_CI_natural),
    names_to = c(".value", "origin"),
    names_pattern = "(mean|lower_CI|upper_CI)_(hatchery|natural)"
  )%>%
  rename(estimate=mean)%>%
  mutate(true_val=ifelse(origin=="hatchery",hatchery,natural))%>%
  select(origin,age,n_iterations,n_scales,estimate,lower_CI,upper_CI,true_val)

pooling_estimates<-pooling_estimates%>%
  mutate(error=estimate-true_val,
         rel_error=((estimate-true_val)/true_val),
         coverage=ifelse(lower_CI<=true_val & true_val<=upper_CI,
                         TRUE,FALSE),
         moe=((upper_CI-lower_CI)/2))


stats<-pooling_estimates%>%
  group_by(origin,age,n_scales)%>%
  summarise(
    mean_rel_error = mean(rel_error, na.rm = TRUE),
    lower_CI_rel_error = quantile(rel_error, probs = (1-CI)/2, 
                                  na.rm = TRUE),
    upper_CI_rel_error = quantile(rel_error, probs = 1-(1-CI)/2, 
                                  na.rm = TRUE),
    mean_moe=mean(moe, na.rm = TRUE),
    lower_CI_moe = quantile(moe, probs = (1-CI)/2, 
                                  na.rm = TRUE),
    upper_CI_moe = quantile(moe, probs = 1-(1-CI)/2, 
                                  na.rm = TRUE),
    pct_coverage = mean(coverage)
  )

p_rel_error<-ggplot(data=stats,
           aes(x=(n_scales),y=mean_rel_error))+
  geom_point( size=3, aes(color=factor(age)),
              position = position_dodge(width = 200))+
  geom_errorbar(aes(ymin = lower_CI_rel_error, 
                    ymax = upper_CI_rel_error,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2,
                position = position_dodge(width = 200))+
  scale_x_continuous(limits=c(800,10200),
                     breaks=seq(from=1000,to=10000,by=1000))+
  scale_y_continuous(labels = scales::percent)+
  theme_bw()+
  labs(x = NULL, y = "relative error",color="Age",
       linetype="Age",
       tag="A")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        legend.position = "top")+
  facet_wrap(.~origin)

p_moe<-ggplot(data=stats,
                    aes(x=(n_scales),y=mean_moe))+
  geom_point( size=3, aes(color=factor(age)),
              position = position_dodge(width = 200))+
  geom_errorbar(aes(ymin = lower_CI_moe, 
                    ymax = upper_CI_moe,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2,
                position = position_dodge(width = 200))+
  scale_x_continuous(limits=c(800,10200),
                     breaks=seq(from=1000,to=10000,by=1000))+
  
  theme_bw()+
  labs(x = NULL, y = "Margin of Error",color="Age",
       linetype="Age",
       tag="B")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.text.x = element_blank(),  
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.text = element_blank())+
  facet_wrap(.~origin)

p_coverage<-ggplot(data=stats,
                  aes(x=(n_scales),y=pct_coverage))+
  geom_point( size=3, aes(color=factor(age)),
              position = position_dodge(width = 200))+
  geom_line(aes(color=factor(age),
                linetype= factor(age)))+
  scale_x_continuous(limits=c(800,10200),
                     breaks=seq(from=1000,to=10000,by=1000))+
  scale_y_continuous(labels = scales::percent)+
  theme_bw()+
  labs(x = "Scale sample size", y = "CI coverage",color="Age",
       linetype="Age",
       tag="C")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none",
        strip.text = element_blank())+
  facet_wrap(.~origin)

p<-grid.arrange(p_rel_error,p_moe,p_coverage,
                ncol = 1,
                heights=c(1.2,.8,1))

ggsave(p,file="outputs/p_largesamples.png",width=7.5,height=6)


################################
#regular pooling small samples
################################
pooling_iteration_results<-readRDS("outputs/methods_comparisons_lowsamples.Rds")
#plotting

pooling_estimates<-data.frame()
for(i in 1:length(pooling_iteration_results)){
  d<-pooling_iteration_results[[i]]
  for(s in 1:length(d)){
    ds<-d[[s]]$comparison
    ds$n_scales<-as.numeric(names(d[s]))
    
    pooling_estimates<-pooling_estimates%>%rbind(ds)
  }
}

pooling_estimates <- pooling_estimates %>%
  pivot_longer(
    cols = c(mean_hatchery, lower_CI_hatchery, upper_CI_hatchery,
             mean_natural, lower_CI_natural, upper_CI_natural),
    names_to = c(".value", "origin"),
    names_pattern = "(mean|lower_CI|upper_CI)_(hatchery|natural)"
  )%>%
  rename(estimate=mean,
         age=tag_age)%>%
  mutate(true_val=ifelse(origin=="hatchery",hatchery,natural))%>%
  select(origin,age,n_iterations,n_scales,estimate,lower_CI,upper_CI,true_val)

pooling_estimates<-pooling_estimates%>%
  mutate(error=estimate-true_val,
         rel_error=((estimate-true_val)/true_val),
         coverage=ifelse(lower_CI<=true_val & true_val<=upper_CI,
                         TRUE,FALSE),
         moe=((upper_CI-lower_CI)/2))


stats<-pooling_estimates%>%
  group_by(origin,age,n_scales)%>%
  summarise(
    mean_rel_error = mean(rel_error, na.rm = TRUE),
    lower_CI_rel_error = quantile(rel_error, probs = (1-CI)/2, 
                                  na.rm = TRUE),
    upper_CI_rel_error = quantile(rel_error, probs = 1-(1-CI)/2, 
                                  na.rm = TRUE),
    mean_moe=mean(moe, na.rm = TRUE),
    lower_CI_moe = quantile(moe, probs = (1-CI)/2, 
                                na.rm = TRUE),
    upper_CI_moe = quantile(moe, probs = 1-(1-CI)/2, 
                                na.rm = TRUE),
    pct_coverage = mean(coverage)
  )

p_rel_error<-ggplot(data=stats,
                    aes(x=(n_scales),y=mean_rel_error))+
  geom_point( size=3, aes(color=factor(age)),
              position = position_dodge(width = 20))+
  geom_errorbar(aes(ymin = lower_CI_rel_error, 
                    ymax = upper_CI_rel_error,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2,
                position = position_dodge(width = 20))+
  scale_x_continuous(limits=c(80,1020),
                     breaks=seq(from=100,to=1000,by=100))+
  scale_y_continuous(labels = scales::percent)+
  theme_bw()+
  labs(x = NULL, y = "relative error",color="Age",
       linetype="Age",
       tag="A")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        legend.position = "top")+
  facet_wrap(.~origin)

p_moe<-ggplot(data=stats,
                  aes(x=(n_scales),y=mean_moe))+
  geom_point( size=3, aes(color=factor(age)),
              position = position_dodge(width = 20))+
  geom_errorbar(aes(ymin = lower_CI_moe, 
                    ymax = upper_CI_moe,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2,
                position = position_dodge(width = 20))+
  scale_x_continuous(limits=c(80,1020),
                     breaks=seq(from=100,to=1000,by=100))+
  theme_bw()+
  labs(x = NULL, y = "Margin of Error",color="Age",
       linetype="Age",
       tag="B")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.text.x = element_blank(),  
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.text = element_blank())+
  facet_wrap(.~origin)
p_moe

p_coverage<-ggplot(data=stats,
                   aes(x=(n_scales),y=pct_coverage))+
  geom_point( size=3, aes(color=factor(age)),
              position = position_dodge(width = 20))+
  geom_line(aes(color=factor(age),
                linetype= factor(age)))+
  scale_x_continuous(limits=c(80,1020),
                     breaks=seq(from=100,to=1000,by=100))+
  scale_y_continuous(labels = scales::percent,
                     limits=c(.5,1))+
  theme_bw()+
  labs(x = "Scale sample size", y = "CI coverage",color="Age",
       linetype="Age",
       tag="C")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none",
        strip.text = element_blank())+
  facet_wrap(.~origin)
p_coverage

p<-grid.arrange(p_rel_error,p_moe,p_coverage,
                ncol = 1,
                heights=c(1.2,.8,1))

ggsave(p,file="outputs/p_smallsamples.png",width=7.5,height=6)



################################
#regular pooling theta impact
################################
pooling_iteration_results<-readRDS("outputs/methods_comparisons_inc_theta.Rds")
#plotting

pooling_estimates<-data.frame()
for(i in 1:length(pooling_iteration_results)){
  d<-pooling_iteration_results[[i]]
  for(s in 1:length(d)){
    ds<-d[[s]]$'500'$comparison
    ds$theta<-as.numeric(names(d[s]))
    ds$iteration=i
    
    pooling_estimates<-pooling_estimates%>%rbind(ds)
  }
}


pooling_estimates<-pooling_estimates%>%
  mutate(error=estimate-true_val,
         rel_error=((estimate-true_val)/true_val),
         coverage=ifelse(lower_CI<=true_val & true_val<=upper_CI,
                         TRUE,FALSE),
         moe=((upper_CI-lower_CI)/2))


stats<-pooling_estimates%>%
  group_by(origin,age,theta)%>%
  summarise(
    mean_rel_error = mean(rel_error, na.rm = TRUE),
    lower_CI_rel_error = quantile(rel_error, probs = (1-CI)/2, 
                                  na.rm = TRUE),
    upper_CI_rel_error = quantile(rel_error, probs = 1-(1-CI)/2, 
                                  na.rm = TRUE),
    mean_moe=mean(moe, na.rm = TRUE),
    lower_CI_moe = quantile(moe, probs = (1-CI)/2, 
                                na.rm = TRUE),
    upper_CI_moe = quantile(moe, probs = 1-(1-CI)/2, 
                                na.rm = TRUE),
    pct_coverage = mean(coverage)
  )

p_rel_error<-ggplot(data=stats,
                    aes(x=(theta),y=mean_rel_error))+
  geom_point( size=3, aes(color=factor(age)),
              position = position_dodge(width = .02))+
  geom_errorbar(aes(ymin = lower_CI_rel_error, 
                    ymax = upper_CI_rel_error,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2,
                position = position_dodge(width = .02))+
  scale_x_continuous(limits=c(-0.1,1.1),
                     breaks=seq(from=-0.1,to=1.1,by=0.1))+
  scale_y_continuous(labels = scales::percent)+
  theme_bw()+
  labs(x = NULL, y = "relative error",color="Age",
       linetype="Age",
       tag="A")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        legend.position = "top")+
  facet_wrap(.~origin)

p_moe<-ggplot(data=stats,
                  aes(x=(theta),y=mean_moe))+
  geom_point( size=3, aes(color=factor(age)),
              position = position_dodge(width = .02))+
  geom_errorbar(aes(ymin = lower_CI_moe, 
                    ymax = upper_CI_moe,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2,
                position = position_dodge(width = .02))+
  scale_x_continuous(limits=c(-0.1,1.1),
                     breaks=seq(from=-0.1,to=1.1,by=0.1))+
  theme_bw()+
  labs(x = "theta", y = "Margin of Error",color="Age",
       linetype="Age",
       tag="B")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),,
        legend.position = "none")+
  facet_wrap(.~origin)

p_coverage<-ggplot(data=stats,
                   aes(x=(theta),y=pct_coverage))+
  geom_point( size=3, aes(color=factor(age)),
              position = position_dodge(width = .02))+
  geom_line(aes(color=factor(age),
                linetype= factor(age)))+
  scale_x_continuous(limits=c(-0.1,1.1),
                     breaks=seq(from=-0.1,to=1.1,by=0.1))+
  scale_y_continuous(labels = scales::percent,
                     limits=c(.50,1))+
  theme_bw()+
  labs(x = "Theta", y = "CI coverage",color="Age",
       linetype="Age",
       tag="C")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none")+
  facet_wrap(.~origin)

p<-grid.arrange(p_rel_error,p_moe,
                ncol = 1,
                heights=c(1,1))

ggsave(p,file="outputs/p_theta.png",width=7.5,height=6)

################################
#regular pooling tag_rate impact
################################
pooling_iteration_results<-readRDS("outputs/methods_comparisons_inc_tagrate.Rds")
#plotting

pooling_estimates<-data.frame()
for(i in 1:length(pooling_iteration_results)){
  d<-pooling_iteration_results[[i]]
  for(s in 1:length(d)){
    ds<-d[[s]]$'500'$comparison
    ds$theta<-as.numeric(names(d[s]))
    
    pooling_estimates<-pooling_estimates%>%rbind(ds)
  }
}

pooling_estimates<-pooling_estimates%>%
  mutate(error=estimate-true_val,
         rel_error=((estimate-true_val)/true_val),
         coverage=ifelse(lower_CI<=true_val & true_val<=upper_CI,
                         TRUE,FALSE),
         moe=((upper_CI-lower_CI)/2))


stats<-pooling_estimates%>%
  group_by(origin,age,tag_rate)%>%
  summarise(
    mean_rel_error = mean(rel_error, na.rm = TRUE),
    lower_CI_rel_error = quantile(rel_error, probs = (1-CI)/2, 
                                  na.rm = TRUE),
    upper_CI_rel_error = quantile(rel_error, probs = 1-(1-CI)/2, 
                                  na.rm = TRUE),
    mean_moe=mean(moe, na.rm = TRUE),
    lower_CI_moe = quantile(moe, probs = (1-CI)/2, 
                                na.rm = TRUE),
    upper_CI_moe = quantile(moe, probs = 1-(1-CI)/2, 
                                na.rm = TRUE),
    pct_coverage = mean(coverage)
  )

p_rel_error<-ggplot(data=stats,
                    aes(x=(tag_rate),y=mean_rel_error))+
  geom_point( size=3, aes(color=factor(age)),
              position = position_dodge(width = .02))+
  geom_errorbar(aes(ymin = lower_CI_rel_error, 
                    ymax = upper_CI_rel_error,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2,
                position = position_dodge(width = .02))+
  scale_x_continuous(limits=c(0,1.1),
                     breaks=seq(from=0,to=1.1,by=0.1))+
  scale_y_continuous(labels = scales::percent)+
  coord_cartesian(ylim = c(-1, 1))+
  theme_bw()+
  labs(x = NULL, y = "relative error",color="Age",
       linetype="Age",
       tag="A")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        legend.position = "top")+
  facet_wrap(.~origin)

p_moe<-ggplot(data=stats,
                  aes(x=(tag_rate),y=mean_moe))+
  geom_point( size=3, aes(color=factor(age)),
              position = position_dodge(width = .02))+
  geom_errorbar(aes(ymin = lower_CI_moe, 
                    ymax = upper_CI_moe,
                    color=factor(age),
                    linetype= factor(age)), 
                width = 0.2,
                position = position_dodge(width = .02))+
  scale_x_continuous(limits=c(0,1.1),
                     breaks=seq(from=0,to=1.1,by=0.1))+
  theme_bw()+
  labs(x = NULL, y = "relative moe",color="Age",
       linetype="Age",
       tag="B")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.text.x = element_blank(),  
        axis.title.x = element_blank(),
        legend.position = "none")+
  facet_wrap(.~origin)

p_coverage<-ggplot(data=stats,
                   aes(x=(tag_rate),y=pct_coverage))+
  geom_point( size=3, aes(color=factor(age)),
              position = position_dodge(width = .02))+
  geom_line(aes(color=factor(age),
                linetype= factor(age)))+
  scale_x_continuous(limits=c(0,1.1),
                     breaks=seq(from=0,to=1.1,by=0.1))+
  scale_y_continuous(labels = scales::percent,
                     limits=c(0.5,1))+
  theme_bw()+
  labs(x = "tag rate", y = "CI coverage",color="Age",
       linetype="Age",
       tag="C")+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.position = "none")+
  facet_wrap(.~origin)

p<-grid.arrange(p_rel_error,p_moe,p_coverage,
                ncol = 1,
                heights=c(1.3,1,1.2))

ggsave(p,file="outputs/p_tagrate.png",width=7.5,height=6)
