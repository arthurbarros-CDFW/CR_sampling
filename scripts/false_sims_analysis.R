rm( list = ls()) #clear env
#looking at false positives and negatives
library(tidyverse)
library(data.table)
library(ggplot2)
library(gtools)

#Load functions
sapply(list.files("scripts/functions", pattern = "\\.R$", full.names = TRUE), source)

sim_h100<-readRDS("outputs/false_positives_test.Rds")
sim_h95<-readRDS("outputs/false_negatives_test.Rds")

#for h100 how many false positives?
false_positives<-sim_h100[[1]]$age_results
false_positives<-false_positives%>%
  group_by(scales_sampled,age)%>%
  summarise(false_positives=mean(fake_natural_estimated)*100)

ggplot(false_positives, aes(x = scales_sampled, y = false_positives, color = factor(age))) +
  geom_line() +
  geom_point() +
  scale_y_continuous(limits=c(90,100),
                     breaks=seq(from=90,to=100,by=2))+
  scale_x_continuous(limits=c(min(false_positives$scales_sampled),max(false_positives$scales_sampled)),
                     breaks=seq(from=min(false_positives$scales_sampled),to=max(false_positives$scales_sampled),by=100))+
  labs(title = "percent iterations with false positives",
       x = "Number of scale samples", 
       y = "percent false positives",
       color = "Age") +
  theme_minimal()

ggsave(file="outputs/false_positives.png",,width=10,height=6)

#false negatives
false_negatives<-sim_h95[[1]]$age_results
false_negatives<-false_negatives%>%
  group_by(scales_sampled,age)%>%
  summarise(false_negatives=mean(true_natural_missed)*100)

ggplot(false_negatives, aes(x = scales_sampled, y = false_negatives, color = factor(age))) +
  geom_line() +
  geom_point() +
  scale_y_continuous(limits=c(0,10),
                     breaks=seq(from=0,to=10,by=2))+
  scale_x_continuous(limits=c(min(false_negatives$scales_sampled),max(false_negatives$scales_sampled)),
                     breaks=seq(from=min(false_negatives$scales_sampled),to=max(false_negatives$scales_sampled),by=100))+
  labs(title = "percent iterations with false negatives",
       x = "Number of scale samples", 
       y = "percent false negatives",
       color = "Age") +
  theme_minimal()

ggsave(file="outputs/false_negatives.png",width=10,height=6)
