#follow up questions
library(tidyverse)
#1) based on historical escapement, how many scales should a given location
#collect using the escapement contribution proportion
trib_data<-readRDS("data/all_returns_clean.Rds")
year_list<-2010:2020
escapement_est<-trib_data%>%
  filter(return_year%in%year_list)%>%
  group_by(iteration,location,return_year)%>%
  summarise(escapement=sum(count))

year_esc<-trib_data%>%
  filter(return_year%in%year_list)%>%
  group_by(iteration,return_year)%>%
  summarise(year_escapement=sum(count))

escapement_est<-escapement_est%>%
  left_join(year_esc)%>%
  mutate(prop_esc=escapement/year_escapement)

ggplot(data=escapement_est)+
  geom_histogram(aes(x=prop_esc))+
  facet_wrap(.~location)

max_proportions<-escapement_est%>%
  group_by(location)%>%
  summarize(max_prop=max(prop_esc))
#instead of max do top 95%
top95_proportions<-escapement_est%>%
  group_by(location)%>%
  summarise(top95_prop=round(quantile(prop_esc,probs=.95),3))

#2)When we set tag_rate=1, I learned that weighting by contribution to total 
# escapement was much less accurate, weighting by contribution to untagged 
# escapement is best.
year_list<-2010:2020
untagged_est<-trib_data%>%
  filter(return_year%in%year_list)%>%
  group_by(iteration,location,return_year)%>%
  summarise(total_esc=sum(count),mean_tag_rate=mean(tag_rate))

untagged_est<-untagged_est%>%
  mutate(untagged_esc=total_esc*(1-mean_tag_rate))

year_untagged<-untagged_est%>%
  group_by(iteration,return_year)%>%
  summarise(year_untagged=sum(untagged_esc))

untagged_est<-untagged_est%>%
  left_join(year_untagged)%>%
  mutate(prop_untagged=untagged_esc/year_untagged)

ggplot(data=untagged_est)+
  geom_histogram(aes(x=prop_untagged))+
  facet_wrap(.~location)+
  theme_bw()

top95_proportions<-untagged_est%>%
  group_by(location)%>%
  summarise(top95_prop=round(quantile(prop_untagged,probs=.95),3))

top95_proportions<-top95_proportions%>%
  mutate(scale_count_500=ceiling(top95_prop*500))
ggsave(file="outputs/untagged_proportions.png",width=6,height=6)
