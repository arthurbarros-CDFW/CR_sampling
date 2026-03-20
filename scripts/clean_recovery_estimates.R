#clean and format location recovery data produced by E. Chen
#data source: https://github.com/echenfishbitch/SRFC-cohort-reconstruction-wNF
rm( list = ls()) #clear env
#data simulations
library(tidyverse)
library(data.table)
library(ggplot2)
library(gtools)
##################################
#get natural-origin spawner estimates

#read natural-origin spawner files
rds_files <- list.files(path = "data/echen_outputs/natural_spawners", 
                        pattern = "\\.Rds$", 
                        full.names = TRUE)

#read all files into a list
all_data <- lapply(rds_files, readRDS)

#name the list elements with filenames
names(all_data) <- gsub("\\.rds$", "", basename(rds_files))

new_names <- names(all_data) %>%
  str_remove("\\.Rds$") %>%  #remove .rds
  word(-1)                    #get the last word

names(all_data) <- new_names

dt_list <- list()

for(location in names(all_data)) {
  #for each location, combine all iterations
  location_dt <- rbindlist(lapply(all_data[[location]], as.data.table), 
                           idcol = "iteration")
  location_dt[, location := location]
  dt_list[[location]] <- location_dt
}

#combine all locations
all_dt <- rbindlist(dt_list)

#reshape from wide to long format
all_dt_long <- melt(all_dt, 
                    id.vars = c("iteration", "brood_year", "location"),
                    measure.vars = patterns("Age[0-9]+Sp"),
                    variable.name = "age_class",
                    value.name = "count")

#get age from the age_class column
all_dt_long[, age := as.numeric(gsub("Age|Sp", "", age_class))]

#calculate return_year
all_dt_long[, return_year := brood_year + age]

nattie_trib_returns<- all_dt_long[, .(count = sum(count)), 
                                  by = .(iteration, location, return_year,age)]
nattie_trib_returns$origin<-"natural"

nattie_trib_returns<-nattie_trib_returns%>%
  filter(location!="All")
nattie_trib_returns$tag_rate=0

##################################
#get hatchery-origin trib spawner estimates

hatch_returns<-readRDS("data/echen_outputs/Escape to Spawning Grounds by River.Rds")

hatch_returns<-rbindlist(hatch_returns,idcol = "iteration")

hatch_trib_returns <- melt(hatch_returns, 
                           id.vars = c("iteration", "run_year", "recovery_location_name"),
                           measure.vars = patterns("Age[0-9]+Sp_Hatchery"),
                           variable.name = "age_class",
                           value.name = "hatch_count")

hatch_cwt_returns <- melt(hatch_returns, 
                          id.vars = c("iteration", "run_year", "recovery_location_name"),
                          measure.vars = patterns("Age[0-9]+Sp_CWT"),
                          variable.name = "age_class",
                          value.name = "cwt_count")



#get age from the age_class column (remove "Age" and "Sp")
hatch_trib_returns[, age := as.numeric(gsub("Age|Sp_Hatchery", "", age_class))]
hatch_cwt_returns[, age := as.numeric(gsub("Age|Sp_CWT", "", age_class))]

hatch_trib_returns<-hatch_trib_returns%>%
  rename(location=recovery_location_name,
         return_year=run_year)%>%
  select(-age_class)

hatch_cwt_returns<-hatch_cwt_returns%>%
  rename(location=recovery_location_name,
         return_year=run_year)%>%
  select(-age_class)

hatch_trib_returns<-hatch_trib_returns%>%
  left_join(hatch_cwt_returns)%>%
  mutate(tag_rate=cwt_count/hatch_count)

hatch_trib_returns<-hatch_trib_returns%>%
  mutate(tag_rate=ifelse(cwt_count==0,.25,tag_rate))

hatch_trib_returns<-hatch_trib_returns%>%
  rename(count=hatch_count)%>%
  select(-cwt_count)

hatch_trib_returns$origin<-"hatchery"

#other tribs
other_tribs<-c("COTTONWOOD CREEK","COW CREEK","MILL CREEK",
               "DEER CREEK","BUTTE CREEK")
battle_trib<-c("BATTLE CREEK BELOW CNFH","COLEMAN NFH")
  
#fix hatch locations
hatch_trib_returns<-hatch_trib_returns%>%
  mutate(location=ifelse(location%in%other_tribs,"Other",location))%>%
  mutate(location=ifelse(location%in%battle_trib,"Battle",location))%>%
  mutate(location= recode(location,
                          "AMERICAN RIVER"="American",
                          "CLEAR CREEK"="Clear",
                          "FEATHER RIVER"="Feather",
                          "SAC R AB RBDD"="Sac",
                          "YUBA RIVER"="Yuba"))
hatch_trib_returns<-hatch_trib_returns%>%
  filter(location%in%unique(nattie_trib_returns$location))

#add FRH and NFH
FRH_hatch_returns<-read.csv("data/echen_outputs/CWT Hatchery FRH.csv")
FRH_hatch_returns$location="FRH"
NFH_hatch_returns<-read.csv("data/echen_outputs/CWT Hatchery NFH.csv")
NFH_hatch_returns$location="NFH"

hatchery_returns<-FRH_hatch_returns%>%rbind(NFH_hatch_returns)

hatchery_returns <- melt(data.table(hatchery_returns), 
                    id.vars = c("brood_year", "location"),
                    measure.vars = patterns("Age[0-9]+Hat"),
                    variable.name = "age_class",
                    value.name = "count")
hatchery_returns[, age := as.numeric(gsub("Age|Hat", "", age_class))]

#calculate return_year
hatchery_returns[, return_year := brood_year + age]

hatchery_returns<- hatchery_returns[, .(count = sum(count)), 
                                  by = .( location, return_year,age)]

hatchery_returns<-hatchery_returns%>%
  select(return_year,age,count,location)

#add in 0s
all_combinations <- hatchery_returns[, CJ(return_year = unique(return_year),
                                          age = unique(age),
                                          location = unique(location),
                                          unique = TRUE)]

hatchery_returns <- merge(all_combinations, 
                                   hatchery_returns, 
                                   by = c("return_year", "age", "location"),
                                   all.x = TRUE)
hatchery_returns[is.na(count), count := 0]
hatchery_returns$origin="hatchery"
hatchery_returns$tag_rate=.25
hatchery_returns<-hatchery_returns%>%
  filter(age>1,
         return_year%in%unique(nattie_trib_returns$return_year))%>%
  mutate(age=ifelse(age==5,4,age))

#add iterations with noise
n_iterations <- 1000
hatchery_returns_iter <- rbindlist(
  lapply(1:n_iterations, function(i) {
    #create a copy of the original data for this iteration
    dt_iter <- copy(hatchery_returns)
    
    #add iteration
    dt_iter[, iteration := i]
    
    #add random noise to count
    #use log-normal noise to keep counts positive
    dt_iter[, count := ifelse(count == 0,
                              #for zeros: generate small positive numbers
                              runif(.N, min = 0, max = 5),
                              #for non-zeros: multiplicative noise
                              count * exp(rnorm(.N, mean = 0, sd = 0.1)))]

    return(dt_iter)
  })
)

hatch_trib_returns<-hatch_trib_returns%>%
  rbind(hatchery_returns_iter)

hatch_trib_returns<-hatch_trib_returns%>%
  mutate(age=ifelse(age==5,4,age))

#all iterative return estimates
all_returns<-hatch_trib_returns%>%
  rbind(nattie_trib_returns)

all_returns<-all_returns%>%
  mutate(count=round(count))

#now get and clean theta estimates based on sampling effort (RMIS estimated_number)
FRH_recoveries<-read.csv("data/echen_outputs/CWTRecoveries FRH.csv")
FRH_recoveries<-FRH_recoveries%>%
  filter(fishery%in%c(54,50))
FRH_recoveries<-unique(select(FRH_recoveries,
                               run_year,
                               recovery_location_name,
                               estimated_number))
FRH_recoveries$location<-FRH_recoveries$recovery_location_name
FRH_recoveries$return_year<-FRH_recoveries$run_year

FRH_recoveries<-FRH_recoveries%>%
  mutate(location=ifelse(location%in%other_tribs,"Other",location))%>%
  mutate(location=ifelse(location%in%battle_trib,"Battle",location))%>%
  mutate(location= recode(location,
                          "AMERICAN RIVER"="American",
                          "CLEAR CREEK"="Clear",
                          "FEATHER RIVER"="Feather",
                          "SAC R AB RBDD"="Sac",
                          "YUBA RIVER"="Yuba",
                          "NIMBUS FISH HATCHERY"= "NFH",
                          "FEATHER R HATCHERY"="FRH"))

FRH_recoveries<-FRH_recoveries%>%
  filter(location%in%unique(all_returns$location))%>%
  group_by(location,return_year)%>%
  summarise(theta=1/sum(estimated_number))

FRH_recoveries<-FRH_recoveries%>%
  filter(as.numeric(return_year)%in%unique(all_returns$return_year))

average_year_theta<-FRH_recoveries%>%
  group_by(return_year)%>%
  summarise(year_theta=mean(theta,na.rm=T))

FRH_recoveries<-FRH_recoveries%>%left_join(average_year_theta)
FRH_recoveries<-FRH_recoveries%>%
  mutate(theta=ifelse(is.na(theta),year_theta,theta))%>%
  select(-year_theta)

location_year_theta<-FRH_recoveries

saveRDS(all_returns,"data/all_returns_clean.Rds")
saveRDS(location_year_theta,"data/location_year_theta.Rds")
