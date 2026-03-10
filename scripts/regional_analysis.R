rm( list = ls()) #clear env
#reading and prepping and analyzing natural origin age data
library(tidyverse)
library(data.table)
library(ggplot2)
library(gtools)
library(vegan)
###################
#Data Prep

#read files
rds_files <- list.files(path = "C:/Users/ABarros/OneDrive - California Department of Fish and Wildlife/Analysis/Cohort Reconstruction/R work/data/echen_outputs", 
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

#calculate means
mean_ages_dt <- all_dt[, lapply(.SD, mean, na.rm = TRUE), 
                    by = .(location, brood_year), 
                    .SDcols = patterns("^Age")]

#all_regions
all_regions<-mean_ages_dt%>%filter(location=="All")
regions<-mean_ages_dt%>%filter(!location%in%c("All","Other"))

###################
#Bray-Curtis
regions_df <- as.data.frame(regions)

community_matrix <- regions_df %>%
  select(location, brood_year, Age2Sp, Age3Sp, Age4Sp) %>%
  # Create a unique ID for each sample (location_broodyear)
  mutate(sample_id = paste(location, brood_year, sep = "_")) %>%
  column_to_rownames("sample_id") %>%
  select(-location, -brood_year)

#create dissimilarity matrix using vegan
bray_dist <- vegdist(community_matrix, method = "bray")

#get metadata for grouping
metadata <- regions_df %>%
  mutate(sample_id = paste(location, brood_year, sep = "_")) %>%
  select(sample_id, location, brood_year)

#set row names to match the dissimilarity matrix
rownames(metadata) <- metadata$sample_id

#check that orders match
all(rownames(metadata) == rownames(community_matrix))  # Should be TRUE

#PERMANOVA location
permanova_location <- adonis2(bray_dist ~ location, 
                              data = metadata, 
                              permutations = 999, 
                              method = "bray")

#PERMANOVA location + brood year
permanova_both <- adonis2(bray_dist ~ location * brood_year, 
                          data = metadata, 
                          permutations = 999, 
                          method = "bray")

#look at dispersion
dispersion <- betadisper(bray_dist, group = metadata$location)
permutest(dispersion) 
plot(dispersion, main = "Multivariate Dispersion by Location")

#vizualize with nMDS (non-Metric Multi Dimensional Scaling)
nmds <- metaMDS(community_matrix, distance = "bray", k = 2, trymax = 100)

plot(nmds, type = "n")
points(nmds, display = "sites", pch = 19, 
       col = as.numeric(as.factor(metadata$location)))
legend("topright", legend = levels(as.factor(metadata$location)), 
       col = 1:length(unique(metadata$location)), pch = 19)

nmds_scores <-as.data.frame(vegan::scores(nmds)$sites)
nmds_scores$location <- metadata$location
nmds_scores$brood_year <- metadata$brood_year

ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = location)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95) +
  theme_minimal() +
  labs(title = "NMDS Ordination of Age Class Composition",
       color = "Location")

#additional analysis for groupings
location_avg <- regions_df %>%
  group_by(location) %>%
  summarise(across(starts_with("Age"), mean, na.rm = TRUE))

dist_mat <- dist(location_avg[,2:4])
hc <- hclust(dist_mat, method = "ward.D2")
plot(hc, labels = location_avg$location, main = "Location Clustering Based on Age Composition")

########################
#proportional analysis
library(compositions)
community_matrix_prop <- community_matrix %>%
  rowwise() %>%
  mutate(total = sum(c(Age2Sp, Age3Sp, Age4Sp))) %>%
  mutate(across(starts_with("Age"), ~ ./total)) %>%
  select(-total) %>%
  ungroup()

community_matrix_pseudo <- community_matrix_prop + 0.001

comp_data <- acomp(community_matrix_pseudo)
aitchison_dist <- dist(comp_data)

permanova_aitchison <- adonis2(aitchison_dist ~ location, 
                               data = metadata, 
                               permutations = 999)

nmds_prop <- metaMDS(community_matrix_prop, distance = "bray", k = 2, trymax = 100)

nmds_scores_prop <- as.data.frame(vegan::scores(nmds_prop)$sites)
nmds_scores_prop$location <- metadata$location
nmds_scores_prop$brood_year <- metadata$brood_year

ggplot(nmds_scores_prop, aes(x = NMDS1, y = NMDS2, color = location)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95) +
  theme_minimal() +
  labs(title = "NMDS of Age Class PROPORTIONS by Location",
       subtitle = "Looking at relative composition, not absolute abundance",
       color = "Location")

permanova_prop <- adonis2(community_matrix_prop ~ location, 
                          data = metadata, 
                          permutations = 999, 
                          method = "bray")

dispersion_prop <- betadisper(vegdist(community_matrix_prop, "bray"), 
                              group = metadata$location)
permutest(dispersion_prop)

regions_prop <- regions_df %>%
  rowwise() %>%
  mutate(total = sum(c(Age2Sp, Age3Sp, Age4Sp), na.rm = TRUE),
         prop_Age2 = Age2Sp/total,
         prop_Age3 = Age3Sp/total,
         prop_Age4 = Age4Sp/total) %>%
  select(-total, -Age2Sp, -Age3Sp, -Age4Sp) %>%
  pivot_longer(cols = starts_with("prop"), 
               names_to = "age_class", 
               values_to = "proportion")

ggplot(regions_prop, aes(x = brood_year, y = proportion, fill = age_class)) +
  geom_area(position = "fill") +
  facet_wrap(~location, scales = "free_y") +
  theme_minimal() +
  labs(title = "Age Class Proportions Over Time by Location",
       x = "Brood Year", y = "Proportion", fill = "Age Class")

##########
#grouping locations
location_prop_avg <- regions_df %>%
  rowwise() %>%
  mutate(total = sum(c(Age2Sp, Age3Sp, Age4Sp))) %>%
  mutate(prop_Age2 = Age2Sp/total,
         prop_Age3 = Age3Sp/total,
         prop_Age4 = Age4Sp/total) %>%
  group_by(location) %>%
  summarise(across(starts_with("prop"), mean, na.rm = TRUE)) %>%
  column_to_rownames("location")

comp_avg <- acomp(as.matrix(location_prop_avg + 0.001))
dist_comp <- dist(comp_avg)

hc_prop <- hclust(dist_comp, method = "ward.D2")
plot(hc_prop, main = "Location Clustering Based on Age PROPORTIONS")

groups_prop <- cutree(hc_prop, k = 3)
print(groups_prop)

group_assignments <- data.frame(
  location = names(groups_prop),
  group = paste("Group", groups_prop)
)

regions_prop_grouped <- regions_df %>%
  rowwise() %>%
  mutate(total = sum(c(Age2Sp, Age3Sp, Age4Sp), na.rm = TRUE),
         prop_Age2 = Age2Sp/total,
         prop_Age3 = Age3Sp/total,
         prop_Age4 = Age4Sp/total) %>%
  select(-total, -Age2Sp, -Age3Sp, -Age4Sp) %>%
  pivot_longer(cols = starts_with("prop"), 
               names_to = "age_class", 
               values_to = "proportion") %>%
  left_join(group_assignments, by = "location")

group_plots <- list()

for(i in seq_along(unique(regions_prop_grouped$group))) {
  g <- unique(regions_prop_grouped$group)[i]
  
  group_data <- regions_prop_grouped %>% 
    filter(group == g)
  
  p <- ggplot(group_data, aes(x = brood_year, y = proportion, fill = age_class)) +
    geom_area(position = "fill") +
    facet_wrap(~location, ncol = length(unique(group_data$location))) +
    theme_minimal() +
    labs(title = g,
         x = if(i==1)"Brood Year", y = "Proportion") +
    scale_fill_manual(values = c("prop_Age2" = "#66c2a5", 
                                 "prop_Age3" = "#fc8d62", 
                                 "prop_Age4" = "#8da0cb"),
                      labels = c("Age 2", "Age 3", "Age 4"),
                      guide = if(i == 3) guide_legend(title = "Age Class") else "none") +  #control legend
    theme(axis.text.x = element_text(hjust = 1),
          strip.text = element_text(size = 10, face = "bold"),
          strip.background = element_rect(fill = "lightgray", color = NA))
  
  group_plots[[g]] <- p
}

#individual plots
group_plots[[1]]  # First group
group_plots[[2]]  # Second group
group_plots[[3]]  # Third group

library(gridExtra)

grid.arrange(grobs = group_plots, ncol = 1,nrow=3, 
             top = "Age Class Proportions by Location Group")
