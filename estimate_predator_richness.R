###################################
#### clean working environment ####
###################################
rm(list = ls())

################################################
#### load packages and data, format dataset ####
################################################
library(dplyr)
library(ggplot2)
library(grid)

predation_allometry <- read.delim("~/My Drive/predation_cooperation/data/final_datasets/predation_allometry.txt")

distribution_overlap <- read.delim("~/My Drive/predation_cooperation/data/final_datasets/distribution_overlap.txt")


#####################################################
#### create range of prey mass for each predator ####
#####################################################

# these colums are already present in the "cooperation_pred_final.txt" final dataset used for analyses, 
# and this R script is only present to show how these were computed.

# The "predation" dataset only contains bird species that eat birds as part of their diet.
# Species of predators that rarely eat birds were already removed.

# estimate min prey mass depending on predator mass: log-log predator-prey mass allometry
modmin <- lm(predation_allometry$min.prey.log ~ predation_allometry$masslog)
summary(modmin)
fmin <- function(x){
  as.numeric(modmin$coefficients[1]) + as.numeric(modmin$coefficients[2]) * x
}

# estimate max prey mass depending on predator mass: log-log predator-prey mass allometry
modmax <- lm(predation_allometry$max.prey.log ~ predation_allometry$masslog)
summary(modmax)
fmax <- function(x){
  as.numeric(modmax$coefficients[1]) + as.numeric(modmax$coefficients[2]) * x
}

ggplot(data = predation_allometry, aes(x=masslog, y=max.prey.log))+
  geom_point(color="black", size = 2)+
  geom_point(aes(y = min.prey.log), color="red", size = 2)+
  geom_smooth(method = "lm", size=1.2, color="black", fullrange=T, fill="black", alpha=0.08)+
  geom_smooth(method = "lm", size=1.2, aes(y=min.prey.log), color="red", fullrange=T, fill="red", alpha=0.08)+
  xlab("Log predator body mass")+
  ylab("Log prey body mass")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size=22), axis.text = element_text(size=18))

length(predation_allometry$min.prey.log[!is.na(predation_allometry$min.prey.log)])
length(predation_allometry$max.prey.log[!is.na(predation_allometry$max.prey.log)])

# add min and max prey mass estimated to the dataset
predation_allometry$estimated_min_prey <- exp(fmin(predation_allometry$masslog))
predation_allometry$estimated_max_prey <- exp(fmax(predation_allometry$masslog))




#################################################################################################
#### use predation allometry and distribution overlap to generate average predation richness ####
#################################################################################################

distribution_overlap <- distribution_overlap[,c(1:307)] # the five last columns are not needed for this step (latitude, and predation richness metrics that we are estimating in this script)


# this loop check if each species fall in the allometric prey range of each predator. If they don't, replace number of shared cell per 0 (i.e not a predator of the species)

for (i in 1:dim(distribution_overlap)[1]) {          # loop through ech focal species
  for (j in 6:(dim(distribution_overlap)[2]-1)) {    # loop through each predator species
    
    if (distribution_overlap$mass[i] <= predation_allometry$estimated_min_prey[predation_allometry$jetz.name == colnames(distribution_overlap[j])]
        || distribution_overlap$mass[i] >= predation_allometry$estimated_max_prey[predation_allometry$jetz.name == colnames(distribution_overlap[j])]
        || is.na(distribution_overlap[i,j])
    ) {
      distribution_overlap[i,j] <- 0
    }
    
  }
}

# create new variable to store average predation richness
distribution_overlap$average_predation_richness <- NA

# calculate average predation richness
for(i in 1:dim(distribution_overlap)[1]){
  distribution_overlap$average_predation_richness[i] <- sum(distribution_overlap[i, 6:(dim(distribution_overlap)[2]-1)]) /
    distribution_overlap$Total_cells[i]
}

distribution_overlap$average_predation_richness


# create new variable to store total predation richness
distribution_overlap$total_predation_richness <- NA

# calculate total predation richness
for(i in 1:dim(distribution_overlap)[1]){
  distribution_overlap$total_predation_richness[i] <- sum(distribution_overlap[i, 6:(dim(distribution_overlap)[2]-2)] > 0)
}

distribution_overlap$total_predation_richness

