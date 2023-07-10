###################################
#### clean working environment ####
###################################
rm(list = ls())

#############################################
#### load packages and formatted dataset ####
#############################################
library(readxl)
library(dplyr)
library(brms)
library(ape)
library(phangorn)
library(bayesplot)
library(ggplot2)
library(grid)
library(cmdstanr)
library(patchwork)

cooperation_pred <- read.csv("~/My Drive/predation_cooperation/data/final_datasets/cooperation_pred.csv", row.names=1)

### If running the model without migratory species ###
# if  not, skip this line
# cooperation_pred <- subset(cooperation_pred, mov_min=="sed")

# create a social system variable, as an ordered factor
cooperation_pred$coop <- factor(cooperation_pred$fam_sys_known50, ordered = T, levels = c("no_fam", "family", "coop_families"))

##########################
#### Import phylogeny ####
##########################

#MyTree <- read.nexus("~/My Drive/predation_cooperation/data/output_ericson.nex") # import distribution of 100 trees
#MyTree <- read.nexus("~/My Drive/predation_cooperation/data/output_hackett.nex") # import distribution of 100 trees
MyTree <- read.tree("~/My Drive/predation_cooperation/data/Prum_Jetz_Cooney_9993.tree") 

tree_phylo <- maxCladeCred(MyTree) # get maximum clade credibility tree

setdiff(cooperation_pred$tip_label, tree_phylo$tip.label) # check that all species are in the phylogeny

tree_phylo <- keep.tip(tree_phylo, cooperation_pred$tip_label) # drop tip from phylo: remove the species not present in the dataset

phylo_mcct <- vcv.phylo(tree_phylo) # transform into a variance covariance matrix

cooperation_pred$phylo <- cooperation_pred$tip_label

#####################################################################
#### Drop the few species with unreliable data for social system ####
#####################################################################

cooperation_pred <- cooperation_pred[cooperation_pred$tip_label != "Melanodryas_vittata",]
cooperation_pred <- cooperation_pred[cooperation_pred$tip_label != "Tangara_cyanoventris",]
cooperation_pred <- cooperation_pred[cooperation_pred$tip_label != "Fulica_ardesiaca",]
cooperation_pred <- cooperation_pred[cooperation_pred$tip_label != "Tangara_larvata",]

####################################
#### Scale continuous variables ####
####################################

sd_predation <- sd(cooperation_pred$average_predation_richness)
mean_predation <- mean(cooperation_pred$average_predation_richness)
cooperation_pred$average_predation_richness <- scale(cooperation_pred$average_predation_richness)

cooperation_pred$log_mass <- log(cooperation_pred$mass)
sd_logmass <- sd(cooperation_pred$log_mass)
mean_logmass <- mean(cooperation_pred$log_mass)
cooperation_pred$log_mass <- scale(cooperation_pred$log_mass)

sd_habitat <- sd(cooperation_pred$habitat)
mean_habitat <- mean(cooperation_pred$habitat)
cooperation_pred$habitat <- scale(cooperation_pred$habitat)

cooperation_pred$latitude_mean <- abs(cooperation_pred$latitude_mean)
sd_latitude <- sd(cooperation_pred$latitude_mean)
mean_latitude <- mean(cooperation_pred$latitude_mean)
cooperation_pred$latitude_mean <- scale(cooperation_pred$latitude_mean)

####################
#### Full model ####
####################

m2 <- brm(
  coop ~ average_predation_richness * habitat +
    log_mass +
    latitude_mean +
    scale(Prcp.P) + scale(Temp.P) + 
    scale(Prcp.Mean) + scale(Temp.Mean) +
    scale(Prcp.Var) + scale(Temp.Var) +
    (1|gr(phylo, cov = phylo_mcct)), data = cooperation_pred, 
  family = cumulative("logit"), data2 = list(phylo_mcct = phylo_mcct),
  warmup = 1000, iter = 2000, thin = 1, chains = 3, seed = 25, cores = 3,
  prior = c(
    set_prior("normal(0,2)", class = "b"),
    set_prior("normal(0,2)", class = "Intercept"),
    set_prior("normal(0,1.5)", class = "sd", coef = "Intercept", group = "phylo")),
  control = list(adapt_delta = 0.99), save_ranef = F, backend = "cmdstanr")

summary(m2)

# posterior predictive checks
pp_check(m2, ndraws = 50)
pp_check(m2, ndraws = 100, type = "stat_2d")
plot(m2)


################################
#### Plot from model output ####
################################


## Plot model output ##

# extract posterior distributions
posterior <- as.matrix(m2)
dat_plot <- as.data.frame(posterior)

# prediction across range of predation richness 
x2.sim <- seq(min(cooperation_pred$average_predation_richness), max(cooperation_pred$average_predation_richness), by = 0.01)

## from model output, it is possible to get the probability of being in each sociality category across values of predation richness
## Given that the ordinal category are 0=non family living, 1=family living, 2 cooperative breeding
## We plot prediction as 0 * p(non-family living) + 1 * p(non-family living) + 2 * p(cooperative breeding)


#### plot across predation richness for close habitat (i.e habitat set to -1SD)

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- (2*(1 - ((exp(mean(dat_plot$`b_Intercept[2]`) - dat_plot$b_habitat * -1 - dat_plot$`b_average_predation_richness:habitat` * -1 * (x2.sim[i]) - dat_plot$b_average_predation_richness * (x2.sim[i])))/
                             (1+exp(mean(dat_plot$`b_Intercept[2]`) - dat_plot$b_habitat * -1 - dat_plot$`b_average_predation_richness:habitat` * -1 * (x2.sim[i]) - dat_plot$b_average_predation_richness * (x2.sim[i])))))) + 
    (1 - ((exp(mean(dat_plot$`b_Intercept[1]`) - dat_plot$b_habitat * -1 - dat_plot$`b_average_predation_richness:habitat` * -1 * (x2.sim[i]) - dat_plot$b_average_predation_richness * (x2.sim[i])))/
            (1+exp(mean(dat_plot$`b_Intercept[1]`) - dat_plot$b_habitat * -1 - dat_plot$`b_average_predation_richness:habitat` * -1 * (x2.sim[i]) - dat_plot$b_average_predation_richness * (x2.sim[i])))) 
     - (1 - ((exp(mean(dat_plot$`b_Intercept[2]`) - dat_plot$b_habitat * -1 - dat_plot$`b_average_predation_richness:habitat` * -1 * (x2.sim[i]) - dat_plot$b_average_predation_richness * (x2.sim[i])))/
               (1+exp(mean(dat_plot$`b_Intercept[2]`) - dat_plot$b_habitat * -1 - dat_plot$`b_average_predation_richness:habitat` * -1 * (x2.sim[i]) - dat_plot$b_average_predation_richness * (x2.sim[i]))))))
}


bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.025)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.975)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean, bayes.c.eff.lower, bayes.c.eff.upper, bayes.c.eff.lower.bis, bayes.c.eff.upper.bis)



mean_pred <- attributes(cooperation_pred$average_predation_richness)$`scaled:center`
sd_pred <- attributes(cooperation_pred$average_predation_richness)$`scaled:scale`

inv_fun <- function(x){(x*sd_pred)+mean_pred}




cooperation_pred$coopbis <- as.numeric(cooperation_pred$coop)
cooperation_pred$coopbis <- cooperation_pred$coopbis-1

p <- ggplot(plot.dat, aes(x = inv_fun(x2.sim), y = bayes.c.eff.mean)) +
  geom_line(color = "black", alpha = 0.8, size = 2)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "black", alpha = 0.15)+
  ggtitle("-1 SD")+
  xlab("Average predation richness")+
  ylab("")+
  scale_y_continuous(breaks = seq(0, 2, by = 1))+
  coord_cartesian(clip = 'off')+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size=16),
        axis.text.y=element_blank(),
        axis.text.x = element_text(size=14),
        plot.title = element_text(hjust = 0.5))+
  theme(plot.margin = unit(c(1,1,1,9), "lines"))

p <- p + geom_point(data=cooperation_pred, aes(inv_fun(average_predation_richness), (coopbis), alpha=coopbis), color = "black", size = 2.5,
                    position = position_jitter(w=0, h = 0.15))+
  scale_alpha(range = c(0.05, 0.2))+
  theme(legend.position = "none")

p <- p + annotation_custom(
  grob = textGrob(label = "Non-family living", hjust = 1, gp = gpar(fontsize = 14),
                  y = 0.12,
                  x = -0.015))
p <- p + annotation_custom(
  grob = textGrob(label = "Family living", hjust = 1, gp = gpar(fontsize = 14),
                  y = 0.51,
                  x = -0.015))
p <- p + annotation_custom(
  grob = textGrob(label = "Cooperative breeding", hjust = 1, gp = gpar(fontsize = 14),
                  y = 0.91,
                  x = -0.015))
p



#### plot across predation richness for open habitat (i.e habitat set to +1SD)

int.sim <- matrix(rep(NA, nrow(dat_plot)*length(x2.sim)), nrow = nrow(dat_plot))
for(i in 1:length(x2.sim)){
  int.sim[, i] <- (2*(1 - ((exp(mean(dat_plot$`b_Intercept[2]`) - dat_plot$b_habitat * 1 - dat_plot$`b_average_predation_richness:habitat` * 1 * (x2.sim[i]) - dat_plot$b_average_predation_richness * (x2.sim[i])))/
                             (1+exp(mean(dat_plot$`b_Intercept[2]`) - dat_plot$b_habitat * 1 - dat_plot$`b_average_predation_richness:habitat` * 1 * (x2.sim[i]) - dat_plot$b_average_predation_richness * (x2.sim[i])))))) + 
    (1 - ((exp(mean(dat_plot$`b_Intercept[1]`) - dat_plot$b_habitat * 1 - dat_plot$`b_average_predation_richness:habitat` * 1 * (x2.sim[i]) - dat_plot$b_average_predation_richness * (x2.sim[i])))/
            (1+exp(mean(dat_plot$`b_Intercept[1]`) - dat_plot$b_habitat * 1 - dat_plot$`b_average_predation_richness:habitat` * 1 * (x2.sim[i]) - dat_plot$b_average_predation_richness * (x2.sim[i])))) 
     - (1 - ((exp(mean(dat_plot$`b_Intercept[2]`) - dat_plot$b_habitat * 1 - dat_plot$`b_average_predation_richness:habitat` * 1 * (x2.sim[i]) - dat_plot$b_average_predation_richness * (x2.sim[i])))/
               (1+exp(mean(dat_plot$`b_Intercept[2]`) - dat_plot$b_habitat * 1 - dat_plot$`b_average_predation_richness:habitat` * 1 * (x2.sim[i]) - dat_plot$b_average_predation_richness * (x2.sim[i]))))))
}


bayes.c.eff.mean <- apply(int.sim, 2, mean)
bayes.c.eff.lower <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.025)))
bayes.c.eff.upper <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.975)))
bayes.c.eff.lower.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.25)))
bayes.c.eff.upper.bis <- apply(int.sim, 2, function(x) quantile(x, probs = c(0.75)))
plot.dat <- data.frame(x2.sim, bayes.c.eff.mean, bayes.c.eff.lower, bayes.c.eff.upper, bayes.c.eff.lower.bis, bayes.c.eff.upper.bis)


mean_pred <- attributes(cooperation_pred$average_predation_richness)$`scaled:center`
sd_pred <- attributes(cooperation_pred$average_predation_richness)$`scaled:scale`

inv_fun <- function(x){(x*sd_pred)+mean_pred}




cooperation_pred$coopbis <- as.numeric(cooperation_pred$coop)
cooperation_pred$coopbis <- cooperation_pred$coopbis-1

q <- ggplot(plot.dat, aes(x = inv_fun(x2.sim), y = bayes.c.eff.mean)) +
  geom_line(color = "black", alpha = 0.8, size = 2)+
  geom_ribbon(aes(ymin = bayes.c.eff.lower, ymax = bayes.c.eff.upper), fill = "black", alpha = 0.15)+
  ggtitle("+1 SD")+
  xlab("Average predation richness")+
  ylab("")+
  scale_y_continuous(breaks = seq(0, 2, by = 1))+
  coord_cartesian(clip = 'off')+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.title = element_text(size=16),
        axis.text.y=element_blank(),
        axis.text.x = element_text(size=14),
        plot.title = element_text(hjust = 0.5))+
  theme(plot.margin = unit(c(1,1,1,1), "lines"))

q <- q + geom_point(data=cooperation_pred, aes(inv_fun(average_predation_richness), (coopbis), alpha=coopbis), color = "black", size = 2.5,
                    position = position_jitter(w=0, h = 0.15))+
  scale_alpha(range = c(0.05, 0.2))+
  theme(legend.position = "none")

q


p + q



###########################################################################
#### Model without Paleartic/Nearctic/Holarctic and Widespread species ####
###########################################################################

cooperation_pred <- read.csv("~/My Drive/predation_cooperation/data/final_datasets/cooperation_pred.csv", row.names=1)

# create a social system variable, as an ordered factor
cooperation_pred$coop <- factor(cooperation_pred$fam_sys_known50, ordered = T, levels = c("no_fam", "family", "coop_families"))


## Subset dataset and phylogeny

cooperation_pred <- cooperation_pred[cooperation_pred$region!="Holarctic",]
cooperation_pred <- cooperation_pred[cooperation_pred$region!="Nearctic",]
cooperation_pred <- cooperation_pred[cooperation_pred$region!="Palearctic",]
cooperation_pred <- cooperation_pred[cooperation_pred$region!="Widespread",]
cooperation_pred <- cooperation_pred[cooperation_pred$region!="Nearctic, Neotropical",]

cooperation_pred <- cooperation_pred[cooperation_pred$tip_label != "Melanodryas_vittata",]
cooperation_pred <- cooperation_pred[cooperation_pred$tip_label != "Tangara_cyanoventris",]
cooperation_pred <- cooperation_pred[cooperation_pred$tip_label != "Fulica_ardesiaca",]
cooperation_pred <- cooperation_pred[cooperation_pred$tip_label != "Tangara_larvata",]

MyTree <- read.tree("~/My Drive/predation_cooperation/data/Prum_Jetz_Cooney_9993.tree") 

tree_phylo <- maxCladeCred(MyTree) # get maximum clade credibility tree

setdiff(cooperation_pred$tip_label, tree_phylo$tip.label) # check that all species are in the phylogeny

tree_phylo <- keep.tip(tree_phylo, cooperation_pred$tip_label) # drop tip from phylo: remove the species not present in the dataset

phylo_mcct <- vcv.phylo(tree_phylo) # transform into a variance covariance matrix

cooperation_pred$phylo <- cooperation_pred$tip_label


sd_predation <- sd(cooperation_pred$average_predation_richness)
mean_predation <- mean(cooperation_pred$average_predation_richness)
cooperation_pred$average_predation_richness <- scale(cooperation_pred$average_predation_richness)

cooperation_pred$log_mass <- log(cooperation_pred$mass)
sd_logmass <- sd(cooperation_pred$log_mass)
mean_logmass <- mean(cooperation_pred$log_mass)
cooperation_pred$log_mass <- scale(cooperation_pred$log_mass)

sd_habitat <- sd(cooperation_pred$habitat)
mean_habitat <- mean(cooperation_pred$habitat)
cooperation_pred$habitat <- scale(cooperation_pred$habitat)

cooperation_pred$latitude_mean <- abs(cooperation_pred$latitude_mean)
sd_latitude <- sd(cooperation_pred$latitude_mean)
mean_latitude <- mean(cooperation_pred$latitude_mean)
cooperation_pred$latitude_mean <- scale(cooperation_pred$latitude_mean)



model_no_holarctic <- brm(
  coop ~ average_predation_richness * habitat +
    log_mass +
    latitude_mean +
    scale(Prcp.P) + scale(Temp.P) + 
    scale(Prcp.Mean) + scale(Temp.Mean) +
    scale(Prcp.Var) + scale(Temp.Var) +
    (1|gr(phylo, cov = phylo_mcct)), data = cooperation_pred, 
  family = cumulative("logit"), data2 = list(phylo_mcct = phylo_mcct),
  warmup = 1000, iter = 2000, thin = 1, chains = 3, seed = 25, cores = 3,
  prior = c(
    set_prior("normal(0,2)", class = "b"),
    set_prior("normal(0,2)", class = "Intercept"),
    set_prior("normal(0,1.5)", class = "sd", coef = "Intercept", group = "phylo")),
  control = list(adapt_delta = 0.99), save_ranef = F, backend = "cmdstanr")

summary(model_no_holarctic)
