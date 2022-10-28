## plot BD rates

library(RevGadgets)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# read the trace files 

# read the trace files from by different priors

BD_lognormal <- readTrace("./output/ALL50K/ALL50K_BD_logNormal.log")[[1]]

BD_Beta_18_2 <- readTrace("./output/ALL50K/ALL50K_BD_T_B18_2.log")[[1]]

#transform the data frame for ploting

BD_lognormal <- BD_lognormal %>% select(c(birth_rate, death_rate, diversification)) %>% mutate(PriorType="logNormal")

BD_Beta_18_2 <- BD_Beta_18_2 %>% select(c(birth_rate, death_rate, diversification)) %>% mutate(PriorType="Beta_18_2")

# bind dataframes

#combined_BD <- rbind(BD_lognormal,BD_Beta_18_2)

#BD_birth_rate <- ggplot(combined_BD, aes(x=birth_rate, fill=PriorType)) + geom_density(alpha=0.5) + xlim(0.0, 0.20) + xlab("Birth Rate") + ylab("Density") + scale_fill_discrete(name = "Prior distribtuions")

#BD_death_rate <-  ggplot(combined_BD, aes(x=death_rate, fill=PriorType)) + geom_density(alpha=0.5) + xlim(0.0, 0.20) + xlab("Death Rate") + ylab("Density") + scale_fill_discrete(name = "Prior distribtuions")

#BD_div_rate <- ggplot(combined_BD, aes(x=diversification, fill=PriorType)) + geom_density(alpha=0.5) + xlim(0.0, 0.20) + xlab("Diversification Rate") + ylab("Density") + scale_fill_discrete(name = "Prior distribtuions")

# summary stats



summary(BD_lognormal[,-4]) 

summary(BD_Beta_18_2[,-4])

# t-test?

t.test(BD_lognormal$birth_rate, BD_Beta_18_2$birth_rate)

t.test(BD_lognormal$death_rate, BD_Beta_18_2$death_rate)

t.test(BD_lognormal$diversification, BD_Beta_18_2$diversification)


# version 2 

# Transform the lognormal and beta prior data  

log_transformed_bd_diversification <- BD_lognormal %>% select(birth_rate) %>% mutate(type="Birth_rate") %>% rename(Value=birth_rate)

log_transformed_bd_death <- BD_lognormal %>% select(death_rate) %>% mutate(type="Death_rate") %>% rename(Value=death_rate)

log_transformed_bd_birth <- BD_lognormal %>% select(diversification) %>% mutate(type="Diversification") %>% rename(Value=diversification)

# be aware of birth rate vs. diversification replacement 

beta_transformed_bd_birth <- BD_Beta_18_2 %>% select(birth_rate) %>% mutate(type="Birth_rate") %>% rename(Value=birth_rate)

beta_transformed_bd_death <- BD_Beta_18_2 %>% select(death_rate) %>% mutate(type="Death_rate") %>% rename(Value=death_rate)

beta_transformed_bd_diversification <- BD_Beta_18_2 %>% select(diversification) %>% mutate(type="Diversification") %>% rename(Value=diversification)


# combine datasets lognormal

log_transformed_all_BD <- rbind(log_transformed_bd_birth,log_transformed_bd_death,log_transformed_bd_diversification)

beta_transformed_all_BD <- rbind(beta_transformed_bd_birth,beta_transformed_bd_death,beta_transformed_bd_diversification)


#plot 

log_BD_plot <- ggplot(log_transformed_all_BD, aes(x=Value, fill=type)) + geom_density(alpha=0.5) + xlim(0,0.15) + xlab("Value") + ylab("Density") + scale_fill_discrete(name = "Parameters") + ggtitle("Estimated parameters under logNormal prior dist.")

beta_BD_plot <- ggplot(beta_transformed_all_BD, aes(x=Value, fill=type)) + geom_density(alpha=0.5) + xlim(0,0.15) + xlab("Value") + ylab("Density") + scale_fill_discrete(name = "Parameters") + ggtitle("Estimated parameters under Beta(18,2) prior dist.")

# plot together 

BD_together <- ggarrange(log_BD_plot, beta_BD_plot, labels = c("A", "B"), ncol = 1, nrow = 2,
  common.legend = TRUE, legend="right")


BD_together

# compute bayesian factors 




candidateModels <- c("ConstBD"=marginalLikelihoodConstBD,
                     "DecrBD"=marginalLikelihoodDecrBD,
                     "EpisodicBD"=marginalLikelihoodEpisodicBD,
                     "MassExtinctionBD"=marginalLikelihoodMassExtinctionBD)
