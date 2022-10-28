################################################################################
#
# Plot estimates of the Yule model
#
# authors: Sebastian Höhna
#
################################################################################

library(RevGadgets)
library(ggplot2)
library(dplyr)
library(data.table)


# read the posterior and prior first_output

ALL50K_yule_posterior <- readTrace("first_output/ALL50K/ALL50KTree_Yule.log")[[1]] %>%
  mutate(birth_rate_prior=readTrace("first_output/ALL50K/ALL50KTree_Yule_prior.log")[[1]]$birth_rate) %>%
  mutate(DataType="ALL50K")



intron_yule_posterior <- readTrace("first_output/intron/intronTree_Yule.log")[[1]] %>%
  mutate(birth_rate_prior=readTrace("first_output/intron/intronTree_Yule_prior.log")[[1]]$birth_rate) %>%
  mutate(DataType="intron")

CDS_yule_posterior <- readTrace("first_output/CDS/CDSTree_Yule.log")[[1]] %>%
  mutate(birth_rate_prior=readTrace("first_output/CDS/CDSTree_Yule_prior.log")[[1]]$birth_rate) %>%
  mutate(DataType="CDS")

lnc_RNA_yule_posterior <- readTrace("first_output/lnc_RNA/lnc_RNATree_Yule.log")[[1]] %>%
  mutate(birth_rate_prior=readTrace("first_output/lnc_RNA/lnc_RNATree_Yule_prior.log")[[1]]$birth_rate) %>%
  mutate(DataType="lncRNA")

other_yule_posterior <- readTrace("first_output/other/otherTree_Yule.log")[[1]] %>%
  mutate(birth_rate_prior=readTrace("first_output/other/otherTree_Yule_prior.log")[[1]]$birth_rate) %>%
  mutate(DataType="other")

pseudogene_yule_posterior <- readTrace("first_output/pseudogene/pseudogeneTree_Yule.log")[[1]] %>%
  mutate(birth_rate_prior=readTrace("first_output/pseudogene/pseudogeneTree_Yule_prior.log")[[1]]$birth_rate) %>%
  mutate(DataType="pseudogene")

unannotated_yule_posterior <- readTrace("first_output/unannotated/unannotatedTree_Yule.log")[[1]] %>%
  mutate(birth_rate_prior=readTrace("first_output/unannotated/unannotatedTree_Yule_prior.log")[[1]]$birth_rate) %>%
  mutate(DataType="unannotated")


##combine data frames

combined_yule <- rbind(ALL50K_yule_posterior, intron_yule_posterior, CDS_yule_posterior,
      lnc_RNA_yule_posterior, other_yule_posterior,  pseudogene_yule_posterior,
      unannotated_yule_posterior)


## plot the yule rates

subset_combined_yule <- combined_yule[,6:8]

melted_yule <- melt(subset_combined_yule)

ggplot(melted_yule, aes(x=value, fill=variable)) + geom_density(adjust=1.5, alpha=.4) + facet_wrap(.~DataType) + theme(legend.position = c(0.80,0.80))


# plot the prior vs the posterior
plot_ALL50K <- plotTrace(list(ALL50K_yule_posterior), vars=c("birth_rate", "birth_rate_prior"))[[1]]  +
  # modify legend location using ggplot2
  

plot

ggsave("ALL50K_HKY85_birth_rate_prior_posterior.png", plot, height=5, width=5)


mean(yule_posterior$birth_rate)
mean(yule_posterior$birth_rate_prior)

z_scores <-  


