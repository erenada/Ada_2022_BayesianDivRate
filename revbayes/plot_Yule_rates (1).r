
library(RevGadgets)
library(ggplot2)
library(dplyr)
library(tidyr)

# read the trace files from by different priors

yule_lognormal <- readTrace("./first_output/ALL50K/All50KTree_Yule.log")[[1]]

#yule_Beta_2_18 <- readTrace("./output/ALL50K/ALL50KTree_YuleBeta2-18.log")[[1]]

yule_Beta_18_2 <- readTrace("./first_output/ALL50K/ALL50K_BD_Beta.log")[[1]]


# summary stats

sd(yule_Beta_18_2$birth_rate)

sd(yule_lognormal$birth_rate)



# t-tests

t.test(yule_lognormal$birth_rate, yule_Beta_18_2$birth_rate)


#transform the data frame for ploting

ALL50K_yule_lognormal <- readTrace("./output/071222/")[[1]] %>% select(birth_rate) %>% mutate(PriorType="logNormal")

CDS_yule_lognormal <- readTrace("./first_output/CDS/CDSTree_Yule.log")[[1]] %>% select(birth_rate) %>% mutate(PriorType="logNormal")

intron_yule_lognormal <- readTrace("./first_output/intron/intronTree_Yule.log")[[1]] %>% select(birth_rate) %>% mutate(PriorType="logNormal")

lncRNA_yule_lognormal <- readTrace("./first_output/lnc_RNA/lnc_RNATree_Yule.log")[[1]] %>% select(birth_rate) %>% mutate(PriorType="logNormal")

other_yule_lognormal <- readTrace("./first_output/other/otherTree_Yule.log")[[1]] %>% select(birth_rate) %>% mutate(PriorType="logNormal")

unannotated_yule_lognormal <- readTrace("./first_output/unannotated/unannotatedTree_Yule.log")[[1]] %>% select(birth_rate) %>% mutate(PriorType="logNormal")

UTR_yule_lognormal <- readTrace("./first_output/UTR/UTRTree_Yule.log")[[1]] %>% select(birth_rate) %>% mutate(PriorType="logNormal")

pseudogene_yule_lognormal <- readTrace("./first_output/pseudogene/pseudogeneTree_Yule.log")[[1]] %>% select(birth_rate) %>% mutate(PriorType="logNormal")


#Beta

ALL50K_yule_Beta_2_18 <- readTrace("./first_output/ALL50K/ALL50K_BD_Beta.log")[[1]] %>% select(birth_rate) %>% mutate(PriorType="Beta_2_18")
CDS_yule_Beta_2_18 <- readTrace("./first_output/CDS/CDS_BD_Beta.log")[[1]] %>% select(birth_rate) %>% mutate(PriorType="Beta_2_18")
intron_yule_Beta_2_18 <- readTrace("./first_output/intron/intron_BD_Beta.log")[[1]] %>% select(birth_rate) %>% mutate(PriorType="Beta_2_18")
lncRNA_yule_Beta_2_18 <- readTrace("./first_output/lnc_RNA/lnc_RNA_BD_Beta.log")[[1]] %>% select(birth_rate) %>% mutate(PriorType="Beta_2_18")
other_yule_Beta_2_18 <- readTrace("./first_output/other/other_BD_Beta.log")[[1]] %>% select(birth_rate) %>% mutate(PriorType="Beta_2_18")
unannotated_yule_Beta_2_18 <- readTrace("./first_output/unannotated/unannotated_BD_Beta.log")[[1]] %>% select(birth_rate) %>% mutate(PriorType="Beta_2_18")
UTR_yule_Beta_2_18 <- readTrace("./first_output/UTR/UTR_BD_Beta.log")[[1]] %>% select(birth_rate) %>% mutate(PriorType="Beta_2_18")
pseudogene_yule_Beta_2_18 <- readTrace("./first_output/pseudogene/pseudogene_BD_Beta.log")[[1]] %>% select(birth_rate) %>% mutate(PriorType="Beta_2_18")



#yule_Beta_18_2 <- readTrace("./first_output/ALL50K/ALL50K_BD_Beta.log")[[1]] %>% select(birth_rate) %>% mutate(PriorType="Beta_18_2")


min(yule_lognormal$birth_rate)

max(yule_lognormal$birth_rate)

median(yule_lognormal$birth_rate)

mean(yule_lognormal$birth_rate)

sd(yule_lognormal$birth_rate)

# bind dataframes
ALL50K_combined_yule <- rbind(ALL50K_yule_lognormal,ALL50K_yule_Beta_2_18)
CDS_combined_yule <- rbind(CDS_yule_lognormal,CDS_yule_Beta_2_18)
intron_combined_yule <- rbind(intron_yule_lognormal,intron_yule_Beta_2_18)
lncRNA_combined_yule <- rbind(lncRNA_yule_lognormal,lncRNA_yule_Beta_2_18)
other_combined_yule <- rbind(other_yule_lognormal,other_yule_Beta_2_18)
unannotated_combined_yule <- rbind(unannotated_yule_lognormal,unannotated_yule_Beta_2_18)
UTR_combined_yule <- rbind(UTR_yule_lognormal,UTR_yule_Beta_2_18)
pseudogene_combined_yule <- rbind(UTR_yule_lognormal,UTR_yule_Beta_2_18)



ALL50K_plot <-ggplot(ALL50K_combined_yule, aes(x=birth_rate, fill=PriorType)) + geom_density(alpha=0.5) + xlim(0.0, 0.15) + xlab("Birth Rate") + ylab("Density") + scale_fill_discrete(name = "Prior distribtuions")

CDS_plot <-ggplot(CDS_combined_yule, aes(x=birth_rate, fill=PriorType)) + geom_density(alpha=0.5) + xlim(0.0, 0.15) + xlab("Birth Rate") + ylab("Density") + scale_fill_discrete(name = "Prior distribtuions")

intron_plot <-ggplot(intron_combined_yule, aes(x=birth_rate, fill=PriorType)) + geom_density(alpha=0.5) + xlim(0.0, 0.15) + xlab("Birth Rate") + ylab("Density") + scale_fill_discrete(name = "Prior distribtuions")

lncRNA_plot <-ggplot(lncRNA_combined_yule, aes(x=birth_rate, fill=PriorType)) + geom_density(alpha=0.5) + xlim(0.0, 0.15) + xlab("Birth Rate") + ylab("Density") + scale_fill_discrete(name = "Prior distribtuions")

other_plot <-ggplot(other_combined_yule, aes(x=birth_rate, fill=PriorType)) + geom_density(alpha=0.5) + xlim(0.0, 0.15) + xlab("Birth Rate") + ylab("Density") + scale_fill_discrete(name = "Prior distribtuions")

unannotated_plot <-ggplot(unannotated_combined_yule, aes(x=birth_rate, fill=PriorType)) + geom_density(alpha=0.5) + xlim(0.0, 0.15) + xlab("Birth Rate") + ylab("Density") + scale_fill_discrete(name = "Prior distribtuions")

UTR_plot <-ggplot(UTR_combined_yule, aes(x=birth_rate, fill=PriorType)) + geom_density(alpha=0.5) + xlim(0.0, 0.15) + xlab("Birth Rate") + ylab("Density") + scale_fill_discrete(name = "Prior distribtuions")

pseudogene_plot <-ggplot(pseudogene_combined_yule, aes(x=birth_rate, fill=PriorType)) + geom_density(alpha=0.5) + xlim(0.0, 0.15) + xlab("Birth Rate") + ylab("Density") + scale_fill_discrete(name = "Prior distribtuions")


#  ggarrange()

library(gridExtra)

grid.arrange(ALL50K_plot,CDS_plot, intron_plot,lncRNA_plot,other_plot,unannotated_plot,UTR_plot,pseudogene_plot, ncol = 4, nrow = 2)
