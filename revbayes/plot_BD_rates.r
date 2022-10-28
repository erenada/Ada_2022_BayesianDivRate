################################################################################
#
# Plot estimates of the birth-death model
#
# authors: Sebastian Höhna
#
################################################################################

#sp_info <- read.csv("/Users/eren/Documents/GitHub/Bird_SISRS/SUPPLEMENTS/Species_Info.csv")

library(RevGadgets)
library(ggplot2)

# read the posterior and prior output
ALL50K_bd_posterior <- readTrace("output/ALL50K/ALL50K_BD.log")
intron_bd_posterior <- readTrace("output/intron/intron_BD.log")
CDS_bd_posterior <- readTrace("output/CDS/CDS_BD.log")
other_bd_posterior <- readTrace("output/other/other_BD.log")
pseudogene_bd_posterior <- readTrace("output/pseudogene/pseudogene_BD.log")
lnc_RNA_bd_posterior <- readTrace("output/lnc_RNA/lnc_RNA_BD.log")
UTR_bd_posterior <- readTrace("output/UTR/UTR_BD.log")




# plot the prior vs the posterior
plot <- plotTrace(ALL50K_bd_posterior, vars=c("birth_rate", "death_rate"))[[1]]  +
     # modify legend location using ggplot2
     theme(legend.position = c(0.80,0.80))

plot

ggsave("ALL50K_birth_death_rate.png", plot, height=5, width=5)

mean(bd_posterior[[1]]$birth_rate)

mean(bd_posterior[[1]]$death_rate)




