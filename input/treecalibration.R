## Tree calibration

## libraries

library(plyr)
library(tidyverse)
library(ape)
library(phytools)
library(jtools)
library(scales)
library(treeio)
library(ggtree)

## time constrains

neoaves_min <- 66.0
neoaves_max <- 99.6 
neoaves_outgroup <- c('StrCam')

#Read Trees

ALL_50K_chronos <- root(read.tree("../RAxML_br/RAxML_result.ALL50K"), outgroup = neoaves_outgroup)

plot(ALL_50K_chronos)

ggplot(ALL_50K_chronos, aes(x, y)) + geom_tree() + theme_tree()

# calibration

ALL_50K_cal <- makeChronosCalib(ALL_50K_chronos,age.min = neoaves_min,age.max = neoaves_max, node = 34)

ALL_50K_edge_list <- compute.brlen(ALL_50K_chronos)$edge.length
ALL_50K_branch_list <- branching.times(compute.brlen(ALL_50K_chronos))  


for(i in 1:10){
  ALL_50K_chronos_iter <- chronos(ALL_50K_chronos,calibration = ALL_50K_cal)
  ALL_50K_edge_list <- rbind(ALL_50K_edge_list, ALL_50K_chronos_iter$edge.length)
  ALL_50K_branch_list <- rbind(ALL_50K_branch_list, branching.times(ALL_50K_chronos_iter))
}




ALL_50K_alex <- chronos(ALL_50K_chronos,calibration = ALL_50K_cal, model = "discrete")


plot(ALL_50K_chronos)
plot(ALL_50K_alex)


#write.tree(ALL_50K_chronos_iter, file = "ALL_50Kiteration.tre")

# extra

ALL_50K_edge_list <- ALL_50K_edge_list[-1,]
ALL_50K_branch_list <- ALL_50K_branch_list[-1,]

ALL_50K_branches <- colnames(ALL_50K_branch_list)

#median ages

ALL_50K_medians <- c()
for(i in 1:ncol(ALL_50K_branch_list)){
  ALL_50K_medians[colnames(ALL_50K_branch_list)[i]] <- median(ALL_50K_branch_list[,i])
}

ALL_50K_timetree2 <- compute.brlen(ALL_50K_chronos)

ALL_50K_alex <- compute.brlen(ALL_50K_chronos)

plot(ALL_50K_alex)

ALL_50K_median_edge <- c()

for(i in 1:ncol(ALL_50K_edge_list)){
  ALL_50K_median_edge<-c(ALL_50K_median_edge,median(ALL_50K_edge_list[,i]))
}

ALL_50K_timetree2$edge.length <- ALL_50K_median_edge

if(!is.ultrametric(ALL_50K_timetree2)){
  ALL_50K_timetree2 <- force.ultrametric(ALL_50K_timetree2,"extend")
}

write.tree(ALL_50K_timetree2, file = "ALL_50KTimeTree.tre")

plot(ALL_50K_chronos)

