## ALL50K ROOT + Neognathae + Galliformes CALIBRATION


library(tidyverse)
library(ape)
library(phytools)
library(treeio)
library(ggtree)


neoaves_min <- 66.0
neoaves_max <- 99.6 
neoaves_outgroup <- c('StrCam')



ALL50K_chronos <- root(read.tree("../RAxML_br/HKY85/RAxML_result.ALL50K"), outgroup = "StrCam", resolve.root = T)

plot(ALL50K_chronos)

#ggplot(ALL50K_chronos, aes(x, y)) + geom_tree() + theme_tree()

#fix the root

ALL50K_chronos$edge.length[which(!(ALL50K_chronos$edge[,1] %in% ALL50K_chronos$edge[,2]))] <- sum(ALL50K_chronos$edge.length[which(!(ALL50K_chronos$edge[,1] %in% ALL50K_chronos$edge[,2]))])/2

#plot.phylo(ALL50K_chronos, use.edge.length = F)

nodelabels()

#ggtree(ALL50K_chronos, ladderize=T, branch.length="none")

#calibration start


ALL50K_cal <- add_row(makeChronosCalib(ALL50K_chronos,age.min = neoaves_min, age.max = neoaves_max, soft.bounds=FALSE), node = c(35,34), age.min = c(66,50.03), age.max=c(70,50.03),soft.bounds = FALSE)


#get edge and branch lists

ALL50K_edge_list <- compute.brlen(ALL50K_chronos)$edge.length
ALL50K_branch_list <- branching.times(compute.brlen(ALL50K_chronos))

#repeat

for(i in 1:100){
  ALL50K_chronos_iter <- chronos(ALL50K_chronos,calibration = ALL50K_cal, model = "clock")
  ALL50K_edge_list <- rbind(ALL50K_edge_list, ALL50K_chronos_iter$edge.length)
  ALL50K_branch_list <- rbind(ALL50K_branch_list, branching.times(ALL50K_chronos_iter))
}

plot(ALL50K_chronos_iter)

#choosing median ages


ALL50K_edge_list <- ALL50K_edge_list[-1,]
ALL50K_branch_list <- ALL50K_branch_list[-1,]

ALL50K_branches <- colnames(ALL50K_branch_list)

ALL50K_medians <- c()
for(i in 1:ncol(ALL50K_branch_list)){
  ALL50K_medians[colnames(ALL50K_branch_list)[i]] <- median(ALL50K_branch_list[,i])
}

# 
ALL50K_timetree <- compute.brlen(ALL50K_chronos)

is.ultrametric(ALL50K_timetree)
ALL50K_median_edge <- c()
for(i in 1:ncol(ALL50K_edge_list)){
  ALL50K_median_edge<-c(ALL50K_median_edge,median(ALL50K_edge_list[,i]))
}
ALL50K_timetree$edge.length <- ALL50K_median_edge

plot(ALL50K_timetree)

write.tree(ALL50K_timetree, file = "../TimeTrees/HKY85/ALL50K_medianTT.tree")
