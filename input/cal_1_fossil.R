##Calib script on 04/12


library(tidyverse)
library(ape)
library(phytools)
library(treeio)
library(ggtree)

## time constrains

neoaves_min <- 66.0
neoaves_max <- 99.6 
neoaves_outgroup <- c('StrCam')

#read tree
intron_chronos <- root(read.tree("../RAxML_br/RAxML_result.intron"), outgroup = neoaves_outgroup, resolve.root = T)

plot(intron_chronos)

#ggplot(intron_chronos, aes(x, y)) + geom_tree() + theme_tree()


#fix the root

#intron_chronos$edge.length[which(!(intron_chronos$edge[,1] %in% intron_chronos$edge[,2]))] <- sum(intron_chronos$edge.length[which(!(intron_chronos$edge[,1] %in% intron_chronos$edge[,2]))])/2

#plot.phylo(intron_chronos, use.edge.length = F)

#nodelabels()

## 

#ggtree(intron_chronos, ladderize=T, branch.length="none")


#calibration start


 intron_cal <- makeChronosCalib(intron_chronos,age.min = neoaves_min, age.max = neoaves_max, soft.bounds=FALSE)


#get edge and branch lists

intron_edge_list <- compute.brlen(intron_chronos)$edge.length
intron_branch_list <- branching.times(compute.brlen(intron_chronos))

#repeat

nodelabels()

for(i in 1:10){
  intron_chronos_iter <- chronos(intron_chronos,calibration = intron_cal, model = "discrete")
  intron_edge_list <- rbind(intron_edge_list, intron_chronos_iter$edge.length)
  intron_branch_list <- rbind(intron_branch_list, branching.times(intron_chronos_iter))
}

#plot
plot(intron_chronos_iter)

#edge list

intron_edge_list <- intron_edge_list[-1,]
intron_branch_list <- intron_branch_list[-1,]

intron_branches <- colnames(intron_branch_list)

# median

intron_medians <- c()
for(i in 1:ncol(intron_branch_list)){
  intron_medians[colnames(intron_branch_list)[i]] <- median(intron_branch_list[,i])
}


intron_timetree <- compute.brlen(intron_chronos)

intron_median_edge <- c()

for(i in 1:ncol(intron_edge_list)){
  intron_median_edge<-c(intron_median_edge,median(intron_edge_list[,i]))
}
intron_timetree$edge.length <- intron_median_edge
if(!is.ultrametric(intron_timetree)){
  intron_timetree <- force.ultrametric(intron_timetree,"extend")
}

plot(intron_timetree, type = "cladogram")


write.tree(intron_timetree, file = "intron_TT_1CAL.tree")
