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
Cal2_intron_chronos <- root(read.tree("../RAxML_br/RAxML_result.intron"), outgroup = neoaves_outgroup, resolve.root = T)

plot(Cal2_intron_chronos)

#ggplot(Cal2_intron_chronos, aes(x, y)) + geom_tree() + theme_tree()


#fix the root

#Cal2_intron_chronos$edge.length[which(!(Cal2_intron_chronos$edge[,1] %in% Cal2_intron_chronos$edge[,2]))] <- sum(Cal2_intron_chronos$edge.length[which(!(Cal2_intron_chronos$edge[,1] %in% Cal2_intron_chronos$edge[,2]))])/2

#plot.phylo(Cal2_intron_chronos, use.edge.length = F)

#nodelabels()

## 

#ggtree(Cal2_intron_chronos, ladderize=T, branch.length="none")


#calibration start


Cal2_intron_cal <- add_row(Cal2_intron_cal <- makeChronosCalib(ALL50K_chronos,age.min = neoaves_min, age.max = neoaves_max), node = c(34,48,40), age.min=c(66,55.8,49), age.max = c(70,55.8,52.1), soft.bounds=FALSE)


#get edge and branch lists

Cal2_intron_edge_list <- compute.brlen(Cal2_intron_chronos)$edge.length
Cal2_intron_branch_list <- branching.times(compute.brlen(Cal2_intron_chronos))

#repeat

for(i in 1:10){
  Cal2_intron_chronos_iter <- chronos(Cal2_intron_chronos,calibration = Cal2_intron_cal, model = "discrete")
  Cal2_intron_edge_list <- rbind(Cal2_intron_edge_list, Cal2_intron_chronos_iter$edge.length)
  Cal2_intron_branch_list <- rbind(Cal2_intron_branch_list, branching.times(Cal2_intron_chronos_iter))
}

#plot
plot(Cal2_intron_chronos_iter)

#edge list

Cal2_intron_edge_list <- Cal2_intron_edge_list[-1,]
Cal2_intron_branch_list <- Cal2_intron_branch_list[-1,]

Cal2_intron_branches <- colnames(Cal2_intron_branch_list)

# median

Cal2_intron_medians <- c()
for(i in 1:ncol(Cal2_intron_branch_list)){
  Cal2_intron_medians[colnames(Cal2_intron_branch_list)[i]] <- median(Cal2_intron_branch_list[,i])
}


Cal2_intron_timetree <- compute.brlen(Cal2_intron_chronos)

Cal2_intron_median_edge <- c()

for(i in 1:ncol(Cal2_intron_edge_list)){
  Cal2_intron_median_edge<-c(Cal2_intron_median_edge,median(Cal2_intron_edge_list[,i]))
}
Cal2_intron_timetree$edge.length <- Cal2_intron_median_edge
if(!is.ultrametric(Cal2_intron_timetree)){
  Cal2_intron_timetree <- force.ultrametric(Cal2_intron_timetree,"extend")
}

plot(Cal2_intron_timetree)


write.tree(Cal2_intron_timetree, file = "Cal2_intron_TT_1CAL.tree")
